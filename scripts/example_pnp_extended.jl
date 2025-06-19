## Example script demonstrating the use of TSSOS for perspective-n-point
# Extended version which shows more variable manipulation, custom solvers, etc.
# Many features are only available in my fork
# 
# Lorenzo Shaikewitz, 6/19/2025

## Imports
using LinearAlgebra, Random
using JuMP, Dualization
using Clarabel
using TSSOS, DynamicPolynomials
using TutorialTSSOS

## Data Generation
Random.seed!(0) # fix the RNG

N = 10           # number of points
σ = ones(N)*0.05 # covariance of each measurement [m]
camK = [572.41  0.0     325.26;
        0.0     573.57  242.05;
        0 0 1]    # camera calibration matrix (from LM-O dataset)

y, b, gt = genPnPproblem(N, σ, camK)


## Use TSSOS to relax the problem to an SDP
# We'll write the problem as a QCQP and use TSSOS to relax it to SDP
# this is a very different workflow than JuMP!

# 1) create variables (DynamicPolynomials)
#    Here we show an example where `t` is analytically eliminated as a function of `R`
@polyvar R[1:3,1:3]
vars = vec(R) # need to vectorize any matrix vars

# Analytically eliminate t as a function of R
e3 = [0;0;1]
U = Array{Any}(undef, N)
for i = 1:N
    U[i] = (I - y[:,i]*e3')*camK
end
H = sum([U[i]'*U[i] / (σ[i]^2) for i = 1:N])
t = -inv(H)*sum([U[i]'*U[i]*R*b[:,i] / (σ[i]^2) for i = 1:N])

# 2) Build objective
# express this in polynomial form
obj = 0.
for i = 1:N
    global obj
    residual_i = U[i]*(R*b[:,i] + t)
    obj += residual_i'*residual_i / σ[i]^2
end

# 3) Build constraints
# We write these as lists of polynomials with float coefficients (Poly{Float64})
# expr ≥ 0
ineq = Vector{TSSOS.Poly{Float64}}()
# expr = 0
eq = Vector{TSSOS.Poly{Float64}}()

# impose chirality constraints
for i = 1:N
    # note the use of brackets around a scalar term with `append!`
    append!(ineq, [e3'*camK*(R*b[:,i] + t)]) # ≥ 0
    # we could also have used `push!` without the brackets:
    # push!(ineq, e3'*camK*(R*b[:,i] + t))
end

# SO(3) equality constraints
# orthogonality
append!(eq, vec(R'*R - I))
# right hand rule
append!(eq, R[1:3,3] .- cross(R[1:3,1],R[1:3,2]))
append!(eq, R[1:3,1] .- cross(R[1:3,2],R[1:3,3]))
append!(eq, R[1:3,2] .- cross(R[1:3,3],R[1:3,1]))

# 4) Solve
pop = [obj; ineq; eq]
order = 2 # relaxation order
# use TSSOS only to reformulate the problem (`solve=false`)
opt, sol, data, gap, model = cs_tssos_first(pop, vars, order, numeq=length(eq), TS="block", CS="MF", QUIET=false, solve=false, solution=false, 
                                            MomentOne=true)
# to tell TSSOS not to solve, we set `solve=false` and `solution=false`.
# `MomentOne=true` is needed to extrac the solution.
# note that we also have `gap` and `model` returned.
# `opt`, `sol`, and `gap` are all `nothing` since there was no solution.
# Warning: other TS/CS settings may not work with returning the model / solve=false.

# now we have a JuMP model, and we can set a custom optimizer.
# Let's use Clarabel
set_optimizer(model, Clarabel.Optimizer)

# before optimizing, we could add additional (convex) terms
# if desired. None are needed for this polynomial problem

# we could also print the model to inspect it. We need to dualize it
# to get a minimization problem
print(dualize(model))

# solve the model with JuMP
optimize!(model)

# since we are solving on our own, we need to check optimality
# and round the solution ourselves. The following code is mostly
# copied from `cs_tssos_first`

# status of SDP solve
SDP_status = termination_status(model)
opt = objective_value(model)

# extract rank 1 solution
measure = -dual(model[:con])
momone = get_moment(measure, data.tsupp, data.cliques, data.cql, data.cliquesize, nb=data.nb)
sol,gap,data.flag = TSSOS.approx_sol(momone, opt, data.n, data.cliques, data.cql, data.cliquesize, data.supp, data.coe, numeq=data.numeq, gtol=data.gtol, ftol=data.ftol, QUIET=false)

# optionally, refine the solution using IPOPT
# if data.flag == 1
#     sol = gap > 0.5 ? randn(data.n) : sol
#     sol,data.flag = TSSOS.refine_sol(opt, sol, data, QUIET=false, gtol=data.gtol)
# end

# Higher order relaxation
# opt, sol, data, gap, model = cs_tssos_higher!(data; QUIET=false, solve=false, solution=false, MomentOne=true)
# can now re-optimize with the new model


## Solution processing
# extract the solution
R_val = (vars=>sol) .|> R
t_val = (vars=>sol) .|> t
# this syntax is equivalent to `t[i](vars=>sol)`
# each entry of the vector `t` is a function.

# inspect the constraints
ineq_val = (vars=>sol) .|> ineq
eq_val = (vars=>sol) .|> eq

println("Solution satisfies $(sum(ineq_val .≥ 0)) inequality constraints and $(sum(abs.(eq_val) .≤ 1e-4)) equality constraints")