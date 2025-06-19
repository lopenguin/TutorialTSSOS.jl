## Example script demonstrating the use of TSSOS for perspective-n-point
# This simple script works with the standard version of TSSOS
# Run this in the REPL with 
#    `include("scripts/example_pnp_simple.jl")`
# 
# Lorenzo Shaikewitz, 6/19/2025

# Going further: 
# 1) The first time you run this it will be quite slow. Try running it again!
#    You should see a significant speedup--this is the magic of Julia precompilation.
# 2) Try changing the relaxation order from 2 to 1. Notice that we lose tightness.
#    (Aside: TSSOS automatically calls IPOPT on the QCQP to refine a non-tight solution;
#     for this set of constraints, IPOPT always fails).
#    - After solving with order 1, uncomment the higher order relaxation. Notice the 
#      higher order solution is tight! This functionality makes it easy to step up along
#      the hierarchy when tightness is not satisfied.
# 3) Take a look at the returned moment (`data.moment`). These is the optimal values of the SDP variables.
#    TSSOS is doing a term sparse relaxation, so these may not be immediately interpretable.
#    For interpretablility, try changing `order=1` and `CS=false`. This turns off the correlative sparsity,
#    now `data.moment` is just a 13 x 13 matrix corresponding to [1;vec(R);t]*[1;vec(R);t]' at optimality.

## Imports
using LinearAlgebra, Random
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
@polyvar R[1:3,1:3]
@polyvar t[1:3]
vars = [vec(R); t] # need to vectorize any matrix vars

# constants
const e3 = [0;0;1]

# 2) Build objective (to be minimized)
# express this in polynomial form
obj = 0.
for i = 1:N
    global obj
    residual_i = (I - y[:,i]*e3')*camK*(R*b[:,i] + t)
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
opt, sol, data = cs_tssos_first(pop, vars, order, numeq=length(eq), TS="block", CS="MF", QUIET=false, solve=true, solution=true)
# for more discussion on term sparsity (TS) and correlative sparsity (CS), see the TSSOS paper.
# `solve=true` tells TSSOS to solve the SDP
# `solution=true` tells TSSOS to round for an explicit solution

# Higher order relaxation
# opt, sol, data = cs_tssos_higher!(data)


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