# TSSOS Tutorial
A quick TSSOS tutorial aimed at practioners in robotics / computer vision. I assume the user is familiar with 
Shor's relaxation / the moment-SOS hierarchy, and focus on how to use TSSOS. Specifically, I focus on:
- The basics of using TSSOS to automatically relax a polynomial optimization problem to a convex semidefinite program
- Extracting the JuMP model from TSSOS to add constraints / use a different solver

Notably, the tutorial does not include:
- Complex, symmetric, matrix, etc. sum of squares
- Using binary variables in TSSOS
- Different types of sparsity

If there is interest, I may expand this tutorial to cover these problems in the future.

Written by [Lorenzo Shaikewitz](lorenzos@mit.edu).

## Quick Start
First, clone the repo:
```shell
git clone https://github.com/lopenguin/TSSOSTutorial.jl
```

In the repo directly, enter Julia with the repo environment:
```shell
julia --project
```

You'll need to manually add TSSOS and my helper package SimpleRotations:
```julia
(TutorialTSSOS) pkg> add https://github.com/lopenguin/SimpleRotations.jl, https://github.com/lopenguin/TSSOS
# you can use the official TSSOS instead for the features in `simple`
```

You are now ready to run either example script! Return to the Julia REPL and type:
```julia
include("scripts/example_pnp_simple.jl")
```

## Certifiable Perspective-n-Point
I briefly review the backprojection form of perspective-n-point (PnP) and state the optimization problems solved in these examples.

The maximum likelihood form of the backprojection PnP problem is given by the following optimization problem:
$$
\min_{\substack{R\in\mathrm{SO}(3)\\t\in\mathbb{R}^3}} \sum_{i=1}^N \frac1{\sigma_i^2} \|(I_3 - y_i \hat{e}_3^T)K(R b_i + t)\|_2^2
$$
$$
\text{s.t. }\ \ \hat{e}_3^T K(R b_i + t) > 0,\ i=1,...,N
$$

This is the optimization problem solved in [pnp_simple](scripts/example_pnp_simple.jl).

To define the varaibles: $y_i$ are the 2D pixel measurements, $b_i$ are the corresponding 3D points, $K$ is the camera calibration matrix, and $\sigma$ is a vector of isotropic variances for each of the $N$ measurements. We also enforce chirality constraints.

With a little algebra (including dropping the chirality constraints), we can solve for the optimal $t$ as a function of $R$. In [pnp_extended](scripts/example_pnp_extended.jl), we solve the previous optimization problem but fix $t$ as below:
$$
t^\star = -\left(\sum_{j=1}^N \frac1{\sigma_j^2} U_j^T U_j\right)^{-1} \sum_{i=1}^N \frac1{\sigma_i^2}U_i^T U_i R b_i
$$
where $U_i\triangleq I_3 - y_i\hat{e}_3^T$.