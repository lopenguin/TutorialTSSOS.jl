# TSSOS Tutorial
A TSSOS tutorial aimed at practioners in robotics / computer vision.

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

## Certifiable Perspective-N-Point
Plan:
- "simple" script demonstrates basic usage
- "extended" demonstrates adding constraints / imposing custom 
- Describe optimization problem, eliminating t, etc. below