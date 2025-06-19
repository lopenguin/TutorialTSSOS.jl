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
```
using Pkg
Pkg.add(url=["https://github.com/lopenguin/SimpleRotations.jl", "https://github.com/lopenguin/TSSOS"])
# you can use the official TSSOS instead for the features in `simple`
```

In the repo directory, run `julia` and type `]` to enter package mode. Then type `instantiate` to download all dependencies. The key dependency is [TSSOS](https://github.com/wangjie212/TSSOS).

TODO: 
- may be an issue with my version of TSSOS.
- Mention simplerotations

## Certifiable Perspective-N-Point
Plan:
- "simple" script demonstrates basic usage
- "extended" demonstrates adding constraints / imposing custom 
- Describe optimization problem, eliminating t, etc. below