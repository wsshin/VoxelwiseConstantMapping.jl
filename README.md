# VoxelwiseConstantMapping.jl

[![CI](https://github.com/wsshin/VoxelwiseConstantMapping.jl/workflows/CI/badge.svg)](https://github.com/wsshin/VoxelwiseConstantMapping.jl/actions)
[![Codecov](http://codecov.io/github/wsshin/VoxelwiseConstantMapping.jl/coverage.svg?branch=main)](http://codecov.io/github/wsshin/VoxelwiseConstantMapping.jl?branch=main)

## Introduction
**VoxelwiseConstantMapping** is a simple Julia package that generalizes piecewise-constant functions to higher dimensions.  In the numerical implementation of a piecewise-constant function, the main subtlety occurs at the transition point˜˜.  For example, the [Heaviside step function](https://en.wikipedia.org/wiki/Heaviside_step_function) is defined as `H(x) = 0` for `x < 0` and `H(x) = 1` for `x > 0`, but its value at the transition point `x = 0` is ambiguous.  One popular choice is to assign the average of the values on the two sides of the transition point (i.e., `H(0) = 0.5` for the Heaviside step function), because such a choice reduces the discontinuity at the transition point.

**VoxelwiseConstantMapping** extends this concept to an arbitrary dimension.  For example, consider a function `f(x,y)` defined on R².  Suppose that we have pixel boundaries at `x = 0` and `y = 0`.  In other words,
- `f(x>0, y>0) = v₁`
- `f(x<0, y>0) = v₂`
- `f(x<0, y<0) = v₃`
- `f(x>0, y<0) = v₄`

Then, **VoxelwiseConstantMapping** defines the values of `f` at the transition lines as follows:
- `f(x=0, y>0) = (v₁+v₂)/2`
- `f(x=0, y<0) = (v₃+v₄)/2`
- `f(x>0, y=0) = (v₁+v₄)/2`
- `f(x<0, y=0) = (v₂+v₃)/2`
- `f(x=0, y=0) = (v₁+v₂+v₃+v₄)/4`

## Example
The main constructor to use is `VoxelwiseConstant(...)`.  The voxels are assumed having rectangular faces, and you can have as many voxels as you want.  The voxel sizes can be nonuniform.

```julia
using VoxelwiseConstantMapping

vals = rand(3,4)  # 3 × 4 pixels (= 2D voxels)
xbounds = [-1, 0, 0.5, 0.6]  # 4 boundaries of 3 pixels in x-direction
ybounds = [-0.1, -0.05, 0, 4, 10]  # 5 boundaries of 4 pixels in y-direction
vc = VoxelwiseConstant(vals, (xbounds, ybounds))

v₀₀  = vc([0,0])  # returns value at (x,y) = (0,0)
```
