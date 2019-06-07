# Examples

This page is a work in progress

## Defining a hexagonal grid


```julia
grid = HASlib.HexGrid(4e-10)
# output
HASlib.HexGrid(4.0e-10, Generic2DLattice([4.0e-10, 0.0], [-2.0e-10, 3.4641e-10], [1.5708e10, 9.069e9], [-0.0, 1.8138e10]))
```

## Orient the the first coordinate along a specific direction
```julia
b1,b2 = HASlib.rotate_reciprocal(grid.l,[1,0])
# output
([1.8138e10, 0.0], [9.069e9, 1.5708e10])
```
```julia
b1,b2 = HASlib.rotate_reciprocal(grid.l,[1,1])
# output
([1.5708e10, -9.069e9], [1.5708e10, 9.069e9])
```
```julia
b1,b2 = HASlib.rotate_reciprocal(grid.l,[0,1])
# output
([9.069e9, -1.5708e10], [1.8138e10, 1.11063e-6])
```
