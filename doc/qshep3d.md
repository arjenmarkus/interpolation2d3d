# QSHEP3D - interpolation in 3D using a once continuously differentiable approximation

The package **qshep3d* handles interpolation via scattered 3D points.
It offers the following methods:

Create the mesh, including the values:

```fortran
    use quadratic_shepard_3d

    type(qshep3d) :: mesh

    real, dimension(n) :: x, y, z, value

    call mesh%create_mesh( x, y, z, value )
```

Note that due to the approximation method the values are required when setting up the interpolation structure.
New data require a new structure, but this is handled within the ``create_mesh`` routine.

Interpolate for a single data point:

```fortran
    type(qshep3d) :: mesh

    real :: px, py, pz

    value = mesh%interpolate_point( px, py, pz )
```

Estimate the gradient at a single data point:

```fortran
    type(qshep3d) :: mesh

    real :: px, py, py
    real :: dxdydz(3)

    dxdydz = mesh%gradient_point( px, py, pz )
```

## Error handling

The method ``error_status`` can be used to retrieve the error code and the error message

```fortran
    type(qshep3d) :: mesh

    integer :: error_code
    character(len=80) :: error_message

    call mesh%error_status( error_code )

    ! Or, to get a descriptive string of the error status:

    call mesh%error_status( error_code, error_message )
```


