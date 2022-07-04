# SPHPAK - interpolation on a sphere

The package **sphpak** handles interpolation on a sphere. The "latitude" and "longitude" coordinates have to be
provided as radian, not degrees. It offers the following methods:

Create the network (lat/lon coordinates - in radians - only, note the order):

```fortran
    use sphpak

    type(network_sphere) :: network

    real, dimension(n) :: lat, lon

    call network%create_mesh( lat, lon )
```

Set the values for the data points:

```fortran
    type(network_sphere) :: network

    real, dimension(n) :: w

    call network%set_values( w )
```

Note: since the mesh is independent of the values for each data point, you can change the values without having
to recalculate the network.

Interpolate for a single data point:

```fortran
    type(network_sphere) :: network

    real :: plat, plon, value

    value = network%interpolate_point( plat, plon )
```

Interpolate for a single data point using estimated gradients (once-differential approximation):

```fortran
    type(network_sphere) :: network

    real :: plat, plon, value

    real :: wxwywz(3,n) ! From gradient_global

    value = network%interpolate_diff( plat, plon )

    !
    ! Or, using globally estimated gradients
    !
    value = network%interpolate_diff( plat, plon, wxwywz )
```

Interpolate on a grid:

```fortran
    type(network_sphere) :: network

    real :: plat(nlat), plon(nlon)
    real :: value(nlat,nlon)

    real :: wxwywz(3,n) ! From gradient_global

    call network%interpolate_grid( plat, plon, value )

    !
    ! Or, using globally estimated gradients
    !
    call network%interpolate_grid( plat, plon, value, zxzy )
```

Estimate the gradient at a network point:

```fortran
    type(network_sphere) :: network

    integer :: indx
    real    :: dx, dy, dz

    call network%gradient_local( indx, dx, dy, dz )
```

Estimate the gradient for all network points, using a global method:

```fortran
    type(network_sphere) :: network

    real, allocatable :: wxwywz(:,:) ! Allocated by the method

    call network%gradient_global( wxwywz )
```

## Error handling

The method ``error_status`` can be used to retrieve the error code and the error message

```fortran
    type(network_sphere) :: network

    integer :: error_code
    character(len=80) :: error_message

    call network%error_status( error_code )

    ! Or, to get a descriptive string of the error status:

    call network%error_status( error_code, error_message )
```


