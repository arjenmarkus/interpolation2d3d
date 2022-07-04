# SRFPAK - interpolation in two dimensions

The package **srfpak** has the following methods:

Create the network (x/y coordinates only):

```fortran
    use srfpak

    type(network_2d) :: network

    real, dimension(n) :: x, y

    call network%create_mesh( x, y )
```

Set the values for the data points:

```fortran
    type(network_2d) :: network

    real, dimension(n) :: z

    call network%set_values( z )
```

Note: since the mesh is independent of the values for each data point, you can change the values without having
to recalculate the network.

Interpolate for a single data point:

```fortran
    type(network_2d) :: network

    real :: px, py, value

    value = network%interpolate_point( px, py )
```

Interpolate for a single data point using estimated gradients (cubic interpolation):

```fortran
    type(network_2d) :: network

    real :: px, py, value

    real :: zxzy(2,n) ! From gradient_global

    value = network%interpolate_cubic( px, py )

    !
    ! Or, using globally estimated gradients
    !
    value = network%interpolate_cubic( px, py, zxzy )
```

Interpolate on a grid:

```fortran
    type(network_2d) :: network

    real :: px(nx), py(ny)
    real :: value(nx,ny)

    real :: zxzy(2,n) ! From gradient_global

    call network%interpolate_grid( px, py, value )

    !
    ! Or, using globally estimated gradients
    !
    call network%interpolate_grid( px, py, value, zxzy )
```

Estimate the gradient at a network point:

```fortran
    type(network_2d) :: network

    integer :: indx
    real    :: dx, dy

    call network%gradient_local( indx, dx, dy )
```

Estimate the gradient for all network points, using a global method:

```fortran
    type(network_2d) :: network

    real, allocatable :: zxzy(:,:) ! Allocated by the method

    call network%gradient_global( zxzy )
```

Volume and area of the network with data (considered as a digital terrain model):

```fortran
    type(network_2d) :: network

    real :: volume, area
    real :: average

    volume = network%volume()
    area   = network%area()

    average = volume / area
```

Convex hull (set of points forming the edges of the network):

```fortran
    type(network_2d) :: network

    real, allocatable :: xb(:), yb(:)

    call network%convex_hull( xb, yb )
```

## Auxiliary routines (independent of the network):

Barycentric coordinates of a point w.r.t. a triangle:

```fortran
    real :: xtriangle(3), ytriangle(3)
    real :: coords(3)

    real :: x, y

    coords = barycentric_coords( x, y, xtriangle, ytriangle )
```

Centre of the circum circle of a triangle:

```fortran
    real :: xtriangle(3), ytriangle(3)
    real :: coords(2)

    real :: x, y

    coords = circum_circle( xtriangle, ytriangle )
```

Range of coordinates (useful for generating a rectilinear grid):

```fortran
    real :: coords(m)

    real :: crdmin, crdmax
    real :: rate      ! Determines the spacing (growth) via a geometric series

    call coord_range( coords, crdmin,crdmax, rate )
```

## Error handling

The method ``error_status`` can be used to retrieve the error code and the error message

```fortran
    type(network_2d) :: network

    integer :: error_code
    character(len=80) :: error_message

    call network%error_status( error_code )

    ! Or, to get a descriptive string of the error status:

    call network%error_status( error_code, error_message )
```


