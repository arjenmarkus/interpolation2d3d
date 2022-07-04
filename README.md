# Interpolation in 2D and 3D

Interpolation in two and three dimensions based on the packages available at [netlib.org](http://netlib.org)

The packages srfpak, sphpak, qshep2d and qpshep3d implement an object-oriented interface to the interpolation
packages by Robert Renka. The code as found on netlib has been modernised using the [spag program](https://www.fortran.uk/plusfortmanual/spag.html).
Some further minor changes to the code have been introduced to avoid error messages about array-bound
violations, due to the old-style FORTRAN 77 declaration "iadj(1)" and the like.

The packages offer the following functionality:

Given a set of data points, estimate the value at an arbitrary point by means of suitable interpolation
methods.

The packages differ in their methods and geometrical features, but all provide an interface that allows:

* Creating a suitable network with the given data points
* Interpolating the data to an arbitrary point
* Estimating the gradient of the interpolation function

In particular:

* **srfpak** interpolates on a two-dimensional network of data points (TOMS algorithm 624)
* **sphpak** interpolates on a sphere covered (partially) with data points (TOMS algorithm 623)
* **qshep2d** interpolates in two dimensions using a modified quadratic method, with the resulting function once continously differentiable (TOMS algorithm 660)
* **qshep3d** interpolates in three dimensions using a modified quadratic method, with the resulting function once continously differentiable (TOMS algorithm 661)

See the documentation for the individual packages for more details.

Simple example:
```fortran
program example
    use srfpak

    implicit none

    type(network_2d) :: network

    real              :: x(4) = [ 0.0, 0.0, 1.0, 1.0]
    real              :: y(4) = [ 0.0, 1.0, 0.0, 1.0]
    real              :: z(4)
    real, allocatable :: xp(:), yp(:), zp(:), zxzy(:,:)
    real, allocatable :: xcrd(:), ycrd(:), zint(:,:)
    real              :: coord(5)

    real              :: result, xc, yc, dx, dy
    integer           :: i

    call network%create_mesh( x, y )

    z = 1.0

    call network%set_values( z )

    write(*,*) 'Point (0.5,0.5):   ', network%interpolate_point( 0.5, 0.5 )

    z = [ 0.0, 0.0, 1.0, 1.0 ]

    call network%set_values( z )

    write(*,*) 'Point (0.5,0.5):   ', network%interpolate_point( 0.5, 0.5 )

    write(*,*) 'Volume:            ', network%volume()

    write(*,*) 'Area:              ', network%area()
end program example
```
