! srfpak.f90 --
!     Wrapper around the SRFPAK/TRIPAK collection of routines for
!     interpolation using a triangular net.
!
!     To do:
!     - wrap GETNP, not the highest priority?
!     - wrap TRFIND, seems useful
!
!     - check for valid network when using it
!
module srfpak
    use interpolate_frame
    use alg624

    implicit none

    type, extends(interpolation_super) :: network_2d
        real, allocatable    :: x(:), y(:), z(:)
        integer, allocatable :: iadj(:), iend(:)
        integer              :: start
        real                 :: enclosed_area
    contains
        procedure :: create_mesh       => create_mesh_2d
        procedure :: set_values        => set_values_2d
        procedure :: interpolate_point => interpolate_point_2d
        procedure :: interpolate_cubic => interpolate_cubic_2d
        procedure :: interpolate_grid  => interpolate_grid_2d
        procedure :: volume            => volume_2d
        procedure :: area              => area_2d
        procedure :: gradient_global   => gradient_global_2d
        procedure :: gradient_local    => gradient_local_2d
    end type network_2d

    private
    public :: network_2d, coord_range, enclosed_area, barycentric_coords, circum_circle

    public :: volume_func_2d ! Temporary

contains

! create_mesh_2d --
!     Create a triangular network from the given points
!
! Arguments:
!     this              Network structure
!     x                 X-coordinates of the points
!     y                 Y-coordinates of the points
!
subroutine create_mesh_2d( this, x, y )
    class(network_2d), intent(inout)  :: this
    real, intent(in)                  :: x(:)
    real, intent(in)                  :: y(:)

    integer                           :: ierr, number

    if ( size(x) /= size(y) ) then
        this%error_code = code_different_sizes
        return
    endif

    if ( size(x) < 3 ) then
        this%error_code = code_too_few_points
        return
    endif

    if ( allocated(this%x) ) then
        deallocate( this%x, this%y, this%z, this%iadj, this%iend )
    endif

    number = size(x)
    allocate( this%x(number), this%y(number), this%z(number) )
    allocate( this%iadj(6*number) )
    allocate( this%iend(number) )

    this%error_code = code_no_errors
    this%x          = x
    this%y          = y
    this%start      = 1

    !
    ! Reordering the points would make the triangulation process
    ! more efficient, but this gives a problem when introducing
    ! the values. So, for the moment, leave it out.
    !
    ! call reordr( number, iflag, this%x this%y, dummy(1), ind )
    !

    !
    ! Determine the triangular network
    !
    call trmesh( number, this%x, this%y, this%iadj, this%iend, ierr )

    if ( ierr /= 0 ) then
        this%error_code = code_all_collinear
        deallocate( this%x, this%y, this%z, this%iadj, this%iend )
        return
    endif

    !
    ! Determine the enclosed area (easier than first querying the
    ! convex hull)
    !
    this%z             = 1.0
    this%enclosed_area = this%volume()

    !
    ! Reset the array of values
    this%z             = 0.0

end subroutine create_mesh_2d

! set_values --
!     Store the values for all points
!
! Arguments:
!     this              Network structure
!     z                 Array of values for the nodal points in the network
!
subroutine set_values_2d( this, z )
    class(network_2d), intent(inout) :: this
    real, intent(in)                 :: z(:)

    if ( size(z) /= size(this%z) ) then
        this%error_code = code_different_sizes
        return
    endif

    this%error_code = code_no_errors
    this%z          = z
end subroutine set_values_2d

! interpolate_point_2d --
!     Interpolate the values for a single point
!
! Arguments:
!     this              Network structure
!     px                X-coordinate of the point
!     py                Y-coordinate of the point
!
real function interpolate_point_2d( this, px, py )
    class(network_2d), intent(inout) :: this
    real, intent(in)                 :: px, py

    real                             :: pz
    integer                          :: ierr

    this%error_code = code_no_errors

    call intrc0( size(this%x), px, py, this%x, this%y, this%z, this%iadj, this%iend, this%start, pz, ierr )

    interpolate_point_2d = pz

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (1)
            this%error_code = code_extrapolation

        case (-1)
            this%error_code = code_programming_error

        case (-2)
            this%error_code = code_all_collinear

        case default
            this%error_code = code_programming_error

    end select

end function interpolate_point_2d

! interpolate_cubic_2d --
!     Interpolate the values for a single point using a cubic interpolant
!
! Arguments:
!     this              Network structure
!     px                X-coordinate of the point
!     py                Y-coordinate of the point
!     zxzy              Array of gradient values (optional)
!
real function interpolate_cubic_2d( this, px, py, zxzy )
    class(network_2d), intent(inout) :: this
    real, intent(in)                 :: px, py
    real, intent(in), optional       :: zxzy(:,:)

    real                             :: zxzy_dummy(1)
    real                             :: pz
    integer                          :: flag
    integer                          :: ierr

    if ( present(zxzy) ) then
        flag = 1
        call intrc1( size(this%x), px, py, this%x, this%y, this%z, this%iadj, this%iend, flag, zxzy, this%start, pz, ierr )
    else
        flag = 0
        call intrc1( size(this%x), px, py, this%x, this%y, this%z, this%iadj, this%iend, flag, zxzy_dummy, this%start, pz, ierr )
    endif

    interpolate_cubic_2d = pz

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (1)
            this%error_code = code_extrapolation

        case (-1)
            this%error_code = code_programming_error

        case (-2)
            this%error_code = code_all_collinear

        case default
            this%error_code = code_programming_error

    end select

end function interpolate_cubic_2d

! interpolate_grid_2d --
!     Interpolate the values on a grid
!
! Arguments:
!     this              Network structure
!     px                X-coordinates of the grid lines
!     py                Y-coordinates of the grid lines
!     z                 Result of the interpolation
!     zxzy              If present, to use these estimated gradients (see interpolate_cubic_2d
!                       and the gradient routines)
!
subroutine interpolate_grid_2d( this, px, py, z, zxzy )
    class(network_2d), intent(inout) :: this
    real, intent(in)                 :: px(:), py(:)
    real, intent(out)                :: z(:,:)
    real, intent(in), optional       :: zxzy(:,:)

    real                             :: zxzy_dummy(1)
    integer                          :: flag
    integer                          :: ierr

    if ( size(px) /= size(z,1) .or. size(py) /= size(z,2) ) then
        this%error_code = code_different_sizes
        return
    endif

    if ( present(zxzy) ) then
        flag = 1
        call unif( size(this%x), this%x, this%y, this%z, this%iadj, this%iend, size(z,1), size(px), size(py), &
                   px, py, flag, zxzy, z, ierr )
    else
        flag = 0
        call unif( size(this%x), this%x, this%y, this%z, this%iadj, this%iend, size(z,1), size(px), size(py), &
                   px, py, flag, zxzy_dummy, z, ierr )
    endif

    if ( ierr > 0 ) then
        this%error_code = code_extrapolation
    endif

end subroutine interpolate_grid_2d

! volume_2d --
!     Calculate the volume spanned by the nodal points and the values
!
! Arguments:
!     this              Network structure
!
real function volume_2d( this )
    class(network_2d), intent(inout) :: this

    volume_2d = volume( size(this%x), this%x, this%y, this%z, this%iadj, this%iend )

end function volume_2d

! area_2d --
!     Return the area spanned by the nodal points
!
! Arguments:
!     this              Network structure
!
real function area_2d( this )
    class(network_2d), intent(inout) :: this

    area_2d = this%enclosed_area
end function area_2d

! volume_func_2d --
!     Calculate the volume spanned by the nodal points and values of the given function
!
! Arguments:
!     this              Network structure
!     func              Function of x and y
!
real function volume_func_2d( this, func )
    class(network_2d), intent(inout) :: this

    interface
        real function func( x, y )
            real, intent(in) :: x, y
        end function func
    end interface

    real, allocatable :: z(:)
    integer           :: i

    allocate( z(size(this%x)) )

    do i = 1,size(this%x)
        z(i) = func( this%x(i), this%y(i) )
    enddo

    volume_func_2d = volume( size(this%x), this%x, this%y, z, this%iadj, this%iend )

end function volume_func_2d

! convex_hull_2d --
!     Return the coordinates of the convex hull
!
! Arguments:
!     this              Network structure
!     xb                Array of x-coordinates of the vertices of the convex hull
!     yb                Array of y-coordinates of the vertices of the convex hull
!
! Note:
!     The arrays xb and yb are allocated by the routine
!
subroutine convex_hull_2d( this, xb, yb )
    class(network_2d), intent(inout) :: this
    real, intent(out), allocatable   :: xb(:)
    real, intent(out), allocatable   :: yb(:)

    integer                          :: i, idx, nb, na, nt
    integer, allocatable             :: nodes(:)

    allocate( nodes(size(this%x)) )

    call bnodes( size(this%x), this%iadj, this%iend, nb, na, nt, nodes )

    allocate( xb(nb), yb(nb) )

    do i = 1,nb
        idx   = nodes(i)
        xb(i) = this%x(idx)
        yb(i) = this%y(idx)
    enddo
end subroutine convex_hull_2d

! gradient_global_2d --
!     Estimathe the global gradients of the given dataset
!
! Arguments:
!     this              Network structure
!     zxzy              Array of gradients in the x- and y-ditections
!
! Note:
!     The array zxzy is allocated to the right size
!
subroutine gradient_global_2d( this, zxzy )
    class(network_2d), intent(inout) :: this
    real, intent(out), allocatable   :: zxzy(:,:)

    integer                          :: nit, ierr
    integer, allocatable             :: nodes

    real                             :: eps

    !
    ! Use the recommended values
    !
    eps = 0.01
    nit = 4

    allocate( zxzy(2,size(this%x)) )
    zxzy = 0.0

    call gradg( size(this%x), this%x, this%y, this%z, this%iadj, this%iend, eps, nit, zxzy, ierr )

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (1)
            this%error_code = code_no_convergence

        case (2)
            this%error_code = code_programming_error

        case default
            this%error_code = code_programming_error

    end select

end subroutine gradient_global_2d

! gradient_local_2d --
!     Estimathe the gradients of the given dataset using a local method
!
! Arguments:
!     this              Network structure
!     idx               Index of the network node for which to determine the gradient
!     dx                Gradient in x-direction
!     dy                Gradient in y-direction
!
subroutine gradient_local_2d( this, idx, dx, dy )
    class(network_2d), intent(inout) :: this
    integer, intent(in)              :: idx
    real, intent(out)                :: dx
    real, intent(out)                :: dy

    integer                          :: ierr

    call gradl( size(this%x), idx, this%x, this%y, this%z, this%iadj, this%iend, dx, dy, ierr )

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (-1)
            this%error_code = code_programming_error

        case (-2)
            this%error_code = code_all_collinear

        case default
            this%error_code = code_programming_error

    end select

end subroutine gradient_local_2d

! ----------------------------------------------------
! Routines that do not depend on the network structure
!


! enclosed_area --
!     Wrapper for the function AREA
!
! Arguments:
!     x                 Array of x-coordinates of the vertices of a polygon
!     y                 Array of y-coordinates of the vertices of a polygon
!
! Note:
!     The area is negative if the vertices are traversed in clockwise
!     direction.
!
real function enclosed_area( x, y )
    real, intent(in)     :: x(:), y(:)

    integer              :: i, nb
    integer, allocatable :: nodes(:)

    nb = size(x)

    if ( size(x) /= size(y) ) then
        write(*,*) 'Enclosed area: coordinate arrays not the same size'
        error stop
    endif

    allocate( nodes(nb) )
    nodes = [ (i, i =1,nb)]

    enclosed_area = area( x, y, nb, nodes )
end function enclosed_area

! barycentric_coords --
!     Wrapper for the subroutine COORDS: determine the barycentric
!     coordinates for a point inside a triangle
!
! Arguments:
!     x                 X-coordinate of the point
!     y                 Y-coordinate of the point
!     xtriangle         Array of x-coordinates of the vertices of the triangle
!     ytriangle         Array of y-coordinates of the vertices of the triangle
!
! Returns:
!     Barycentric coordinates as an array of three elements
!
function barycentric_coords( x, y, xtriangle, ytriangle )
    real, intent(in)     :: x, y
    real, intent(in)     :: xtriangle(:), ytriangle(:)

    real                 :: barycentric_coords(3)

    integer              :: ierr

    call coords( x, y, xtriangle(1), xtriangle(2), xtriangle(3), &
                       ytriangle(1), ytriangle(2), ytriangle(3), &
                       barycentric_coords, ierr )

    if ( ierr /= 0 ) then
        write(*,*) 'Vertices of the triangle are collinear'
        error stop
    endif

end function barycentric_coords

! circum_circle --
!     Calculate the circum circle of a triangle
!
! Arguments:
!     xtriangle         Array of x-coordinates of the vertices of the triangle
!     ytriangle         Array of y-coordinates of the vertices of the triangle
!
! Returns:
!     The coordinates of the centre of the circle
!
function circum_circle( xtriangle, ytriangle )
    real, intent(in)     :: xtriangle(:), ytriangle(:)

    real                 :: circum_circle(2)

    integer              :: ierr

    call circum( xtriangle(1), xtriangle(2), xtriangle(3), &
                 ytriangle(1), ytriangle(2), ytriangle(3), &
                 circum_circle(1), circum_circle(2), ierr )

    if ( ierr /= 0 ) then
        write(*,*) 'Vertices of the triangle are collinear'
        error stop
    endif

end function circum_circle

! coord_range --
!     Simple procedure to set a coordinate range
!
! Arguments:
!     coord             Array of coordinates to be filled
!     crdmin            Minimum coordinate
!     crdmax            Maximum coordinate
!     rate              Rate at which the intervals should grow
!                       (optional, defaults to 1)
!
subroutine coord_range( coord, crdmin, crdmax, rate )
    real, intent(out)          :: coord(:)
    real, intent(in)           :: crdmin
    real, intent(in)           :: crdmax
    real, intent(in), optional :: rate

    real                       :: c, dc, crate
    integer                    :: i, n

    if ( present(rate) ) then
        crate = rate
    else
        crate = 1.0
    endif

    n = size(coord)

    if ( crate == 1.0 ) then
        dc = (crdmax - crdmin) / (n - 1)
        do i = 1,n
            coord(i) = crdmin + (i-1) * dc
        enddo
    else
        dc = (crdmax - crdmin) * (1.0 - crate) / (1.0 - crate**(n+1))
        do i = 1,n
            coord(i) = crdmin + (1.0 - crate**(i+1)) / (1.0 - crate) * dc
        enddo
    endif
end subroutine coord_range

!!include 'alg624.f90'

end module srfpak
