! sphpak.f90 --
!     Wrapper around the SPHPAK collection of routines for
!     interpolation on a sphere using a triangular net.
!
!     To do:
!     - wrap GETNP, not the highest priority?
!     - wrap TRFIND, seems useful
!
!     - check for valid network when using it
!
module sphpak
    use interpolate_frame
    use alg623

    implicit none

    type, extends(interpolation_super) :: network_sphere
        real, allocatable    :: x(:), y(:), z(:), w(:), lat(:), lon(:)
        integer, allocatable :: iadj(:), iend(:)
        integer              :: start
        real                 :: enclosed_area
    contains
        procedure :: create_mesh       => create_mesh_sphere
        procedure :: set_values        => set_values_sphere
        procedure :: interpolate_point => interpolate_point_sphere
        procedure :: interpolate_diff  => interpolate_diff_sphere
        procedure :: interpolate_grid  => interpolate_grid_sphere
        procedure :: gradient_global   => gradient_global_sphere
        procedure :: gradient_local    => gradient_local_sphere
    end type network_sphere

    private
    public :: network_sphere

contains

! create_mesh_sphere --
!     Create a triangular network from the given points
!
! Arguments:
!     this              Network structure
!     lat               Latitudes of the points (radians!)
!     lon               Longitudes of the points (radians!)
!
subroutine create_mesh_sphere( this, lat, lon )
    class(network_sphere), intent(inout)  :: this
    real, intent(in)                      :: lat(:)
    real, intent(in)                      :: lon(:)

    integer                               :: ierr, number

    if ( size(lat) /= size(lon) ) then
        this%error_code = code_different_sizes
        return
    endif

    if ( size(lat) < 3 ) then
        this%error_code = code_too_few_points
        return
    endif

    if ( allocated(this%lat) ) then
        deallocate( this%lat, this%lon, this%z, this%iadj, this%iend )
    endif

    number = size(lat)
    allocate( this%lat(number), this%lon(number), this%x(number), this%y(number), this%z(number), this%w(number) )
    allocate( this%iadj(6*number) )
    allocate( this%iend(number) )

    this%error_code = code_no_errors
    this%lat        = lat
    this%lon        = lon
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
    call trans( number, lat, lon, this%x, this%y, this%z )

    call trmesh( number, this%x, this%y, this%z, this%iadj, this%iend, ierr )

    if ( ierr /= 0 ) then
        this%error_code = code_all_collinear
        deallocate( this%lat, this%lon, this%x, this%y, this%z, this%w, this%iadj, this%iend )
        return
    endif

end subroutine create_mesh_sphere

! set_values --
!     Store the values for all points
!
! Arguments:
!     this              Network structure
!     w                 Array of values for the nodal points in the network
!
subroutine set_values_sphere( this, w )
    class(network_sphere), intent(inout) :: this
    real, intent(in)                 :: w(:)

    if ( size(w) /= size(this%w) ) then
        this%error_code = code_different_sizes
        return
    endif

    this%error_code = code_no_errors
    this%w          = w
end subroutine set_values_sphere

! interpolate_point_sphere --
!     Interpolate the values for a single point
!
! Arguments:
!     this              Network structure
!     plat              Latitude of the point (radians!)
!     plon              Longitude of the point (radians!)
!
real function interpolate_point_sphere( this, plat, plon )
    class(network_sphere), intent(inout) :: this
    real, intent(in)                     :: plat, plon

    real                                 :: pw
    integer                              :: ierr

    this%error_code = code_no_errors

    call intrc0( size(this%x), plat, plon, this%x, this%y, this%z, this%w, this%iadj, this%iend, this%start, pw, ierr )

    interpolate_point_sphere = pw

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (1)
            this%error_code = code_extrapolation

        case (-1)
            this%error_code = code_programming_error

        case (-2)
            this%error_code = code_all_collinear

        case (-3)
            this%error_code = code_outside_spanned

        case default
            this%error_code = code_programming_error

    end select

end function interpolate_point_sphere

! interpolate_diff_sphere --
!     Interpolate the values for a single point using a once differentiable approximation
!
! Arguments:
!     this              Network structure
!     plat              Latitude of the point
!     plon              Longitude of the point
!     wxzy              Array of gradient values (optional) -- TODO
!
real function interpolate_diff_sphere( this, plat, plon, wxwywz )
    class(network_sphere), intent(inout) :: this
    real, intent(in)                     :: plat, plon
    real, intent(in), optional           :: wxwywz(:,:)

    real                                 :: wxwywz_dummy(3,1)
    real                                 :: pw
    integer                              :: flag
    integer                              :: ierr

    if ( present(wxwywz) ) then
        flag = 1
        call intrc1( size(this%x), plat, plon, this%x, this%y, this%z, this%w, this%iadj, this%iend, flag, &
                     wxwywz, this%start, pw, ierr )
    else
        flag = 0
        call intrc1( size(this%x), plat, plon, this%x, this%y, this%z, this%w, this%iadj, this%iend, flag, &
                     wxwywz_dummy, this%start, pw, ierr )
    endif

    interpolate_diff_sphere = pw

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (1)
            this%error_code = code_extrapolation

        case (-1)
            this%error_code = code_programming_error

        case (-2)
            this%error_code = code_all_collinear

        case (-3)
            this%error_code = code_outside_spanned

        case default
            this%error_code = code_programming_error

    end select

end function interpolate_diff_sphere

! interpolate_grid_sphere --
!     Interpolate the values on a grid
!
! Arguments:
!     this              Network structure
!     plat              Latitudes of the grid lines
!     plon              Longitudes of the grid lines
!     z                 Result of the interpolation
!     wxwy              If present, to use these estimated gradients (see interpolate_diff_sphere
!                       and the gradient routines)
!
subroutine interpolate_grid_sphere( this, plat, plon, w,  wxwywz )
    class(network_sphere), intent(inout) :: this
    real, intent(in)                     :: plat(:), plon(:)
    real, intent(out)                    :: w(:,:)
    real, intent(in), optional           :: wxwywz(:,:)

    real                                 :: wxwywz_dummy(3,1)
    integer                              :: flag
    integer                              :: ierr

    ! TODO!
    if ( size(plat) /= size(w,1) .or. size(plon) /= size(w,2) ) then
        this%error_code = code_different_sizes
        return
    endif

    if ( present(wxwywz) ) then
        flag = 1
        call unif( size(this%x), this%x, this%y, this%z, this%w, this%iadj, this%iend, size(w,1), size(plat), size(plon), &
                   plat, plon, flag, wxwywz, w, ierr )
    else
        flag = 0
        call unif( size(this%x), this%x, this%y, this%z, this%w, this%iadj, this%iend, size(w,1), size(plat), size(plon), &
                   plat, plon, flag, wxwywz_dummy, w, ierr )
    endif

    if ( ierr > 0 ) then
        this%error_code = code_extrapolation
    endif

end subroutine interpolate_grid_sphere

! gradient_global_sphere --
!     Estimate the global gradients of the given dataset
!
! Arguments:
!     this              Network structure
!     wxwy              Array of gradients in the x-, y- and z-directions (!)
!
! Note:
!     The array zxzy is allocated to the right size
!
subroutine gradient_global_sphere( this, wxwywz )
    class(network_sphere), intent(inout) :: this
    real, intent(out), allocatable       :: wxwywz(:,:)

    integer                              :: nit, ierr
    integer, allocatable                 :: nodes

    real                                 :: eps

    !
    ! Use the recommended values
    !
    eps = 0.001
    nit = 5

    allocate( wxwywz(3,size(this%x)) )
    wxwywz = 0.0

    call gradg( size(this%x), this%x, this%y, this%z, this%w, this%iadj, this%iend, eps, nit, wxwywz, ierr )

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

end subroutine gradient_global_sphere

! gradient_local_sphere --
!     Estimathe the gradients of the given dataset using a local method
!
! Arguments:
!     this              Network structure
!     idx               Index of the network node for which to determine the gradient
!     dx                Gradient in x-direction
!     dy                Gradient in y-direction
!     dz                Gradient in y-direction
!
subroutine gradient_local_sphere( this, idx, dx, dy, dz )
    class(network_sphere), intent(inout) :: this
    integer, intent(in)              :: idx
    real, intent(out)                :: dx
    real, intent(out)                :: dy
    real, intent(out)                :: dz

    real                             :: g(3)

    integer                          :: ierr

    call gradl( size(this%x), idx, this%x, this%y, this%z, this%w, this%iadj, this%iend, g, ierr )

    dx = g(1)
    dy = g(2)
    dz = g(3)

    select case ( ierr )
        case (0)
            this%error_code = code_no_errors

        case (-1)
            this%error_code = code_too_few_points  ! At least seven required

        case (-2)
            this%error_code = code_all_collinear

        case default
            this%error_code = code_no_errors ! The value represents the number of nodes used, not passed on

    end select

end subroutine gradient_local_sphere

end module sphpak
