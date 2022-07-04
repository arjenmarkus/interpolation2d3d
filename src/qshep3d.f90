! qshep3d.f90 --
!     Wrapper around the QSHEP3 collection of routines:
!     interpolation of scattered 3D data, using a modified
!     quadratic Shepard method. The result is a once continuously
!     differentiable function.
!
!     To do:
!     - qs3val should indicate an error, not simply return
!       In the test program it returns 0 or 4.5 - at random it seems!!
!
module quadratic_shepard_3d
    use interpolate_frame

    use alg661
    implicit none

    type, extends(interpolation_super) :: qshep_3d
        private
        integer              :: nr
        real                 :: xyzmin(3), xyzdel(3), rmax
        real, allocatable    :: x(:)
        real, allocatable    :: y(:)
        real, allocatable    :: z(:)
        real, allocatable    :: f(:)
        real, allocatable    :: rsq(:)
        real, allocatable    :: a(:,:)
        integer, allocatable :: lcell(:,:,:)
        integer, allocatable :: lnext(:)
    contains
        procedure :: create_mesh       => create_mesh_3d
        procedure :: interpolate_point => interpolate_point_3d
        procedure :: gradient_point    => gradient_point_3d
    end type qshep_3d

contains

! cleanup_mesh_3d --
!     Clean up the mesh
!
! Arguments:
!     this              Interpolation structure
!
subroutine cleanup_mesh_3d( this )
    class(qshep_3d), intent(inout) :: this

    if ( allocated(this%x)     ) deallocate( this%x )
    if ( allocated(this%y)     ) deallocate( this%y )
    if ( allocated(this%z)     ) deallocate( this%z )
    if ( allocated(this%f)     ) deallocate( this%f )
    if ( allocated(this%lcell) ) deallocate( this%lcell )
    if ( allocated(this%lnext) ) deallocate( this%lnext )
    if ( allocated(this%a)     ) deallocate( this%a )
    if ( allocated(this%rsq)   ) deallocate( this%rsq )

    this%error_code = code_no_errors

end subroutine cleanup_mesh_3d

! create_mesh_3d --
!     Create an interpolation mesh based on the scattered data
!
! Arguments:
!     this              Interpolation structure
!     x                 X-coordinate of the data points
!     y                 Y-coordinate of the data points
!     f                 Values at the data points
!
! Note:
!     Use the recommended parameter values. No choices here.
!
subroutine create_mesh_3d( this, x, y, z, f )
    class(qshep_3d), intent(inout) :: this
    real, intent(in)               :: x(:)
    real, intent(in)               :: y(:)
    real, intent(in)               :: z(:)
    real, intent(in)               :: f(:)

    integer                        :: n, nq, nr, nw
    integer                        :: ierr

    !
    ! Simple checks
    !
    if ( size(x) < 10 ) then
        this%error_code = code_too_few_points
        return
    endif

    if ( size(x) /= size(y) .or. size(x) /= size(f) ) then
        this%error_code = code_different_sizes
        return
    endif

    !
    ! Copy the information
    !
    call cleanup_mesh_3d( this )

    this%error_code = code_no_errors

    this%x  = x
    this%y  = y
    this%z  = z
    this%f  = f

    n       = size(x)
    nr      = min( 1.0, ( n / 3.0 ) ** (1.0/3.0) )
    this%nr = nr

    allocate( this%lcell(nr,nr,nr) )
    allocate( this%lnext(n) )
    allocate( this%rsq(n) )
    allocate( this%a(9,n) )

    !
    ! Now calculate the interpolation parameters
    !
    nq = min( 17, n - 1 )
    nw = min( 32, n - 1 )

    call qshep3( n, this%x, this%y, this%z, this%f, nq, nw, nr, this%lcell, this%lnext, &
                 this%xyzmin, this%xyzdel, this%rmax, this%rsq, this%a, ierr )

    select case ( ierr )
        case ( 0 )
            this%error_code = code_no_errors

        case ( 2 )
            this%error_code = code_duplicate_nodes

        case ( 3 )
            this%error_code = code_all_collinear

        case default
            this%error_code = code_programming_error

    end select

    write(*,*) 'n:', n
    write(*,*) 'nr:', nr
    write(*,*) 'rmax:', this%rmax
    write(*,*) 'xyzdel:', this%xyzdel
    write(*,*) 'ierr:', ierr
end subroutine create_mesh_3d

! interpolate_point_3d --
!     Return the interpolation value at point (px,py)
!
! Arguments:
!     this              Interpolation structure
!     x                 X-coordinate of the point
!     y                 Y-coordinate of the point
!
real function interpolate_point_3d( this, px, py, pz )
    class(qshep_3d), intent(inout) :: this
    real, intent(in)               :: px
    real, intent(in)               :: py
    real, intent(in)               :: pz

    integer                        :: ierr

   ! interface
   !     real function qs2val( px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,a)
   !         real, intent(in)    :: px, py
   !         integer, intent(in) :: n
   !         real, intent(in)    :: x(*), y(*), f(*)
   !         integer, intent(in) :: nr
   !         integer, intent(in) :: lcell(nr,*), lnext(*)
   !         real, intent(in)    :: xmin, ymin, dx, dy, rmax
   !         real, intent(in)    :: rsq(*)
   !         real, intent(in)    :: a(5,*)
   !     end function qs2val
   ! end interface

    interpolate_point_3d = qs3val( px, py, pz, size(this%x), this%x, this%y, this%z, this%f, &
                                   this%nr, this%lcell, this%lnext, this%xyzmin,             &
                                   this%xyzdel, this%rmax, this%rsq, this%a, ierr            )

    select case ( ierr )
        case ( 1 )
            this%error_code = code_programming_error ! This should not occur - construction

        case ( 2 )
            this%error_code = code_outside_spanned

        case default
            this%error_code = code_no_errors
    end select

end function interpolate_point_3d

! gradient_point_3d --
!     Return the estimated gradient at point (px,py)
!
! Arguments:
!     this              Interpolation structure
!     x                 X-coordinate of the point
!     y                 Y-coordinate of the point
!
function gradient_point_3d( this, px, py, pz )
    class(qshep_3d), intent(inout) :: this
    real, intent(in)               :: px
    real, intent(in)               :: py
    real, intent(in)               :: pz
    real                           :: gradient_point_3d(3)

   ! interface
   !     subroutine qs2grd( px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,rsq,a,q,qx,qy,ierr)
   !         real, intent(in)     :: px, py
   !         integer, intent(in)  :: n
   !         real, intent(in)     :: x(*), y(*), f(*)
   !         integer, intent(in)  :: nr
   !         integer, intent(in)  :: lcell(nr,*), lnext(*)
   !         real, intent(in)     :: xmin, ymin, dx, dy, rmax
   !         real, intent(in)     :: rsq(*)
   !         real, intent(in)     :: a(5,*)
   !         real, intent(out)    :: q, qx, qy
   !         integer, intent(out) :: ierr
   !     end subroutine qs2grd
   ! end interface

    real    :: q, qx, qy, qz
    integer :: ierr

    call qs3grd( px, py, pz, size(this%x), this%x, this%y, this%z, this%f, &
                 this%nr, this%lcell, this%lnext, this%xyzmin,             &
                 this%xyzdel, this%rmax, this%rsq, this%a,                 &
                 q, qx, qy, qz, ierr                                       )

    gradient_point_3d = [qx, qy, qz]

    select case ( ierr )
        case ( 0 )
            this%error_code = code_no_errors

        case ( 2 )
            this%error_code = code_extrapolation

        case default
            this%error_code = code_programming_error

    end select
end function gradient_point_3d

end module quadratic_shepard_3d
