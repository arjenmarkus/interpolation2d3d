! interpolate_frame.f90 --
!     Overall class to store the information common to all interpolation classes
!     For the moment: only the error messages
!
module interpolate_frame
    implicit none

    integer, parameter :: code_no_errors          = 0
    integer, parameter :: code_too_few_points     = 1
    integer, parameter :: code_all_collinear      = 2
    integer, parameter :: code_no_convergence     = 3
    integer, parameter :: code_invalid_tolerance  = 4
    integer, parameter :: code_index_out_of_range = 5
    integer, parameter :: code_different_sizes    = 6
    integer, parameter :: code_extrapolation      = 7
    integer, parameter :: code_programming_error  = 8
    integer, parameter :: code_outside_spanned    = 9
    integer, parameter :: code_duplicate_nodes    = 10

    character(len=60), parameter :: status_message(0:10) = [             &
         'No errors encountered                                       ', &    ! 0
         'Too few points (at least three required)                    ', &    ! 1
         'All points collinear                                        ', &    ! 2
         'Convergence not reached within given limit                  ', &    ! 3
         'Invalid value of tolerance (must be positive)               ', &    ! 4
         'Index out of range                                          ', &    ! 5
         'Arrays have different sizes                                 ', &    ! 6
         'Extrapolation performed                                     ', &    ! 7
         'Programming error - impossible case encountered             ', &    ! 8
         'Point outside spanned data area                             ', &    ! 9
         'Duplicate nodes in the data                                 '  ]    ! 10

    type :: interpolation_super
        integer :: error_code = 0
    contains
        procedure :: error_status => error_status_super
    end type interpolation_super

contains

! error_status_super --
!     Return error status
!
! Arguments:
!     this              Network structure
!     code              Error code
!     message           Status message (optional)
!
subroutine error_status_super( this, code, message )
    class(interpolation_super), intent(in)  :: this
    integer, intent(out)                    :: code
    character(len=*), intent(out), optional :: message

    code = this%error_code
    if ( present(message) ) then
        message = status_message(this%error_code)
    endif
end subroutine error_status_super

end module interpolate_frame
