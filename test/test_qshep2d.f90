! test_qshep2d.f90
!     Test program for the quadratic_shepard_2d module
!
program test_qshep2d
    use quadratic_shepard_2d

    implicit none

    type(qshep_2d)    :: shepard
    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: f(:)
    real              :: px, py
    integer           :: i, error_code
    character(len=80) :: error_message

    !
    ! Construct a simple set of data
    !
    allocate( x(10), y(10), f(10) )

    call random_number( x )
    call random_number( y )

    f = x + y

    write(*,*) 'create'
    call shepard%create_mesh( x, y ,f )

    px = 0.5
    py = 0.5

    write(*,*) 'interpolate'
    write(*,*) shepard%interpolate_point( px, py ), ' -- expected: ', px+py
    write(*,*) shepard%gradient_point( px, py ),    ' -- expected: ', 1.0, 1.0

    px = 1.5
    py = 1.5

    write(*,*) px, py, shepard%interpolate_point( px, py )
    call shepard%error_status( error_code, error_message )
    write(*,*) 'Error: ', error_code, ' -- ', trim(error_message)

end program test_qshep2d
