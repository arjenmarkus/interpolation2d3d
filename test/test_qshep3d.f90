! test_qshep3d.f90
!     Test program for the quadratic_shepard_3d module
!
program test_qshep3d
    use quadratic_shepard_3d

    implicit none

    type(qshep_3d)    :: shepard
    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: z(:)
    real, allocatable :: f(:)
    real              :: px, py, pz
    integer           :: ierr
    character(len=80) :: message

    !
    ! Construct a simple set of data
    !
    allocate( x(20), y(20), z(20), f(20) )

    call random_number( x )
    call random_number( y )
    call random_number( z )

    f = x + y + z

    write(*,*) 'create'
    call shepard%create_mesh( x, y, z, f )

    px = 0.5
    py = 0.5
    pz = 0.5

    write(*,*) 'interpolate'
    write(*,*) shepard%interpolate_point( px, py, pz ), ' -- expected: ', px+py+pz
    write(*,*) shepard%gradient_point( px, py, pz ),    ' -- expected: ', 1.0, 1.0, 1.0

    px = 1.5
    py = 1.5
    pz = 1.5

    write(*,*) px, py, pz, shepard%interpolate_point( px, py, pz )

    call shepard%error_status( ierr, message )
    write(*,*) ierr, ' - ', trim(message)

end program test_qshep3d
