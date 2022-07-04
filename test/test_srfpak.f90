! test_srfpak.f90 --
!     Straightforward test program for SRFPAK
!
program test_srfpak
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

    write(*,*) 'Volume (function): ', volume_func_2d(network, f)


    !
    ! Interpolation for a pyramid:
    !
    !              +
    !            / | \  <-- this plane: z = 1-(x+y)
    !           /  |  \
    !          +---+---+    zero at the corners, one in the centre
    !           \  |  /
    !            \ | /
    !              +
    xp = [-1.0, 0.0, 0.0, 1.0, 0.0]
    yp = [ 0.0,-1.0, 1.0, 0.0, 0.0]
    zp = [ 0.0, 0.0, 0.0, 0.0, 1.0]

  ! When I used this, it immediately hung
  !  call network%create_mesh( x, y )
  !  call network%set_values( z )

  ! With this at least the first interpolation gives a result
  !
    call network%create_mesh( xp, yp )
    call network%set_values( zp )

    result = network%interpolate_point( 0.0,  0.0  )
    write(*,*) 'Point ( 0.0,  0.0 ):   ', result, ' - expected: 1.0'
    result = network%interpolate_point( 0.5,  0.5  )
    write(*,*) 'Point ( 0.5,  0.5 ):   ', result, ' - expected: 0.0'


    write(*,*) 'Point ( 0.0,  0.0 ):   ', network%interpolate_point( 0.0,  0.0  ), ' - expected: 1.0'
   ! write(*,*) 'Point ( 0.5,  0.5 ):   ', network%interpolate_point( 0.5,  0.5  ), ' - expected: 0.0' -- extrapolation
    write(*,*) 'Point ( 0.25, 0.25):   ', network%interpolate_point( 0.25, 0.25 ), ' - expected: 0.5'
    write(*,*) 'Point (-0.25,-0.25):   ', network%interpolate_point(-0.25,-0.25 ), ' - expected: 0.5'
    write(*,*) 'Point (-0.25, 0.25):   ', network%interpolate_point(-0.25, 0.25 ), ' - expected: 0.5'
    write(*,*) 'Point ( 0.25,-0.25):   ', network%interpolate_point( 0.25,-0.25 ), ' - expected: 0.5'
    result = network%interpolate_point( 1.00, 1.00 )
    write(*,*) 'Point ( 1.00, 1.00):   ', result, ' - extrapolation - -1?'
  ! write(*,*) 'Point ( 1.00, 1.00):   ', network%interpolate_point( 0.25,-0.25 ), ' - extrapolation - -1?'

    zp = [ 0.0, 0.0, 0.0, 1.0, 1.0]
    call network%set_values( zp )

    do i = 1,15
        xc = 0.1 * i
        yc = 0.1 * i
        result = network%interpolate_point( xc, yc )
        write(*,*) i, xc, yc, result
    enddo

    write(*,*) 'Check border'
    xc = 0.5 * (1.0 - epsilon(xc))
    yc = 0.5 * (1.0 - epsilon(yc))
    result = network%interpolate_point( xc, yc )
    write(*,*) i, xc, yc, result
    xc = 0.5 * (1.0 + epsilon(xc))
    yc = 0.5 * (1.0 + epsilon(yc))
    result = network%interpolate_point( xc, yc )
    write(*,*) i, xc, yc, result

    ! Gradient
    write(*,*) 'Local gradients:'
    do i = 1,5
        call network%gradient_local( i, dx, dy )
        write(*,*) i, xp(i), yp(i), dx, dy
    enddo

    write(*,*) 'Global gradients:'
    call network%gradient_global( zxzy )
    do i = 1,5
        write(*,*) i, xp(i), yp(i), zxzy(:,i)
    enddo

!   write(*,*) network%error_code

    ! Gradient
    zp = xp + yp
    call network%set_values( zp )

    write(*,*) 'Local gradients:'
    do i = 1,5
        call network%gradient_local( i, dx, dy )
        write(*,*) i, xp(i), yp(i), dx, dy
    enddo

    write(*,*) 'Global gradients:'
    call network%gradient_global( zxzy )
    do i = 1,5
        write(*,*) i, xp(i), yp(i), zxzy(:,i)
    enddo

    !
    ! Grids
    !
    call coord_range( coord, 1.0, 2.0 )
    write(*,*) coord
    call coord_range( coord, 1.0, 2.0, 1.25 )
    write(*,*) coord
    call coord_range( coord, 1.0, 2.0, 0.25 )
    write(*,*) coord

    allocate( xcrd(4), ycrd(4), zint(4,4) )
    call coord_range( xcrd, -1.0, 1.0 )
    call coord_range( ycrd, -1.0, 1.0 )

    call network%interpolate_grid( xcrd, ycrd, zint )
    write(*,'(10x,4f10.3)') xcrd
    write(*,'(5f10.3)') (ycrd(i), zint(:,i) , i = 1,size(ycrd))

contains
real function f(x,y)
    real, intent(in) :: x, y

    f = x**2 + y**2
end function

end program test_srfpak
