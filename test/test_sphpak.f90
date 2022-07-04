! test_sphpak.f90 --
!     Straightforward test program for SPHPAK
!
!     TODO:
!     - Conversion to radians in methods (extra argument)?
!     - Better understand the sort of function - in what way linear? Distance or angle?
!       (via a more complicated function or via a line of points along the meridian)
!
program test_sphpak
    use sphpak

    implicit none

    type(network_sphere) :: network

    real, parameter   :: torad = 3.1415926 / 180.0

    real              :: lat(9)
    real              :: lon(9)
    real              :: w(9)

    real, allocatable :: wxwywz(:,:)
    real, allocatable :: latc(:), lonc(:), wc(:,:)

    real              :: result, xc, yc, dx, dy, dz
    integer           :: i

    !
    ! North pole and great circle at 60 degrees - 8 points along the circle
    !
    lat(1)  = 90.0 ; lon(1)  = 0.0
    lat(2:) = 60.0 ; lon(2:) = [ (360.0 / 8.0 * (i-1), i = 1,8) ]

    call network%create_mesh( torad*lat, torad*lon )

    w    = 1.0
    w(1) = 0.0

    call network%set_values( w )

    write(*,*) 'Point (60,0):    ', network%interpolate_point( torad*60.0, torad*0.0 )
    write(*,*) 'Point (60,15):   ', network%interpolate_point( torad*60.0, torad*15.0 )
    write(*,*) 'Point (75,15):   ', network%interpolate_point( torad*75.0, torad*15.0 )
    write(*,*) 'Point (75,-15):  ', network%interpolate_point( torad*75.0,-torad*15.0 )
    write(*,*) 'Point (75,60):   ', network%interpolate_point( torad*75.0, torad*(15.0 + 360.0/8.0) )
    write(*,*) 'Point (90,60):   ', network%interpolate_point( torad*90.0, torad*(15.0 + 360.0/8.0) )
    write(*,*) 'Point (60,90):   ', network%interpolate_point( torad*60.0, torad*90.0 )
    write(*,*) 'Point (90,90):   ', network%interpolate_point( torad*90.0, torad*90.0 )
    write(*,*) 'Point (90,0):    ', network%interpolate_point( torad*90.0, torad*0.0 )

    ! Gradient
    write(*,*) 'Local gradients:'
    do i = 1,5
        call network%gradient_local( i, dx, dy, dz )
        write(*,*) i, lat(i), lon(i), dx, dy, dz
    enddo

    write(*,*) 'Global gradients:'
    call network%gradient_global( wxwywz )
    do i = 1,5
        write(*,*) i, lat(i), lon(i), wxwywz(:,i)
    enddo

    !
    ! Grids
    !
    allocate( latc(4), lonc(4), wc(4,4) )
    call coord_range( latc, 60.0,  90.0 )
    call coord_range( lonc,  0.0, 180.0 )

    call network%interpolate_grid( torad*latc, torad*lonc, wc )
    write(*,'(a10,4f10.3)') 'lon/lat', latc
    write(*,'(5f10.3)') (lonc(i), wc(:,i) , i = 1,size(lonc))

    !
    ! Check with uniform function
    w    = 1.0

    call network%set_values( w )

    write(*,*) 'NOTE: all values should be 1'

    write(*,*) 'Point (60,0):    ', network%interpolate_point( torad*60.0, torad*0.0 )
    write(*,*) 'Point (60,15):   ', network%interpolate_point( torad*60.0, torad*15.0 )
    write(*,*) 'Point (75,15):   ', network%interpolate_point( torad*75.0, torad*15.0 )

contains

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

end program test_sphpak
