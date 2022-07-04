! srfpak_c.f90 --
!     C interface for the Fortran routines and methods of SRFPAK
!
!     Note:
!     The header file srfpak.h defines the network type and the
!     interfaces for the C side.
!
module srfpak_c
    use iso_c_binding
    use srfpak

    implicit none

    type, bind(C) :: network_2d_c
        type(c_ptr) :: network_f = c_null_ptr
    end type network_2d_c

contains
subroutine create_mesh_2d_c( network, n, x, y ) bind(C, name = "createMesh" )
    type(network_2d_c), intent(inout) :: network
    integer, value                    :: n
    real, intent(in)                  :: x(n), y(n)

    type(network_2d), pointer         :: network_f

    allocate( network_f )
    network%network_f = c_loc( network_f )

    call network_f%create_mesh( x, y )
end subroutine create_mesh_2d_c

subroutine delete_mesh_2d_c( network ) bind(C, name = "deleteMesh" )
    type(network_2d_c), intent(inout) :: network

    type(network_2d), pointer         :: network_f

    if ( c_associated( network%network_f ) ) then
        call c_f_pointer( network%network_f, network_f )

        if ( allocated( network_f%x ) ) then
            deallocate( network_f%x )
            deallocate( network_f%y )
            deallocate( network_f%z )
            deallocate( network_f%iadj )
            deallocate( network_f%iend )
        endif

        network%network_f = c_null_ptr
    endif
end subroutine delete_mesh_2d_c

subroutine set_values_2d_c( network, z ) bind(C, name = "setMeshValues" )
    type(network_2d_c), intent(inout) :: network
    real, intent(in), target          :: z(*)

    type(network_2d), pointer         :: network_f
    real, pointer                     :: zp(:)
    integer                           :: n

    call c_f_pointer( network%network_f, network_f )

    n = size(network_f%x)

    zp => z(1:n)

    call network_f%set_values( zp )
end subroutine set_values_2d_c

real function interpolate_point_2d_c( network, px, py ) bind(C, name = "interpolatePoint" )
    type(network_2d_c), intent(inout) :: network
    real, value                       :: px, py

    type(network_2d), pointer         :: network_f

    call c_f_pointer( network%network_f, network_f )

    interpolate_point_2d_c = network_f%interpolate_point( px, py )
end function interpolate_point_2d_c

real function volume_2d_c( network ) bind(C, name = "meshVolume" )
    type(network_2d_c), intent(inout) :: network

    type(network_2d), pointer         :: network_f

    call c_f_pointer( network%network_f, network_f )

    volume_2d_c = network_f%volume()
end function volume_2d_c

end module srfpak_c
