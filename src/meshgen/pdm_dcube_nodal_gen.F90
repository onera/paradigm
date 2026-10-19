#include "pdm_configf.h"

module pdm_dcube_nodal_gen

  use iso_c_binding
  use pdm

  implicit none

  contains

  subroutine PDM_dcube_nodal_gen_create(dcube,     &
                                        f_comm,    &
                                        n_x,       &
                                        n_y,       &
                                        n_z,       &
                                        length,    &
                                        xmin,      &
                                        ymin,      &
                                        zmin,      &
                                        elt_type,  &
                                        order,     &
                                        ownership)
    ! Create a distributed nodal cube mesh
    implicit none

    type (c_ptr), intent(out) :: dcube     ! Pointer to PDM_dcube_nodal_t instance
    integer,      intent(in)  :: f_comm    ! Communicator
    integer,      intent(in)  :: n_x       ! Number of elements in x-direction
    integer,      intent(in)  :: n_y       ! Number of elements in y-direction
    integer,      intent(in)  :: n_z       ! Number of elements in z-direction
    real(8),      intent(in)  :: length    ! Length of cube side
    real(8),      intent(in)  :: xmin      ! Minimal x-coordinate
    real(8),      intent(in)  :: ymin      ! Minimal y-coordinate
    real(8),      intent(in)  :: zmin      ! Minimal z-coordinate
    integer,      intent(in)  :: elt_type  ! Element type
    integer,      intent(in)  :: order     ! Element order
    integer,      intent(in)  :: ownership ! instance ownership

    type(c_ptr)                :: c_comm
    integer(c_int)             :: c_n_x, c_n_y, c_n_z
    integer(c_int)             :: c_elt_type, c_order, c_ownership
    real(c_double)             :: c_length
    real(c_double)             :: c_xmin, c_ymin, c_zmin


    interface
      function PDM_dcube_nodal_gen_create_c(comm,      &
                                            n_x,       &
                                            n_y,       &
                                            n_z,       &
                                            length,    &
                                            xmin,      &
                                            ymin,      &
                                            zmin,      &
                                            elt_type,  &
                                            order,     &
                                            ownership) &
        result(dcube) &
        bind (c, name = 'PDM_dcube_nodal_gen_create')

        use iso_c_binding
        implicit none

        type(c_ptr), value    :: comm
        integer(c_int), value :: n_x, n_y, n_z
        integer(c_int), value :: elt_type, order, ownership

        real(c_double), value :: length
        real(c_double), value :: xmin, ymin, zmin

        type (c_ptr) :: dcube

      end function PDM_dcube_nodal_gen_create_c
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    c_n_x       = n_x
    c_n_y       = n_y
    c_n_z       = n_z
    c_length    = length
    c_xmin      = xmin
    c_ymin      = ymin
    c_zmin      = zmin
    c_elt_type  = elt_type
    c_order     = order
    c_ownership = ownership

    dcube = PDM_dcube_nodal_gen_create_c(c_comm,     &
                                         c_n_x,      &
                                         c_n_y,      &
                                         c_n_z,      &
                                         c_length,   &
                                         c_xmin,     &
                                         c_ymin,     &
                                         c_zmin,     &
                                         c_elt_type, &
                                         c_order,    &
                                         c_ownership)

  end subroutine PDM_dcube_nodal_gen_create



  subroutine PDM_dcube_nodal_gen_random_factor_set(dcube,         &
                                                   random_factor)
    ! Set randomization factor
    implicit none

    type(c_ptr), intent(in) :: dcube         ! Pointer to PDM_dcube_nodal_t instance
    real(8),     intent(in) :: random_factor ! Randomization factor (between 0 and 1)

    interface
      subroutine PDM_dcube_nodal_gen_random_factor_set_c(dcube,         &
                                                         random_factor) &
      bind (c, name="PDM_dcube_nodal_gen_random_factor_set")
        use iso_c_binding
        implicit none

        type (c_ptr),   value :: dcube
        real(c_double), value :: random_factor

      end subroutine PDM_dcube_nodal_gen_random_factor_set_c
    end interface

    call PDM_dcube_nodal_gen_random_factor_set_c(dcube, &
                                                 random_factor)

  end subroutine PDM_dcube_nodal_gen_random_factor_set



  subroutine PDM_dcube_nodal_gen_build(dcube, &
                                       dmn)

    ! Generate the mesh
    implicit none

    type(c_ptr), intent(in)  :: dcube ! Pointer to PDM_dcube_nodal_t instance
    type(c_ptr), intent(out) :: dmn   ! Pointer to PDM_mesh_nodal_t instance

    interface
      function PDM_dcube_nodal_gen_build_c(dcube) result(dmn) &
        bind (c, name = 'PDM_dcube_nodal_gen_build')

        use iso_c_binding
        implicit none

        type (c_ptr), value  :: dcube
        type (c_ptr)         :: dmn

      end function PDM_dcube_nodal_gen_build_c
    end interface

    dmn = PDM_dcube_nodal_gen_build_c(dcube)

  end subroutine PDM_dcube_nodal_gen_build



  subroutine PDM_dcube_nodal_gen_dmesh_nodal_get(dcube, &
                                                 dmn)

    ! Get the PDM_dmesh_nodal_t associated to a PDM_dcube_nodal_t
    implicit none

    type(c_ptr), intent(in)  :: dcube ! Pointer to PDM_dcube_nodal_t instance
    type(c_ptr), intent(out) :: dmn   ! Pointer to PDM_mesh_nodal_t instance

    interface
      function PDM_dcube_nodal_gen_dmesh_nodal_get_c(dcube) result(dmn) &
        bind (c, name = 'PDM_dcube_nodal_gen_dmesh_nodal_get')

        use iso_c_binding
        implicit none

        type (c_ptr), value  :: dcube
        type (c_ptr)         :: dmn

      end function PDM_dcube_nodal_gen_dmesh_nodal_get_c
    end interface

    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get_c(dcube)

  end subroutine PDM_dcube_nodal_gen_dmesh_nodal_get


  
  subroutine PDM_dcube_nodal_gen_free(dcube)
    ! Free a PDM_dcube_nodal_gen_t instance
    implicit none

    type (c_ptr), intent(inout) :: dcube

    interface
      subroutine PDM_dcube_nodal_gen_free_c(dcube) &
        bind(c, name='PDM_dcube_nodal_gen_free')

        use iso_c_binding
        implicit none

        type (c_ptr), value :: dcube
      end subroutine PDM_dcube_nodal_gen_free_c
    end interface

    call PDM_dcube_nodal_gen_free_c(dcube)

  end subroutine PDM_dcube_nodal_gen_free

end module pdm_dcube_nodal_gen
