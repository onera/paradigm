#include "pdm_configf.h"

program testf

  use iso_c_binding
  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_io
  use pdm_writer
  use pdm_sphere_surf_gen
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer, parameter                :: comm        = MPI_COMM_WORLD
  character(len=99)                 :: arg
  integer                           :: i
  integer                           :: n_part      = 1
  integer                           :: part_method = 3

  integer(pdm_g_num_s)              :: n        = 5
  double precision                  :: x_center = 1.d0
  double precision                  :: y_center = 2.d0
  double precision                  :: z_center = 3.d0
  double precision                  :: radius   = 4.d0

  integer(pdm_l_num_s), pointer     :: pn_vtx(:)  => null()
  type(PDM_pointer_array_t), pointer:: pvtx_coord => null()
  type(PDM_pointer_array_t), pointer:: pvtx_ln_to_gn => null()
  integer(pdm_l_num_s), pointer     :: pn_face(:) => null()
  type(PDM_pointer_array_t), pointer:: pface_vtx_idx => null()
  type(PDM_pointer_array_t), pointer:: pface_vtx => null()
  type(PDM_pointer_array_t), pointer:: pface_ln_to_gn => null()

  ! double precision,     pointer     :: vtx_coord(:)     => null()
  double precision,     pointer     :: vtx_coord2(:,:)  => null()
  integer(pdm_g_num_s), pointer     :: vtx_ln_to_gn(:)  => null()
  integer(pdm_l_num_s), pointer     :: face_vtx(:)      => null()
  integer(pdm_g_num_s), pointer     :: face_ln_to_gn(:) => null()

  ! Writer
  type(c_ptr)                       :: wrt
  integer                           :: id_geom
  integer                           :: id_block

  ! MPI
  integer                           :: code
  integer                           :: i_rank
  integer                           :: n_rank

  integer                           :: ipart

  !-----------------------------------------------------------
  !                Read command line arguments
  !-----------------------------------------------------------

  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    select case(arg)
      case ("-n")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) n
      case ("-n_part")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) n_part
      case ("-parmetis")
        part_method = PDM_SPLIT_DUAL_WITH_PARMETIS
      case ("-pt-scotch")
        part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH
      case ("-hilbert")
        part_method = PDM_SPLIT_DUAL_WITH_HILBERT
      case default
        print *, "Invalid command argument ", arg
        stop
    end select

    i = i + 1
  enddo

  !-----------------------------------------------------------
  !                  MPI init
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  !-----------------------------------------------------------
  !                Generate cube
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  !                Generate sphere
  !-----------------------------------------------------------

 
 
  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
