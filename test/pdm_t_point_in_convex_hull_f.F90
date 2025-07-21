#include "pdm_configf.h"

program point_in_convex_hull_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_linear_programming
  use pdm_vtk
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !---------------------------------------------------------------
  integer                       :: dim            ! Spatial dimension
  integer                       :: n_src          ! Number of source points
  integer                       :: n_tgt          ! Number of target points
  integer                       :: seed           ! Random seed
  logical                       :: visu           ! Enable export for visualization
  real(8),              pointer :: src_coord(:,:) ! Coordinates of source points
  real(8),              pointer :: tgt_coord(:,:) ! Coordinates of target points
  integer(pdm_l_num_s), pointer :: tgt_status(:)  ! Status of target points
  integer(pdm_g_num_s), pointer :: ln_to_gn(:)    ! Not used
  integer(pdm_l_num_s), pointer :: color(:)       ! Not used
  integer                       :: rnd_n
  integer, allocatable          :: rnd_seed(:)
  !---------------------------------------------------------------
  
  ! Default values
  dim   = 2
  n_src = 10
  n_tgt = 1000
  seed  = 0
  visu  = .false.

  ! Parse command line arguments
  call read_args(dim,   &
                 n_src, &
                 n_tgt, &
                 seed,  &
                 visu)

  ! Initialize random generator
  call random_seed(size=rnd_n)
  allocate(rnd_seed(rnd_n))
  rnd_seed(:) = seed
  call random_seed(put=rnd_seed)
  deallocate(rnd_seed)

  ! Generate source points
  call random_points(dim,       &
                     n_src,     &
                     -1.d0,     &
                      1.d0,     &
                     src_coord)

  ! Generate target points
  call random_points(dim,       &
                     n_tgt,     &
                     -2.d0,     &
                      2.d0,     &
                     tgt_coord)

  ! Classify target points w.r.t. convex hull of source points
  allocate(tgt_status(n_tgt))
  call PDM_lp_pts_inside_convex_hull(dim,        &
                                     n_src,      &
                                     src_coord,  &
                                     n_tgt,      &
                                     tgt_coord,  &
                                     tgt_status)

  if (visu) then
    ! Export for visualization
    ln_to_gn => null()
    color    => null()
    call PDM_vtk_write_point_cloud("point_in_convex_hull_f_src.vtk", &
                                   n_src,                            &
                                   src_coord,                        &
                                   ln_to_gn,                         &
                                   color)

    call PDM_vtk_write_point_cloud("point_in_convex_hull_f_tgt.vtk", &
                                   n_tgt,                            &
                                   tgt_coord,                        &
                                   ln_to_gn,                         &
                                   tgt_status)
  endif

  ! Free memory
  deallocate(src_coord, &
             tgt_coord, &
             tgt_status)

  contains

  subroutine read_args(dim,   &
                       n_src, &
                       n_tgt, &
                       seed,  &
                       visu)
    ! Parse command line arguments
    implicit none
    integer, intent(inout) :: dim
    integer, intent(inout) :: n_src
    integer, intent(inout) :: n_tgt
    integer, intent(inout) :: seed
    logical, intent(inout) :: visu
    integer                :: i_arg
    character(len=99)      :: arg

    i_arg = 1
    do while (i_arg <= command_argument_count())
      call get_command_argument(i_arg, arg)
      select case(arg)

        case ("-dim")
          i_arg = i_arg + 1
          call get_command_argument(i_arg, arg)
          read(arg, *) dim

        case ("-n_src")
          i_arg = i_arg + 1
          call get_command_argument(i_arg, arg)
          read(arg, *) n_src

        case ("-n_tgt")
          i_arg = i_arg + 1
          call get_command_argument(i_arg, arg)
          read(arg, *) n_tgt

        case ("-seed")
          i_arg = i_arg + 1
          call get_command_argument(i_arg, arg)
          read(arg, *) seed

        case ("-visu")
          visu = .true.

      endselect

      i_arg = i_arg + 1
    enddo
  end subroutine read_args


  subroutine random_points(dim,   &
                           n_pts, &
                           lo,    &
                           hi,    &
                           coord)
    ! Generate random points uniformly distributed in the cuboid [lo,hi]^dim
    implicit none
    integer,          intent(in)  :: dim        ! Spatial dimension
    integer,          intent(in)  :: n_pts      ! Number of points
    real(8),          intent(in)  :: lo         ! Lower bound
    real(8),          intent(in)  :: hi         ! Upper bound
    real(8), pointer, intent(out) :: coord(:,:) ! Coordinates of points

    allocate(coord(3, n_pts))
    call random_number(coord(1:dim,1:n_pts))
    coord = lo + (hi-lo)*coord
    coord(dim+1:3,1:n_pts) = 0.d0

  end subroutine random_points

end program point_in_convex_hull_f