#include "pdm_configf.h"

module pdm_convex

  use pdm
  use iso_c_binding

  implicit none

  contains

    subroutine PDM_points_inside_convex_hull(dim,        &
                                             n_src,      &
                                             src_coord,  &
                                             n_tgt,      &
                                             tgt_coord,  &
                                             tgt_status)
    ! Classify target points with respect to the convex hull of given source points.
    !
    ! ..warning:: Coordinates must always be defined in dimension 3, even if ``dim`` is lower
    implicit none

    integer,                       intent(in)    :: dim            ! Spatial dimension (<= 3)
    integer,                       intent(in)    :: n_src          ! Number of source points
    real(8),              pointer, intent(in)    :: src_coord(:,:) ! Coordinates of source points (shape = [3, n_src])
    integer,                       intent(in)    :: n_tgt          ! Number of target points
    real(8),              pointer, intent(in)    :: tgt_coord(:,:) ! Coordinates of target points (shape = [3, n_tgt])
    integer(pdm_l_num_s), pointer, intent(inout) :: tgt_status(:)  ! Status of each target point (0 = outside, 1 = inside) (size = n_tgt)

    interface
      subroutine PDM_points_inside_convex_hull_c(dim,        &
                                                 n_src,      &
                                                 src_coord,  &
                                                 n_tgt,      &
                                                 tgt_coord,  &
                                                 tgt_status) &
      bind (c, name="PDM_points_inside_convex_hull")
        use iso_c_binding
        integer(c_int), value :: dim
        integer(c_int), value :: n_src
        type(c_ptr),    value :: src_coord
        integer(c_int), value :: n_tgt
        type(c_ptr),    value :: tgt_coord
        type(c_ptr),    value :: tgt_status
      end subroutine PDM_points_inside_convex_hull_c
    end interface

    if (n_src == 0 .or. n_tgt == 0) then
      return
    endif

    if (.not.associated(src_coord)) then
      print *, "PDM_points_inside_convex_hull : src_coord is not associated"
      stop
    endif

    if (.not.associated(tgt_coord)) then
      print *, "PDM_points_inside_convex_hull : tgt_coord is not associated"
      stop
    endif

    if (.not.associated(tgt_status)) then
      print *, "PDM_points_inside_convex_hull : tgt_status is not associated"
      stop
    endif

    call PDM_points_inside_convex_hull_c(dim,               &
                                         n_src,             &
                                         c_loc(src_coord),  &
                                         n_tgt,             &
                                         c_loc(tgt_coord),  &
                                         c_loc(tgt_status))

  end subroutine PDM_points_inside_convex_hull

end module pdm_convex
