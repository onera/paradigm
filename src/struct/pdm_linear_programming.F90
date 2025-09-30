!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2025  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_linear_programming

  use pdm
  use iso_c_binding

  implicit none

  !!
  !! Enum type PDM_lp_status_t
  !!
  integer(c_int), parameter :: PDM_LP_FEASIBLE   = 0
  integer(c_int), parameter :: PDM_LP_UNFEASIBLE = 1
  integer(c_int), parameter :: PDM_LP_UNBOUNDED  = 2

  contains

  subroutine PDM_lp_solve_nd(dim,  &
                             n,    &
                             a,    &
                             b,    &
                             l,    &
                             u,    &
                             c,    &
                             x,    &
                             stat)
    ! Solve the d-dimensional linear optimization problem
    !   maximize c.x
    !   subject to constraints ai.x <= bi
    !                          l <= x <= u
    !
    ! .. note:: The matrix ``a`` is defined in row-major (C) order, i.e. ``a_{i,j} = a(dim*(i-1)+j)`` (1 <= i <= n, 1 <= j <= dim)
    implicit none

    integer,          intent(in)    :: dim  ! Dimension
    integer,          intent(in)    :: n    ! Number of inequality constraints
    real(8), pointer, intent(in)    :: a(:) ! a in ax <= b (size = n * dim)
    real(8), pointer, intent(in)    :: b(:) ! b in ax <= b (size = n)
    real(8), pointer, intent(in)    :: l(:) ! Lower bounds l <= x (size = dim)
    real(8), pointer, intent(in)    :: u(:) ! Upper bounds x <= u (size = dim)
    real(8), pointer, intent(in)    :: c(:) ! Constant in the objective function (size = dim)
    real(8), pointer, intent(inout) :: x(:) ! b in ax <= b (size = dim)
    integer,          intent(out)   :: stat ! Problem status

    interface
      function PDM_lp_solve_nd_c(dim,  &
                                 n,    &
                                 a,    &
                                 b,    &
                                 l,    &
                                 u,    &
                                 c,    &
                                 x)    &
      result (stat)                    &
      bind (c, name="PDM_lp_solve_nd")
        use iso_c_binding
        implicit none
        integer(c_int),   value :: dim
        integer(c_int),   value :: n
        type(c_ptr),      value :: a
        type(c_ptr),      value :: b
        type(c_ptr),      value :: l
        type(c_ptr),      value :: u
        type(c_ptr),      value :: c
        type(c_ptr),      value :: x
        integer(c_int)          :: stat
      end function PDM_lp_solve_nd_c

    end interface

    stat = PDM_lp_solve_nd_c(dim,      &
                             n,        &
                             c_loc(a), &
                             c_loc(b), &
                             c_loc(l), &
                             c_loc(u), &
                             c_loc(c), &
                             c_loc(x))

  end subroutine PDM_lp_solve_nd

end module pdm_linear_programming