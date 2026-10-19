#include "pdm_configf.h"

module pdm_fortran

  use pdm

  implicit none

  integer, parameter :: PDM_TYPE_INT    = 0
#ifdef PDM_LONG_G_NUM
  integer, parameter :: PDM_TYPE_G_NUM  = 1
#else
  integer, parameter :: PDM_TYPE_G_NUM  = 0
#endif
  integer, parameter :: PDM_TYPE_DOUBLE   = 2
  integer, parameter :: PDM_TYPE_COMPLEX8 = 3
  integer, parameter :: PDM_TYPE_COMPLEX4 = 4
  integer, parameter :: PDM_TYPE_REAL4    = 5
  integer, parameter :: PDM_TYPE_CPTR     = 6

  interface

    subroutine pdm_fortran_free_c (ptrC) &
      bind (c, name = 'free')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: ptrC


    end subroutine pdm_fortran_free_c

  end interface

end module pdm_fortran
