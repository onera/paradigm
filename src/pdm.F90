module pdm
#include "pdmf.h"

integer, parameter :: PDM_FALSE = 0
integer, parameter :: PDM_TRUE  = 1

interface

function PDM_MPI_Comm_f2c (f_comm) &
result (c_comm)                    &
bind (c, name = 'PDM_MPI_Comm_f2c')

  use iso_c_binding
  implicit none

  integer(c_int), value :: f_comm
  type(c_ptr)           :: c_comm

end function PDM_MPI_Comm_f2c

function PDM_MPI_Comm_c2f (c_comm) &
result (f_comm)                    &
bind (c, name = 'PDM_MPI_Comm_c2f')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: c_comm
  integer(c_int)     :: f_comm

end function PDM_MPI_Comm_c2f

function PDM_MPI_Type_f2c (f_datatype) &
result (datatype)                      &
bind (c, name = 'PDM_MPI_Type_f2c')

    use iso_c_binding
    implicit none

    integer(c_int), value :: f_datatype
    type(c_ptr)           :: datatype

end function PDM_MPI_Type_f2c

function PDM_MPI_Op_f2c (f_op) &
  result (op)                  &
  bind (c, name = 'PDM_MPI_Op_f2c')

    use iso_c_binding
    implicit none

    integer(c_int), value :: f_op
    type(c_ptr)           :: op

end function PDM_MPI_Op_f2c

end interface

end module pdm
