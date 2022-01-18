module pdm
#include "pdmf.h"

interface

function PDM_MPI_Comm_f2c (f_comm) &
 result (c_comm)                   &

 bind (c, name = 'PDM_MPI_Comm_f2c')

 use iso_c_binding
 implicit none

 integer(c_int), value :: f_comm
 integer(c_int)        :: c_comm

 end function PDM_MPI_Comm_f2c

 end interface

end module pdm
