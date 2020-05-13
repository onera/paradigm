/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_cuda_print.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif


/*============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/



/*============================================================================
 * Global variable
 *============================================================================*/


/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
__global__
void
_print_from_gpu
(
 void
 )
{
  printf("Hello World! from thread [%d,%d] From device\n", threadIdx.x,blockIdx.x);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


void
CUDA_print
(
 void
 )
{
  _print_from_gpu<<<1,10>>>();
  cudaDeviceSynchronize();
}


#ifdef	__cplusplus
}
#endif

