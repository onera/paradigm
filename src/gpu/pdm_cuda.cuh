#ifndef __PDM_CUDA_CUH__
#define __PDM_CUDA_CUH__

/*============================================================================
 * CUDA functions
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2017 ONERA

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

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */

#include <stdarg.h>

/* BFT library headers */

#include "pdm_config.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/


//Kernel reduction
#define Reduce_kernel(n_threads, n_blocks, shared_size, blockSize, length, data_in, data_out) {reduce6<<<n_blocks,n_threads,shared_size>>>(data_in, data_out, blockSize, length);}
inline __device__ void warpReduce
(
  volatile int *sdata, 
  unsigned int tid,
  int blockSize
  ) 
{
  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
  if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
  if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
  if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
  if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}

inline __global__ void reduce6
(
  int *g_idata, 
  int *g_odata,
  int blockSize,
  int length
  ) 
{
  extern __shared__ int sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  if ((i >= length))  
  {
    sdata[tid] = 0;
  }
  else if ((i + blockDim.x) >= length)
  {
    sdata[tid] = g_idata[i];
  }
  else
  {
    sdata[tid] = g_idata[i] + g_idata[i + blockDim.x];
  }
  __syncthreads();


  if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }

  if (tid < 32) warpReduce(sdata, tid, blockSize);

  if (tid == 0) {
    g_odata[blockIdx.x] = sdata[0];
    printf("sendcount = %d\n", g_odata[blockIdx.x]);
  }
}


/*============================================================================
 * Public types
 *============================================================================*/

/* Function pointers for printf() and fflush(stdout) type functions */

typedef int (PDM_printf_proxy_t) (const char     *const format,
                                  va_list               arg_ptr);

typedef int (PDM_printf_flush_proxy_t) (void);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Replacement for printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by PDM_printf_proxy_set().
 *
 * parameters:
 *   format: <-- format string, as printf() and family.
 *   ... :   <-- variable arguments based on format string.
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0' used
 *   to end output strings
 */
 
__device__
int
PDM_printf_GPU(const char  *const format,
           ...);

/*
 * Flush for output of PDM_printf() with modifiable behavior.
 *
 * This function calls fflush(stdout) if PDM_printf()'s default behavior is
 * used. If PDM_printf's behavior is modified with PDM_printf_proxy_set(),
 * PDM_printf_flush()'s behavior may have to be also adjusted with
 * PDM_printf_flush_proxy_set().
 *
 * returns:
 *   using the default behavior, the return value is that of
 *   fflush(stdout): O upon successful completion, EOF otherwise
 *   (with errno set to indicate the error).
 */

__device__
int
PDM_printf_flush_GPU(void);

/*
 * Returns function associated with the PDM_printf() function.
 *
 * returns:
 *   pointer to the vprintf() or replacement function.
 */

__device__
PDM_printf_proxy_t *
PDM_printf_proxy_get_GPU(void);

/*
 * Associates a vprintf() type function with the PDM_printf() function.
 *
 * parameters:
 *   fct: <-- pointer to a vprintf() type function.
 */

__device__
void
PDM_printf_proxy_set_GPU(PDM_printf_proxy_t  *const fct);

/*
 * Returns function associated with PDM_printf_flush().
 *
 * returns:
 *   pointer to the PDM_printf_flush() proxy.
 */

__device__
PDM_printf_flush_proxy_t *
PDM_printf_flush_proxy_get_GPU(void);

/*
 * Associates a proxy function with PDM_printf_flush().
 *
 * warning:
 *   PDM_printf() is called by the default PDM_error() error handler
 *   (so as to ensure that the error text appears at the end of the
 *   program output), so a PDM_print_flush replacement must not itself
 *   call (directly or indirectly) PDM_error() if the default error
 *   handler is used.
 *
 * parameter:
 *   fct <-- pointer to a function similar to {return fflush(stdout)}.
 */

__device__
void
PDM_printf_flush_proxy_set_GPU(PDM_printf_flush_proxy_t  *const fct);

/*
 * Realloc implementation for cuda allocation
 *
 * parameters:
 *  ptr:        pointer to realloc
 *  oldLength:  length of the pointer to realloc
 *  newLength:  new length to allocate
 * 
 * returns:
 *  reallocated pointer
 */

inline
#ifdef __CUDACC__
__host__ __device__
#endif
void*
cudaRealloc(void* ptr, 
            size_t oldLength, 
            size_t newLength);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_CUDA_CUH_ */