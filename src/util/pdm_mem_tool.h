/*
 * \file
 */


#ifndef __PDM_MEM_TOOL_H__
#define __PDM_MEM_TOOL_H__

/*--------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Headers for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"


/*----------------------------------------------------------------------------
 *  Macro definition
 *----------------------------------------------------------------------------*/

//
// First stage of integration of the dev_memory_analysis branch (Langloff internship)
//
// creation of a pdm_realloc macro to fix the problem 
// about the return return of a pointer allocated at 0 to NULL
//


//#define MAX_CALL_TRACE 20
//#define BUFFER_SIZE 1048576

//#define _MEM_ANALYZER_

//#define _OPTIMIZE_
//#define _SUMMARY_
//#define SANITIZE
//#define NO_MPI // incompatible avec _SUMMARY_


// #ifdef _MEM_ANALYZER_
// 	#define PDM_malloc(_ptr,_ni,_type) 
// 		_ptr = (_type*) PDM_MEM_MALLOC(_ni,sizeof(_type),__FILE__,__LINE__,#_ptr)

// 	#define PDM_realloc(_old_ptr,_new_ptr,_ni,_type) 
// 		_new_ptr = (_type*) PDM_MEM_REALLOC(_old_ptr,_new_ptr,_ni,sizeof(_type),__FILE__,__LINE__,#_new_ptr)

// 	#define PDM_calloc(_ptr,_ni,_type) 
// 		_ptr = (_type*) PDM_MEM_CALLOC(_ni,sizeof(_type),__FILE__,__LINE__,#_ptr)
	
// 	#define PDM_free(_ptr) 
// 		PDM_MEM_FREE(_ptr,__FILE__,__LINE__,#_ptr)

//  	#ifdef _OPTIMIZE_
// 		#define backtrace(_arr,MAX_CALL_TRACE) 0
// 		#define backtrace_symbols(_arr,_size) NULL
// 	#endif
// #else
	#define PDM_malloc(_ptr,_ni,_type) \
		_ptr = (_type*) PDM_MALLOC(_ni,sizeof(_type),__func__,__FILE__,__LINE__,#_ptr)
	
	#define PDM_realloc(_old_ptr,_new_ptr,_ni,_type) \
		_new_ptr = (_type*) PDM_REALLOC(_old_ptr,_ni,sizeof(_type),__func__,__FILE__,__LINE__,#_new_ptr)
	
//		_new_ptr = (_type*) PDM_REALLOC(_old_ptr,_new_ptr,_ni,sizeof(_type),__func__,__FILE__,__LINE__,#_new_ptr)
	#define PDM_calloc(_ptr,_ni,_type) \
		_ptr = (_type*) PDM_CALLOC(_ni,sizeof(_type),__func__,__FILE__,__LINE__,#_ptr)
	
	#define PDM_free(_ptr) \
		PDM_FREE((void ** ) &(_ptr))

// #endif



/*============================================================================
 * Types
 *============================================================================*/

// #ifndef NO_MPI
// typedef struct PDM_mempeak_t {

// 	unsigned long long int peak;
//   double time_peak;
//   char cs[MAX_CALL_TRACE][9];
//   int size_cs;


// 	unsigned long int nb_alloc;
// 	unsigned long int nb_free;
// 	unsigned long long int total_alloc;
// 	unsigned long long int total_free;
// 	double total_time;
// } PDM_mempeak_t;
// #endif

// typedef struct PDM_mem_tool_t {
// 	#ifndef _SUMMARY_ 
// 	// local file
// 	const char *local_proc_file;
// 	int fd_local_proc;
// 	char buffer[BUFFER_SIZE];
// 	int index_buffer;
// 	#endif

// 	// boolean
// 	int is_running;
	
// 	#ifndef NO_MPI
// 	int is_mpi_finished;

// 	// MPI
// 	int rank;
// 	PDM_MPI_Comm comm;
// 	int size_mpi;

// 	// memory
// 	PDM_mempeak_t mempeak;
// 	unsigned long long int real_time_mem;
// 	char **cs_peak;
// 	#endif
// } PDM_mem_tool_t;


/*=============================================================================
 * Private function prototypes
 *============================================================================*/

// /**
//  *
//  * \brief Memory error handler
//  *
//  * \param [in]  ptr    		Pointer to check
//  * \param [in]  file_name	Source file of the pointer __FILE__
//  * \param [in]  line		Declaration line of the pointer __LINE__
//  *
//  */

// void PDM_mem_error
// (
//  void *ptr,
//  const char *file_name,
//  unsigned int line 
// );


// /**
//  *
//  * \brief Build a line for allocated pointer
//  *
//  * \param [in]  ptr		Pointer allocated
//  * \param [in]  var_name	Name of the pointer
//  * \param [in]  alloc_used	Requested allocate size
//  * \param [in]  chunk_size	Size given by glibc 
//  * \param [in] 	file_name	Source file of the pointer __FILE__
//  * \param [in]  line		Source line of the pointer __LINE__
//  * \param [in]  size_cs		Size of the call stack  
//  * \param [in]  cs		Array of string for the call stack 
//  * \param [in]  time		Time at the allocation  
//  * \param [in]  alloc_type	Allocation type M(malloc),R(realloc),C(calloc) 
//  *
//  */


// char *
// _buffer_line_malloc
// (
//  void *ptr,
//  const char *var_name,
//  unsigned long int alloc_used,
//  unsigned long int chunk_size,
//  const char *file_name,
//  unsigned int line,
//  int size_cs,
//  char **cs,
//  double time,
//  const char *alloc_type
// );

// /**
//  *
//  * \brief Build a line for freed pointer
//  *
//  * \param [in]  ptr		Pointer allocated
//  * \param [in]  var_name	Name of the pointer
//  * \param [in]  chunk_size	Size given by glibc 
//  * \param [in] 	file_name	Source file of the pointer __FILE__
//  * \param [in]  line		Source line of the pointer __LINE__
//  * \param [in]  size_cs		Size of the call stack  
//  * \param [in]  cs		Array of string for the call stack 
//  * \param [in]  time		Time at the allocation  
//  *
//  */


// char *
// _buffer_line_free
// (
//  void *ptr,
//  const char *var_name,
//  unsigned long int chunk_size,
//  const char *file_name,
//  unsigned int line,
//  int size_cs,
//  char **cs,
//  double time
// );

// /**
//  *
//  * \brief Write line to the buffer
//  *
//  * \param [in]  line_buffer	Line for the buffer
//  *
//  */

// void 
// _write_to_buffer
// (
//  char *line_buffer
// );

// /**
//  *
//  * \brief Pick the address function 
//  *
//  * \param [in]  path	Line output from backtrace
//  *
//  */

// char *_get_id_fct
// (
//  char *path
// );


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

// /**
//  *
//  * \brief Inititialize parameter of the pdm_mem_tool
//  *
//  * \param [in]  local_proc_file	Name of the output memory file for 1 MPI rank
//  * \param [in]  rank		MPI rank
//  * \param [in]	size_world	MPI world size
//  * \param [in]  comm		MPI communicator 
//  *
//  */

// void PDM_MEM_INIT(
// #ifdef _MEM_ANALYZER_
//  	const char *local_proc_file
// 	#ifndef NO_MPI
//  		, int rank, 
//    		int size_world,
//    		PDM_MPI_Comm comm
// 	#endif
// #else
// 	void
// #endif
// );

// /**
//  *
//  * \brief Start timer and recording
//  *
//  */

// void PDM_MEM_START(void);

// /**
//  *
//  * \brief Gather rank peak data before MPI communication end
//  *
//  */

// void PDM_MEM_FINALIZE(void);

// /**
//  *
//  * \brief End recording and write buffer to file if buffer is not empty
//  *
//  */

// void PDM_MEM_END(void);

// /**
//  *
//  * \brief Stop recording, timer continue
//  *
//  */

// void PDM_MEM_STOP(void);

// /**
//  *
//  * \brief Resume recording
//  *
//  */

// void PDM_MEM_RESUME(void);

// /**
//  *
//  * \brief Allocate a pointer while recording
//  *
//  * \param [in]  nb		Number of elements
//  * \param [in]  size_type	Size type of the elements
//  * \param [in]  file_name 	Source file of the pointer __FILE__
//  * \param [in]  line		Source line __LINE__ 
//  * \param [out] var_name	Pointer name
//  *
//  */

// void *PDM_MEM_MALLOC
// (
//  size_t nb,
//  size_t size_type,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// );

// /**
//  *
//  * \brief Reallocate a pointer while recording
//  *
//  * \param [in]  old_ptr		Old pointer to reallocate
//  * \param [in]  new_ptr		New pointer to be reallocated
//  * \param [in]  nb		Number of elements
//  * \param [in]  size_type	Size type of the elements
//  * \param [in]  file_name 	Source file of the pointer __FILE__
//  * \param [in]  line		Source line __LINE__ 
//  * \param [out] var_name	Pointer name
//  *
//  */


// void *PDM_MEM_REALLOC
// (
//  void *old_ptr,
//  void *new_ptr,
//  size_t nb,
//  size_t size_type,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// );

// /**
//  *
//  * \brief Allocate a pointer while recording
//  *
//  * \param [in]  nb		Number of elements
//  * \param [in]  size_type	Size type of the elements
//  * \param [in]  file_name 	Source file of the pointer __FILE__
//  * \param [in]  line		Source line __LINE__ 
//  * \param [out] var_name	Pointer name
//  *
//  */


// void *PDM_MEM_CALLOC
// (
//  size_t nb,
//  size_t size_type,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// );

// /**
//  *
//  * \brief Free a pointer while recording
//  *
//  * \param [in]  ptr		Pointer to be freed
//  * \param [in]  file_name 	Source file of the pointer __FILE__
//  * \param [in]  line		Source line __LINE__ 
//  * \param [out] var_name	Pointer name
//  *
//  */


// void PDM_MEM_FREE
// (
//  void *ptr,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// );

/**
 *
 * \brief Allocate a pointer whithout recording
 *
 * \param [in]  nb		Number of elements
 * \param [in]  size_type	Size type of the elements
 * \param [in]  func_name	Source function __func__
 * \param [in]  file_name 	Source file of the pointer __FILE__
 * \param [in]  line		Source line __LINE__ 
 * \param [out] var_name	Pointer name
 *
 */


void *PDM_MALLOC
(
 size_t nb,
 size_t size_type,
 const char *func_name,
 const char *file_name,
 unsigned int line,
 const char *var_name
);

/**
 *
 * \brief Reallocate a pointer whithout recording
 *
 * \param [in]  old_ptr		Old pointer to reallocate
 * \param [in]  new_ptr		New pointer to be reallocated
 * \param [in]  nb		Number of elements
 * \param [in]  size_type	Size type of the elements
 * \param [in]  func_name	Source function __func__
 * \param [in]  file_name 	Source file of the pointer __FILE__
 * \param [in]  line		Source line __LINE__ 
 * \param [out] var_name	Pointer name
 *
 */


void *PDM_REALLOC
(
 void *old_ptr,
// void *new_ptr,
 size_t nb,
 size_t size_type,
 const char *func_name,
 const char *file_name,
 unsigned int line,
 const char *var_name
);

/**
 *
 * \brief Allocate a pointer whithout recording
 *
 * \param [in]  nb		Number of elements
 * \param [in]  size_type	Size type of the elements
 * \param [in]  func_name	Source function __func__
 * \param [in]  file_name 	Source file of the pointer __FILE__
 * \param [in]  line		Source line __LINE__ 
 * \param [out] var_name	Pointer name
 *
 */


void *PDM_CALLOC
(
 size_t nb,
 size_t size_type,
 const char *func_name,
 const char *file_name,
 unsigned int line,
 const char *var_name
);

/**
 *
 * \brief Free a pointer whithout recording
 *
 * \param [in]  ptr		Pointer to be freed
 *
 */


void PDM_FREE
(
 void **ptr
);


#endif /*__PDM_MEM_TOOL_H__*/
