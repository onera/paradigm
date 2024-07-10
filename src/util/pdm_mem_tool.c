/*
 * \file
 */


/*--------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Headers for the current file
 *----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <execinfo.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <mpi.h>
#include <time.h>
#include <malloc.h>

#include "pdm_mem_tool.h"
#include "pdm_timer.h"


/*============================================================================
 * Global variables
 *============================================================================*/

// static PDM_mem_tool_t pdm_mem_tool;
// struct _pdm_timer_t *timer_prog;

// #ifndef NO_MPI
// MPI_Datatype MPI_mempeak_t;
// #endif

/*=============================================================================
 * Private function prototypes
 *============================================================================*/

// /**
//  *
//  * \brief Memory error handler
//  *
//  * \param [in]  ptr             Pointer to check
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Declaration line of the pointer __LINE__
//  *
//  */


// void PDM_mem_error
// (
//  void *ptr,
//  const char *file_name,
//  unsigned int line
// )
// {
// 	if (ptr == NULL) {
// 		printf("allocation fail in %s:%d\n",file_name,line);
// 		#ifndef _SUMMARY_
// 		if (pdm_mem_tool.index_buffer > 0) {
//                 	if (pdm_mem_tool.fd_local_proc < 0) {
//                         	pdm_mem_tool.fd_local_proc = open(pdm_mem_tool.local_proc_file,O_APPEND |  O_RDWR ,0666);
//                 	}
//                 	write(pdm_mem_tool.fd_local_proc,&(pdm_mem_tool.buffer),pdm_mem_tool.index_buffer);
//                 	close(pdm_mem_tool.fd_local_proc);
//         	}
// 		#endif
// 		if (timer_prog != NULL) {
//                 	PDM_timer_hang_on(timer_prog);
//                 	PDM_timer_free(timer_prog);
// 		}

// 		exit(EXIT_FAILURE);
// 	}
// }

// /**
//  *
//  * \brief Pick the address function 
//  *
//  * \param [in]  path    Line output from backtrace
//  *
//  */


// char *
// _get_id_fct
// (
//  char *path
// )
// {
//         char *start_name = strchr(path,'[') + 1;
//         char *end_name = strchr(path,']');

//         int len = (int) (end_name - start_name);

//         char *name = (char *) malloc(len * sizeof(char));
// 	PDM_mem_error(name,__FILE__,__LINE__);

//         if (name != NULL) {
// 		memcpy(name,start_name,len);
//                 name[len] = '\0';

//                 return name;
//         } else {
//                 return NULL;
//         }
// }

// /**
//  *
//  * \brief Build a line for allocated pointer
//  *
//  * \param [in]  ptr             Pointer allocated
//  * \param [in]  var_name        Name of the pointer
//  * \param [in]  alloc_used      Requested allocate size
//  * \param [in]  chunk_size      Size given by glibc 
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Source line of the pointer __LINE__
//  * \param [in]  size_cs         Size of the call stack  
//  * \param [in]  cs              Array of string for the call stack 
//  * \param [in]  time            Time at the allocation  
//  * \param [in]  alloc_type      Allocation type M(malloc),R(realloc),C(calloc) 
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
// )
// {
// 	// pointer's address
// 	char addr_to_str[21] = {0};
// 	sprintf(addr_to_str,"%14ld",(uintptr_t)ptr);
// 	unsigned char size_str_addr = strlen(addr_to_str);

// 	// malloc size used
// 	char alloc_size_used_to_str[21] = {0};
// 	sprintf(alloc_size_used_to_str,"%14ld",alloc_used);
// 	unsigned char size_str_alloc_used_size = strlen(alloc_size_used_to_str);

// 	// chunk size given
// 	char chunk_size_to_str[21] = {0};
// 	sprintf(chunk_size_to_str,"%14ld",chunk_size);
// 	unsigned char size_str_chunk_size = strlen(chunk_size_to_str);

// 	// variable name
// 	unsigned char size_str_varname = strlen(var_name);

// 	// file name
// 	unsigned short int size_str_file_name = strlen(file_name);
	
// 	// line
// 	char line_to_str[21] = {0};
// 	sprintf(line_to_str,"%11d",line);
// 	unsigned char size_str_line = strlen(line_to_str);

// 	// call stack
// 	unsigned short int size_str_call_stack = 9 * (size_cs-2);

// 	char *call_stack_to_str = (char*) calloc(size_str_call_stack,sizeof(char));
// 	PDM_mem_error(call_stack_to_str,__FILE__,__LINE__);

// 	unsigned char nb_fct = 0;
// 	for (int i=size_cs-1;i>0;i--) {
// 		char *function = _get_id_fct(cs[i]);
// 		unsigned char size_str_fct = strlen(function);
		
// 		// suppose function's address use at most 8 bits (on arch 64-bits), 
// 		// ignore shared libraries
// 		if (size_str_fct <= 8) { 			
// 			memcpy(call_stack_to_str + 9 * nb_fct,
// 				function,
// 				size_str_fct
// 				);
// 			if (nb_fct > 0) {
// 				memcpy(call_stack_to_str + 9 * nb_fct -1,":",1);
// 			}
// 			nb_fct++;
// 		}
// 		free(function);
// 	}
// 	size_str_call_stack = strlen(call_stack_to_str);

// 	// time
// 	char time_to_str[30] = {0};
// 	sprintf(time_to_str,"%8.15e",time);
// 	unsigned char size_str_time = strlen(time_to_str);

// 	unsigned int size_line = size_str_addr +
// 				size_str_alloc_used_size +
// 				size_str_chunk_size +
// 				size_str_varname +
// 				size_str_file_name +
// 				size_str_line +
// 				size_str_call_stack +
// 				size_str_time +
// 				11;

// 	char *line_buffer = (char*) calloc(size_line,sizeof(char));
// 	PDM_mem_error(line_buffer,__FILE__,__LINE__);

// 	// buffer line assembly
// 	unsigned short int index = 0;
// 	memcpy(line_buffer+index,addr_to_str,size_str_addr);
// 	index += size_str_addr;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,alloc_type,1);
// 	index += 1;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,chunk_size_to_str,size_str_chunk_size);
// 	index += size_str_chunk_size;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,alloc_size_used_to_str,size_str_alloc_used_size);
// 	index += size_str_alloc_used_size;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,var_name,size_str_varname);
// 	index += size_str_varname;
// 	memcpy(line_buffer+index,"\t",1);
// 	index +=1;
// 	memcpy(line_buffer+index,file_name,size_str_file_name);
// 	index += size_str_file_name;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,line_to_str,size_str_line);
// 	index += size_str_line;
// 	memcpy(line_buffer+index,"\t",1);
// 	index +=1;
// 	memcpy(line_buffer+index,call_stack_to_str,size_str_call_stack);
// 	index += size_str_call_stack;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,time_to_str,size_str_time);
// 	index += size_str_time;
// 	memcpy(line_buffer+index,"\n\0",2);

// 	free(call_stack_to_str);

// 	return line_buffer;
// }


// /**
//  *
//  * \brief Build a line for freed pointer
//  *
//  * \param [in]  ptr             Pointer allocated
//  * \param [in]  var_name        Name of the pointer
//  * \param [in]  chunk_size      Size given by glibc 
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Source line of the pointer __LINE__
//  * \param [in]  size_cs         Size of the call stack  
//  * \param [in]  cs              Array of string for the call stack 
//  * \param [in]  time            Time at the allocation  
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
// )
// {
// 	char addr_to_str[21] = {0};
// 	char chunk_size_to_str[21] = {0};
// 	unsigned char size_str_addr;
// 	unsigned char size_str_chunk_size;

// 	if (ptr != NULL) {
// 		// pointer's address
// 		sprintf(addr_to_str,"%14ld",(uintptr_t)ptr);
// 		size_str_addr = strlen(addr_to_str);
	
// 		// chunk size
// 		sprintf(chunk_size_to_str,"%14ld",chunk_size);
// 		size_str_chunk_size = strlen(chunk_size_to_str);
// 	} else {
// 		sprintf(addr_to_str,"%s","NULL");
// 		size_str_addr = 4;

// 	       	sprintf(chunk_size_to_str,"%s","0");
// 		size_str_chunk_size = 1;
// 	}

// 	// variable name
// 	unsigned char size_str_varname = strlen(var_name);

// 	// file name
// 	unsigned short int size_str_file_name = strlen(file_name);
	
// 	// line
// 	char line_to_str[21] = {0};
// 	sprintf(line_to_str,"%11d",line);
// 	unsigned char size_str_line = strlen(line_to_str);

// 	// call stack
// 	unsigned short int size_str_call_stack = 9 * (size_cs-2);
	
// 	char *call_stack_to_str = (char*) calloc(size_str_call_stack,sizeof(char));	
// 	PDM_mem_error(call_stack_to_str,__FILE__,__LINE__);

// 	unsigned char nb_fct = 0;
// 	for (int i=size_cs-1;i>0;i--) {
// 		char *function = _get_id_fct(cs[i]);
// 		unsigned char size_str_fct = strlen(function);
		 
// 		// suppose function's address use at most 8 bits (on arch 64-bits) 
// 		// ignore shared libraries
// 		if (size_str_fct <= 8) {			
// 			memcpy(call_stack_to_str + 9 * nb_fct,
// 				function,
// 				size_str_fct
// 				);
// 			if (nb_fct > 0) {
// 				memcpy(call_stack_to_str + 9 * nb_fct -1,":",1);
// 			}
// 			nb_fct++;
// 		}
// 		free(function);
// 	}
// 	size_str_call_stack = strlen(call_stack_to_str);

// 	// time
// 	char time_to_str[30] = {0};
// 	sprintf(time_to_str,"%8.15e",time);
// 	unsigned char size_str_time = strlen(time_to_str);

// 	unsigned int size_line = size_str_addr +
// 				size_str_chunk_size +
// 				size_str_varname +
// 				size_str_file_name +
// 				size_str_line +
// 				size_str_call_stack +
// 				size_str_time +
// 				13;

// 	char *line_buffer = (char*) calloc(size_line,sizeof(char));
// 	PDM_mem_error(line_buffer,__FILE__,__LINE__);

// 	// buffer line assembly
// 	unsigned short int index = 0;
// 	memcpy(line_buffer+index,addr_to_str,size_str_addr);
// 	index += size_str_addr;
// 	memcpy(line_buffer+index,"\tF\t",3);
// 	index += 3;
// 	memcpy(line_buffer+index,chunk_size_to_str,size_str_chunk_size);
// 	index += size_str_chunk_size;
// 	memcpy(line_buffer+index,"\t/\t",3);
// 	index += 3;
// 	memcpy(line_buffer+index,var_name,size_str_varname);
// 	index += size_str_varname;
// 	memcpy(line_buffer+index,"\t",1);
// 	index +=1;
// 	memcpy(line_buffer+index,file_name,size_str_file_name);
// 	index += size_str_file_name;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,line_to_str,size_str_line);
// 	index += size_str_line;
// 	memcpy(line_buffer+index,"\t",1);
// 	index +=1;
// 	memcpy(line_buffer+index,call_stack_to_str,size_str_call_stack);
// 	index += size_str_call_stack;
// 	memcpy(line_buffer+index,"\t",1);
// 	index += 1;
// 	memcpy(line_buffer+index,time_to_str,size_str_time);
// 	index += size_str_time;
// 	memcpy(line_buffer+index,"\n\0",2);

// 	if (call_stack_to_str != NULL) {
// 		free(call_stack_to_str);
// 	}

// 	return line_buffer;
// }


// /**
//  *
//  * \brief Write line to the buffer
//  *
//  * \param [in]  line_buffer     Line for the buffer
//  *
//  */


// #ifndef _SUMMARY_
// void _write_to_buffer(char *line_buffer) {
// 	if (pdm_mem_tool.fd_local_proc < 0) {
// 		printf("erreur fichier non ouvert PDM_MEM_START oublié ?\n");
// 		exit(1);
// 	}

// 	unsigned short int size_line = strlen(line_buffer);
// 	if (pdm_mem_tool.index_buffer + size_line >= BUFFER_SIZE) {
// 		write(pdm_mem_tool.fd_local_proc,&(pdm_mem_tool.buffer),pdm_mem_tool.index_buffer);
// 		pdm_mem_tool.index_buffer = 0;
// 		memset(pdm_mem_tool.buffer,0,BUFFER_SIZE);
// 	}
// 	memcpy(pdm_mem_tool.buffer+pdm_mem_tool.index_buffer,line_buffer,size_line);
// 	pdm_mem_tool.index_buffer += size_line;
// }
// #endif



/*=============================================================================
 * Public function prototypes
 *============================================================================*/

// /**
//  *
//  * \brief Inititialize parameter of the pdm_mem_tool
//  *
//  * \param [in]  local_proc_file Name of the output memory file for 1 MPI rank
//  * \param [in]  rank            MPI rank
//  * \param [in]  size_world      MPI world size
//  * \param [in]  comm            MPI communicator 
//  *
//  */


// void PDM_MEM_INIT(
// #ifdef _MEM_ANALYZER_
//  	const char *file_name 
//  	#ifndef NO_MPI
//  		, int rank, 
//    		int size_world,
//    		PDM_MPI_Comm comm
// 	#endif
// #else
//    	void
// #endif
// ) {
// #ifdef _MEM_ANALYZER_

// 	#ifndef _SUMMARY_
// 		#ifndef NO_MPI
// 		char rank_to_str[10];
// 		sprintf(rank_to_str,"%d",rank);
	
// 		int file_name_size = strlen(file_name)+strlen(rank_to_str)+1;
	
// 		char *mem_file = (char *) malloc(file_name_size * sizeof(char));
// 		PDM_mem_error(mem_file,__FILE__,__LINE__);

// 		strcpy(mem_file,file_name);
// 		strcat(mem_file,rank_to_str);
// 		pdm_mem_tool.local_proc_file = mem_file;
// 		#else
// 		pdm_mem_tool.local_proc_file = file_name;
// 		#endif
			
// 	memset(pdm_mem_tool.buffer,0,BUFFER_SIZE);
// 	pdm_mem_tool.index_buffer = 0;
// 	#endif

// 	pdm_mem_tool.is_running = -1;

// 	#ifndef NO_MPI
// 	pdm_mem_tool.is_mpi_finished = 0;
// 	pdm_mem_tool.rank = rank;
// 	pdm_mem_tool.size_mpi = size_world;
// 	pdm_mem_tool.comm = comm;
// 	pdm_mem_tool.real_time_mem = 0;

// 	pdm_mem_tool.mempeak.nb_alloc = 0;
// 	pdm_mem_tool.mempeak.nb_free = 0;
// 	pdm_mem_tool.mempeak.total_alloc = 0;
// 	pdm_mem_tool.mempeak.total_free = 0;
// 	pdm_mem_tool.mempeak.total_time = 0.;


// 	// MPI Datatype PDM_mempeak_t creation
// 	int lg_bloc[9] = {1,1,180,1,1,1,1,1,1};
// 	MPI_Datatype types[9] = {MPI_LONG_LONG_INT,
// 					MPI_DOUBLE,
// 					MPI_CHAR,
// 					MPI_INT,
// 					MPI_LONG_INT,
// 					MPI_LONG_INT,
// 					MPI_LONG_LONG_INT,
// 					MPI_LONG_LONG_INT,
// 					MPI_DOUBLE};

// 	PDM_mempeak_t peak_t; // dummy peak
        
// 	MPI_Aint addr[10],offset[10];

//         MPI_Get_address(&(peak_t),              &(addr[0]));
//         MPI_Get_address(&(peak_t.peak),         &(addr[1]));
//         MPI_Get_address(&(peak_t.time_peak),         &(addr[2]));
//         MPI_Get_address(&(peak_t.cs),           &(addr[3]));
//         MPI_Get_address(&(peak_t.size_cs),      &(addr[4]));
//  	MPI_Get_address(&(peak_t.nb_alloc),    	&(addr[5]));
//         MPI_Get_address(&(peak_t.nb_free),     	&(addr[6]));
//  	MPI_Get_address(&(peak_t.total_alloc), 	&(addr[7]));
//         MPI_Get_address(&(peak_t.total_free),  	&(addr[8]));
//         MPI_Get_address(&(peak_t.total_time),  	&(addr[9]));



//         for (int i=0;i<9;i++) {
//                 offset[i] = MPI_Aint_diff(addr[i+1],addr[0]);
//         }

//         MPI_Type_create_struct(9,lg_bloc,offset,types,&MPI_mempeak_t);

//         MPI_Type_commit(&MPI_mempeak_t);
// 	#endif
// #endif
// }

// /**
//  *
//  * \brief Start timer and recording
//  *
//  */



// void PDM_MEM_START(void) {
// #ifdef _MEM_ANALYZER_
// 	#ifndef _SUMMARY_
// 	if (pdm_mem_tool.fd_local_proc > 0) {
// 		printf("erreur fichier déjà ouvert\n");
// 		exit(1);
// 	}

// 	pdm_mem_tool.fd_local_proc = open(pdm_mem_tool.local_proc_file,O_CREAT |  O_RDWR ,0666);
	
// 	if (pdm_mem_tool.fd_local_proc < 0) {
// 		printf("erreur ouverture fichier\n");
// 		exit(1);
// 	}
// 	#endif

// 	pdm_mem_tool.is_running = 1;

	
// 	timer_prog = PDM_timer_create();
// 	PDM_timer_resume(timer_prog);
// #endif
// } 


// /**
//  *
//  * \brief Gather rank peak data before MPI communication end
//  *
//  */


// void PDM_MEM_FINALIZE(void) {
// #ifdef _MEM_ANALYZER_
// 	#ifndef NO_MPI
// 	if (!pdm_mem_tool.is_mpi_finished) {
// 		pdm_mem_tool.is_mpi_finished = 1;

// 		// write call_stack in mem_peak struct
//                 for (int i=0;i<pdm_mem_tool.mempeak.size_cs;i++) {
//                         char *fct = _get_id_fct(pdm_mem_tool.cs_peak[i]);
//                         char len_fct = strlen(fct);
//                         memcpy(pdm_mem_tool.mempeak.cs[i],fct,len_fct);
//                         pdm_mem_tool.mempeak.cs[i][len_fct] = '\0';
//                         free(fct);
//                 }
// 		free(pdm_mem_tool.cs_peak);


// 		// total time
// 	        if (timer_prog != NULL) {
//         	        PDM_timer_hang_on(timer_prog);
//                 	pdm_mem_tool.mempeak.total_time = PDM_timer_elapsed(timer_prog);
// 			PDM_timer_resume(timer_prog);
// 		} else {
//                 	pdm_mem_tool.mempeak.total_time = 0;
//         	}

// 		// receive data variable
// 		PDM_mempeak_t *recvpeaks = NULL;

// 		if (pdm_mem_tool.rank == 0) {
// 			recvpeaks = (PDM_mempeak_t *) malloc(pdm_mem_tool.size_mpi * sizeof(PDM_mempeak_t));
// 			PDM_mem_error(recvpeaks,__FILE__,__LINE__);
// 		}
	
// 		// gather proc's data
// 		MPI_Gather(&pdm_mem_tool.mempeak,1,MPI_mempeak_t,recvpeaks,1,MPI_mempeak_t,0,*((MPI_Comm *)PDM_MPI_2_mpi_comm(pdm_mem_tool.comm)));

// 		// write to file
// 		if (pdm_mem_tool.rank == 0) {
// 			char buffer[4096] = {0};
// 			int index_buffer = 0;
// 			int fd_summary = open("Mem_summary.txt",O_CREAT |  O_RDWR ,0666);
	
// 			if (fd_summary < 0) {
// 				printf("erreur ouverture fichier\n");
// 				exit(1);
// 			}

// 			for (int i=0;i<pdm_mem_tool.size_mpi;i++) {
// 				// rank
// 				char rank_to_str[12] = {0};
// 				sprintf(rank_to_str,"%5d",i);
// 				char size_str_rank = strlen(rank_to_str);

// 				// peak
// 				char peak_to_str[21] = {0};
// 				sprintf(peak_to_str,"%20lld",recvpeaks[i].peak);
// 				unsigned char size_str_peak = strlen(peak_to_str);
			
// 				// time_peak
// 				char time_peak_to_str[30] = {0};
// 				sprintf(time_peak_to_str,"%8.15e",recvpeaks[i].time_peak);
// 				char size_str_time_peak = strlen(time_peak_to_str);
			
// 				// call stack
// 				unsigned short int size_str_call_stack = 9 * (recvpeaks[i].size_cs-2);
// 				char *call_stack_to_str = (char*) calloc(size_str_call_stack,sizeof(char));
// 				PDM_mem_error(call_stack_to_str,__FILE__,__LINE__);
			
// 				unsigned char nb_fct = 0;
// 				for (int j=recvpeaks[i].size_cs-1;j>0;j--) {
// 					unsigned char size_str_fct = strlen(recvpeaks[i].cs[j]);
					
// 					// suppose function's address use at most 8 bits (on arch 64-bits) 
// 					// ignore shared libraries
// 					if (size_str_fct <= 8) { 
// 						memcpy(call_stack_to_str + 9 * nb_fct,
// 							recvpeaks[i].cs[j],
// 							size_str_fct
// 							);
// 						if (nb_fct > 0) {
// 							memcpy(call_stack_to_str + 9 * nb_fct -1,":",1);
// 						}
// 						nb_fct++;
// 					}
// 				}

// 				size_str_call_stack = strlen(call_stack_to_str);

// 				// nb alloc
// 				char nb_alloc_to_str[18] = {0};
// 				sprintf(nb_alloc_to_str,"%14ld",recvpeaks[i].nb_alloc);
// 				char size_str_nb_alloc = strlen(nb_alloc_to_str);
				
// 				// nb free
// 				char nb_free_to_str[18] = {0};
// 				sprintf(nb_free_to_str,"%14ld",recvpeaks[i].nb_free);
// 				char size_str_nb_free = strlen(nb_free_to_str);
				
			
// 				// total alloc
// 				char total_alloc_to_str[21] = {0};
// 				sprintf(total_alloc_to_str,"%20lld",recvpeaks[i].total_alloc);
// 				char size_str_total_alloc = strlen(total_alloc_to_str);
				

// 				// total free
// 				char total_free_to_str[21] = {0};
// 				sprintf(total_free_to_str,"%20lld",recvpeaks[i].total_free);
// 				char size_str_total_free = strlen(total_free_to_str);
				
// 				// total time
// 				char total_time_to_str[30] = {0};
//                                 sprintf(total_time_to_str,"%8.15e",recvpeaks[i].total_time);
//                                 char size_str_total_time = strlen(total_time_to_str);

// 				unsigned int size_line = size_str_rank +
// 							size_str_peak +
// 							size_str_time_peak +
// 							size_str_call_stack +
// 							size_str_nb_alloc + 
// 							size_str_nb_free +
// 							size_str_total_alloc +
// 							size_str_total_free +
// 							size_str_total_time +
// 							15;
			
// 				char *line_buffer = (char *) calloc(size_line,sizeof(char));
// 				PDM_mem_error(line_buffer,__FILE__,__LINE__);
		
// 				// buffer line assembly
// 				unsigned short int index = 0;
// 				memcpy(line_buffer+index,rank_to_str,size_str_rank);
// 				index += size_str_rank;
// 				memcpy(line_buffer+index,"\t",1);
// 				index += 1;
// 				memcpy(line_buffer+index,peak_to_str,size_str_peak);
// 				index += size_str_peak;
// 				memcpy(line_buffer+index,"\t",1);
//                        		index += 1;
// 				memcpy(line_buffer+index,time_peak_to_str,size_str_time_peak);
// 				index += size_str_time_peak;
// 				memcpy(line_buffer+index,"\t",1);
//                        		index += 1;
// 				memcpy(line_buffer+index,call_stack_to_str,size_str_call_stack);
// 				index += size_str_call_stack;
// 				memcpy(line_buffer+index,"\t",1);
//                        		index += 1;
// 				memcpy(line_buffer+index,nb_alloc_to_str,size_str_nb_alloc);
// 				index += size_str_nb_alloc;
// 				memcpy(line_buffer+index,"\t",1);
//                        		index += 1;
// 				memcpy(line_buffer+index,nb_free_to_str,size_str_nb_free);
// 				index += size_str_nb_free;
// 				memcpy(line_buffer+index,"\t",1);
//                        		index += 1;
// 				memcpy(line_buffer+index,total_alloc_to_str,size_str_total_alloc);
// 				index += size_str_total_alloc;
// 				memcpy(line_buffer+index,"\t",1);
//                        		index += 1;
// 				memcpy(line_buffer+index,total_free_to_str,size_str_total_free);
// 				index += size_str_total_free;
// 				memcpy(line_buffer+index,"\t",1);
//                                 index += 1;
//                                 memcpy(line_buffer+index,total_time_to_str,size_str_total_time);
//                                 index += size_str_total_time;
// 				memcpy(line_buffer+index,"\n\0",2);
	
// 				free(call_stack_to_str);

// 				// write to buffer
// 				size_line = strlen(line_buffer);
// 				if (index_buffer + size_line >= 4096) {
// 					write(fd_summary,&buffer,index_buffer);
// 					index_buffer = 0;
// 					memset(buffer,0,4096);
// 				}
// 				memcpy(buffer+index_buffer,line_buffer,size_line);
// 				index_buffer += size_line;

// 				free(line_buffer);
// 			}
	
// 			write(fd_summary,&buffer,index_buffer);	
// 			free(recvpeaks);
// 			close(fd_summary);
// 		}

// 		MPI_Type_free(&MPI_mempeak_t);

// 	}
// 	#endif
// #endif
// }

// /**
//  *
//  * \brief End recording and write buffer to file if buffer is not empty
//  *
//  */


// void PDM_MEM_END(void) {
// #ifdef _MEM_ANALYZER_
// 	#ifndef _SUMMARY_
// 	if (pdm_mem_tool.index_buffer > 0) {
// 		if (pdm_mem_tool.fd_local_proc < 0) {
//         		pdm_mem_tool.fd_local_proc = open(pdm_mem_tool.local_proc_file,O_APPEND |  O_RDWR ,0666);
// 		} 
// 		write(pdm_mem_tool.fd_local_proc,&(pdm_mem_tool.buffer),pdm_mem_tool.index_buffer);
//                 pdm_mem_tool.index_buffer = 0;
//                 memset(pdm_mem_tool.buffer,0,BUFFER_SIZE);
// 		close(pdm_mem_tool.fd_local_proc);
// 	}
// 	#endif
	
// 	double elapse_time;
	
// 	if (timer_prog != NULL) {
// 		PDM_timer_hang_on(timer_prog);
// 		elapse_time = PDM_timer_elapsed(timer_prog);
// 		PDM_timer_free(timer_prog);
// 	} else {
// 		elapse_time = 0;
// 	}
// 	pdm_mem_tool.is_running = 0;
// 	printf("fin mem analyze, temps : %.5f sec\n",elapse_time);
// #endif
// }

// /**
//  *
//  * \brief Stop recording, timer continue
//  *
//  */


// void PDM_MEM_STOP(void) {
// #ifdef _MEM_ANALYZER_
// 	pdm_mem_tool.is_running = 0;
// #endif
// }

// /**
//  *
//  * \brief Resume recording
//  *
//  */


// void PDM_MEM_RESUME(void) {
// #ifdef _MEM_ANALYZER_
// 	pdm_mem_tool.is_running = 1;
// #endif
// }


// /**
//  *
//  * \brief Allocate a pointer while recording
//  *
//  * \param [in]  nb              Number of elements
//  * \param [in]  size_type       Size type of the elements
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Source line __LINE__ 
//  * \param [out] var_name        Pointer name
//  *
//  */


// void *PDM_MEM_MALLOC
// (
//  size_t nb,
//  size_t size_type,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// )
// {

// 	unsigned long int alloc_used = nb * size_type;
// 	void *ptr = malloc(alloc_used);
// 	PDM_mem_error(ptr,file_name,line);

// 	if (pdm_mem_tool.is_running) {
		
// 		PDM_timer_hang_on(timer_prog);

// 		void *array[MAX_CALL_TRACE];
// 		char **functions;
// 		int size_cs;

// 		size_cs = backtrace(array,MAX_CALL_TRACE);
// 		functions = backtrace_symbols(array,size_cs);

// 		double time = PDM_timer_elapsed(timer_prog);
		
// 		#ifndef SANITIZE
// 		unsigned long int chunk_size = (unsigned long int) (*((size_t*)ptr-1)) & ~7ULL;
// 		#else
// 		unsigned long int chunk_size = malloc_usable_size(ptr);
// 		#endif

// 		#ifndef _SUMMARY_
// 		char *line_buffer = _buffer_line_malloc(ptr,
//  					var_name,
// 	 				alloc_used,
// 					chunk_size,
// 					file_name,
//  					line,
//  					size_cs,
// 					functions,
//  					time,
// 					"M");

// 		_write_to_buffer(line_buffer);
// 		free(line_buffer);
// 		#endif	             	
		
// 		#ifndef NO_MPI
// 		pdm_mem_tool.real_time_mem += chunk_size;
//        		pdm_mem_tool.mempeak.nb_alloc += 1;
// 	       	pdm_mem_tool.mempeak.total_alloc += chunk_size;

// 		if (pdm_mem_tool.real_time_mem > pdm_mem_tool.mempeak.peak) {
//                 	pdm_mem_tool.mempeak.peak = pdm_mem_tool.real_time_mem;
//                 	pdm_mem_tool.mempeak.time_peak = time;
//                 	pdm_mem_tool.mempeak.size_cs = size_cs;
// 			free(pdm_mem_tool.cs_peak);
// 			pdm_mem_tool.cs_peak = functions;
// 		} else {
// 			free(functions);
// 		}
// 		#else
// 		free(functions);
// 		#endif
	
// 		PDM_timer_resume(timer_prog);
// 	}

// 	return ptr;
// }


// /**
//  *
//  * \brief Reallocate a pointer while recording
//  *
//  * \param [in]  old_ptr         Old pointer to reallocate
//  * \param [in]  new_ptr         New pointer to be reallocated
//  * \param [in]  nb              Number of elements
//  * \param [in]  size_type       Size type of the elements
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Source line __LINE__
//  * \param [out] var_name        Pointer name
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
// )
// {

// 	if (pdm_mem_tool.is_running){

// 		void *array[MAX_CALL_TRACE];
// 		char **functions;
// 		int size_cs;
	
// 		PDM_timer_hang_on(timer_prog);	

// 		double time = PDM_timer_elapsed(timer_prog);

// 		size_cs = backtrace(array,MAX_CALL_TRACE);
// 		functions = backtrace_symbols(array,size_cs);
		
// 		unsigned long int chunk_size;
			
// 		if (old_ptr != NULL) {
// 			#ifndef SANITIZE
// 			chunk_size = (unsigned long int) (*((size_t*)old_ptr-1)) & ~7ULL;
// 			#else
// 			chunk_size = malloc_usable_size(old_ptr);
// 			#endif

// 			#ifndef NO_MPI	
// 			pdm_mem_tool.real_time_mem -= chunk_size;
// 			pdm_mem_tool.mempeak.total_free += chunk_size;
// 			pdm_mem_tool.mempeak.nb_free += 1;
// 			#endif
// 		} else {
// 			chunk_size = 0;
// 		}
			
// 		#ifndef _SUMMARY_
// 		char *line_buffer_free = _buffer_line_free(old_ptr,
// 						var_name,
// 						chunk_size,
// 						file_name,
// 						line,
// 						size_cs,
// 						functions,
// 						time);

// 		_write_to_buffer(line_buffer_free);
// 		free(line_buffer_free);
// 		#endif

// 		PDM_timer_resume(timer_prog);

// 		unsigned long int alloc_used = nb * size_type;

// 		if (alloc_used != 0) {
// 			new_ptr = realloc(old_ptr,alloc_used);
// 			PDM_mem_error(new_ptr,file_name,line);
// 		}

// 		PDM_timer_hang_on(timer_prog);
	
// 		if (new_ptr != NULL) {
// 			#ifndef SANITIZE
// 			chunk_size = (unsigned long int) (*((size_t*)new_ptr-1)) & ~7ULL;
// 			#else
// 			chunk_size = malloc_usable_size(new_ptr);
// 			#endif
// 		} else {
// 			chunk_size = 0;
// 		}

// 		time = PDM_timer_elapsed(timer_prog);

// 		#ifndef _SUMMARY_
// 		char *line_buffer_malloc = _buffer_line_malloc(new_ptr,
// 						var_name,
// 						alloc_used,
// 						chunk_size,
// 						file_name,
// 						line,
// 						size_cs,
// 						functions,
// 						time,
// 						"R");
	
// 		_write_to_buffer(line_buffer_malloc);
// 		free(line_buffer_malloc);
// 		#endif

// 		#ifndef NO_MPI
// 		pdm_mem_tool.real_time_mem += chunk_size;
// 		pdm_mem_tool.mempeak.nb_alloc += 1;
// 		pdm_mem_tool.mempeak.total_alloc += chunk_size;
		
//                 if (pdm_mem_tool.real_time_mem > pdm_mem_tool.mempeak.peak) {
//                         pdm_mem_tool.mempeak.peak = pdm_mem_tool.real_time_mem;
//                         pdm_mem_tool.mempeak.time_peak = time;
//                         pdm_mem_tool.mempeak.size_cs = size_cs;
//                         free(pdm_mem_tool.cs_peak);
//                         pdm_mem_tool.cs_peak = functions;
//                 } else {
//                         free(functions);
// 		}
// 		#else
// 		free(functions);
// 		#endif

//  		PDM_timer_resume(timer_prog);
	
// 		return new_ptr;

// 	} else {
// 		unsigned long int alloc_used = nb * size_type;

// 		if (alloc_used != 0) {
// 			new_ptr = realloc(old_ptr,alloc_used);
// 			PDM_mem_error(new_ptr,file_name,line);
// 		}

// 		return new_ptr;
// 	}
// }


// /**
//  *
//  * \brief Allocate a pointer while recording
//  *
//  * \param [in]  nb              Number of elements
//  * \param [in]  size_type       Size type of the elements
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Source line __LINE__
//  * \param [out] var_name        Pointer name
//  *
//  */


// void *PDM_MEM_CALLOC
// (
//  size_t nb,
//  size_t size_type,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// )
// {

// 	void *ptr = calloc(nb,size_type);
// 	PDM_mem_error(ptr,file_name,line);

// 	if (pdm_mem_tool.is_running) {
		
// 		PDM_timer_hang_on(timer_prog);

// 		void *array[MAX_CALL_TRACE];
// 		char **functions;
// 		int size_cs;
	

// 		size_cs = backtrace(array,MAX_CALL_TRACE);
// 		functions = backtrace_symbols(array,size_cs);

// 		double time = PDM_timer_elapsed(timer_prog);
// 		#ifndef SANITIZE
// 		unsigned long int chunk_size = (unsigned long int) (*((size_t*)ptr-1)) & ~7ULL;
// 		#else
// 		unsigned long int chunk_size = malloc_usable_size(ptr);
// 		#endif
		
// 		#ifndef _SUMMARY_
// 		unsigned long int alloc_used = nb * size_type;
// 		char *line_buffer = _buffer_line_malloc(ptr,
//  					var_name,
// 	 				alloc_used,
// 					chunk_size,
// 					file_name,
//  					line,
//  					size_cs,
// 					functions,
//  					time,
// 					"C");
		
// 		_write_to_buffer(line_buffer);
// 		free(line_buffer);	
// 		#endif	

// 		#ifndef NO_MPI
// 		pdm_mem_tool.real_time_mem += chunk_size;
//         	pdm_mem_tool.mempeak.nb_alloc += 1;
// 		pdm_mem_tool.mempeak.total_alloc += chunk_size;

// 		if (pdm_mem_tool.real_time_mem > pdm_mem_tool.mempeak.peak) {
//                 	pdm_mem_tool.mempeak.peak = pdm_mem_tool.real_time_mem;
//                 	pdm_mem_tool.mempeak.time_peak = time;
//                 	pdm_mem_tool.mempeak.size_cs = size_cs;
// 			free(pdm_mem_tool.cs_peak);
// 			pdm_mem_tool.cs_peak = functions;
// 		} else {
// 			free(functions);
// 		}	
// 		#else 
// 		free(functions);
// 		#endif

// 		PDM_timer_resume(timer_prog);
// 	}

// 	return ptr;
// }

// /**
//  *
//  * \brief Free a pointer while recording
//  *
//  * \param [in]  ptr             Pointer to be freed
//  * \param [in]  file_name       Source file of the pointer __FILE__
//  * \param [in]  line            Source line __LINE__ 
//  * \param [out] var_name        Pointer name
//  *
//  */


// void PDM_MEM_FREE
// (
//  void *ptr,
//  const char *file_name,
//  unsigned int line,
//  const char *var_name
// )
// {
//         if (pdm_mem_tool.is_running) {

// 		unsigned long int chunk_size;

// 		if (ptr != NULL) {
// 			#ifndef SANITIZE
// 			chunk_size = (unsigned long int) (*((size_t*)ptr-1)) & ~7ULL;
// 			#else
// 		       	chunk_size = malloc_usable_size(ptr);	
// 			#endif

// 			#ifndef NO_MPI
// 			pdm_mem_tool.real_time_mem -= chunk_size;
// 			#endif
// 		} else {
// 			chunk_size = 0;
// 		}


// 		#ifndef _SUMMARY_
// 		PDM_timer_hang_on(timer_prog);

// 		double time = PDM_timer_elapsed(timer_prog);

//                 void *array[MAX_CALL_TRACE];
//                 char **functions;
//                 int size_cs;

//                 size_cs = backtrace(array,MAX_CALL_TRACE);
//                 functions = backtrace_symbols(array,size_cs);

			
// 		char *line_buffer = _buffer_line_free(ptr,
// 						var_name,
// 						chunk_size,
// 						file_name,
// 						line,
// 						size_cs,
// 						functions,
// 						time);

// 		free(functions);

// 		_write_to_buffer(line_buffer);
// 		free(line_buffer);
                
// 		PDM_timer_resume(timer_prog);
// 		#endif	

// 		#ifndef NO_MPI
// 		pdm_mem_tool.mempeak.nb_free += 1;
// 		pdm_mem_tool.mempeak.total_free += chunk_size;
// 		#endif
// 	}

// 	free(ptr);
// }


/**
 *
 * \brief Allocate a pointer whithout recording
 *
 * \param [in]  nb              Number of elements
 * \param [in]  size_type       Size type of the elements
 * \param [in]  func_name       Source function __func__
 * \param [in]  file_name       Source file of the pointer __FILE__
 * \param [in]  line            Source line __LINE__ 
 * \param [out] var_name        Pointer name
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
)
{
	size_t s_alloc = nb * size_type;
	void *ptr = malloc(s_alloc);

	if (ptr == NULL && s_alloc != 0) {
		printf("Error - PDM_malloc : Null pointer\n");
		printf("%s %s : %d into %s\n",var_name,file_name,line,func_name);
		abort();
	}
	
	return ptr;
}

/**
 *
 * \brief Reallocate a pointer whithout recording
 *
 * \param [in]  old_ptr         Old pointer to reallocate
 * \param [in]  new_ptr         New pointer to be reallocated
 * \param [in]  nb              Number of elements
 * \param [in]  size_type       Size type of the elements
 * \param [in]  func_name       Source function __func__
 * \param [in]  file_name       Source file of the pointer __FILE__
 * \param [in]  line            Source line __LINE__
 * \param [out] var_name        Pointer name
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
)
{

  void *new_ptr = NULL;
	size_t s_alloc = nb * size_type;

  if (s_alloc == 0) {
    if (old_ptr == NULL) {
      new_ptr = NULL;
    }
    else  {
  	  free (old_ptr);
		  new_ptr = malloc (s_alloc);
    }
	}

  else {
		new_ptr = realloc(old_ptr, nb * size_type);

		if (new_ptr == NULL) {
			printf("Error - PDM_realloc : Null pointer\n");
			printf("%s %s : %d into %s\n",var_name,file_name,line,func_name);
			abort();
		}
	}

	return new_ptr;
}


/**
 *
 * \brief Allocate a pointer whithout recording
 *
 * \param [in]  nb              Number of elements
 * \param [in]  size_type       Size type of the elements
 * \param [in]  func_name       Source function __func__
 * \param [in]  file_name       Source file of the pointer __FILE__
 * \param [in]  line            Source line __LINE__
 * \param [out] var_name        Pointer name
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
)
{
	void *ptr = calloc(nb,size_type);
	if (ptr == NULL && (nb * size_type) != 0) {
		printf("Error - PDM_calloc : Null pointer\n");
		printf("%s %s : %d into %s\n",var_name,file_name,line,func_name);
		abort();
	}
	
	return ptr;
}


/**
 *
 * \brief Free a pointer whithout recording
 *
 * \param [in]  ptr             Pointer to be freed
 *
 */


void PDM_FREE
(
 void **ptr
)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
}

