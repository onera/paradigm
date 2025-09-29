/*
 * \file
 */

#ifndef __PDM_MPI_H__
#define __PDM_MPI_H__

/*============================================================================
 * Bibliotheque de messagerie
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/
#include "pdm.h"
#include <mpi.h>

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

typedef MPI_Request  PDM_MPI_Request;
typedef MPI_Win      PDM_MPI_Win;
typedef MPI_Comm     PDM_MPI_Comm;
typedef MPI_Datatype PDM_MPI_Datatype;
typedef MPI_File     PDM_MPI_File;
typedef MPI_Group    PDM_MPI_Group;
typedef MPI_Op       PDM_MPI_Op;

typedef MPI_Offset PDM_MPI_Offset;

typedef MPI_Aint PDM_MPI_Aint;
typedef MPI_Fint PDM_MPI_Fint;

#define PDM_MPI_IN_PLACE  MPI_IN_PLACE

#define PDM_MPI_MODE_NOPRECEDE MPI_MODE_NOPRECEDE
#define PDM_MPI_UNDEFINED      MPI_UNDEFINED

enum {
  PDM_MPI_COMM_UNDEFINED, /* dist graph topology */
  PDM_MPI_DIST_GRAPH,     /* dist graph topology */
  PDM_MPI_CART,           /* cartesian topology  */
  PDM_MPI_GRAPH           /* graph topology      */
};

#define PDM_MPI_SUCCESS                   MPI_SUCCESS
#define PDM_MPI_ERR_BUFFER                MPI_ERR_BUFFER
#define PDM_MPI_ERR_COUNT                 MPI_ERR_COUNT
#define PDM_MPI_ERR_TYPE                  MPI_ERR_TYPE
#define PDM_MPI_ERR_TAG                   MPI_ERR_TAG
#define PDM_MPI_ERR_COMM                  MPI_ERR_COMM
#define PDM_MPI_ERR_RANK                  MPI_ERR_RANK
#define PDM_MPI_ERR_ROOT                  MPI_ERR_ROOT
#define PDM_MPI_ERR_TRUNCATE              MPI_ERR_TRUNCATE
#define PDM_MPI_ERR_GROUP                 MPI_ERR_GROUP
#define PDM_MPI_ERR_OP                    MPI_ERR_OP
#define PDM_MPI_ERR_REQUEST               MPI_ERR_REQUEST
#define PDM_MPI_ERR_TOPOLOGY              MPI_ERR_TOPOLOGY
#define PDM_MPI_ERR_DIMS                  MPI_ERR_DIMS
#define PDM_MPI_ERR_ARG                   MPI_ERR_ARG
#define PDM_MPI_ERR_UNKNOWN               MPI_ERR_UNKNOWN
#define PDM_MPI_ERR_OTHER                 MPI_ERR_OTHER
#define PDM_MPI_ERR_INTERN                MPI_ERR_INTERN
#define PDM_MPI_ERR_IN_STATUS             MPI_ERR_IN_STATUS
#define PDM_MPI_ERR_PENDING               MPI_ERR_PENDING
#define PDM_MPI_MAX_ERROR_STRING          MPI_MAX_ERROR_STRING
#define PDM_MPI_ERR_ACCESS                MPI_ERR_ACCESS
#define PDM_MPI_ERR_AMODE                 MPI_ERR_AMODE
#define PDM_MPI_ERR_BAD_FILE              MPI_ERR_BAD_FILE
#define PDM_MPI_ERR_CONVERSION            MPI_ERR_CONVERSION
#define PDM_MPI_ERR_DUP_DATAREP           MPI_ERR_DUP_DATAREP
#define PDM_MPI_ERR_FILE_EXISTS           MPI_ERR_FILE_EXISTS
#define PDM_MPI_ERR_FILE_IN_USE           MPI_ERR_FILE_IN_USE
#define PDM_MPI_ERR_FILE                  MPI_ERR_FILE
#define PDM_MPI_ERR_INFO_KEY              MPI_ERR_INFO_KEY
#define PDM_MPI_ERR_INFO_NOKEY            MPI_ERR_INFO_NOKEY
#define PDM_MPI_ERR_INFO_VALUE            MPI_ERR_INFO_VALUE
#define PDM_MPI_ERR_IO                    MPI_ERR_IO
#define PDM_MPI_ERR_NO_MEM                MPI_ERR_NO_MEM
#define PDM_MPI_ERR_NOT_SAME              MPI_ERR_NOT_SAME
#define PDM_MPI_ERR_NO_SPACE              MPI_ERR_NO_SPACE
#define PDM_MPI_ERR_NO_SUCH_FILE          MPI_ERR_NO_SUCH_FILE
#define PDM_MPI_ERR_QUOTA                 MPI_ERR_QUOTA
#define PDM_MPI_ERR_READ_ONLY             MPI_ERR_READ_ONLY
#define PDM_MPI_ERR_UNSUPPORTED_DATAREP   MPI_ERR_UNSUPPORTED_DATAREP
#define PDM_MPI_ERR_UNSUPPORTED_OPERATION MPI_ERR_UNSUPPORTED_OPERATION
#define PDM_MPI_ERR_WIN                   MPI_ERR_WIN
#define PDM_MPI_ERR_LASTCODE              MPI_ERR_LASTCODE
#define PDM_MPI_ERR_ASSERT                MPI_ERR_ASSERT
#define PDM_MPI_ERR_BASE                  MPI_ERR_BASE
#define PDM_MPI_ERR_DISP                  MPI_ERR_DISP
#define PDM_MPI_ERR_KEYVAL                MPI_ERR_KEYVAL
#define PDM_MPI_ERR_LOCKTYPE              MPI_ERR_LOCKTYPE
#define PDM_MPI_ERR_RMA_CONFLICT          MPI_ERR_RMA_CONFLICT
#define PDM_MPI_ERR_RMA_SYNC              MPI_ERR_RMA_SYNC
#define PDM_MPI_ERR_SIZE                  MPI_ERR_SIZE

// enum  {

//   PDM_MPI_SUCCESS,
//   PDM_MPI_ERR_BUFFER,
//   PDM_MPI_ERR_COUNT,
//   PDM_MPI_ERR_TYPE,
//   PDM_MPI_ERR_TAG,
//   PDM_MPI_ERR_COMM,
//   PDM_MPI_ERR_RANK,
//   PDM_MPI_ERR_ROOT,
//   PDM_MPI_ERR_TRUNCATE,
//   PDM_MPI_ERR_GROUP,
//   PDM_MPI_ERR_OP,
//   PDM_MPI_ERR_REQUEST,
//   PDM_MPI_ERR_TOPOLOGY,
//   PDM_MPI_ERR_DIMS,
//   PDM_MPI_ERR_ARG,
//   PDM_MPI_ERR_UNKNOWN,
//   PDM_MPI_ERR_OTHER,
//   PDM_MPI_ERR_INTERN,
//   PDM_MPI_ERR_IN_STATUS,
//   PDM_MPI_ERR_PENDING,
//   PDM_MPI_MAX_ERROR_STRING,
//   PDM_MPI_ERR_ACCESS,
//   PDM_MPI_ERR_AMODE,
//   PDM_MPI_ERR_BAD_FILE,
//   PDM_MPI_ERR_CONVERSION,
//   PDM_MPI_ERR_DUP_DATAREP,
//   PDM_MPI_ERR_FILE_EXISTS,
//   PDM_MPI_ERR_FILE_IN_USE,
//   PDM_MPI_ERR_FILE,
//   PDM_MPI_ERR_INFO_KEY,
//   PDM_MPI_ERR_INFO_NOKEY,
//   PDM_MPI_ERR_INFO_VALUE,
//   PDM_MPI_ERR_IO,
//   PDM_MPI_ERR_NO_MEM,
//   PDM_MPI_ERR_NOT_SAME,
//   PDM_MPI_ERR_NO_SPACE,
//   PDM_MPI_ERR_NO_SUCH_FILE,
//   PDM_MPI_ERR_QUOTA,
//   PDM_MPI_ERR_READ_ONLY,
//   PDM_MPI_ERR_UNSUPPORTED_DATAREP,
//   PDM_MPI_ERR_UNSUPPORTED_OPERATION,
//   PDM_MPI_ERR_WIN,
//   PDM_MPI_ERR_LASTCODE,
//   PDM_MPI_ERR_ASSERT,
//   PDM_MPI_ERR_BASE,
//   PDM_MPI_ERR_DISP,
//   PDM_MPI_ERR_KEYVAL,
//   PDM_MPI_ERR_LOCKTYPE,
//   PDM_MPI_ERR_RMA_CONFLICT,
//   PDM_MPI_ERR_RMA_SYNC,
//   PDM_MPI_ERR_SIZE
// };

#define PDM_MPI_MODE_CREATE          MPI_MODE_CREATE
#define PDM_MPI_MODE_RDONLY          MPI_MODE_RDONLY
#define PDM_MPI_MODE_WRONLY          MPI_MODE_WRONLY
#define PDM_MPI_MODE_RDWR            MPI_MODE_RDWR
#define PDM_MPI_MODE_DELETE_ON_CLOSE MPI_MODE_DELETE_ON_CLOSE
#define PDM_MPI_MODE_UNIQUE_OPEN     MPI_MODE_UNIQUE_OPEN
#define PDM_MPI_MODE_EXCL            MPI_MODE_EXCL
#define PDM_MPI_MODE_APPEND          MPI_MODE_APPEND
#define PDM_MPI_MODE_SEQUENTIAL      MPI_MODE_SEQUENTIAL
#define PDM_MPI_DISPLACEMENT_CURRENT MPI_DISPLACEMENT_CURRENT
#define PDM_MPI_SEEK_SET             MPI_SEEK_SET
#define PDM_MPI_SEEK_CUR             MPI_SEEK_CUR
#define PDM_MPI_SEEK_END             MPI_SEEK_END
#define PDM_MPI_MODE_WRONLY_APPEND   MPI_MODE_WRONLY | PDM_MPI_MODE_APPEND
#define PDM_MPI_MODE_WRONLY_CREATE   MPI_MODE_WRONLY | PDM_MPI_MODE_CREATE

// enum {

//   PDM_MPI_MODE_CREATE,
//   PDM_MPI_MODE_RDONLY,
//   PDM_MPI_MODE_WRONLY,
//   PDM_MPI_MODE_RDWR ,
//   PDM_MPI_MODE_DELETE_ON_CLOSE,
//   PDM_MPI_MODE_UNIQUE_OPEN,
//   PDM_MPI_MODE_EXCL,
//   PDM_MPI_MODE_APPEND,
//   PDM_MPI_MODE_SEQUENTIAL,
//   PDM_MPI_DISPLACEMENT_CURRENT,
//   PDM_MPI_SEEK_SET,
//   PDM_MPI_SEEK_CUR,
//   PDM_MPI_SEEK_END,
//   PDM_MPI_MODE_WRONLY_APPEND,
//   PDM_MPI_MODE_WRONLY_CREATE

// };

#define PDM_MPI_MAX MPI_MAX
#define PDM_MPI_MIN MPI_MIN
#define PDM_MPI_SUM MPI_SUM
#define PDM_MPI_MINLOC MPI_MINLOC
#define PDM_MPI_MAXLOC MPI_MAXLOC
#define PDM_MPI_OP_NULL MPI_OP_NULL

// typedef enum {
//   PDM_MPI_MAX,
//   PDM_MPI_MIN,
//   PDM_MPI_SUM,
//   PDM_MPI_MINLOC,
//   PDM_MPI_MAXLOC,
//   PDM_MPI_OP_NULL
// } PDM_MPI_Op;

#define PDM_MPI_BYTE               MPI_BYTE
#define PDM_MPI_PACKED             MPI_PACKED
#define PDM_MPI_CHAR               MPI_CHAR
#define PDM_MPI_SHORT              MPI_SHORT
#define PDM_MPI_INT                MPI_INT
#define PDM_MPI_LONG               MPI_LONG
#define PDM_MPI_FLOAT              MPI_FLOAT
#define PDM_MPI_DOUBLE             MPI_DOUBLE
#define PDM_MPI_LONG_DOUBLE        MPI_LONG_DOUBLE
#define PDM_MPI_UNSIGNED_CHAR      MPI_UNSIGNED_CHAR
#define PDM_MPI_UNSIGNED_SHORT     MPI_UNSIGNED_SHORT
#define PDM_MPI_UNSIGNED_LONG      MPI_UNSIGNED_LONG
#define PDM_MPI_UNSIGNED           MPI_UNSIGNED
#define PDM_MPI_FLOAT_INT          MPI_FLOAT_INT
#define PDM_MPI_DOUBLE_INT         MPI_DOUBLE_INT
#define PDM_MPI_LONG_DOUBLE_INT    MPI_LONG_DOUBLE_INT
#define PDM_MPI_LONG_INT           MPI_LONG_INT
#define PDM_MPI_SHORT_INT          MPI_SHORT_INT
#define PDM_MPI_2INT               MPI_2INT
#define PDM_MPI_CHARACTER          MPI_CHARACTER
#define PDM_MPI_INTEGER            MPI_INTEGER
#define PDM_MPI_REAL               MPI_REAL
#define PDM_MPI_DOUBLE_PRECISION   MPI_DOUBLE_PRECISION
#define PDM_MPI_DATATYPE_NULL      MPI_DATATYPE_NULL
#define PDM_MPI_INT8_T             MPI_INT8_T
#define PDM_MPI_INT16_T            MPI_INT16_T
#define PDM_MPI_INT32_T            MPI_INT32_T
#define PDM_MPI_INT64_T            MPI_INT64_T
#define PDM_MPI_UINT8_T            MPI_UINT8_T
#define PDM_MPI_UINT16_T           MPI_UINT16_T
#define PDM_MPI_UINT32_T           MPI_UINT32_T
#define PDM_MPI_UINT64_T           MPI_UINT64_T
#define PDM_MPI_UNSIGNED_LONG_LONG MPI_UNSIGNED_LONG_LONG

// enum {
//   PDM_MPI_BYTE                   = -1,
//   PDM_MPI_PACKED                 = -2,
//   PDM_MPI_CHAR                   = -3,
//   PDM_MPI_SHORT                  = -4,
//   PDM_MPI_INT                    = -5,
//   PDM_MPI_LONG                   = -6,
//   PDM_MPI_FLOAT                  = -7,
//   PDM_MPI_DOUBLE                 = -8,
//   PDM_MPI_LONG_DOUBLE            = -9,
//   PDM_MPI_UNSIGNED_CHAR          = -10,
//   PDM_MPI_UNSIGNED_SHORT         = -11,
//   PDM_MPI_UNSIGNED_LONG          = -12,
//   PDM_MPI_UNSIGNED               = -13,
//   PDM_MPI_FLOAT_INT              = -14,
//   PDM_MPI_DOUBLE_INT             = -15,
//   PDM_MPI_LONG_DOUBLE_INT        = -16,
//   PDM_MPI_LONG_INT               = -17,
//   PDM_MPI_SHORT_INT              = -18,
//   PDM_MPI_2INT                   = -19,
//   PDM_MPI_CHARACTER              = -20,
//   PDM_MPI_INTEGER                = -21,
//   PDM_MPI_REAL                   = -22,
//   PDM_MPI_DOUBLE_PRECISION       = -23,
//   PDM_MPI_DATATYPE_NULL          = -24,
//   PDM_MPI_INT8_T                 = -25,
//   PDM_MPI_INT16_T                = -26,
//   PDM_MPI_INT32_T                = -27,
//   PDM_MPI_INT64_T                = -28,
//   PDM_MPI_UINT8_T                = -29,
//   PDM_MPI_UINT16_T               = -30,
//   PDM_MPI_UINT32_T               = -31,
//   PDM_MPI_UINT64_T               = -32,
//   PDM_MPI_UNSIGNED_LONG_LONG     = -33,
// };

#define PDM_MPI_COMM_NULL MPI_COMM_NULL
#define PDM_MPI_COMM_WORLD MPI_COMM_WORLD

// enum {
//   PDM_MPI_COMM_NULL  = -1,
//   PDM_MPI_COMM_WORLD = -2
// };

#define PDM_MPI_FILE_NULL MPI_FILE_NULL
// enum {
//   PDM_MPI_FILE_NULL = -1
// };

#define PDM_MPI_REQUEST_NULL MPI_REQUEST_NULL
// enum {
//   PDM_MPI_REQUEST_NULL  = -123456789
// };

#define PDM_MPI_WIN_NULL MPI_WIN_NULL
// enum {
//   PDM_MPI_WIN_NULL  = -1
// };


#define PDM_MPI_GROUP_NULL MPI_GROUP_NULL
// enum {
//   PDM_MPI_GROUP_NULL  = -1
// };

enum {
  PDM_MPI_SPLIT_SHARED  = 1,
  PDM_MPI_SPLIT_NUMA    = 2
};

typedef struct _pdm_mpi_win_shared_t PDM_mpi_win_shared_t;

/*============================================================================
 * Prototype des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *----------------------------------------------------------------------------*/
int PDM_MPI_Init(int *argc, char ***argv);

/*----------------------------------------------------------------------------
 * PDM_MPI_Finalize -> MPI_Finalize
 *----------------------------------------------------------------------------*/
int PDM_MPI_Finalize(void);

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/
MPI_Comm PDM_MPI_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/
PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm(void *);

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/
void *PDM_MPI_free_mpi_comm(void *mpi_comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_open (wrapping de la fonction MPI_File_open)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_open
(
  PDM_MPI_Comm  comm,
  char         *filename,
  int           amode,
  PDM_MPI_File *fh
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_close (wrapping de la fonction MPI_File_close)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_close(PDM_MPI_File *fh);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_seek (wrapping de la fonction MPI_File_seek)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_seek(PDM_MPI_File, PDM_MPI_Offset, int);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_size (wrapping de la fonction MPI_File_get_size)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_get_size(PDM_MPI_File, PDM_MPI_Offset *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_position (wrapping de la fonction MPI_File_get_position)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_get_position(PDM_MPI_File, PDM_MPI_Offset *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_set_view (wrapping de la fonction MPI_File_set_view)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_set_view(PDM_MPI_File, PDM_MPI_Offset, PDM_MPI_Datatype, PDM_MPI_Datatype, const char *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_view (wrapping de la fonction MPI_File_get_view)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_get_view(PDM_MPI_File, PDM_MPI_Offset *, PDM_MPI_Datatype *, PDM_MPI_Datatype *, char *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at (wrapping de la fonction MPI_File_read_at)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_read_at
(
  PDM_MPI_File      fh,
  PDM_MPI_Offset    offset,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_lus
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at_all (wrapping de la fonction MPI_File_read_at_all)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_read_at_all
(
  PDM_MPI_File      fh,
  PDM_MPI_Offset    offset,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_lus
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at (wrapping de la fonction MPI_File_write_at)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_write_at
(
  PDM_MPI_File      fh,
  PDM_MPI_Offset    offset,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_ecrits
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at_all (wrapping de la fonction MPI_File_write_at_all)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_write_at_all
(
  PDM_MPI_File      fh,
  PDM_MPI_Offset    offset,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_ecrits
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read (wrapping de la fonction MPI_File_read)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_read
(
  PDM_MPI_File      fh,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_lus
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_all (wrapping de la fonction MPI_File_read_all)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_read_all
(
  PDM_MPI_File      fh,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_lus
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write (wrapping de la fonction MPI_File_write)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_write
(
  PDM_MPI_File      fh,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_ecrits
);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_all (wrapping de la fonction MPI_File_write_all)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_File_write_all
(
  PDM_MPI_File      fh,
  void             *buf,
  int               count,
  PDM_MPI_Datatype  datatype,
  int              *n_octet_ecrits
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Gather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Igather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                    void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
                    int root, PDM_MPI_Comm comm,
                    PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Gatherv (wrapping de la fonction MPI_Gatherv)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Gatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                    void *recvbuf, int *recvcounts, int *displs,
                    PDM_MPI_Datatype recvtype, int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Recv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
                 int tag, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Irecv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
                  int tag, PDM_MPI_Comm comm, PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Send (wrapping de la fonction MPI_Send)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Send(void *buf, int count, PDM_MPI_Datatype datatype, int dest,
                 int tag, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Isend (wrapping de la fonction MPI_Issend)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Isend(const void *buf, int count, PDM_MPI_Datatype datatype, int dest, int tag,
                  PDM_MPI_Comm comm, PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Send_init (wrapping de la fonction MPI_Send_init)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Send_init
(
  const void             *buf,
        int               count,
        PDM_MPI_Datatype  datatype,
        int               dest,
        int               tag,
        PDM_MPI_Comm      comm,
        PDM_MPI_Request  *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv_init (wrapping de la fonction MPI_Recv_init)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Recv_init
(
        void             *buf,
        int               count,
        PDM_MPI_Datatype  datatype,
        int               dest,
        int               tag,
        PDM_MPI_Comm      comm,
        PDM_MPI_Request  *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Wait (wrapping de la fonction MPI_Wait)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Wait(PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Waitall (wrapping de la fonction MPI_Wait)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Waitall(int count, PDM_MPI_Request array_of_requests[]);

/*----------------------------------------------------------------------------
 * PDM_MPI_Request_free (wrapping de la fonction MPI_Request_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Request_free(PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Test (wrapping de la fonction MPI_Test)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Test(PDM_MPI_Request *request, int *flag);

/*----------------------------------------------------------------------------
 * PDM_MPI_Startall (wrapping de la fonction MPI_Startall)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Startall(int count, PDM_MPI_Request array_of_requests[]);

/*----------------------------------------------------------------------------
 * PDM_MPI_Start (wrapping de la fonction MPI_Start)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Start(PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Type_create_hindexed
(
        int               count,
  const int               array_of_blocklengths[],
  const PDM_MPI_Aint      array_of_displacements[],
        PDM_MPI_Datatype  oldtype,
        PDM_MPI_Datatype *newtype
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_create_contiguous (wrapping de la fonction MPI_Type_create_contiguous)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Type_create_contiguous
(
  int               count,
  PDM_MPI_Datatype  oldtype,
  PDM_MPI_Datatype *newtype
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_commit (wrapping de la fonction MPI_Type_commit)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_commit(PDM_MPI_Datatype *datatype);

/*----------------------------------------------------------------------------
 * MPI_Type_size (wrapping de la fonction MPI_Type_commit)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_size(PDM_MPI_Datatype datatype, int *size);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_free (wrapping de la fonction MPI_Type_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_free(PDM_MPI_Datatype *datatype);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_f2c (wrapping de la fonction MPI_comm_f2c)
 *----------------------------------------------------------------------------*/
PDM_MPI_Comm PDM_MPI_Comm_f2c(PDM_MPI_Fint comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_c2f (wrapping de la fonction MPI_comm_c2f)
 *----------------------------------------------------------------------------*/
PDM_MPI_Fint PDM_MPI_Comm_c2f(PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Scatter (wrapping de la fonction MPI_Scatter)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Scatter
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  int               root,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Barrier (wrapping de la fonction MPI_Barrier)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Barrier(PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Wtime (wrapping de la fonction MPI_Wtime)
 *----------------------------------------------------------------------------*/
double PDM_MPI_Wtime(void);

/*----------------------------------------------------------------------------
 * PDM_MPI_Alloc_mem (wrapping de la fonction MPI_Alloc_mem)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Alloc_mem(PDM_MPI_Aint size, void *baseptr);

/*----------------------------------------------------------------------------
 * PDM_MPI_Free_mem (wrapping de la fonction MPI_Free_mem)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Free_mem(void *baseptr);

/*----------------------------------------------------------------------------
 * PDM_MPI_Bcast (wrapping de la fonction MPI_Bcast)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Bcast
(
  void             *buffer,
  int               count,
  PDM_MPI_Datatype  datatype,
  int               root,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_IBcast (wrapping de la fonction MPI_IBcast)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Ibcast
(
  void              *buffer,
  int                count,
  PDM_MPI_Datatype   datatype,
  int                root,
  PDM_MPI_Comm       comm,
  PDM_MPI_Request   *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Allgather
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgatherv (wrapping de la fonction MPI_Allgatherv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Allgatherv
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *displs,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce (wrapping de la fonction MPI_Reduce)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Reduce
(
  void             *sendbuf,
  void             *recvbuf,
  int               count,
  PDM_MPI_Datatype  datatype,
  PDM_MPI_Op        op,
  int               root,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce_scatter (wrapping de la fonction MPI_Reduce_scatter)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Reduce_scatter
(
  void             *sendbuf,
  void             *recvbuf,
  int              *counts,
  PDM_MPI_Datatype  datatype,
  PDM_MPI_Op        op,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Allreduce (wrapping de la fonction MPI_Allreduce)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Allreduce
(
   void             *sendbuf,
   void             *recvbuf,
   int               count,
   PDM_MPI_Datatype  datatype,
   PDM_MPI_Op        op,
   PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Scan (wrapping de la fonction MPI_Scan)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Scan
(
  const void             *sendbuf,
        void             *recvbuf,
        int               count,
        PDM_MPI_Datatype  datatype,
        PDM_MPI_Op        op,
        PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Exscan (wrapping de la fonction MPI_Exscan)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Exscan
(
  const void             *sendbuf,
        void             *recvbuf,
        int               count,
        PDM_MPI_Datatype  datatype,
        PDM_MPI_Op        op,
        PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Iscan (wrapping de la fonction MPI_Iscan)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Iscan
(
  const void             *sendbuf,
        void             *recvbuf,
        int               count,
        PDM_MPI_Datatype  datatype,
        PDM_MPI_Op        op,
        PDM_MPI_Comm      comm,
        PDM_MPI_Request  *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoall (wrapping de la fonction MPI_Alltoall)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Alltoall
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Alltoall)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Ialltoall
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm,
  PDM_MPI_Request  *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv (wrapping de la fonction MPI_Alltoallv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Alltoallv
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoallv (wrapping de la fonction MPI_Ialltoallv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Ialltoallv
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm,
  PDM_MPI_Request  *request
);


/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv_Init (wrapping de la fonction MPI_Alltoallv_Init)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Alltoallv_init
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm,
  PDM_MPI_Request  *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoallv_init (wrapping de la fonction MPI_Neighbor_alltoallv_init)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Neighbor_alltoallv_init
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm,
  PDM_MPI_Request  *request
);


/*----------------------------------------------------------------------------
 * PDM_MPI_Win_create (wrapping de la fonction MPI_Win_create)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_create(void         *baseptr,
                       PDM_MPI_Aint  size,
                       int           disp_unit,
                       PDM_MPI_Comm  comm,
                       PDM_MPI_Win  *win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_allocate (wrapping de la fonction MPI_Win_allocate)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_allocate(PDM_MPI_Aint  size,
                         int           disp_unit,
                         PDM_MPI_Comm  comm,
                         void         *baseptr,
                         PDM_MPI_Win  *win);
/*----------------------------------------------------------------------------
 * PDM_Win_start (wrapping de la fonction Win_start)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_start(PDM_MPI_Group group, int mpi_assert, PDM_MPI_Win win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_post (wrapping de la fonction MPI_Win_post)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_post(PDM_MPI_Group group, int mpi_assert, PDM_MPI_Win win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_complete (wrapping de la fonction MPI_Win_complete)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_complete(PDM_MPI_Win win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_wait (wrapping de la fonction MPI_Win_wait)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_wait(PDM_MPI_Win win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_free (wrapping de la fonction MPI_Win_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_free(PDM_MPI_Win  *win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_fence (wrapping de la fonction MPI_Win_fence)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_fence(int assert, PDM_MPI_Win win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Group_free (wrapping de la fonction MPI_Group_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Group_free(PDM_MPI_Group *group);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_group (wrapping de la fonction MPI_Comm_group)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_group(PDM_MPI_Comm comm, PDM_MPI_Group *group);

/*----------------------------------------------------------------------------
 * PDM_MPI_Group_incl (wrapping de la fonction MPI_Group_incl)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Group_incl(PDM_MPI_Group group, int n, const int ranks[],
                       PDM_MPI_Group *newgroup);

/*----------------------------------------------------------------------------
 * PDM_MPI_Topo_test (wrapping de la fonction MPI_Topo_test)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Topo_test(PDM_MPI_Comm comm, int *status);

/*----------------------------------------------------------------------------
 * PDM_MPI_Error_string (wrapping de la fonction MPI_Error_string)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Error_string(int errorcode, char *string, int *resultlen);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_rank(PDM_MPI_Comm comm, int *rank);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_size(PDM_MPI_Comm comm, int *size);

/*----------------------------------------------------------------------------
 * PDM_MPI_get_max_error_string
 *----------------------------------------------------------------------------*/
int PDM_MPI_get_max_error_string(void);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_free
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_free(PDM_MPI_Comm *comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_split(PDM_MPI_Comm comm, int color, int key, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_dup
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_dup(PDM_MPI_Comm comm, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_split_type_numa(PDM_MPI_Comm comm, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_split_type(PDM_MPI_Comm comm, int split_type, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_rand_tag_get
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_Comm comm, void *attribute_val, int *flag);
int PDM_MPI_Rand_tag            (PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Dist_graph_create_adjacent
 *----------------------------------------------------------------------------*/
int PDM_MPI_Dist_graph_create_adjacent(PDM_MPI_Comm  comm_old,
                                             int     indegree,
                                       const int     sources[],
                                             int     outdegree,
                                       const int     destinations[],
                                       int           reorder,
                                       PDM_MPI_Comm *comm_dist_graph);


/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Neighbor_allgather
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_allgatherv (wrapping de la fonction MPI_Neighbor_allgatherv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Neighbor_allgatherv
(
 void             *sendbuf,
 int               sendcount,
 PDM_MPI_Datatype  sendtype,
 void             *recvbuf,
 int              *recvcounts,
 int              *displs,
 PDM_MPI_Datatype  recvtype,
 PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoall (wrapping de la fonction MPI_Neighbor_alltoall)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Neighbor_alltoall
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Ineighbor_alltoall
(
  void             *sendbuf,
  int               sendcount,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int               recvcount,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm,
  PDM_MPI_Request  *request
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoallv (wrapping de la fonction MPI_Neighbor_alltoallv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Neighbor_alltoallv
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ineighbor_alltoallv (wrapping de la fonction MPI_Ineighbor_alltoallv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Ineighbor_alltoallv
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm,
  PDM_MPI_Request  *request
);

void
PDM_MPI_setup_hybrid_dist_comm_graph
(
  PDM_MPI_Comm   comm,
  PDM_MPI_Comm  *comm_shared_out,
  PDM_MPI_Comm  *comm_dist_graph_out,
  int           *n_degree,
  int          **neighbor
);

int
PDM_MPI_Dist_graph_neighbors_count
(
  PDM_MPI_Comm  comm,
  int          *n_degree_in,
  int          *n_degree_out,
  int          *is_weighted
);

int
PDM_MPI_Dist_graph_neighbors
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *sources,
  int            n_degree_out,
  int           *destinations
);

void
PDM_MPI_setup_dist_graph_from_neighbor_in
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *neighbor_in,
  PDM_MPI_Comm  *comm_dist_graph_out
);


/*----------------------------------------------------------------------------
 * PDM_mpi_Win_allocate_shared_get
 *----------------------------------------------------------------------------*/
PDM_mpi_win_shared_t*
PDM_mpi_win_shared_create(PDM_MPI_Aint          size,
                          int                   disp_unit,
                          PDM_MPI_Comm          comm);

void* PDM_mpi_win_shared_get(PDM_mpi_win_shared_t *wins);

void PDM_mpi_win_shared_free(PDM_mpi_win_shared_t *wins);

PDM_MPI_Comm PDM_MPI_get_group_of_master(PDM_MPI_Comm comm, PDM_MPI_Comm sub_comm);

int PDM_mpi_win_shared_lock_all(int assert, PDM_mpi_win_shared_t* win);
int PDM_mpi_win_shared_unlock_all(PDM_mpi_win_shared_t* win);
int PDM_mpi_win_shared_sync(PDM_mpi_win_shared_t* win);


/*----------------------------------------------------------------------------
 * MPI Standard extension :
 *   - Add method to ease MPI
 *   - Futur method but not yet supported (ex: MPI_Neigbor_init)
 *----------------------------------------------------------------------------*/

void PDM_mpi_purge(PDM_ownership_t comm_ownership);

/**
 * \brief Initializes multiple non-blocking send operations to a selected set of processes.
 *
 * This function is a wrapper for MPI_Send_init, designed to simplify the
 * initialization of multiple non-blocking sends to a predefined list of
 * destination ranks. It dynamically allocates an array of MPI requests and
 * initializes a separate send request for each active destination.
 *
 * This function does not perform any communication; it only prepares the
 * requests. The communication must be started later using functions like
 * MPI_Start, and completed with functions like MPI_Waitall or MPI_Testall.
 *
 * \param[in]  sendbuf Pointer to the start of the send buffer.
 * \param[in]  sendcounts Array of size n_active_send, where sendcounts[i] is
 *             the number of elements to send to the i-th active
 *             destination.
 * \param[in]  sdispls Array of size n_active_send, where sdispls[i] is the
 *             displacement from sendbuf to the starting element for
 *             the i-th active destination.
 * \param[in]  datatype The type of the data elements in the send buffer.
 * \param[in]  n_active_send The number of active destination processes.
 * \param[in]  active_send Array of size n_active_send, containing the ranks
 *             of the destination processes.
 * \param[in]  tag The message tag to be used for the send operations.
 * \param[in]  comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 * dynamically allocated array of MPI requests.
 *
 * \return An MPI error code. PDM_MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays sendcounts, sdispls, and active_send must be allocated
 * and correctly sized (at least n_active_send).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests,
 * which must be freed by the user after communication is complete.
 *
 * \see PDM_MPI_Recvs_init for the corresponding receive operation.
 * \see MPI_Send_init, MPI_Start, MPI_Waitall
 */
int
PDM_MPI_Sends_init
(
  const void              *sendbuf,
        int               *sendcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_send,
        int               *active_send,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);

/**
 * \brief Initializes and starts multiple non-blocking send operations to selected processes.
 *
 * This function is a convenience wrapper that combines the initialization of
 * multiple non-blocking send requests and their immediate start. It acts as
 * a single call to perform a series of non-blocking sends to a predefined
 * list of destination ranks.
 *
 * The function allocates an array of MPI requests and initiates a separate
 * `MPI_Isend` for each active destination, but unlike `PDM_MPI_Sends_init`,
 * it does not create persistent requests.
 *
 * \param[in] sendbuf Pointer to the start of the send buffer.
 * \param[in] sendcounts Array of size n_active_send, where sendcounts[i] is
 * the number of elements to send to the i-th active destination.
 * \param[in] sdispls Array of size n_active_send, where sdispls[i] is the
 * displacement from sendbuf to the starting element for the i-th active
 * destination.
 * \param[in] datatype The type of the data elements in the send buffer.
 * \param[in] n_active_send The number of active destination processes.
 * \param[in] active_send Array of size n_active_send, containing the ranks
 * of the destination processes.
 * \param[in] tag The message tag to be used for the send operations.
 * \param[in] comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 * dynamically allocated array of non-blocking requests.
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays sendcounts, sdispls, and active_send must be allocated
 * and correctly sized (at least n_active_send).
 * \post A non-null pointer to an array of MPI requests is returned via
 * out_requests, which must be freed by the user after communication is
 * completed (e.g., with MPI_Waitall).
 *
 * \see PDM_MPI_Sends_init for persistent requests.
 * \see PDM_MPI_Irecvs for the corresponding receive operation.
 * \see MPI_Isend, MPI_Waitall, MPI_Testall
 */
int
PDM_MPI_Isends
(
  const void              *sendbuf,
        int               *sendcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_send,
        int               *active_send,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);

/**
 * \brief Initializes and starts multiple non-blocking receive operations from selected processes.
 *
 * This function is a convenience wrapper that combines the initialization of
 * multiple non-blocking receive requests and their immediate start. It acts as
 * a single call to perform a series of non-blocking receives from a predefined
 * list of source ranks.
 *
 * The function allocates an array of MPI requests and initiates a separate
 * `MPI_Irecv` for each active source.
 *
 * \param[out] recvbuf Pointer to the start of the receive buffer.
 * \param[in] recvcounts Array of size n_active_recv, where recvcounts[i] is
 * the number of elements to receive from the i-th active source.
 * \param[in] sdispls Array of size n_active_recv, where sdispls[i] is the
 * displacement from recvbuf to the starting element for the i-th active source.
 * \param[in] datatype The type of the data elements in the receive buffer.
 * \param[in] n_active_recv The number of active source processes.
 * \param[in] active_recv Array of size n_active_recv, containing the ranks
 * of the source processes.
 * \param[in] tag The message tag to be used for the receive operations.
 * \param[in] comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 * dynamically allocated array of non-blocking requests.
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays recvcounts, sdispls, and active_recv must be allocated
 * and correctly sized (at least n_active_recv).
 * \post A non-null pointer to an array of MPI requests is returned via
 * out_requests, which must be freed by the user after communication is
 * completed (e.g., with MPI_Waitall).
 *
 * \see PDM_MPI_Recvs_init for persistent requests.
 * \see PDM_MPI_Isends for the corresponding send operation.
 * \see MPI_Irecv, MPI_Waitall, MPI_Testall
 */
int
PDM_MPI_Irecvs
(
  const void              *recvbuf,
        int               *recvcounts,
        int               *sdispls,
        PDM_MPI_Datatype   datatype,
        int                n_active_recv,
        int               *active_recv,
        int                tag,
        PDM_MPI_Comm       comm,
        PDM_MPI_Request  **out_requests
);


/**
 * \brief Initializes multiple non-blocking receive operations from a selected set of processes.
 *
 * This function is a wrapper for MPI_Recv_init, designed to simplify the
 * initialization of multiple non-blocking receives from a predefined list of
 * source ranks. It dynamically allocates an array of MPI requests and
 * initializes a separate receive request for each active source.
 *
 * This function does not perform any communication; it only prepares the
 * requests. The communication must be started later using functions like
 * MPI_Start, and completed with functions like MPI_Waitall or MPI_Testall.
 *
 * \param[in]  recvbuf Pointer to the start of the receive buffer.
 * \param[in]  recvcounts Array of size n_active_recv, where recvcounts[i] is
 *             the number of elements to receive from the i-th active
 *             source.
 * \param[in]  rdispls Array of size n_active_recv, where rdispls[i] is the
 *             displacement from recvbuf to the starting element for
 *             the i-th active source.
 * \param[in]  datatype The type of the data elements in the receive buffer.
 * \param[in]  n_active_recv The number of active source processes.
 * \param[in]  active_recv Array of size n_active_recv, containing the ranks
 *             of the source processes.
 * \param[in]  tag The message tag to be used for the receive operations.
 * \param[in]  comm The MPI communicator to be used.
 * \param[out] out_requests Pointer to a pointer that will be set to the
 *             dynamically allocated array of MPI requests.
 *
 * \return An MPI error code. MPI_SUCCESS on success, otherwise an error code.
 *
 * \pre The arrays recvcounts, rdispls, and active_recv must be allocated
 *      and correctly sized (at least n_active_recv).
 * \post A non-null pointer to an array of MPI requests is returned via out_requests,
 *       which must be freed by the user after communication is complete.
 *
 * \see PDM_MPI_Sends_init for the corresponding send operation.
 * \see MPI_Recv_init, MPI_Start, MPI_Waitall
 */
int
PDM_MPI_Recvs_init
(
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   datatype,
  int                n_active_recv,
  int               *active_recv,
  int                tag,
  PDM_MPI_Comm       comm,
  PDM_MPI_Request  **out_requests
);

/*----------------------------------------------------------------------------
 *
 * PDM_MPI_Partofactiverank()
 *
 *----------------------------------------------------------------------------*/
void
PDM_MPI_Partofactiverank
(
  int          *sendcounts,
  int          *recvcounts,
  PDM_MPI_Comm  comm,
  double       *part_active_rank
);

/*----------------------------------------------------------------------------
 *
 * PDM_MPI_Ialltoallv_p2p (MPI_Ialltoallv in p2p)
 *
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Ialltoallv_p2p_l
(
  void              *sendbuf,
  int               *sendcounts,
  size_t            *sdispls,
  PDM_MPI_Datatype   sendtype,
  void              *recvbuf,
  int               *recvcounts,
  size_t            *rdispls,
  PDM_MPI_Datatype   recvtype,
  PDM_MPI_Comm       comm,
  PDM_MPI_Request  **request_s,
  PDM_MPI_Request  **request_r,
  int               *n_request_s,
  int               *n_request_r
);

/**
 * @brief Emulates MPI_Ialltoallv using non-blocking P2P messages.
 *
 * This function provides a flexible implementation of a non-blocking all-to-all communication.
 * It dynamically posts non-blocking P2P sends and receives based on non-zero entries in the
 * sendcounts and recvcounts arrays. If explicit rank lists are provided (@a send_rank/@a recv_rank),
 * it delegates to the more specialized PDM_MPI_Ialltoallv_select_p2p function. Otherwise, it
 * assumes a dense communication pattern and iterates through all ranks of the communicator.
 * This function is useful for scenarios where a native collective may not be optimal
 * (e.g., for sparse communication patterns or debugging).
 *
 * \param sendbuf       The send buffer.
 * \param sendcounts    An array of integers specifying the number of elements to send to each rank.
 * \param sdispls       An array of integers specifying the displacement in @a sendbuf for each message.
 * \param sendtype      The datatype of send buffer elements.
 * \param n_send_rank   (Optional) The number of ranks to send to. Used only if @a send_rank is not NULL.
 * \param send_rank     (Optional) An array of integers containing the ranks to send to. If NULL, a dense check on sendcounts is performed.
 * \param recvbuf       The receive buffer.
 * \param recvcounts    An array of integers specifying the number of elements to receive from each rank.
 * \param rdispls       An array of integers specifying the displacement in @a recvbuf for each message.
 * \param recvtype      The datatype of receive buffer elements.
 * \param n_recv_rank   (Optional) The number of ranks to receive from. Used only if @a recv_rank is not NULL.
 * \param recv_rank     (Optional) An array of integers containing the ranks to receive from. If NULL, a dense check on recvcounts is performed.
 * \param tag           The message tag for P2P communication.
 * \param comm          The communicator.
 * \param n_send_recv_request A pointer to an integer that will be filled with the total number of requests.
 * \param out_requests  A pointer to a PDM_MPI_Request array that will be allocated and filled with the requests.
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
int
PDM_MPI_Ialltoallv_p2p
(
  void              *sendbuf,
  int               *sendcounts,
  int               *sdispls,
  PDM_MPI_Datatype   sendtype,
  int                n_send_rank,
  int               *send_rank,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  int                n_recv_rank,
  int               *recv_rank,
  int                tag,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
);


/**
 * @brief Performs a non-blocking all-to-all communication for a selected, sparse set of ranks using P2P messages.
 *
 * This function emulates MPI_Ialltoallv for a sparse communication pattern by posting non-blocking
 * P2P sends (MPI_Issend) and receives (MPI_Irecv) to a user-defined list of ranks. It is particularly
 * useful when the communication partners are known explicitly and represent a small subset of the
 * communicator's ranks. The function returns an array of requests that must be waited on
 * (e.g., with PDM_MPI_Waitall) to ensure the completion of the communication.
 *
 * \param sendbuf       The send buffer.
 * \param sendcounts    An array of integers specifying the number of elements to send to each rank in @a send_rank.
 * \param sdispls       An array of integers specifying the displacement in @a sendbuf for each message.
 * \param sendtype      The datatype of send buffer elements.
 * \param n_send_rank   The number of ranks to send to.
 * \param send_rank     An array of integers containing the ranks to send to.
 * \param recvbuf       The receive buffer.
 * \param recvcounts    An array of integers specifying the number of elements to receive from each rank in @a recv_rank.
 * \param rdispls       An array of integers specifying the displacement in @a recvbuf for each message.
 * \param recvtype      The datatype of receive buffer elements.
 * \param n_recv_rank   The number of ranks to receive from.
 * \param recv_rank     An array of integers containing the ranks to receive from.
 * \param tag           The message tag for P2P communication.
 * \param comm          The communicator.
 * \param n_send_recv_request A pointer to an integer that will be filled with the total number of requests.
 * \param out_requests  A pointer to a PDM_MPI_Request array that will be allocated and filled with the requests.
 * \return PDM_SUCCESS or an error code from the underlying MPI calls.
 */
int
PDM_MPI_Ialltoallv_select_p2p
(
  void              *sendbuf,
  int               *sendcounts,
  int               *sdispls,
  PDM_MPI_Datatype   sendtype,
  int                n_send_rank,
  int               *send_rank,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  int                n_recv_rank,
  int               *recv_rank,
  int                tag,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
);


int
PDM_MPI_Ialltoallv_p2p_rma
(
  PDM_MPI_Win        send_win,
  int               *target_disp,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **out_requests
);


int
PDM_MPI_Alltoallv_l
(
  void             *sendbuf,
  int              *sendcounts,
  size_t           *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  size_t           *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

int
PDM_MPI_Alltoallv_p2p
(
  void             *sendbuf,
  int              *sendcounts,
  int              *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  int              *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);

int
PDM_MPI_Alltoallv_p2p_l
(
  void             *sendbuf,
  int              *sendcounts,
  size_t           *sdispls,
  PDM_MPI_Datatype  sendtype,
  void             *recvbuf,
  int              *recvcounts,
  size_t           *rdispls,
  PDM_MPI_Datatype  recvtype,
  PDM_MPI_Comm      comm
);


int
PDM_MPI_Alltoallv_p2p_init
(
  void              *sendbuf,
  int               *sendcounts,
  int               *sdispls,
  PDM_MPI_Datatype   sendtype,
  void              *recvbuf,
  int               *recvcounts,
  int               *rdispls,
  PDM_MPI_Datatype   recvtype,
  int                tag,
  PDM_MPI_Comm       comm,
  int               *n_send_recv_request,
  PDM_MPI_Request  **requests
);

/*----------------------------------------------------------------------------
 * PDM_MPI_Get_ialltoallv (Implemtation of alltoall like with window )
 *----------------------------------------------------------------------------*/
int PDM_MPI_Get_ialltoallv(PDM_MPI_Win       win_send,
                           PDM_MPI_Win       win_recv,
                           void             *sendbuf,
                           int              *sendcounts,
                           int              *sdispls,
                           PDM_MPI_Datatype  sendtype,
                           void             *recvbuf,
                           int              *recvcounts,
                           int              *rdispls,
                           PDM_MPI_Datatype  recvtype,
                           PDM_MPI_Comm      comm);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MPI_H__ */
