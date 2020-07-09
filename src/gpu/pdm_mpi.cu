// /*============================================================================
//  * Encapsulation de MPI
//  *============================================================================*/

// /*----------------------------------------------------------------------------
//  * Standard C library headers
//  *----------------------------------------------------------------------------*/


// #include <stdlib.h>
// #include <stdio.h>
// #include <string.h>
// #include <limits.h>

// /*----------------------------------------------------------------------------
//  *  Header for the current file
//  *----------------------------------------------------------------------------*/

// #include "pdm_mpi.h"
// #include "pdm_cuda_error.cuh"
// #include "pdm_mpi.cuh"

// #include <mpi.h>

// #ifdef __cplusplus
// extern "C" {
// #if 0
// } /* Fake brace to force back Emacs auto-indentation back to column 0 */
// #endif
// #endif /* __cplusplus */


// /*============================================================================
//  * Definition des types
//  *============================================================================*/


// /*============================================================================
//  * Definition des variables globales
//  *============================================================================*/

// /*----------------------------------------------------------------------------
//  * Indirection sur le code d'erreur
//  *----------------------------------------------------------------------------*/

// __device__ static const int mpi_err[] = {

//   MPI_SUCCESS,
//   MPI_ERR_BUFFER,
//   MPI_ERR_COUNT,
//   MPI_ERR_TYPE,
//   MPI_ERR_TAG,
//   MPI_ERR_COMM,
//   MPI_ERR_RANK,
//   MPI_ERR_ROOT,
//   MPI_ERR_TRUNCATE,
//   MPI_ERR_GROUP,
//   MPI_ERR_OP,
//   MPI_ERR_REQUEST,
//   MPI_ERR_TOPOLOGY,
//   MPI_ERR_DIMS,
//   MPI_ERR_ARG,
//   MPI_ERR_UNKNOWN,
//   MPI_ERR_OTHER,
//   MPI_ERR_INTERN,
//   MPI_ERR_IN_STATUS,
//   MPI_ERR_PENDING,
//   MPI_MAX_ERROR_STRING,



//   MPI_ERR_ACCESS,
//   MPI_ERR_AMODE,
//   MPI_ERR_BAD_FILE,
//   MPI_ERR_CONVERSION,
//   MPI_ERR_DUP_DATAREP,
//   MPI_ERR_FILE_EXISTS,
//   MPI_ERR_FILE_IN_USE,
//   MPI_ERR_FILE,
//   MPI_ERR_INFO_KEY,
//   MPI_ERR_INFO_NOKEY,
//   MPI_ERR_INFO_VALUE,
//   MPI_ERR_IO,
//   MPI_ERR_NO_MEM,
//   MPI_ERR_NOT_SAME,
//   MPI_ERR_NO_SPACE,
//   MPI_ERR_NO_SUCH_FILE,
//   MPI_ERR_QUOTA,
//   MPI_ERR_READ_ONLY,
//   MPI_ERR_UNSUPPORTED_DATAREP,
//   MPI_ERR_UNSUPPORTED_OPERATION,
//   MPI_ERR_WIN,
//   MPI_ERR_LASTCODE,
//   MPI_ERR_ASSERT,
//   MPI_ERR_BASE,
//   MPI_ERR_DISP,
//   MPI_ERR_KEYVAL,
//   MPI_ERR_LOCKTYPE,
//   MPI_ERR_RMA_CONFLICT,
//   MPI_ERR_RMA_SYNC,
//   MPI_ERR_SIZE

// };

// /*----------------------------------------------------------------------------
//  * Indirection sur le mode du fichier
//  *----------------------------------------------------------------------------*/

// __device__ static const int mpi_file_mode[] = {

// MPI_MODE_CREATE,
// MPI_MODE_RDONLY,
// MPI_MODE_WRONLY,
// MPI_MODE_RDWR,
// MPI_MODE_DELETE_ON_CLOSE,
// MPI_MODE_UNIQUE_OPEN,
// MPI_MODE_EXCL,
// MPI_MODE_APPEND,
// MPI_MODE_SEQUENTIAL,
// MPI_DISPLACEMENT_CURRENT,
// MPI_SEEK_SET,
// MPI_SEEK_CUR,
// MPI_SEEK_END,
// MPI_MODE_WRONLY | MPI_MODE_APPEND,
// MPI_MODE_WRONLY | MPI_MODE_CREATE

// };

// /*----------------------------------------------------------------------------
//  * Indirection PDM_MPI_Datatype -> MPI_Datatype
//  *----------------------------------------------------------------------------*/

// __device__ static const MPI_Datatype mpi_datatype_cste[] = {

//   MPI_BYTE,
//   MPI_PACKED,
//   MPI_CHAR,
//   MPI_SHORT,
//   MPI_INT,
//   MPI_LONG,
//   MPI_FLOAT,
//   MPI_DOUBLE,
//   MPI_LONG_DOUBLE,
//   MPI_UNSIGNED_CHAR,
//   MPI_UNSIGNED_SHORT,
//   MPI_UNSIGNED_LONG,
//   MPI_UNSIGNED,
//   MPI_FLOAT_INT,
//   MPI_DOUBLE_INT,
//   MPI_LONG_DOUBLE_INT,
//   MPI_LONG_INT,
//   MPI_SHORT_INT,
//   MPI_2INT,
//   MPI_CHARACTER,
//   MPI_INTEGER,
//   MPI_REAL,
//   MPI_DOUBLE_PRECISION,
//   MPI_DATATYPE_NULL,
//   MPI_INT8_T,
//   MPI_INT16_T,
//   MPI_INT32_T,
//   MPI_INT64_T,
//   MPI_UINT8_T,
//   MPI_UINT16_T,
//   MPI_UINT32_T,
//   MPI_UINT64_T
// };

// /*----------------------------------------------------------------------------
//  * Indirection PDM_MPI_Op -> MPI_Op
//  *----------------------------------------------------------------------------*/

// __device__ static const MPI_Op mpi_op[] = {

//   MPI_MAX,
//   MPI_MIN,
//   MPI_SUM,
//   MPI_OP_NULL

// };



// /*----------------------------------------------------------------------------
//  * Indirection constantes PDM_MPI_File ->constantes MPI_File
//  *----------------------------------------------------------------------------*/

// __device__ static const MPI_File mpi_file_cste[] = {

//   MPI_FILE_NULL

// };

// /*----------------------------------------------------------------------------
//  * Indirection constantes PDM_MPI_Comm ->constantes MPI_Comm
//  *----------------------------------------------------------------------------*/

// __device__ static const MPI_Comm mpi_comm_cste[] = {

//   MPI_COMM_NULL,
//   MPI_COMM_WORLD

// };

// /*----------------------------------------------------------------------------
//  * Indirection constantes PDM_MPI_Request ->constantes MPI_Request
//  *----------------------------------------------------------------------------*/

// __device__ static const MPI_Request mpi_request_cste[] = {

//   MPI_REQUEST_NULL,

// };

// /*----------------------------------------------------------------------------
//  * Indirection PDM_MPI_File -> MPI_File
//  * stockage dans un tableau
//  *----------------------------------------------------------------------------*/

// __device__ static MPI_File **mpi_file   = NULL; /* Tableau de stockage */
// __device__ static int       l_mpi_file = 0;     /* Taille du tableau */
// __device__ static int       n_mpi_file = 0;     /* Nombre de fichiers stockes */

// /*----------------------------------------------------------------------------
//  * Indirection PDM_MPI_Comm -> MPI_Comm
//  * stockage dans un tableau
//  *----------------------------------------------------------------------------*/

// __device__ static MPI_Comm **mpi_comm   = NULL; /* Tableau de stockage */
// __device__ static int       l_mpi_comm = 0;     /* Taille du tableau */
// __device__ static int       n_mpi_comm = 0;     /* Nombre de communicateurs stockes */

// /*----------------------------------------------------------------------------
//  * Indirection PDM_MPI_Request -> MPI_Request
//  * stockage dans un tableau
//  *----------------------------------------------------------------------------*/

// __device__ static MPI_Request **mpi_request = NULL; /* Tableau de stockage */
// __device__ static int       l_mpi_request = 0;   /* Taille du tableau */
// __device__ static int       n_mpi_request = 0;   /* Nombre de communicateurs stockes */

// /*----------------------------------------------------------------------------
//  * Indirection PDM_MPI_Datatype -> MPI_Datatype
//  * stockage dans un tableau des types utilisateurs
//  *----------------------------------------------------------------------------*/

// __device__ static MPI_Datatype **mpi_datatype   = NULL; /* Tableau de stockage */
// __device__ static int           l_mpi_datatype = 0;     /* Taille du tableau */
// __device__ static int           n_mpi_datatype = 0;     /* Nombre de communicateurs stockes */


// /*============================================================================
//  * Defintion des fonctions privees
//  *============================================================================*/


// /*----------------------------------------------------------------------------
//  * mpi_err -> pdm_mpi_err
//  *
//  * PDM_MPI_Comm -> MPI_Comm
//  *----------------------------------------------------------------------------*/

// __device__ static int _mpi_2_pdm_mpi_err(int code_mpi)
// {
//   int code = PDM_MPI_ERR_OTHER;
//   switch(code_mpi) {
//   case MPI_SUCCESS:
//     code = PDM_MPI_SUCCESS;
//   break;
//   case MPI_ERR_BUFFER:
//     code = PDM_MPI_ERR_BUFFER;
//       break;
//   case MPI_ERR_COUNT:
//     code = PDM_MPI_ERR_COUNT;
//       break;
//   case MPI_ERR_TYPE:
//     code = PDM_MPI_ERR_TYPE;
//     break;
//   case MPI_ERR_TAG:
//     code = PDM_MPI_ERR_TAG;
//     break;
//   case MPI_ERR_COMM:
//     code = PDM_MPI_ERR_COMM;
//     break;
//   case MPI_ERR_RANK:
//     code = PDM_MPI_ERR_RANK;
//     break;
//   case MPI_ERR_ROOT:
//     code = PDM_MPI_ERR_ROOT;
//     break;
//   case MPI_ERR_TRUNCATE:
//     code = PDM_MPI_ERR_TRUNCATE;
//     break;
//   case MPI_ERR_GROUP:
//     code = PDM_MPI_ERR_GROUP;
//     break;
//   case MPI_ERR_OP:
//     code = PDM_MPI_ERR_OP;
//     break;
//   case MPI_ERR_REQUEST:
//     code = PDM_MPI_ERR_REQUEST;
//     break;
//   case MPI_ERR_TOPOLOGY:
//     code = PDM_MPI_ERR_TOPOLOGY;
//     break;
//   case MPI_ERR_DIMS:
//     code = PDM_MPI_ERR_DIMS;
//     break;
//   case MPI_ERR_ARG:
//     code = PDM_MPI_ERR_ARG;
//     break;
//   case MPI_ERR_UNKNOWN:
//     code = PDM_MPI_ERR_UNKNOWN;
//     break;
//   case MPI_ERR_OTHER:
//     code = PDM_MPI_ERR_OTHER;
//     break;
//   case MPI_ERR_INTERN:
//     code = PDM_MPI_ERR_INTERN;
//     break;
//   case MPI_ERR_IN_STATUS:
//     code = PDM_MPI_ERR_IN_STATUS;
//     break;
//   case MPI_ERR_PENDING:
//     code = PDM_MPI_ERR_PENDING;
//     break;



//   case MPI_ERR_ACCESS:
//     code = PDM_MPI_ERR_ACCESS;
//     break;
//   case MPI_ERR_AMODE:
//     code = PDM_MPI_ERR_AMODE;
//     break;
//   case MPI_ERR_BAD_FILE:
//     code = PDM_MPI_ERR_BAD_FILE;
//     break;
//   case MPI_ERR_CONVERSION:
//     code = PDM_MPI_ERR_CONVERSION;
//     break;
//   case MPI_ERR_DUP_DATAREP:
//     code = PDM_MPI_ERR_DUP_DATAREP;
//     break;
//   case MPI_ERR_FILE_EXISTS:
//     code = PDM_MPI_ERR_FILE_EXISTS;
//     break;
//   case MPI_ERR_FILE_IN_USE:
//     code = PDM_MPI_ERR_FILE_IN_USE;
//     break;
//   case MPI_ERR_FILE:
//     code = PDM_MPI_ERR_FILE;
//     break;
//   case MPI_ERR_INFO_KEY:
//     code = PDM_MPI_ERR_INFO_KEY;
//     break;
//   case MPI_ERR_INFO_NOKEY:
//     code = PDM_MPI_ERR_INFO_NOKEY;
//     break;
//   case MPI_ERR_INFO_VALUE:
//     code = PDM_MPI_ERR_INFO_VALUE;
//     break;
//   case MPI_ERR_IO:
//     code = PDM_MPI_ERR_IO;
//     break;
//   case MPI_ERR_NO_MEM:
//     code = PDM_MPI_ERR_NO_MEM;
//     break;
//   case MPI_ERR_NOT_SAME:
//     code = PDM_MPI_ERR_NOT_SAME;
//     break;
//   case MPI_ERR_NO_SPACE:
//     code = PDM_MPI_ERR_NO_SPACE;
//     break;
//   case MPI_ERR_NO_SUCH_FILE:
//     code = PDM_MPI_ERR_NO_SUCH_FILE;
//     break;
//   case MPI_ERR_QUOTA:
//     code = PDM_MPI_ERR_QUOTA;
//     break;
//   case MPI_ERR_READ_ONLY:
//     code = PDM_MPI_ERR_READ_ONLY;
//     break;
//   case MPI_ERR_UNSUPPORTED_DATAREP:
//     code = PDM_MPI_ERR_UNSUPPORTED_DATAREP;
//     break;
//   case MPI_ERR_UNSUPPORTED_OPERATION:
//     code = PDM_MPI_ERR_UNSUPPORTED_OPERATION;
//     break;
//   case MPI_ERR_WIN:
//     code = PDM_MPI_ERR_WIN;
//     break;
//   case MPI_ERR_LASTCODE:
//     code = PDM_MPI_ERR_LASTCODE;
//     break;



//   case MPI_ERR_ASSERT:
//     code = PDM_MPI_ERR_ASSERT;
//     break;
//   case MPI_ERR_BASE:
//     code = PDM_MPI_ERR_BASE;
//     break;
//   case MPI_ERR_DISP:
//     code = PDM_MPI_ERR_DISP;
//     break;
//   case MPI_ERR_KEYVAL:
//     code = PDM_MPI_ERR_KEYVAL;
//     break;
//   case MPI_ERR_LOCKTYPE:
//     code = PDM_MPI_ERR_LOCKTYPE;
//     break;
//   case MPI_ERR_RMA_CONFLICT:
//     code = PDM_MPI_ERR_RMA_CONFLICT;
//     break;
//   case MPI_ERR_RMA_SYNC:
//     code = PDM_MPI_ERR_RMA_SYNC;
//     break;
//   case MPI_ERR_SIZE:
//     code = PDM_MPI_ERR_SIZE;

//   }
//   return code;
// }

// /*----------------------------------------------------------------------------
//  * _pdm_mpi_2_mpi_comm
//  *
//  * PDM_MPI_Comm -> MPI_Comm
//  *----------------------------------------------------------------------------*/

// __device__ static MPI_Comm _pdm_mpi_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
// {

//   /* Traitement des communicateurs predefinis */

//   if (pdm_mpi_comm < 0)
//     return mpi_comm_cste[-pdm_mpi_comm - 1];

//   /* Traitement des communicateurs utilisateurs */

//   else {
//     if (pdm_mpi_comm < l_mpi_comm)
//       return *(mpi_comm[pdm_mpi_comm]);
//     else {
//       PDM_error_GPU(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_comm :"
//             " pdm_mpi_comm '%d' non valide\n", pdm_mpi_comm);
//       __trap();
//     }
//   }
// }



// /*============================================================================
//  * Defintion des fonctions publiques
//  *============================================================================*/
// /*----------------------------------------------------------------------------
//  * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
//  *
//  *----------------------------------------------------------------------------*/

// __device__ int PDM_MPI_Comm_rank_GPU(PDM_MPI_Comm comm, int *rank)
// {
//   int code = MPI_Comm_rank(_pdm_mpi_2_mpi_comm(comm), rank);
//   return _mpi_2_pdm_mpi_err(code);
// }

// /*----------------------------------------------------------------------------
//  * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
//  *
//  *----------------------------------------------------------------------------*/

// __device__ int PDM_MPI_Comm_size_GPU(PDM_MPI_Comm comm, int *size)
// {
//   int code = MPI_Comm_size(_pdm_mpi_2_mpi_comm(comm), size);
//   return _mpi_2_pdm_mpi_err(code);
// }









// #ifdef __cplusplus
// }
// #endif /* __cplusplus */