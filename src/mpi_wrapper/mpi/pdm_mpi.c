/*============================================================================
 * Encapsulation de MPI
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#ifdef __linux__
#include <sys/syscall.h> //Non portable mettre un ifdef
#endif
#include <unistd.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_mpi_priv.h"
#include "pdm_printf.h"
#include "pdm_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*============================================================================
 * Definition des variables globales
 *============================================================================*/

/*============================================================================
 * Defintion des fonctions pprivees
 *============================================================================*/

/*============================================================================
 * Defintion des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *----------------------------------------------------------------------------*/
int PDM_MPI_Init(int *argc, char ***argv)
{
  return MPI_Init(argc, argv);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *----------------------------------------------------------------------------*/
int PDM_MPI_Finalize (void)
{
  return MPI_Finalize();
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *----------------------------------------------------------------------------*/
MPI_Comm PDM_MPI_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
{
  return pdm_mpi_comm;
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *----------------------------------------------------------------------------*/
void *PDM_MPI_free_mpi_comm(void *pt_mpi_comm)
{
  MPI_Comm *comm = (MPI_Comm *) pt_mpi_comm;
  MPI_Comm_free (comm);
  return NULL;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_mpi_2_pdm_mpi_comm
 *----------------------------------------------------------------------------*/
PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm(void *pt_mpi_comm)
{
  MPI_Comm _mpi_comm = *((MPI_Comm *) pt_mpi_comm);
  return _mpi_comm;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_open (wrapping de la fonction MPI_File_open)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_open(PDM_MPI_Comm comm, char *filename, int amode, PDM_MPI_File *fh)
{
  char *hints = getenv("PDM_IO_HINTS");

  MPI_Info hints_mpi = MPI_INFO_NULL;

  if (hints != NULL) {

    MPI_Info_create (&hints_mpi);

    char *cp_hints;
    char *name;
    char *value;
    PDM_malloc(cp_hints, (strlen(hints) + 1), char);
    PDM_malloc(name    , (strlen(hints) + 1), char);
    PDM_malloc(value   , (strlen(hints) + 1), char);
    strcpy (cp_hints, hints);

    char *pch;
    char *str2 = cp_hints;

    do {
      pch = strtok (str2,"=");
      str2 = NULL;
      if (pch != NULL) {
        strcpy(name, pch);
        pch = strtok (str2, ":");
        if (pch == NULL) {
          PDM_printf ("Error PDM_MPI_File_open : No value for hint \"%s\"."
                  " Check \"PDM_IO_HINTS\" environment variable\n", name);
          exit(1);
        }
        else {
          strcpy(value, pch);
          MPI_Info_set (hints_mpi, name, value);
          PDM_printf ("MPI/IO hint \"%s\" = \"%s\"\n", name, value);
        }
      }
    } while (pch != NULL);

    PDM_free(cp_hints);
    PDM_free(name);
    PDM_free(value);

  }

  int code = MPI_File_open(comm,
                           filename,
                           amode,
                           hints_mpi,
                           fh);

  if (hints != NULL) {
    MPI_Info_free(&hints_mpi);
  }

  if (code != MPI_SUCCESS) {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_close (wrapping de la fonction MPI_File_close)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_close(PDM_MPI_File *fh)
{
  int code =  MPI_File_close(fh);
  if (code != MPI_SUCCESS) {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;
    MPI_Error_string(code, buffer, &buffer_len);
    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);
    abort();
  }
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_seek (wrapping de la fonction MPI_File_seek)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_seek(PDM_MPI_File fh, PDM_MPI_Offset offset, int whence)
{
  int code = MPI_File_seek(fh,
                           (MPI_Offset) offset,
                           whence);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_size (wrapping de la fonction MPI_File_get_size)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_get_size(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  MPI_Offset _tmp_offset;
  int code = MPI_File_get_size(fh,
                           (MPI_Offset*) &_tmp_offset);
  *offset = _tmp_offset;
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_position (wrapping de la fonction MPI_File_get_position)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_get_position(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  MPI_Offset _tmp_offset;
  int code = MPI_File_get_position(fh,
                                    &_tmp_offset);
  *offset = _tmp_offset;
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_set_view (wrapping de la fonction MPI_File_set_view)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_set_view(PDM_MPI_File fh, PDM_MPI_Offset disp, PDM_MPI_Datatype etype,
                          PDM_MPI_Datatype filetype, const char *datarep)
{
  int code = MPI_File_set_view(fh,
                               (MPI_Offset) disp,
                               etype,
                               filetype,
                               datarep,
                               MPI_INFO_NULL);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_view (wrapping de la fonction MPI_File_get_view)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_get_view(PDM_MPI_File fh, PDM_MPI_Offset *disp,
                          PDM_MPI_Datatype *etype, PDM_MPI_Datatype *filetype, char *datarep)
{
  MPI_Datatype mpi_etype;
  MPI_Datatype mpi_filetype;
  MPI_Offset _disp = (MPI_Offset) *disp;

  int code = MPI_File_get_view(fh,
                               &_disp,
                               &mpi_etype,
                               &mpi_filetype,
                               datarep);

  *etype    = mpi_etype;
  *filetype = mpi_filetype;
  *disp     = (PDM_MPI_Offset) _disp;

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at (wrapping de la fonction MPI_File_read_at)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_read_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                     int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code = MPI_File_read_at(fh,
                              (MPI_Offset) offset,
                              buf,
                              count,
                              datatype,
                              &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at_all (wrapping de la fonction MPI_File_read_at_all)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_read_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code = MPI_File_read_at_all(fh,
                                  (MPI_Offset) offset,
                                  buf,
                                  count,
                                  datatype,
                                  &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at (wrapping de la fonction MPI_File_write_at)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_write_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                      int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  MPI_Offset _offset = (MPI_Offset) offset;
  int code = MPI_File_write_at(fh,
                               _offset,
                               buf,
                               count,
                               datatype,
                               &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at_all (wrapping de la fonction MPI_File_write_at_all)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_write_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  MPI_Offset _offset = (MPI_Offset) offset;
  int code = MPI_File_write_at_all(fh,
                                   (MPI_Offset) _offset,
                                   buf,
                                   count,
                                   datatype,
                                   &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read (wrapping de la fonction MPI_File_read)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_read(PDM_MPI_File fh, void *buf, int count,
                  PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code =  MPI_File_read(fh,
                            buf,
                            count,
                            datatype,
                            &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_all (wrapping de la fonction MPI_File_read_all)
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_all(PDM_MPI_File fh, void *buf, int count,
                      PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  int code = MPI_File_read_all(fh,
                                buf,
                                count,
                                datatype,
                                &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write (wrapping de la fonction MPI_File_write)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_write(PDM_MPI_File fh, void *buf, int count,
                   PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  int code =  MPI_File_write(fh,
                             buf,
                             count,
                             datatype,
                             &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_all (wrapping de la fonction MPI_File_write_all)
 *----------------------------------------------------------------------------*/
int PDM_MPI_File_write_all(PDM_MPI_File fh, void *buf, int count,
                       PDM_MPI_Datatype datatype, int *n_octet_lus)

{
  MPI_Status status;

  int code =  MPI_File_write_all(fh,
                                 buf,
                                 count,
                                 datatype,
                                 &status);

  if (code == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  } else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Gather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm)
{
  int code = MPI_Gather(sendbuf, sendcount, sendtype,
                        recvbuf, recvcount, recvtype,
                        root, comm);
  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Igather (wrapping de la fonction MPI_Igather)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Igather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
                int root, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  int code = MPI_Igather(sendbuf, sendcount, sendtype,
                        recvbuf, recvcount, recvtype,
                        root, comm, request);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gatherv (wrapping de la fonction MPI_Gatherv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Gatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int *recvcounts, int *displs,
                PDM_MPI_Datatype recvtype, int root, PDM_MPI_Comm comm)
{
  int code = MPI_Gatherv(sendbuf,
                         sendcount,
                         sendtype,
                         recvbuf,
                         recvcounts,
                         displs,
                         recvtype,
                         root,
                         comm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Recv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
             int tag, PDM_MPI_Comm comm)
{
  int code =  MPI_Recv(buf, count, datatype, source,
                       tag, comm, MPI_STATUS_IGNORE);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Irecv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
              int tag, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  int code =  MPI_Irecv(buf, count, datatype, source,
                       tag, comm, request);
  assert(code == 0);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Send (wrapping de la fonction MPI_Send)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Send(void *buf, int count, PDM_MPI_Datatype datatype, int dest,
             int tag, PDM_MPI_Comm comm)
{
  int code = MPI_Send(buf, count, datatype, dest,
                      tag, comm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Isend (wrapping de la fonction MPI_Isend)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Isend(const void *buf, int count, PDM_MPI_Datatype datatype, int dest, int tag,
              PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  int code = MPI_Isend(buf, count, datatype, dest,
                       tag, comm, request);
  assert(code == 0);
  return code;
}


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
)
{
  int code = MPI_Send_init(buf,
                           count,
                           datatype,
                           dest,
                           tag,
                           comm,
                           request);

  assert(code == 0);
  return code;
}

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
)
{
  int code = MPI_Recv_init(buf,
                           count,
                           datatype,
                           dest,
                           tag,
                           comm,
                           request);
  assert(code == 0);
  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Send_init (wrapping de la fonction MPI_Send_init)
 *
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
)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code = MPI_Send_init(buf,
                           count,
                           _pdm_mpi_2_mpi_datatype(datatype),
                           dest,
                           tag,
                           _pdm_mpi_2_mpi_comm(comm),
                           &_mpi_request);

  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  assert(code == 0);
  return _mpi_2_pdm_mpi_err(code);
}


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
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_send, PDM_MPI_Request);

  int size_send_type;
  MPI_Type_size(_pdm_mpi_2_mpi_datatype(datatype), &size_send_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_send; i++) {
    void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
    MPI_Request _mpi_request = MPI_REQUEST_NULL;
    int t_rank = active_send[i];
    code = MPI_Send_init(buf,
                         sendcounts[i],
                         _pdm_mpi_2_mpi_datatype(datatype),
                         t_rank,
                         tag,
                         _pdm_mpi_2_mpi_comm(comm),
                         &_mpi_request);
    requests[i] = _mpi_2_pdm_mpi_request_add(_mpi_request);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;

  return _mpi_2_pdm_mpi_err(code);
}


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
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_send, PDM_MPI_Request);

  int size_send_type;
  MPI_Type_size(_pdm_mpi_2_mpi_datatype(datatype), &size_send_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_send; i++) {
    void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_send_type);
    MPI_Request _mpi_request = MPI_REQUEST_NULL;
    int t_rank = active_send[i];
    code = MPI_Isend(buf,
                     sendcounts[i],
                     _pdm_mpi_2_mpi_datatype(datatype),
                     t_rank,
                     tag,
                     _pdm_mpi_2_mpi_comm(comm),
                     &_mpi_request);
    requests[i] = _mpi_2_pdm_mpi_request_add(_mpi_request);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv_init (wrapping de la fonction MPI_Recv_init)
 *
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
)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code = MPI_Recv_init(buf,
                           count,
                           _pdm_mpi_2_mpi_datatype(datatype),
                           dest,
                           tag,
                           _pdm_mpi_2_mpi_comm(comm),
                           &_mpi_request);

  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  assert(code == 0);
  return _mpi_2_pdm_mpi_err(code);
}


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
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_recv, PDM_MPI_Request);

  int size_recv_type;
  MPI_Type_size(_pdm_mpi_2_mpi_datatype(datatype), &size_recv_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_recv; i++) {
    void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recv_type);
    MPI_Request _mpi_request = MPI_REQUEST_NULL;
    int t_rank = active_recv[i];
    code = MPI_Recv_init(buf,
                         recvcounts[i],
                         _pdm_mpi_2_mpi_datatype(datatype),
                         t_rank,
                         tag,
                         _pdm_mpi_2_mpi_comm(comm),
                         &_mpi_request);
    requests[i] = _mpi_2_pdm_mpi_request_add(_mpi_request);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;
  return _mpi_2_pdm_mpi_err(code);
}


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
)
{
  PDM_MPI_Request *requests = NULL;
  PDM_malloc(requests, n_active_recv, PDM_MPI_Request);

  int size_recv_type;
  MPI_Type_size(_pdm_mpi_2_mpi_datatype(datatype), &size_recv_type);

  int code = MPI_SUCCESS;
  for (int i = 0; i < n_active_recv; i++) {
    void *buf = (void *) ((unsigned char*) recvbuf + sdispls[i] * size_recv_type);
    MPI_Request _mpi_request = MPI_REQUEST_NULL;
    int t_rank = active_recv[i];
    code = MPI_Irecv(buf,
                     recvcounts[i],
                     _pdm_mpi_2_mpi_datatype(datatype),
                     t_rank,
                     tag,
                     _pdm_mpi_2_mpi_comm(comm),
                     &_mpi_request);
    requests[i] = _mpi_2_pdm_mpi_request_add(_mpi_request);
    if (code != MPI_SUCCESS) {
      break;
    }
  }
  *out_requests = requests;

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Wait (wrapping de la fonction MPI_Wait)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Wait(PDM_MPI_Request *request)
{
  int code = MPI_Wait(request, MPI_STATUS_IGNORE);
  assert(code == 0);

  // If we use persistent comm, the request after the wait is not MPI_REQUEST_NULL
  // The request will be free by user with : MPI_Request_free
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Waitall (wrapping de la fonction PDM_MPI_Waitall)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Waitall(int count, PDM_MPI_Request array_of_requests[]) {
  int code = MPI_Waitall(count, array_of_requests, MPI_STATUSES_IGNORE);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Request_free (wrapping de la fonction MPI_Request_free)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Request_free
(
  PDM_MPI_Request *request
)
{
  if(*request == PDM_MPI_REQUEST_NULL) {
    return PDM_MPI_SUCCESS;
  }
  int code = MPI_Request_free(request);
  return code;;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Test (wrapping de la fonction MPI_Test)
 *----------------------------------------------------------------------------*/

int PDM_MPI_Test(PDM_MPI_Request *request, int *flag)
{
  // Test was already done
  if(*request == MPI_REQUEST_NULL) {
    *flag = 1;
    return MPI_SUCCESS;
  }
  int code = MPI_Test(request, flag, MPI_STATUS_IGNORE);

  if(*flag == 0) {
    return code; // Message was not ready
  }

  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Startall (wrapping de la fonction MPI_Startall)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Startall(int count, PDM_MPI_Request array_of_requests[])
{
  int code = MPI_SUCCESS;
  for(int i = 0; i < count; ++i) {
    code = MPI_Start(&array_of_requests[i]);
    if(code != MPI_SUCCESS) {
      break;
    }
  }
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Start (wrapping de la fonction MPI_Start)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Start(PDM_MPI_Request *request)
{
  int code = MPI_SUCCESS;
  code = MPI_Start(request);
  return code;
}

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
)
{
  MPI_Aint *_array_of_displacements;
  PDM_malloc(_array_of_displacements,count,MPI_Aint);

  for (int i = 0; i < count; i++) {
    _array_of_displacements[i] = array_of_displacements[i];
  }

  int code = MPI_Type_create_hindexed(count,
                                      array_of_blocklengths,
                                      _array_of_displacements,
                                      oldtype,
                                      newtype);
  PDM_free(_array_of_displacements);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_create_contiguous(int               count,
                                   PDM_MPI_Datatype  old_datatype,
                                   PDM_MPI_Datatype *newtype)
{
  int code = MPI_Type_contiguous(count,
                                 old_datatype,
                                 newtype);
  assert(code == 0);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_commit (wrapping de la fonction MPI_Type_commit)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_commit(PDM_MPI_Datatype *datatype)
{
  int code =  MPI_Type_commit(datatype);
  return code;
}

/*----------------------------------------------------------------------------
 * MPI_Type_size (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_size(PDM_MPI_Datatype datatype, int *size)
{
  return MPI_Type_size(datatype, size);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_free (wrapping de la fonction MPI_Type_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_free(PDM_MPI_Datatype *datatype)
{
  if(*datatype == PDM_MPI_DATATYPE_NULL) {
    return PDM_MPI_SUCCESS;
  }
  MPI_Datatype mpi_type = *datatype;
  int code = MPI_Type_free(&mpi_type);
  *datatype = PDM_MPI_DATATYPE_NULL;
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_f2c (wrapping de la fonction MPI_comm_f2c)
 *----------------------------------------------------------------------------*/
PDM_MPI_Comm PDM_MPI_Comm_f2c(PDM_MPI_Fint comm)
{
  /* Conversion Fortran vers C */
  MPI_Comm _mpi_comm = MPI_Comm_f2c(comm);
  return _mpi_comm;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_c2f (wrapping de la fonction MPI_comm_c2f
 *----------------------------------------------------------------------------*/
PDM_MPI_Fint PDM_MPI_Comm_c2f(PDM_MPI_Comm comm)
{
  /* Conversion Fortran vers C */
  MPI_Comm _mpi_comm = comm;
  PDM_MPI_Fint f_comm = (PDM_MPI_Fint) MPI_Comm_c2f(_mpi_comm);
  return f_comm;
}

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
)
{
  int code = MPI_Scatter(sendbuf, sendcount, sendtype,
                         recvbuf, recvcount, recvtype,
                         root, comm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Barrier (wrapping de la fonction MPI_Barrier)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Barrier(PDM_MPI_Comm comm)
{
  int code =  MPI_Barrier(comm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Wtime (wrapping de la fonction MPI_Wtime)
 *----------------------------------------------------------------------------*/
double PDM_MPI_Wtime(void)
{
  return MPI_Wtime();
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Alloc_mem (wrapping de la fonction MPI_Alloc_mem)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Alloc_mem(PDM_MPI_Aint size, void *baseptr) {
  return MPI_Alloc_mem(size, MPI_INFO_NULL, baseptr);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Free_mem (wrapping de la fonction MPI_Free_mem)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Free_mem(void *baseptr) {
  return MPI_Free_mem(baseptr);
}

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
)
{
  int code = MPI_Bcast(buffer,
                       count,
                       datatype,
                       root,
                       comm);
  return code;
}

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
)
{

  int code = MPI_Ibcast(buffer,
                        count,
                        datatype,
                        root,
                        comm, request);
  return code;
}

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
)
{
  int code =  MPI_Allgather(sendbuf, sendcount, sendtype,
                            recvbuf, recvcount,
                            recvtype,
                            comm);
  assert(code == MPI_SUCCESS);
  return code;
}

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
)
{
  int code = MPI_Allgatherv(sendbuf, sendcount, sendtype,
                            recvbuf, recvcounts, displs,
                            recvtype, comm);
  return code;
}

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
)
{
  int code = MPI_Reduce(sendbuf, recvbuf, count,
                           datatype,
                           op, root,
                           comm);
  return code;
}

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
)
{
  int code = MPI_Reduce_scatter(sendbuf, recvbuf, counts,
                                datatype,
                                op,
                                comm);
  return code;
}

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
)
{
  int code = MPI_Allreduce(sendbuf,
                           recvbuf,
                           count,
                           datatype,
                           op,
                           comm);
  return code;
}

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
)
{
  int code = MPI_Scan(sendbuf, recvbuf, count,
                      datatype,
                      op,
                      comm);
  return code;
}

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
)
{
  int code = MPI_Exscan(sendbuf, recvbuf, count,
                        datatype,
                        op,
                        comm);
  return code;
}


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
)
{
  int code = MPI_Iscan(sendbuf, recvbuf, count,
                      datatype,
                      op,
                      comm, request);
  return code;
}


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
)
{
  int code = MPI_Alltoall(sendbuf, sendcount,
                          sendtype,
                          recvbuf, recvcount,
                          recvtype,
                          comm);
  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
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
)
{
  int code = MPI_Ialltoall(sendbuf, sendcount,
                          sendtype,
                          recvbuf, recvcount,
                          recvtype,
                          comm, request);
  return code;
}

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
)
{
  int code = MPI_Alltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           sendtype,
                           recvbuf,
                           recvcounts,
                           rdispls,
                           recvtype,
                           comm);

  return code;
}

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
)
{
  int code = MPI_Ialltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           sendtype,
                           recvbuf,
                           recvcounts,
                           rdispls,
                           recvtype,
                           comm, request);

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Get_ialltoallv (Implemtation of alltoall like with window )
 *
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
                           PDM_MPI_Comm      comm)
{
  int code = 0;

  /*
   * Exchange in target view the correct displacement for MPI_Get
   */
  int n_rank, i_rank;
  MPI_Comm_size(comm, &n_rank);
  MPI_Comm_rank(comm, &i_rank);

  int *target_disp;
  PDM_malloc(target_disp, n_rank, int);

  MPI_Alltoall(sdispls    , 1, MPI_INT,
               target_disp, 1, MPI_INT, comm);

  // double t1 = MPI_Wtime();
  for(int i = 0; i < n_rank; ++i) {

    int   origin_data_size = -1;
    MPI_Type_size(recvtype, &origin_data_size);

    int            origin_displ = rdispls[i]; // + recvcounts[i] *
    unsigned char *origin_addr  = (unsigned char *) recvbuf + origin_displ * origin_data_size;
    int            origin_count = recvcounts[i];

    if(origin_count > 0 && i != i_rank) {
      MPI_Get(origin_addr,
              origin_count,
              recvtype,
              i,
              (MPI_Aint) target_disp[i],
              origin_count,
              sendtype,
              win_send);
    }
  }

  // double dt = MPI_Wtime() - t1;

  // t1 = MPI_Wtime();
  int   origin_data_size = -1;
  MPI_Type_size(recvtype, &origin_data_size);

  int            origin_displ = rdispls[i_rank]; // + recvcounts[i_rank] *
  unsigned char *origin_addr  = (unsigned char *) recvbuf + origin_displ * origin_data_size;
  int            origin_count = recvcounts[i_rank];
  MPI_Get(origin_addr,
          origin_count,
          recvtype,
          i_rank,
          (MPI_Aint) target_disp[i_rank],
          origin_count,
          sendtype,
          win_send);

  // dt = MPI_Wtime() - t1;

  PDM_UNUSED(win_recv  );
  PDM_UNUSED(sendcounts);
  PDM_UNUSED(sendbuf   ); // Implicit in win_send

  PDM_free(target_disp);
  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv_init (wrapping de la fonction MPI_Alltoallv_init)
 *
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
)
{
#ifdef HAVE_MPI_COLLECTIVE_INIT_FUNC
  int code = MPI_Alltoallv_init(sendbuf,
                                sendcounts,
                                sdispls,
                                sendtype,
                                recvbuf,
                                recvcounts,
                                rdispls,
                                recvtype,
                                comm,
                                request);
  return code;
#else
  PDM_UNUSED(sendbuf);
  PDM_UNUSED(sendcounts);
  PDM_UNUSED(sdispls);
  PDM_UNUSED(sendtype);
  PDM_UNUSED(recvbuf);
  PDM_UNUSED(recvcounts);
  PDM_UNUSED(rdispls);
  PDM_UNUSED(recvtype);
  PDM_UNUSED(comm);
  PDM_UNUSED(request);
  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_Alltoallv_Init : Persistent collective communication not available !");
  return -1;
#endif
}

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
)
{
#ifdef HAVE_MPI_COLLECTIVE_INIT_FUNC
  int code = MPI_Neighbor_alltoallv_init(sendbuf,
                                         sendcounts,
                                         sdispls,
                                         sendtype,
                                         recvbuf,
                                         recvcounts,
                                         rdispls,
                                         recvtype,
                                         comm,
                                         request);
  return code;
#else
  PDM_UNUSED(sendbuf);
  PDM_UNUSED(sendcounts);
  PDM_UNUSED(sdispls);
  PDM_UNUSED(sendtype);
  PDM_UNUSED(recvbuf);
  PDM_UNUSED(recvcounts);
  PDM_UNUSED(rdispls);
  PDM_UNUSED(recvtype);
  PDM_UNUSED(comm);
  PDM_UNUSED(request);
  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_Neighbor_alltoallv_init : Persistent collective communication not available !");
  return -1;
#endif
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_create (wrapping de la fonction MPI_Win_create)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_create(void         *baseptr,
                       PDM_MPI_Aint  size,
                       int           disp_unit,
                       PDM_MPI_Comm  comm,
                       PDM_MPI_Win  *win)
{
  int code = MPI_Win_create(baseptr,
                            size,
                            disp_unit,
                            MPI_INFO_NULL,
                            comm,
                            win);
  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Win_allocate (wrapping de la fonction MPI_Win_allocate)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_allocate(PDM_MPI_Aint  size,
                         int           disp_unit,
                         PDM_MPI_Comm  comm,
                         void         *baseptr,
                         PDM_MPI_Win  *win)
{
  int code = MPI_Win_allocate(size,
                              disp_unit,
                              MPI_INFO_NULL,
                              comm,
                              baseptr,
                              win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_free (wrapping de la fonction MPI_Win_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_free(PDM_MPI_Win *win)

{
  if(*win == PDM_MPI_WIN_NULL) {
    return PDM_MPI_SUCCESS;
  }
  int code = MPI_Win_free(win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_start (wrapping de la fonction MPI_Win_start)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_start(PDM_MPI_Group group, int mpi_assert, PDM_MPI_Win win)
{
  int code = MPI_Win_start(group, mpi_assert, win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_post (wrapping de la fonction MPI_Win_post)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_post(PDM_MPI_Group group, int mpi_assert, PDM_MPI_Win win)
{
  int code = MPI_Win_post(group, mpi_assert, win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_complete (wrapping de la fonction MPI_Win_complete)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_complete(PDM_MPI_Win win)
{
  int code = MPI_Win_complete(win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_wait (wrapping de la fonction MPI_Win_wait)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_wait(PDM_MPI_Win win)
{
  int code = MPI_Win_wait(win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_fence (wrapping de la fonction MPI_Win_fence)
 *----------------------------------------------------------------------------*/

int PDM_MPI_Win_fence(int assert, PDM_MPI_Win win)
{
  int code = MPI_Win_fence(assert, win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Group_free (wrapping de la fonction MPI_Group_free)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Group_free(PDM_MPI_Group *group)
{
  if(*group == PDM_MPI_GROUP_NULL) {
    return PDM_MPI_SUCCESS;
  }
  int code = MPI_Group_free(group);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_group (wrapping de la fonction MPI_Comm_group)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_group(PDM_MPI_Comm comm, PDM_MPI_Group *group)
{
  int code = MPI_Comm_group(comm, group);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Group_incl (wrapping de la fonction MPI_Group_incl)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Group_incl(PDM_MPI_Group group, int n, const int ranks[],
                       PDM_MPI_Group *newgroup)
{
  int code = MPI_Group_incl(group, n, ranks, newgroup);
  return code;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Topo_test (wrapping de la fonction MPI_Topo_test)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Topo_test(PDM_MPI_Comm comm, int *status) {
  int _status = -1000;
  int code = MPI_Topo_test(comm, &_status);

  if(_status == MPI_UNDEFINED) {
    *status = PDM_MPI_COMM_UNDEFINED;
  } else if(_status == MPI_DIST_GRAPH) {
    *status = PDM_MPI_DIST_GRAPH;
  } else if(_status == PDM_MPI_CART) {
    *status = PDM_MPI_CART;
  } else if(_status == PDM_MPI_GRAPH) {
    *status = PDM_MPI_GRAPH;
  } else {
    PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_Topo_test :"
            " _status '%d' non valide\n", _status);
  }

  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Error_string (wrapping de la fonction MPI_Error_string)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Error_string(int errorcode, char *string, int *resultlen)
{
  int code = MPI_Error_string(errorcode, string, resultlen);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_rank(PDM_MPI_Comm comm, int *rank)
{
  int code = MPI_Comm_rank(comm, rank);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_size(PDM_MPI_Comm comm, int *size)
{
  int code = MPI_Comm_size(comm, size);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_get_max_error_string
 *----------------------------------------------------------------------------*/
int PDM_MPI_get_max_error_string(void)
{
  return MPI_MAX_ERROR_STRING;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_free
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_free(PDM_MPI_Comm *comm)
{
  int code = 0;
  code = MPI_Comm_free(comm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_split(PDM_MPI_Comm comm, int color, int key, PDM_MPI_Comm *newcomm)
{
  int code = MPI_Comm_split(comm, color, key, newcomm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_dup
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_dup(PDM_MPI_Comm comm, PDM_MPI_Comm *newcomm)
{
  int code = MPI_Comm_dup(comm, newcomm);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split_type_numa // Non portable mettre un ifdef
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Comm_split_type_numa
(
 PDM_MPI_Comm comm,
 PDM_MPI_Comm *comm_numa
)
{
  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_SHARED, &comm_node);

  int i_rank_node;
  PDM_MPI_Comm_rank(comm_node, &i_rank_node);

  int i_cpu;
  int i_numa = 0;
#ifdef __linux__
  syscall(SYS_getcpu, &i_cpu, &i_numa, NULL);
#elif __APPLE__
#else
  printf("PDM_MPI_Comm_split_type_numa : appel a SYS_getcpu commente car non portable : a reintroduire après tests dans CMake\n");
  abort();
#endif

  /* Sur le shared on split par numa */
  int code = PDM_MPI_Comm_split(comm_node, i_numa, i_rank_node, comm_numa);

  // *comm_numa = comm_node;
  // PDM_MPI_Comm_free(&comm_node);
  // int code = 0;
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split_type
 *----------------------------------------------------------------------------*/
int PDM_MPI_Comm_split_type(PDM_MPI_Comm comm, int split_type, PDM_MPI_Comm *newcomm)
{
  int i_rank;
  MPI_Comm_rank(comm, &i_rank);

  // PDM_MPI_Comm _newcomm;
  int code = 0;
  if(split_type == PDM_MPI_SPLIT_SHARED) {
    code = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, i_rank /* Key */,
                               MPI_INFO_NULL,
                               newcomm);
  } else if(split_type == PDM_MPI_SPLIT_NUMA) {
    PDM_MPI_Comm_split_type_numa(comm, newcomm);
  } else {
    PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_Comm_split_type :"
            " split_type '%d' non valide\n", split_type);
    abort();
  }
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_allocate_shared_get *
 *----------------------------------------------------------------------------*/
PDM_mpi_win_shared_t*
PDM_mpi_win_shared_create(PDM_MPI_Aint size,
                          int          disp_unit,
                          PDM_MPI_Comm comm)
{
  PDM_mpi_win_shared_t *wins;
  PDM_malloc(wins,1,PDM_mpi_win_shared_t);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  MPI_Info info;
  MPI_Info_create( &info );
  // MPI_Info_set(info, "no_locks", "true");
  // MPI_Info_set( info, "alloc_shared_noncontig", "true" );

  wins->win = MPI_WIN_NULL;
  wins->ptr = NULL;
  int res = 0;
  if(i_rank == 0) {
    res = MPI_Win_allocate_shared(size * disp_unit, disp_unit, info, comm, &wins->ptr , &wins->win);
  } else {
    res = MPI_Win_allocate_shared(0, disp_unit, info , comm, &wins->ptr , &wins->win );
    MPI_Aint size_0;
    int disp_0;
    MPI_Win_shared_query(wins->win, 0, &size_0, &disp_0, &wins->ptr);
  }
  MPI_Info_free(&info);
  assert(res == PDM_MPI_SUCCESS);
  return wins;
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_get
 *----------------------------------------------------------------------------*/
void* PDM_mpi_win_shared_get(PDM_mpi_win_shared_t *wins){
  return wins->ptr;
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_free
 *----------------------------------------------------------------------------*/
void PDM_mpi_win_shared_free(PDM_mpi_win_shared_t *wins){
  MPI_Win_free(&wins->win);
  wins->ptr = NULL;
  PDM_free(wins);
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_lock_all
 *----------------------------------------------------------------------------*/
int PDM_mpi_win_shared_lock_all(int assert, PDM_mpi_win_shared_t* win)
{
  int code = MPI_Win_lock_all(assert, win->win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_unlock_all
 *----------------------------------------------------------------------------*/
int PDM_mpi_win_shared_unlock_all(PDM_mpi_win_shared_t* win)
{
  int code = MPI_Win_unlock_all(win->win);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_sync
 *----------------------------------------------------------------------------*/
int PDM_mpi_win_shared_sync(PDM_mpi_win_shared_t* win)
{
  int code = MPI_Win_sync(win->win);
  return code;
}

// ------------------------------------------------------------------
PDM_MPI_Comm PDM_MPI_get_group_of_master(PDM_MPI_Comm comm, PDM_MPI_Comm sub_comm)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int i_rank_sub;
  PDM_MPI_Comm_rank(sub_comm, &i_rank_sub);

  PDM_MPI_Comm master_of_sub_comm;
  int res = 0;
  if(i_rank_sub == 0){
    res = PDM_MPI_Comm_split(comm, 0, i_rank, &master_of_sub_comm);
  } else {
    res = PDM_MPI_Comm_split(comm, PDM_MPI_UNDEFINED, i_rank, &master_of_sub_comm);
    assert(master_of_sub_comm == PDM_MPI_COMM_NULL);
  }
  assert(res == PDM_MPI_SUCCESS);
  return master_of_sub_comm;
}

// ------------------------------------------------------------------
int PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_Comm comm, void *attribute_val, int *flag)
{
  int code = MPI_Comm_get_attr(comm, MPI_TAG_UB, attribute_val, flag);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_rand_tag_get
 *----------------------------------------------------------------------------*/
int PDM_MPI_Rand_tag (PDM_MPI_Comm comm)
{
  struct timeval t;
  gettimeofday(&t, NULL);

  long ltag = t.tv_usec + 1000000 * t.tv_sec;

  MPI_Bcast (&ltag, 1, MPI_LONG, 0, comm);

  void  *max_tag_tmp;
  int flag;

  // Mandatory to call with PDM_MPI_COMM_WORLD becuase only this one keep attributes (openMPI implemntation for exemple)
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &max_tag_tmp, &flag);
  long max_tag = (long) (*((int *) max_tag_tmp));

  // printf("max_tag = %li | ltag = %li \n", max_tag, ltag);

  return (int) (ltag % max_tag);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Dist_graph_create_adjacent
 *----------------------------------------------------------------------------*/
int PDM_MPI_Dist_graph_create_adjacent(PDM_MPI_Comm  comm_old,
                                             int     indegree,
                                       const int     sources[],
                                             int     outdegree,
                                       const int     destinations[],
                                       int           reorder,
                                       PDM_MPI_Comm *newcomm)
{
  const int *weight_in  = MPI_UNWEIGHTED;
  const int *weight_out = MPI_UNWEIGHTED;
  int code = MPI_Dist_graph_create_adjacent(comm_old,
                                            indegree,
                                            sources,
                                            weight_in,
                                            outdegree,
                                            destinations,
                                            weight_out,
                                            MPI_INFO_NULL,
                                            reorder,
                                            newcomm);
  return code;
}

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
)
{
  int code =  MPI_Neighbor_allgather(sendbuf,
                                     sendcount,
                                     sendtype,
                                     recvbuf,
                                     recvcount,
                                     recvtype,
                                     comm);
  return code;
}

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
)
{
  int code = MPI_Neighbor_allgatherv(sendbuf,
                                     sendcount,
                                     sendtype,
                                     recvbuf,
                                     recvcounts,
                                     displs,
                                     recvtype,
                                     comm);
  return code;
}


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
)
{
  int code = MPI_Neighbor_alltoall(sendbuf, sendcount,
                                   sendtype,
                                   recvbuf, recvcount,
                                   recvtype,
                                   comm);
  return code;
}

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
)
{
  int code = MPI_Ineighbor_alltoall(sendbuf, sendcount,
                                    sendtype,
                                    recvbuf, recvcount,
                                    recvtype,
                                    comm, request);
  return code;
}

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
)
{
  int code = MPI_Neighbor_alltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           sendtype,
                           recvbuf,
                           recvcounts,
                           rdispls,
                           recvtype,
                           comm);

  return code;
}

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
)
{
  int code = MPI_Ineighbor_alltoallv(sendbuf,
                                     sendcounts,
                                     sdispls,
                                     sendtype,
                                     recvbuf,
                                     recvcounts,
                                     rdispls,
                                     recvtype,
                                     comm, request);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Dist_graph_neighbors_count (wrapping de la fonction MPI_Dist_graph_neighbors_count)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Dist_graph_neighbors_count
(
  PDM_MPI_Comm  comm,
  int          *n_degree_in,
  int          *n_degree_out,
  int          *is_weighted
)
{
  int code = MPI_Dist_graph_neighbors_count(comm,
                                            n_degree_in,
                                            n_degree_out,
                                            is_weighted);
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Dist_graph_neighbors (wrapping de la fonction MPI_Dist_graph_neighbors)
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Dist_graph_neighbors
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *sources,
  int            n_degree_out,
  int           *destinations
)
{

  int *weight_in  = NULL;
  int *weight_out = NULL;
  int code = MPI_Dist_graph_neighbors(comm,
                                      n_degree_in,
                                      sources,
                                      weight_in,
                                      n_degree_out,
                                      destinations,
                                      weight_out);
  return code;
}

/*----------------------------------------------------------------------------
 * MPI Standard extension :
 *   - Add method to ease MPI
 *   - Futur method but not yet supported (ex: MPI_Neigbor_init)
 *----------------------------------------------------------------------------*/

void
PDM_MPI_setup_hybrid_dist_comm_graph
(
  PDM_MPI_Comm   comm,
  PDM_MPI_Comm  *comm_shared_out,
  PDM_MPI_Comm  *comm_dist_graph_out,
  int           *n_degree,
  int          **neighbor
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_MPI_Comm comm_master_of_shm = PDM_MPI_get_group_of_master(comm, comm_shared);

  int i_rank_master_of_shm = -1;
  int n_rank_master_of_shm;
  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Comm_rank(comm_master_of_shm, &i_rank_master_of_shm);
    PDM_MPI_Comm_size(comm_master_of_shm, &n_rank_master_of_shm);
  }
  PDM_MPI_Bcast(&n_rank_master_of_shm, 1, PDM_MPI_INT, 0, comm_shared);
  PDM_MPI_Bcast(&i_rank_master_of_shm, 1, PDM_MPI_INT, 0, comm_shared);

  PDM_mpi_win_shared_t* wnuma_by_numa_n = PDM_mpi_win_shared_create(n_rank_master_of_shm, sizeof(int), comm_shared);
  int *numa_by_numa_n  = PDM_mpi_win_shared_get(wnuma_by_numa_n);
  PDM_mpi_win_shared_lock_all (0, wnuma_by_numa_n);


  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Allgather(&n_rank_in_shm, 1, PDM_MPI_INT,
                      numa_by_numa_n, 1, PDM_MPI_INT, comm_master_of_shm);
  }
  PDM_mpi_win_shared_sync(wnuma_by_numa_n);
  PDM_MPI_Barrier(comm_shared);

  int n_tot_numa = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    n_tot_numa += numa_by_numa_n[i];
  }
  /*
   * Create idx  and  gid of each numa
   */
  PDM_mpi_win_shared_t* wnuma_core_gid    = PDM_mpi_win_shared_create(n_tot_numa               , sizeof(int), comm_shared);
  PDM_mpi_win_shared_t* wnuma_by_numa_idx = PDM_mpi_win_shared_create(n_rank_master_of_shm+1, sizeof(int), comm_shared);
  int *numa_core_gid    = PDM_mpi_win_shared_get(wnuma_core_gid);
  int *numa_by_numa_idx = PDM_mpi_win_shared_get(wnuma_by_numa_idx);
  PDM_mpi_win_shared_lock_all (0, wnuma_core_gid);
  PDM_mpi_win_shared_lock_all (0, wnuma_by_numa_idx);


  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    numa_by_numa_idx[0] = 0;
    for(int i = 0; i < n_rank_master_of_shm; ++i) {
      numa_by_numa_idx[i+1] = numa_by_numa_idx[i] + numa_by_numa_n[i];
    }
  }
  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_by_numa_idx);

  numa_core_gid[numa_by_numa_idx[i_rank_master_of_shm]+i_rank_in_shm] = i_rank;

  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_core_gid);

  /*
   *  Exchange of the global numbering of rank for each NUMA
   */
  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    int *lnuma_core_gid;
    PDM_malloc(lnuma_core_gid,n_rank_in_shm ,int);
    for(int i = 0; i < n_rank_in_shm; ++i) {
      lnuma_core_gid[i] = numa_core_gid[numa_by_numa_idx[i_rank_master_of_shm]+i];
    }
    PDM_MPI_Allgatherv(lnuma_core_gid, n_rank_in_shm, PDM_MPI_INT,
                       numa_core_gid , numa_by_numa_n, numa_by_numa_idx, PDM_MPI_INT, comm_master_of_shm);
    PDM_free(lnuma_core_gid);
  }
  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_core_gid);

  /*
   * Computation of degree_in
   */
  int *send_n;
  int *recv_n;
  int *send_idx;
  int *recv_idx;
  PDM_malloc(send_n  , n_rank  , int);
  PDM_malloc(recv_n  , n_rank  , int);
  PDM_malloc(send_idx, n_rank+1, int);
  PDM_malloc(recv_idx, n_rank+1, int);

  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
    recv_n[i] = 0;
  }

  int n_degrees_in = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int lid_rank = (j - numa_by_numa_idx[i]) % n_rank_in_shm; // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        n_degrees_in++;
      }
    }
  }

  int *neighbor_in;
  PDM_malloc(neighbor_in, (n_degrees_in ) ,int);
  n_degrees_in = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int gid_rank = numa_core_gid[j];
      int lid_rank = (j - numa_by_numa_idx[i]) % n_rank_in_shm;  // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        neighbor_in[n_degrees_in++] = gid_rank;
      }
    }
  }

  for(int i = 0; i < n_degrees_in; ++i) {
    send_n[neighbor_in[i]]++;
  }

  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    send_n[i] = 0;
  }

  int *send_cur_i_rank;
  PDM_malloc(send_cur_i_rank,send_idx[n_rank] ,int);

  for(int i = 0; i < n_degrees_in; ++i) {
    int idx_write = send_idx[neighbor_in[i]] + send_n[neighbor_in[i]]++;
    send_cur_i_rank[idx_write] = i_rank;
  }

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  int *recv_opp_i_rank;
  PDM_malloc(recv_opp_i_rank,recv_idx[n_rank] ,int);

  PDM_MPI_Alltoallv(send_cur_i_rank, send_n, send_idx, PDM_MPI_INT,
                    recv_opp_i_rank, recv_n, recv_idx, PDM_MPI_INT, comm);


  int n_degrees_out = recv_idx[n_rank];
  int *neighbor_out = recv_opp_i_rank; // Already sort normaly

  PDM_free(send_n);
  PDM_free(recv_n);
  PDM_free(send_idx);
  PDM_free(recv_idx);
  PDM_free(send_cur_i_rank);

  PDM_MPI_Comm comm_dist_graph;
  PDM_MPI_Dist_graph_create_adjacent(comm,
                                     n_degrees_in,
                                     neighbor_in,
                                     n_degrees_out,
                                     neighbor_out,
                                     0,
                                     &comm_dist_graph);

  PDM_mpi_win_shared_unlock_all(wnuma_by_numa_n);
  PDM_mpi_win_shared_unlock_all(wnuma_core_gid);
  PDM_mpi_win_shared_unlock_all(wnuma_by_numa_idx);
  PDM_mpi_win_shared_free(wnuma_by_numa_n);
  PDM_mpi_win_shared_free(wnuma_core_gid);
  PDM_mpi_win_shared_free(wnuma_by_numa_idx);

  PDM_free(recv_opp_i_rank);

  *comm_shared_out     = comm_shared;
  *comm_dist_graph_out = comm_dist_graph;

  *n_degree = n_degrees_in;
  *neighbor = neighbor_in;
}



void
PDM_MPI_setup_dist_graph_from_neighbor_in
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *neighbor_in,
  PDM_MPI_Comm  *comm_dist_graph_out
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int *send_n = NULL;
  int *recv_n = NULL;
  PDM_malloc(send_n, n_rank, int);
  PDM_malloc(recv_n, n_rank, int);

  int *send_idx = NULL;
  int *recv_idx = NULL;
  PDM_malloc(send_idx, n_rank+1, int);
  PDM_malloc(recv_idx, n_rank+1, int);

  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
    recv_n[i] = 0;
  }

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  for(int i = 0; i < n_degree_in; ++i) {
    send_n[neighbor_in[i]]++;
  }

  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    send_n[i] = 0;
  }

  int *send_cur_i_rank;
  PDM_malloc(send_cur_i_rank,send_idx[n_rank] ,int);

  for(int i = 0; i < n_degree_in; ++i) {
    int idx_write = send_idx[neighbor_in[i]] + send_n[neighbor_in[i]]++;
    send_cur_i_rank[idx_write] = i_rank;
  }

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  int *recv_opp_i_rank;
  PDM_malloc(recv_opp_i_rank,recv_idx[n_rank] ,int);

  PDM_MPI_Alltoallv(send_cur_i_rank, send_n, send_idx, PDM_MPI_INT,
                    recv_opp_i_rank, recv_n, recv_idx, PDM_MPI_INT, comm);


  int n_degrees_out = recv_idx[n_rank];
  int *neighbor_out = recv_opp_i_rank; // Already sort normaly

  PDM_free(send_n);
  PDM_free(recv_n);
  PDM_free(send_idx);
  PDM_free(recv_idx);
  PDM_free(send_cur_i_rank);

  PDM_MPI_Dist_graph_create_adjacent(comm,
                                     n_degree_in,
                                     neighbor_in,
                                     n_degrees_out,
                                     neighbor_out,
                                     0,
                                     comm_dist_graph_out);
}

int
PDM_MPI_Comm_compare
(
  PDM_MPI_Comm  comm1,
  PDM_MPI_Comm  comm2,
  int          *result
)
{
  return MPI_Comm_compare(comm1, comm2, result);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
