/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/* #if !defined(_XOPEN_SOURCE) || !defined(_BSD_SOURCE) */
/* #define _XOPEN_SOURCE 500 */
/* #endif */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <stddef.h>
#include <stdint.h>

#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Définitions des macro locales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Default hostname size
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Max
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*============================================================================
 * Variables globales
 *============================================================================*/

static const int  local_hostname_default_len = 64;

/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction de hachage adler32
 *
 * parameters :
 *   buf              <-- Buffer
 *   buflength        <-- Taille du buffer
 * return
 *   Cle de hachage
 *----------------------------------------------------------------------------*/

static uint32_t
 _adler32
(
const void *buf,
size_t buflength
)
{

  const uint8_t * buffer = (const uint8_t *)buf;

  uint32_t s1 = 1;
  uint32_t s2 = 0;

  for (size_t n = 0; n < buflength; n++) {
    s1 = (s1 + buffer[n]) % 65521;
    s2 = (s2 + s1) % 65521;
  }

  return (s2 << 16) | s1;
}


/*----------------------------------------------------------------------------
 * Rang du processus courant sur son noeud
 *
 * parameters :
 *   comm             <-- Communicateur MPI global
 * return
 *   Rank
 *----------------------------------------------------------------------------*/

static int
_node_rank_by_hash
(
PDM_MPI_Comm comm
)
{

  char * hostname = NULL;
  size_t hostname_length = 0;

  // TODO: Use MPI_Get_processor_name instead of PDM_io_get_hostname
  //  int MPI_Get_processor_name(char *name, int *resultlen)
  PDM_io_get_hostname(&hostname, &hostname_length);

  /* Détermine la clé de hachage */

  uint32_t check_sum = _adler32(hostname, hostname_length);

  if (0 == 1)
    PDM_printf( "-- hostname adlerCode %s %ud\n", hostname, check_sum);

  int comm_rank = -1;
  PDM_MPI_Comm_rank(comm, &comm_rank);

  PDM_MPI_Comm node_comm = PDM_MPI_COMM_NULL;

  /* Détermine le communicateur local au noeud */

  PDM_MPI_Comm_split(comm, check_sum, comm_rank, &node_comm);

  /* Détermine le rang dans communicateur du noeud */

  int node_rank;
  PDM_MPI_Comm_rank(node_comm, &node_rank);

  int node_size;
  PDM_MPI_Comm_size(node_comm, &node_size);

  /* Vérifie si 2 noms de noeud ne rentrent pas en collision :
     Si c'est le cas, les processus de ces 2 noeuds sont
     regroupés dans un même communicateur */

  int n_send = (int) hostname_length;

  PDM_MPI_Allreduce((void *) &hostname_length, (void *) &n_send, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, node_comm);

  n_send += 1; /* \0 */

  char *send;
  PDM_malloc(send, n_send, char);

  for (int i = 0; i < n_send; i++)
    send[i] = 0x00;

  strncpy(send, hostname, n_send);

  char *recv;
  PDM_malloc(recv, n_send * node_size, char);
  PDM_MPI_Allgather((void * )send, n_send, PDM_MPI_CHAR, (void *)recv, n_send, PDM_MPI_CHAR, node_comm);

  char *neighbor = recv;
  int local_node_rank = 0;

  for (int i = 0; i < node_size; ++i) {
    if (strcmp(send, neighbor) == 0) {
      if (i < node_rank) {
        ++local_node_rank;
      }
      else {
        break;
      }
    }
    neighbor += n_send;
  }

  int is_same_node = (node_rank != local_node_rank);
  int is_same_node_max;

  PDM_MPI_Allreduce((void *) &is_same_node, (void *) &is_same_node_max, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, node_comm);

  if (is_same_node_max == 1) {

    /* Traitement des collisions */

    neighbor = recv;
    char *host_in_comm;
    PDM_malloc(host_in_comm, n_send * node_size, char);
    strncpy(host_in_comm, recv, n_send);
    int n_host_in_comm = 1;

    for (int j = 1; j < node_size; j++) {
      char *pt_host_in_comm = host_in_comm;
      for (int i = 0; i < n_host_in_comm; i++) {
        if (strcmp(pt_host_in_comm, neighbor) != 0) {
          strncpy(host_in_comm + n_host_in_comm * n_send, neighbor, n_send);
          n_host_in_comm += 1;
          break;
        }
        pt_host_in_comm += n_send;
      }
      neighbor += n_send;
    }

    int tag = 0;
    for (int i = 0; i < n_host_in_comm; i++) {
      char *pt_host_in_comm = host_in_comm;
      if (strcmp(pt_host_in_comm, send) == 0) {
        tag = i;
        break;
      }
      pt_host_in_comm += n_send;
    }

    PDM_MPI_Comm node_comm2;
    PDM_MPI_Comm_split(node_comm, tag, node_rank, &node_comm2);

    int node_rank2;
    PDM_MPI_Comm_rank(node_comm2, &node_rank2);

    node_rank = node_rank2;

    PDM_MPI_Comm_free(&node_comm2);
    PDM_free(host_in_comm);

  }

  // Clean up.

  PDM_free(send);
  send = NULL;

  PDM_free(recv);
  recv = NULL;

  PDM_MPI_Comm_free(&node_comm);

  PDM_free(hostname);
  hostname = NULL;

  // Affichage

  if (0 == 1) {

    int master_node_rank = 0;
    int nb_io_rank = 0;
    if (node_rank == 0)
      master_node_rank = 1;


    PDM_MPI_Allreduce((void *) &master_node_rank, (void *) &nb_io_rank, 1,
                      PDM_MPI_INT, PDM_MPI_SUM, comm);

    int comm_size;
    PDM_MPI_Comm_size(comm, &comm_size);
    if (comm_rank == 0)
      PDM_printf("Part of ranks for parallel IO : %d/%d\n", nb_io_rank, comm_size);

  }

  return node_rank;
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Retourne le rang du procesus courant dans le noeud
 *
 * parameters :
 *   comm            <-- Communicateur MPI
 *
 * return :
 *   Rank
 *
 *----------------------------------------------------------------------------*/

int
PDM_io_mpi_node_rank
(
PDM_MPI_Comm comm
)
{

  int node_rank = -1;

  node_rank = _node_rank_by_hash(comm);

  char * hostname;
  size_t hostname_length;

  PDM_io_get_hostname(&hostname, &hostname_length);
  PDM_free(hostname);
  hostname = NULL;

  return node_rank;
}

/*----------------------------------------------------------------------------
 * Retourne le nom du noeud
 *
 * parameters :
 *   hostname_ptr       <-> Pointeur sur la chaine contenant le nom
 *   hostname_length    <-> Longueur de la chaine
 *
 * return :
 *   Code d'erreur
 *
 *----------------------------------------------------------------------------*/

void
PDM_io_get_hostname
(
char  **hostname_ptr,
size_t *hostname_length
)
{

  char * hostname = NULL;
  size_t nHostname = 0;

  int  allocate_more = 0;

  *hostname_ptr = NULL;
  *hostname_length = 0;

  do {
    nHostname += local_hostname_default_len;

    PDM_malloc(hostname, nHostname, char);

    if (hostname == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "allocating %lu bytes of memory for hostname failed: %s.\n",
              sizeof(char) * nHostname, strerror(errno));
      abort();
    }

    int error;

    error = gethostname(hostname, nHostname);

    if (error == -1) {
      if (errno == ENAMETOOLONG) {
        allocate_more = 1;
        PDM_free(hostname); hostname = NULL;
      }
      else {
        PDM_free(hostname);
        hostname = NULL;

        PDM_error(__FILE__, __LINE__, 0, "gethostname failed with error %d: %s.\n", errno, strerror(errno));
        abort();
      }

    }
    else {
      allocate_more = 0;
    }

  } while (allocate_more == 1);

  hostname[nHostname - 1] = 0x00;

  *hostname_length = strlen(hostname);
  *hostname_ptr = hostname;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
