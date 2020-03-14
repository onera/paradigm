/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_distant_neighbor.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_mpi.h"
#include "pdm_order.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/**
 * \struct _distant_neighbor_t
 * \brief  Define a point merge structures
 *
 */
typedef struct  {
  PDM_MPI_Comm comm;                  /*!< MPI communicator */
  int          n_part;                /*!< Number of partitions */
  const int   *n_entity;              /*!< Number of entities for each partition */
  int        **neighbor_idx;          /*!< Indexes of candidate for each current part point
                                       *   (size = number of entities in the current part + 1) */
  int        **neighbor_desc;         /*!< Candidates description (process,
                                       *                           part in the process,
                                       *                           entitiy number in the part) */
  int**        order;
  int**        order_unique;
  int         *requested_data_n;      /*!< Numer of requested data for each process index
                                       * (size : s_comm) */
  int         *requested_data_idx;    /*!< Requested data for each process index
                                       * (size : s_comm) */
  int         *distributed_data_n;    /*!< Numer of distributed data for each process index
                                       * (size : s_comm) */
  int         *distributed_data_idx;  /*!< Distributed data for each process index
                                       * (size : s_comm) */
  int         *distributed_data;      /*!< Distributed data for each process
                                       * (size : requestd_data_idx[s_comm - 1] */
  int         *distributed_part_n;
  int         *distributed_part_idx; /*!< For each part the shift to apply on recv buffer
                                       * (size : n_partà )*/

} _distant_neighbor_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

static PDM_Handles_t *_pdns   = NULL;

/*=============================================================================
 * Static function definitions
 *============================================================================*/


static inline
int
is_same_triplet
(
int iproc1, int ipart1, int ielt1,
int iproc2, int ipart2, int ielt2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        return 1;
      }
    }
  }
  return 0;
}

static
void
compute_unique_idx
(
 int order[],
 int order_unique[],
 int connect_triplet[],
 const size_t nb_ent
)
{
  int idx_unique = -1;
  int last_proc  = -1;
  int last_part  = -1;
  int last_elmt  = -1;

  for(int i = 0; i < nb_ent; i++){

    int old_order = order[i];
    int curr_proc = connect_triplet[3*old_order  ];
    int curr_part = connect_triplet[3*old_order+1];
    int curr_elmt = connect_triplet[3*old_order+2];
    int is_same = is_same_triplet(last_proc, last_part, last_elmt,
                                  curr_proc, curr_part, curr_elmt);
    printf(" curr:: ( %d / %d / %d ) | last:: ( %d / %d / %d ) \n",
           curr_proc, curr_part, curr_elmt,
           last_proc, last_part, last_elmt);

    if(is_same == 0){ // N'est pas le meme
      idx_unique++;
      last_proc = curr_proc;
      last_part = curr_part;
      last_elmt = curr_elmt;
    }
    order_unique[i] = idx_unique;
    printf("[%d] = %d --> %d \n", i, is_same, idx_unique);
  }


  if(1 == 1){
    printf("order_unique:: \n");
    for(int i = 0; i < nb_ent; i++){
      printf(" -------------------------- \n");
      // int pos_unique = order_unique_j1[i];
      // int curr_idx   = order_j1[pos_unique];
      int pos_unique = order_unique[i];
      int curr_idx   = order[i];

      int curr_proc  = connect_triplet[3*curr_idx  ];
      int curr_part  = connect_triplet[3*curr_idx+1];
      int curr_elmt  = connect_triplet[3*curr_idx+2];

      printf("\t pos_unique :: %d \n", pos_unique);
      printf("\t curr_idx   :: %d \n", curr_idx  );
      printf("\t triplet    :: ( %d / %d / %d ) \n", curr_proc, curr_part, curr_elmt);

    }
    printf("\n");
  }
}

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */
static _distant_neighbor_t *
_get_from_id
(
 int  id
)
{
  _distant_neighbor_t *pdn = (_distant_neighbor_t *) PDM_Handles_get (_pdns, id);

  if (pdn == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_distant_neighbor error : Bad identifier\n");
  }

  return pdn;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Return an initialized \ref _distant_neighbor_t structure
 *
 * This function returns an initialized \ref _distant_neighbor_t structure
 *
 * \param [in]   comm          MPI communicator
 * \param [in]   n_part        Number of partitions
 * \param [in]   n_entity      Number of entities for each partition
 * \param [out]  neighbor_idx  Indexes of candidate for each current part point
 *                              (size = number of entities in the current part + 1)
 * \param [out]  neighbor_desc Candidates description (process,
 *                                                     part in the process,
 *                                                     entitiy in the part)
 *
 * \return      A new initialized \ref PDM_distant_neighbor structure
 *
 */
int
PDM_distant_neighbor_create
(
const PDM_MPI_Comm   comm,
const int            n_part,
const int           *n_entity,
      int          **neighbor_idx,
      int          **neighbor_desc
)
{
  printf(" PDM_distant_neighbor_create \n");

  if (_pdns == NULL) {
    _pdns = PDM_Handles_create (4);
  }

  _distant_neighbor_t *pdn = (_distant_neighbor_t *) malloc(sizeof(_distant_neighbor_t));

  int id = PDM_Handles_store (_pdns, pdn);

  pdn->comm          = comm;

  int iRank;
  int nRank;
  PDM_MPI_Comm_rank(pdn->comm, &iRank);
  PDM_MPI_Comm_size(pdn->comm, &nRank);

  char filename[50];
  sprintf(filename, "pdm_logging_%d.log", iRank);
  FILE* fp = fopen(filename, "w");
  log_set_fp(fp);
  log_set_quiet(1);
  log_set_level(-1);
  // log_trace("PDM_distant_neighbor_create::SuperTest");

  pdn->n_part                 = n_part;
  pdn->n_entity               = n_entity;
  pdn->neighbor_idx           = neighbor_idx;
  pdn->neighbor_desc          = neighbor_desc;
  pdn->order                  = (int **) malloc(   pdn->n_part       * sizeof(int **));
  pdn->order_unique           = (int **) malloc(   pdn->n_part       * sizeof(int **));
  pdn->requested_data_n       = (int * ) malloc( ( nRank           ) * sizeof(int * ));
  pdn->requested_data_idx     = (int * ) malloc( ( nRank + 1       ) * sizeof(int * ));
  pdn->distributed_part_n     = (int * ) malloc( ( pdn->n_part     ) * sizeof(int * ));
  pdn->distributed_part_idx   = (int * ) malloc( ( pdn->n_part + 1 ) * sizeof(int * ));

  /*
   * Init the requested counter
   */
  for (int i = 0; i < nRank; i++) {
    pdn->requested_data_idx[i] = 0;
    pdn->requested_data_n[i] = 0;
  }

  /*
   * Sort/unique the triplet (iproc, ipart, ientity)
   * Attention on doit faire le tri pour chaque partttion
   *    --> Dans le cas des vertex on aurai 2 joins qui demande 2 fois le même vertex et le tri
   *       --> Il faudrait copier dans un tableau 1D les pattions les une à la suite des autres !
   *           On change l'interface ?
   *       --> En concatemant avant le MPI est plus simple car on a pas à creer un
   *           tableau d'index sur chaque proc de deplacement entre les différentes partitions
   *       Les 2 sont pratiques car on pourrait imaginer avoir des traitement locaux differnts non ?
   */

  for(int ipart = 0; ipart < pdn->n_part; ipart++){

    printf("[%d] ------------------------------------------ %d \n", ipart, pdn->requested_data_n[0]);

    int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];
    int *_part_neighbor_desc = pdn->neighbor_desc[ipart];

    // printf("[%i] - n_entity:: %d\n", ipart, n_entity[ipart]);

    pdn->order       [ipart] = (int *) malloc( _part_neighbor_idx[n_entity[ipart]] * sizeof(int *));
    pdn->order_unique[ipart] = (int *) malloc( _part_neighbor_idx[n_entity[ipart]] * sizeof(int *));

    // Sort
    PDM_order_lnum_s(pdn->neighbor_desc[ipart],
                     3,
                     pdn->order[ipart],
                     _part_neighbor_idx[n_entity[ipart]]);

    // Compute the unique idx from sort
    compute_unique_idx(pdn->order[ipart],
                       pdn->order_unique[ipart],
                       pdn->neighbor_desc[ipart],
                       _part_neighbor_idx[n_entity[ipart]]);

    // Il faut connaitre le nombre d'occurence une fois trié --> Taille du buffer d'envoie
    // Mais par proc / part
    // Il faut pas parcourir ce pacquet mais le unique !
    int lastidx = -1;
    for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[ipart]]; i_entity++){
      // Go to sorted order
      int u_entity = pdn->order_unique[ipart][i_entity];
      int s_entity = pdn->order[ipart][i_entity];
      // printf("[%d] - order:: %d | unique:: %d | lastidx:: %d \n", i_entity, s_entity, u_entity, lastidx);
      if(lastidx != u_entity){
        int opp_proc = _part_neighbor_desc[3*s_entity  ];
        pdn->requested_data_n[opp_proc]++;
        lastidx = u_entity;
      }
    }
  }

  for (int i = 0; i < nRank; i++) {
    pdn->requested_data_idx[i+1] = pdn->requested_data_idx[i] +
                                   pdn->requested_data_n[i];
  }

  /*
   * Compute size and reset counter
   */
  int s_requested_data = pdn->requested_data_idx[nRank-1] + pdn->requested_data_n[nRank-1];

  printf("s_requested_data:: %d \n", s_requested_data);
  for (int i = 0; i < nRank; i++) {
    pdn->requested_data_n[i] = 0;
  }

  int *requested_data = malloc (sizeof(int) *  2 * s_requested_data); // Store ipart/ientity

  for(int ipart = 0; ipart < pdn->n_part; ipart++){

    int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];
    int *_part_neighbor_desc = pdn->neighbor_desc[ipart];
    pdn->distributed_part_n  [ipart] = 0;
    pdn->distributed_part_idx[ipart] = 0;

    int lastidx = -1;
    for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[ipart]]; i_entity++){
      // Go to sorted order
      int u_entity = pdn->order_unique[ipart][i_entity];
      int s_entity = pdn->order[ipart][i_entity];
      printf("[%d] - order:: %d | unique:: %d | lastidx:: %d \n", i_entity, s_entity, u_entity, lastidx);
      if(lastidx != u_entity){

        int opp_proc = _part_neighbor_desc[3*s_entity  ];
        int opp_part = _part_neighbor_desc[3*s_entity+1];
        int opp_etty = _part_neighbor_desc[3*s_entity+2];

        int idx = pdn->requested_data_idx[opp_proc] + pdn->requested_data_n[opp_proc]++;

        requested_data[2*idx  ] = opp_part;
        requested_data[2*idx+1] = opp_etty;

        pdn->distributed_part_n[ipart]++;

        lastidx = u_entity;
      }
    }
  }

  // Each pacquet have 2 value so we multiply by 2 temporary
  for (int i = 0; i < nRank; i++) {
    pdn->requested_data_n[i] = 2*pdn->requested_data_n[i];
  }

  for (int i = 0; i < pdn->n_part; i++) {
    pdn->distributed_part_idx[i+1] = pdn->distributed_part_n[i] +
                                     pdn->distributed_part_idx[i];
  }

  if(1 == 1){
    log_trace("PDM_distant_neighbor_create::requested_data :: --> ");
    for(int i = 0; i < s_requested_data; ++i){
      log_trace("[%d/%d] ", requested_data[2*i], requested_data[2*i+1]);
    }
    log_trace("\n");
  }

  if(1 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_part_n :: --> ");
    for(int i = 0; i < pdn->n_part; ++i){
      log_trace("%d ", pdn->distributed_part_n[i]);
    }
    log_trace("\n");
    log_trace("PDM_distant_neighbor_create::distributed_part_idx :: --> ");
    for(int i = 0; i < pdn->n_part+1; ++i){
      log_trace("%d ", pdn->distributed_part_idx[i]);
    }
    log_trace("\n");
  }

  /*
   * Exchange the requested data
   */
  pdn->distributed_data_n = malloc (sizeof(int) * nRank);

  PDM_MPI_Alltoall (pdn->requested_data_n,   1, PDM_MPI_INT,
                    pdn->distributed_data_n, 1, PDM_MPI_INT,
                    pdn->comm);

  if(1 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_data_n :: --> ");
    for(int i = 0; i < nRank; ++i){
      log_trace("%d ",  pdn->distributed_data_n[i]);
    }
    log_trace("\n");
  }


  pdn->distributed_data_idx = malloc (sizeof(int) * (nRank + 1));
  pdn->distributed_data_idx[0] = 0;

  for (int i = 0; i < nRank; i++) {
    pdn->distributed_data_idx[i+1] = pdn->distributed_data_n[i] +
                                     pdn->distributed_data_idx[i];
  }

  if(1 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_data_idx :: --> ");
    for(int i = 0; i < nRank+1; ++i){
      log_trace("%d ",  pdn->distributed_data_idx[i]);
    }
    log_trace("\n");
  }


  pdn->distributed_data = malloc (sizeof(int) * 2 * pdn->distributed_data_idx[nRank]);
  // pdn->distributed_data = malloc (sizeof(int) * pdn->distributed_data_idx[nRank]);

  PDM_MPI_Alltoallv (requested_data,
                     pdn->requested_data_n,
                     pdn->requested_data_idx,
                     PDM_MPI_INT,
                     pdn->distributed_data,
                     pdn->distributed_data_n,
                     pdn->distributed_data_idx,
                     PDM_MPI_INT,
                     pdn->comm);

  if(1 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_data :: --> ");
    for(int i = 0; i < pdn->distributed_data_idx[nRank]/2; ++i){
    // for(int i = 0; i < pdn->distributed_data_idx[nRank]; ++i){
      log_trace("[%d/%d] ", pdn->distributed_data[2*i], pdn->distributed_data[2*i+1]);
    }
    log_trace("\n");
  }

  /*
   * Store in structure the ordering in the recv buffer
   */

  /*
   * Re setup size of each member of struct
   */
  pdn->requested_data_idx[0] = 0;
  pdn->distributed_data_idx[0] = 0;
  for (int i = 0; i < nRank; i++) {
    pdn->requested_data_n[i]   = pdn->requested_data_n  [i]/2;
    pdn->distributed_data_n[i] = pdn->distributed_data_n[i]/2;

    pdn->requested_data_idx  [i+1] = pdn->requested_data_idx[i] + pdn->requested_data_n    [i];
    pdn->distributed_data_idx[i+1] = pdn->distributed_data_n[i] + pdn->distributed_data_idx[i];
  }

  if(1 == 1){
    log_trace("Re-Setup --- ");
    log_trace("PDM_distant_neighbor_create::distributed_data :: --> ");
    // for(int i = 0; i < pdn->distributed_data_idx[nRank]/2; ++i){
    for(int i = 0; i < pdn->distributed_data_idx[nRank]; ++i){
      log_trace("[%d/%d] ", pdn->distributed_data[2*i], pdn->distributed_data[2*i+1]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::distributed_data_idx :: --> ");
    for(int i = 0; i < nRank+1; ++i){
      log_trace("%d ",  pdn->distributed_data_idx[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::distributed_data_n :: --> ");
    for(int i = 0; i < nRank; ++i){
      log_trace("%d ",  pdn->distributed_data_n[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::requested_data_idx :: --> ");
    for(int i = 0; i < nRank+1; ++i){
      log_trace("%d ",  pdn->requested_data_idx[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::requested_data_n :: --> ");
    for(int i = 0; i < nRank; ++i){
      log_trace("%d ",  pdn->requested_data_n[i]);
    }
    log_trace("\n");
  }

  /*
   * Free
   */
  free(requested_data);

  assert(fclose(fp) == 0);

  return id;
}




/**
 * \brief Exchange data between \ref _distant_neighbor_t structure
 * \param [in]   id          identifier of internal structre
 *  NB : On va commencer par des entiers en stride constantes
 *
 */
void
PDM_distant_neighbor_exch
(
 const int      id,
 size_t         s_data,
 PDM_stride_t   t_stride,
 int            cst_stride,
 int          **send_entity_stride,
 int          **send_entity_data,
 int         ***recv_entity_stride,
 int         ***recv_entity_data
)
{
  printf(" PDM_distant_neighbor_exchange \n");
  _distant_neighbor_t *pdn = _get_from_id (id);

  int iRank;
  int nRank;
  PDM_MPI_Comm_rank(pdn->comm, &iRank);
  PDM_MPI_Comm_size(pdn->comm, &nRank);

  char filename[50];
  sprintf(filename, "pdm_logging_exch_%d.log", iRank);
  FILE* fp = fopen(filename, "w");
  log_set_fp(fp);
  log_set_quiet(1);
  log_set_level(-1);
  // log_trace("PDM_distant_neighbor_create::Su

  // int s_distributed_data = pdn->distributed_data_idx[nRank]/2;
  int s_distributed_data = pdn->distributed_data_idx[nRank];
  int s_requested_data   = pdn->requested_data_idx[nRank];

  size_t *i_sendBuffer = (size_t *) malloc (sizeof(size_t) * nRank);
  size_t *i_recvBuffer = (size_t *) malloc (sizeof(size_t) * nRank);
  int *n_sendBuffer = (int *) malloc (sizeof(int) * nRank);
  int *n_recvBuffer = (int *) malloc (sizeof(int) * nRank);

  size_t s_sendBuffer = 0;
  size_t s_recvBuffer = 0;

  for (int i = 0; i < nRank; i++) {
    n_sendBuffer[i] = 0;
    n_recvBuffer[i] = 0;
    i_sendBuffer[i] = 0;
    i_recvBuffer[i] = 0;
  }

  // unsigned char *sendBuffer = NULL;
  // unsigned char *recvBuffer = NULL;
  int *sendBuffer = NULL;
  int *recvBuffer = NULL;

  // if(t_stride !=  PDM_STRIDE_CST) {
  //   PDM_error(__FILE__, __LINE__, 0,"PDM_distant_neighbor_exch : STRIDE_CST is only availble \n");
  //   abort ();
  // }

  /*
   * On doit echanger de la même manière que les requested du create mais en multipliant par la stride
   */
  int *recvStride = NULL;
  if (t_stride == PDM_STRIDE_VAR) {

    int *sendStride = (int *) malloc (sizeof(int) * s_distributed_data);
    recvStride = (int *) malloc (sizeof(int) * s_requested_data);

    /*
     * Prepare send stride
     */
    int idx_send = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int ipart = pdn->distributed_data[2*i  ];
      int ienty = pdn->distributed_data[2*i+1];
      printf("sendStride[%d/%d] --> [%d,%d] -> %d \n", idx_send, s_distributed_data, ipart, ienty, send_entity_stride[ipart][ienty]);
      log_trace("sendStride[%d/%d] --> [%d,%d] -> %d \n", idx_send, s_distributed_data, ipart, ienty, send_entity_stride[ipart][ienty]);
      sendStride[idx_send++] = send_entity_stride[ipart][ienty];
    }

    printf("Exch stride %d\n", s_requested_data);
    PDM_MPI_Alltoallv (sendStride,
                       pdn->distributed_data_n,
                       pdn->distributed_data_idx,
                       PDM_MPI_INT,
                       recvStride,
                       pdn->requested_data_n,
                       pdn->requested_data_idx,
                       PDM_MPI_INT,
                       pdn->comm);
    printf("Exch stride end \n");

    /*
     * Fill the recv stride for all parts / entity
     */
    *recv_entity_stride = malloc( pdn->n_part * sizeof(int *) );
    int **_recv_entity_stride = (*(int ***) recv_entity_stride);

    for(int ipart = 0; ipart < pdn->n_part; ipart++){
      int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];
      _recv_entity_stride[ipart] = (int *) malloc( _part_neighbor_idx[pdn->n_entity[ipart]] * sizeof(int *));

      for(int i_entity = 0; i_entity < _part_neighbor_idx[pdn->n_entity[ipart]]; i_entity++){
        int idx = pdn->distributed_part_idx[ipart] + pdn->order_unique[ipart][i_entity];
        printf("recv strid ::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[pdn->n_entity[ipart]], ipart, i_entity, recvStride[idx]);
        log_trace("recv strid::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[pdn->n_entity[ipart]], ipart, i_entity, recvStride[idx]);
        _recv_entity_stride[ipart][i_entity] = recvStride[idx];
      }
    }

    /*
     * Build buffer
     */
    for (int i = 0; i < nRank; i++) {

      /*
       * Setup send data
       */

      int iBeg = pdn->distributed_data_idx[i];
      int iEnd = pdn->distributed_data_idx[i] + pdn->distributed_data_n[i];

      n_sendBuffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++)  {
        n_sendBuffer[i] += sendStride[k];
      }

      // n_sendBuffer[i] *= (int) s_data;

      if (i > 0) {
        i_sendBuffer[i] = i_sendBuffer[i-1] + n_sendBuffer[i-1];
      }
      else {
        i_sendBuffer[i] = 0;
      }

      /*
       * Setup recv data
       */
      iBeg = pdn->requested_data_idx[i];
      iEnd = pdn->requested_data_idx[i] + pdn->requested_data_n[i];

      n_recvBuffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++) {
        n_recvBuffer[i] += recvStride[k];
      }

      // n_recvBuffer[i] *= (int) s_data;

      if (i > 0) {
        i_recvBuffer[i] = i_recvBuffer[i-1] + n_recvBuffer[i-1];
      }
      else {
        i_recvBuffer[i] = 0;
      }
    } /* End iRank loop */

    s_sendBuffer = i_sendBuffer[nRank-1] + n_sendBuffer[nRank-1];
    s_recvBuffer = i_recvBuffer[nRank-1] + n_recvBuffer[nRank-1];

    // sendBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_sendBuffer);
    // recvBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recvBuffer);
    log_trace("PDM_distant_neighbor_exch::s_sendBuffer :: %d --> \n ", s_sendBuffer);
    log_trace("PDM_distant_neighbor_exch::s_recvBuffer :: %d --> \n ", s_recvBuffer);

    sendBuffer = (int *) malloc(sizeof(int) * s_sendBuffer);
    recvBuffer = (int *) malloc(sizeof(int) * s_recvBuffer);

    int *sendStride_idx = (int *) malloc(sizeof(int) * (s_distributed_data+1));
    sendStride_idx[0] = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      sendStride_idx[i+1] = sendStride_idx[i] + sendStride[i];
    }

    /*
     * Compute stride for each part (in the order of send buffer )
     */
    int** stride_idx = (int**) malloc( pdn->n_part * sizeof(int **));
    for(int ipart = 0; ipart < pdn->n_part; ipart++){
      stride_idx[ipart] = (int*) malloc( (pdn->n_entity[ipart] + 1)  * sizeof(int*));
      stride_idx[ipart][0] = 0;
      for(int i_entity = 0; i_entity < pdn->n_entity[ipart]; i_entity++){
        stride_idx[ipart][i_entity+1] = stride_idx[ipart][i_entity] + send_entity_stride[ipart][i_entity];
      }
    }

    /*
     * Fill the sendBuffer
     */
    log_trace("PDM_distant_neighbor_exch::sendBuffer :: --> \n ");
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int ipart = pdn->distributed_data[2*i  ];
      int ienty = pdn->distributed_data[2*i+1];
      for(int idata = 0; idata < send_entity_stride[ipart][ienty]; idata++){
        int idxdata = stride_idx[ipart][ienty] + idata;
        printf("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, ipart, ienty, send_entity_data[ipart][idxdata]);
        log_trace("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, ipart, ienty, send_entity_data[ipart][idxdata]);
        sendBuffer[idx1++] = send_entity_data[ipart][idxdata];
      }
    }
    log_trace("PDM_distant_neighbor_exch::sendBuffer END \n ");

    for(int ipart = 0; ipart < pdn->n_part; ipart++){
      free (stride_idx[ipart]);
    }
    free (stride_idx);
    free (recvStride);
    free (sendStride);
    free (sendStride_idx);

  } else if (t_stride == PDM_STRIDE_CST) {

    int s_block_unit = cst_stride * (int) s_data;

    log_trace("PDM_distant_neighbor_exch::requested_data :: --> \n ");

    for (int i = 0; i < nRank; i++) {

      // i_sendBuffer[i] = pdn->distributed_data_idx[i] * cst_stride * (int) s_data;
      // i_recvBuffer[i] = pdn->requested_data_idx[i] * cst_stride * (int) s_data;

      // n_sendBuffer[i] = pdn->distributed_data_n[i] * cst_stride * (int) s_data;
      // n_recvBuffer[i] = pdn->requested_data_n[i] * cst_stride * (int) s_data;

      i_sendBuffer[i] = pdn->distributed_data_idx[i] * cst_stride;
      i_recvBuffer[i] = pdn->requested_data_idx[i] * cst_stride;

      n_sendBuffer[i] = pdn->distributed_data_n[i] * cst_stride;
      n_recvBuffer[i] = pdn->requested_data_n[i] * cst_stride;

      log_trace("[%d] send::[%d/%d] | recv::[%d/%d] \n ", i, i_sendBuffer[i], n_sendBuffer[i],
                                                             i_recvBuffer[i], n_recvBuffer[i]);

    }

    s_sendBuffer = i_sendBuffer[nRank-1] + n_sendBuffer[nRank-1];
    s_recvBuffer = i_recvBuffer[nRank-1] + n_recvBuffer[nRank-1];

    // sendBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_sendBuffer);
    // recvBuffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recvBuffer);
    log_trace("PDM_distant_neighbor_exch::s_sendBuffer :: %d --> \n ", s_sendBuffer);
    log_trace("PDM_distant_neighbor_exch::s_recvBuffer :: %d --> \n ", s_recvBuffer);
    sendBuffer = (int *) malloc(sizeof(int) * s_sendBuffer);
    recvBuffer = (int *) malloc(sizeof(int) * s_recvBuffer);


    log_trace("PDM_distant_neighbor_exch::sendBuffer :: --> \n ");
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int ipart = pdn->distributed_data[2*i  ];
      int ienty = pdn->distributed_data[2*i+1];
      printf("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, ipart, ienty, send_entity_data[ipart][ienty]);
      log_trace("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, ipart, ienty, send_entity_data[ipart][ienty]);
      sendBuffer[idx1++] = send_entity_data[ipart][ienty];
    }
    log_trace("PDM_distant_neighbor_exch::sendBuffer END \n ");
  }


  /*
   * Exchnage
   */
  // PDM_MPI_Alltoallv_l(sendBuffer,
  //                     n_sendBuffer,
  //                     i_sendBuffer,
  //                     PDM_MPI_BYTE,
  //                     recvBuffer,
  //                     n_recvBuffer,
  //                     i_recvBuffer,
  //                     PDM_MPI_BYTE,
  //                     pdn->comm);
  PDM_MPI_Alltoallv_l(sendBuffer,
                      n_sendBuffer,
                      i_sendBuffer,
                      PDM_MPI_INT,
                      recvBuffer,
                      n_recvBuffer,
                      i_recvBuffer,
                      PDM_MPI_INT,
                      pdn->comm);


  free(sendBuffer);
  free(n_sendBuffer);
  free(n_recvBuffer);
  free(i_recvBuffer);

  /*
   * Une seule valeur est echangé mais plusieurs occurence peuvent exister donc on passe du buffer MPI
   * au donné sans le sort/unique
   */
  *recv_entity_data = malloc( pdn->n_part * sizeof(int *) );
  int **_recv_entity_data = (*(int ***) recv_entity_data);

  if (t_stride == PDM_STRIDE_VAR) {

    /*
     * Compute the recv stride index for each entity
     */
    int** _recv_entity_stride = (*(int ***) recv_entity_stride);
    int** recv_stride_idx = (int**) malloc( pdn->n_part * sizeof(int **));
    for(int ipart = 0; ipart < pdn->n_part; ipart++){
      recv_stride_idx[ipart] = (int*) malloc( (pdn->n_entity[ipart] + 1)  * sizeof(int*));
      recv_stride_idx[ipart][0] = 0;
      for(int i_entity = 0; i_entity < pdn->n_entity[ipart]; i_entity++){
        recv_stride_idx[ipart][i_entity+1] = recv_stride_idx[ipart][i_entity] + _recv_entity_stride[ipart][i_entity];
      }
    }

    /*
     * Copy buffer into the recv_data
     */
    for(int ipart = 0; ipart < pdn->n_part; ipart++){
      int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];

      int recv_part_size = _part_neighbor_idx[pdn->n_entity[ipart]] * recv_stride_idx[ipart][pdn->n_entity[ipart]];

      log_trace("PDM_distant_neighbor_exch::recv_part_size :: --> %d \n ", recv_part_size);
      _recv_entity_data[ipart] = (int *) malloc( recv_part_size * sizeof(int *));

      for(int i_entity = 0; i_entity < _part_neighbor_idx[pdn->n_entity[ipart]]; i_entity++){

        int idx = pdn->distributed_part_idx[ipart] + pdn->order_unique[ipart][i_entity];

        for(int idata = 0; idata < _recv_entity_stride[ipart][i_entity]; idata++){
          int idxdata = recv_stride_idx[ipart][i_entity] + idata;
          // printf("recv ::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[pdn->n_entity[ipart]], ipart, i_entity, recvBuffer[idx]);
          // log_trace("recv ::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[pdn->n_entity[ipart]], ipart, i_entity, recvBuffer[idx]);
          _recv_entity_data[ipart][idxdata] = recvBuffer[idx+idata];
        }
      } /* End ientity */
    } /* End part */




    /*
     * Free
     */

    for(int ipart = 0; ipart < pdn->n_part; ipart++){
     free(recv_stride_idx[ipart]);
    }
    free(recv_stride_idx);


  } else if (t_stride == PDM_STRIDE_CST) {

    // Shift is not good because the buffer contains only one occurence of each elements !!!
    for(int ipart = 0; ipart < pdn->n_part; ipart++){
      int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];
      _recv_entity_data[ipart] = (int *) malloc( _part_neighbor_idx[pdn->n_entity[ipart]] * sizeof(int *));

      log_trace("PDM_distant_neighbor_exch::recvBuffer :: --> \n ");
      for(int i_entity = 0; i_entity < _part_neighbor_idx[pdn->n_entity[ipart]]; i_entity++){
        // int idx = i_sendBuffer[] + order_unique[i_entity]
        int idx = pdn->distributed_part_idx[ipart] + pdn->order_unique[ipart][i_entity];
        printf("recv ::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[pdn->n_entity[ipart]], ipart, i_entity, recvBuffer[idx]);
        log_trace("recv ::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[pdn->n_entity[ipart]], ipart, i_entity, recvBuffer[idx]);
        _recv_entity_data[ipart][i_entity] = recvBuffer[idx];
      }
    }
  }


  free(i_sendBuffer);
  free(recvBuffer);

}


/**
 *
 * \brief Free an distant negihtbor structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_distant_neighbor_free
(
 const int          id
)
{
  _distant_neighbor_t *pdn = _get_from_id (id);

  for(int ipart = 0; ipart < pdn->n_part; ipart++){
    free(pdn->order[ipart]);
    free(pdn->order_unique[ipart]);
  }
  free(pdn->order);
  free(pdn->order_unique);

  free(pdn->requested_data_n);
  free(pdn->requested_data_idx);
  free(pdn->distributed_part_idx);
  free(pdn->distributed_part_n);

  free(pdn->distributed_data);
  free(pdn->distributed_data_n);
  free(pdn->distributed_data_idx);

  free (pdn);

  PDM_Handles_handle_free (_pdns, id, PDM_FALSE);

  const int n_ppm = PDM_Handles_n_get (_pdns);

  if (n_ppm == 0) {
    _pdns = PDM_Handles_free (_pdns);
  }

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
