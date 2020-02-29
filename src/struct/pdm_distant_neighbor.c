/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_distant_neighbor.h"
#include "pdm_printf.h"
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
  PDM_MPI_Comm comm;             /*!< MPI communicator */
  int          n_part;           /*!< Number of partitions */
  const int   *n_entity;         /*!< Number of entities for each partition */
  int        **neighbor_idx;     /*!< Indexes of candidate for each current part point
                                 *   (size = number of entities in the current part + 1) */
  int        **neighbor_desc;    /*!< Candidates description (process,
                                 *                           part in the process,
                                 *                           entitiy number in the part) */
  int          n_rank_exch;      /*!< Number of rank concerned by exchange */
  int*         exch_rank;
  int**        order;
  int**        order_unique;
  int         *requested_data_n; /*!< Numer of requested data for each process index
                                    * (size : s_comm) */
  int         *requested_data_idx;/*!< Requested data for each process index
                                    * (size : s_comm) */
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

  pdn->n_part             = n_part;
  pdn->n_entity           = n_entity;
  pdn->neighbor_idx       = neighbor_idx;
  pdn->neighbor_desc      = neighbor_desc;
  pdn->n_rank_exch        = -1;
  pdn->exch_rank          = NULL;
  pdn->order              = (int **) malloc( pdn->n_part * sizeof(int **));
  pdn->order_unique       = (int **) malloc( pdn->n_part * sizeof(int **));
  pdn->requested_data_n   = (int * ) malloc( (nRank    ) * sizeof(int * ));
  pdn->requested_data_idx = (int * ) malloc( (nRank + 1) * sizeof(int * ));


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

    int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];
    int *_part_neighbor_desc = pdn->neighbor_desc[ipart];

    // For now we order/unique only the part !
    printf("[%i] - n_entity:: %d\n", ipart, n_entity[ipart]);

    pdn->order       [ipart] = (int *) malloc( n_entity[ipart] * sizeof(int *));
    pdn->order_unique[ipart] = (int *) malloc( n_entity[ipart] * sizeof(int *));

    // Sort
    PDM_order_lnum_s(pdn->neighbor_desc[ipart], 3, pdn->order[ipart], n_entity[ipart]);

    // Compute the unique idx from sort
    compute_unique_idx(pdn->order[ipart],
                       pdn->order_unique[ipart],
                       pdn->neighbor_desc[ipart],
                       pdn->n_entity[ipart]);

    // Il faut connaitre le nombre d'occurence une fois trié --> Taille du buffer d'envoie
    // Mais par proc / part
    // Il faut pas parcourir ce pacquet mais le unique !
    for(int i_entity = 0; i_entity < n_entity[ipart]; i_entity++){
      for(int j = _part_neighbor_idx[i_entity]; j < _part_neighbor_idx[i_entity+1]; j++){
        int opp_proc = _part_neighbor_desc[3*j  ];
        // int opp_part = _part_neighbor_desc[3*j+1];
        // int opp_etty = _part_neighbor_desc[3*j+2];
        pdn->requested_data_n[opp_proc]++;
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

    for(int i_entity = 0; i_entity < n_entity[ipart]; i_entity++){
      for(int j = _part_neighbor_idx[i_entity]; j < _part_neighbor_idx[i_entity+1]; j++){
        int opp_proc = _part_neighbor_desc[3*j  ];
        int opp_part = _part_neighbor_desc[3*j+1];
        int opp_etty = _part_neighbor_desc[3*j+2];

        int idx = pdn->requested_data_idx[opp_proc] + pdn->requested_data_n[opp_proc]++;

        requested_data[2*idx  ] = opp_part;
        requested_data[2*idx+1] = opp_etty;


      }
    }
  }


  /*
   * Store in structure the ordering in the recv buffer
   */


  /*
   * Free
   */
  free(requested_data);

  return id;
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

  if (pdn->exch_rank != NULL)
    free (pdn->exch_rank);

  free (pdn);

  PDM_Handles_handle_free (_pdns, id, PDM_FALSE);

  const int n_ppm = PDM_Handles_n_get (_pdns);


  if (n_ppm == 0) {
    _pdns = PDM_Handles_free (_pdns);
  }

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
 int          **recv_entity_stride,
 int          **recv_entity_data
)
{
  printf(" PDM_distant_neighbor_exchange \n");
  _distant_neighbor_t *pdn = _get_from_id (id);

  if(t_stride !=  PDM_STRIDE_CST) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_distant_neighbor_exch : STRIDE_CST is only availble \n");
    abort ();
  }

}

// void
// PDM_distant_neighbor_irecv()
// PDM_distant_neighbor_iisend()
// PDM_distant_neighbor_wait()

// void
// PDM_distant_neighbor_alltoall()
// {
// }


  // /*
  //  * Setup exchange protocol
  //  *   I/ Compute the total size of connected elmts by each proc with others
  //  */
  // int *offered_elmts_rank_idx = (int *) malloc ((nRank + 1) * sizeof(int));
  // for(int i = 0; i < nRank+1; i++){
  //   offered_elmts_rank_idx[i] = 0;
  // }

  // for(int ipart = 0; ipart < pdn->n_part; ipart++){
  //   int *_part_neighbor_idx  = pdn->neighbor_idx[ipart];
  //   int *_part_neighbor_desc = pdn->neighbor_desc[ipart];

  //   for(int i_entity = 0; i_entity < n_entity[ipart]; i_entity++){
  //     for(int j = _part_neighbor_idx[i_entity]; j < _part_neighbor_idx[i_entity+1]; j++){
  //       int opp_proc = _part_neighbor_desc[3*j  ];
  //       // int opp_part = _part_neighbor_desc[3*j+1];
  //       // int opp_etty = _part_neighbor_desc[3*j+2];
  //       offered_elmts_rank_idx[opp_proc+1]++;
  //     }
  //   }
  // }

  // /*
  //  * Panic verbose
  //  */
  // if(0 == 0){
  //   printf("offered_elmts_rank_idx:: ");
  //   for(int i = 0; i < nRank+1; i++){
  //     printf("%d ", offered_elmts_rank_idx[i]);
  //   }
  //   printf("\n");
  // }


  //  * Build list of exchanged ranks

  // pdn->n_rank_exch = 0;
  // for(int i = 0; i < nRank+1; i++){
  //   if((i != iRank) && (offered_elmts_rank_idx[i+1] > 0)) {
  //     pdn->n_rank_exch += 1;
  //   }
  //   offered_elmts_rank_idx[i+1] = offered_elmts_rank_idx[i+1] + offered_elmts_rank_idx[i];
  // }

  // /*
  //  * Allocate and reset
  //  */
  // pdn->exch_rank = (int *) malloc( pdn->n_rank_exch * sizeof(int) );
  // pdn->n_rank_exch = 0;

  // for(int i = 0; i < nRank+1; i++){
  //   if((i != iRank) && (offered_elmts_rank_idx[i+1] > offered_elmts_rank_idx[i])) {
  //     pdn->exch_rank[pdn->n_rank_exch++] = i; // Communication graph
  //   }
  // }

  // /*
  //  * Panic verbose
  //  */
  // if(0 == 0){
  //   printf("offered_elmts_rank_idx shift:: ");
  //   for(int i = 0; i < nRank+1; i++){
  //     printf("%d ", offered_elmts_rank_idx[i]);
  //   }
  //   printf("\n");
  // }

  // /*
  //  * Exchange number of offered elements for each connected element
  //  * to connected ranks
  //  */
  // int *n_recv_elmt    = (int *) malloc (nRank * sizeof(int));
  // for (int i = 0; i < nRank; i++) {
  //   n_recv_elmt[i] = 0;
  // }

  // /*
  //  * Free
  //  */
  // free(offered_elmts_rank_idx);

#ifdef __cplusplus
}
#endif /* __cplusplus */
