/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_distant_neighbor.h"
#include "pdm_logging.h"
#include "pdm_part_extension.h"
#include "pdm_part_extension_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
int
_setup_neighbor
(
 int           i_proc_cur,
 int           i_part_cur,
 int           ntot_part,
 int          *entity_part_bound_part_idx,
 int          *entity_part_bound,
 int          *entity_cell,
 PDM_g_num_t  *cell_ln_to_gn,
 int          *cell_list,
 int          *graph_comm_cell,
 PDM_g_num_t **cell_border_ln_to_gn,
 int         **neighbor_idx,
 int         **neighbor_desc
)
{
  *neighbor_idx         = malloc(     ( entity_part_bound_part_idx[ntot_part] + 1 ) * sizeof(int        ));
  *neighbor_desc        = malloc(  3 *  entity_part_bound_part_idx[ntot_part]       * sizeof(int        )); // i_proc, i_part, i_entity
  *cell_border_ln_to_gn = malloc(       entity_part_bound_part_idx[ntot_part]       * sizeof(PDM_g_num_t));

  int* _neighbor_idx  = *neighbor_idx;
  int* _neighbor_desc = *neighbor_desc;

  PDM_g_num_t* _cell_border_ln_to_gn = *cell_border_ln_to_gn;

  int n_entity_bound = 0;
  for(int i_part = 0; i_part < ntot_part; ++i_part ){
    _neighbor_idx[0] = 0;
    int idx_opp = 0;
    for(int idx_entity = entity_part_bound_part_idx[i_part]; idx_entity < entity_part_bound_part_idx[i_part+1]; ++idx_entity) {

      int i_entity     = entity_part_bound[4*idx_entity  ]-1;
      int i_proc_opp   = entity_part_bound[4*idx_entity+1];
      int i_part_opp   = entity_part_bound[4*idx_entity+2]-1;
      // int i_entity_opp = entity_part_bound[4*idx_entity+3];

      _neighbor_idx [idx_entity+1  ] = _neighbor_idx[idx_entity] + 1;
      _neighbor_desc[3*idx_entity  ] = i_proc_opp;
      _neighbor_desc[3*idx_entity+1] = i_part_opp;
      // _neighbor_desc[3*idx_entity+2] = i_entity_opp;
      _neighbor_desc[3*idx_entity+2] = idx_opp; // Car symétrique

      int icell1 = entity_cell[2*i_entity];

      _cell_border_ln_to_gn[idx_entity] = cell_ln_to_gn[icell1-1];

      cell_list[idx_entity] = icell1;
      graph_comm_cell[3*idx_entity  ] = idx_opp++;
      graph_comm_cell[3*idx_entity+1] = i_proc_cur;
      graph_comm_cell[3*idx_entity+2] = i_part_cur;

      n_entity_bound++;
    }
  }

  return n_entity_bound;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a part extension structure
 *
 * \param [in]   comm           Communicator
 *
 */

PDM_part_extension_t*
PDM_part_extension_create
(
 const int              n_domain,
 const int             *n_part,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
)
{
 PDM_part_extension_t *part_ext = (PDM_part_extension_t *) malloc(sizeof(PDM_part_extension_t));

  part_ext->n_domain = n_domain;
  part_ext->n_part   = n_part;
  part_ext->comm     = comm;
  part_ext->owner    = owner;

  part_ext->parts = malloc(n_domain * sizeof(_part_t *));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    part_ext->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
  }

  return part_ext;
}

void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn
)
{
  part_ext->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  part_ext->parts[i_domain][i_part].cell_face     = cell_face;
  part_ext->parts[i_domain][i_part].face_cell     = face_cell;

  part_ext->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  part_ext->parts[i_domain][i_part].face_edge     = face_edge;

  part_ext->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  part_ext->parts[i_domain][i_part].face_vtx      = face_vtx;

  part_ext->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  part_ext->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  part_ext->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  part_ext->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  part_ext->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = face_part_bound_proc_idx;
  part_ext->parts[i_domain][i_part].face_part_bound_part_idx = face_part_bound_part_idx;
  part_ext->parts[i_domain][i_part].face_part_bound          = face_part_bound;

  part_ext->parts[i_domain][i_part].face_bound_idx      = face_bound_idx;
  part_ext->parts[i_domain][i_part].face_bound          = face_bound;
  part_ext->parts[i_domain][i_part].face_join_idx       = face_join_idx;
  part_ext->parts[i_domain][i_part].face_join           = face_join;
  part_ext->parts[i_domain][i_part].face_bound_ln_to_gn = face_group_ln_to_gn;
}


/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext,
  int                   depth,
  PDM_extend_type_t     extend_type
)
{

  /*
   *  A prevoir : reconstruire les "ghost" sans toute la topologie
   *              par exemple pour l'algébre linéaire
   */

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  assert(part_ext != NULL);
  printf(" PDM_part_extension_compute : depth = %i | extend_type = %i \n", depth, extend_type);

  assert(extend_type == PDM_EXTEND_FROM_FACE);

  int*** neighbor_idx      = (int ***) malloc( part_ext->n_domain * sizeof(int ** ) );
  int*** neighbor_desc     = (int ***) malloc( part_ext->n_domain * sizeof(int ** ) );
  int**  n_face_part_bound = (int  **) malloc( part_ext->n_domain * sizeof(int  * ) );

  int ***cell_list       = (int ***) malloc( part_ext->n_domain * sizeof(int **));
  int ***graph_comm_cell = (int ***) malloc( part_ext->n_domain * sizeof(int **));


  // Begin with exchange by the connectivity the cell opposite number
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    int n_part_total = -1;
    int n_part_loc = part_ext->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &n_part_total, 1, PDM_MPI_INT, PDM_MPI_SUM, part_ext->comm);

    neighbor_idx     [i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *) );
    neighbor_desc    [i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *) );
    n_face_part_bound[i_domain] = (int  *) malloc( part_ext->n_part[i_domain] * sizeof(int  ) );

    printf(" n_part_total = %i \n", n_part_total);
    PDM_g_num_t** cell_border_ln_to_gn = malloc( part_ext->n_part [i_domain] * sizeof(PDM_g_num_t * ));

    cell_list      [i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
    graph_comm_cell[i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_part_face_bound_tot = part_ext->parts[i_domain][i_part].face_part_bound_part_idx[n_part_total];
      cell_list      [i_domain][i_part] = (int *) malloc(     n_part_face_bound_tot * sizeof(int));
      graph_comm_cell[i_domain][i_part] = (int *) malloc( 3 * n_part_face_bound_tot * sizeof(int));

      // Pour chaque partition on hook le graph de comm
      n_face_part_bound[i_domain][i_part] = _setup_neighbor(i_rank,
                                                            i_part,
                                                            n_part_total,
                                                            part_ext->parts[i_domain][i_part].face_part_bound_part_idx,
                                                            part_ext->parts[i_domain][i_part].face_part_bound,
                                                            part_ext->parts[i_domain][i_part].face_cell,
                                                            part_ext->parts[i_domain][i_part].cell_ln_to_gn,
                                                            cell_list[i_domain][i_part],
                                                            graph_comm_cell[i_domain][i_part],
                                                            &cell_border_ln_to_gn[i_part],
                                                            &neighbor_idx [i_domain][i_part],
                                                            &neighbor_desc[i_domain][i_part]);

      printf(" n_face_part_bound[%i][%i] = %i \n", i_domain, i_part, n_face_part_bound[i_domain][i_part]);
    }

    /*
     * Create distant neighbor
     */
    PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                             part_ext->n_part [i_domain],
                                                             n_face_part_bound[i_domain],
                                                             neighbor_idx     [i_domain],
                                                             neighbor_desc    [i_domain]);

    /*
     * Echange du numero :
     *     - i_domain (useless now)
     *     - cell_ln_to_gn
     *     - cell_face (in gn )
     *     - face_ln_to_gn
     */
    PDM_g_num_t **cell_ln_to_gn_opp;
    PDM_distant_neighbor_exch(dn,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_CST,
                              1,
                              NULL,
                    (void **) cell_border_ln_to_gn,
                              NULL,
                   (void ***)&cell_ln_to_gn_opp);

    int **graph_comm_cell_opp;
    PDM_distant_neighbor_exch(dn,
                              sizeof(int),
                              PDM_STRIDE_CST,
                              3,
                              NULL,
                    (void **) graph_comm_cell[i_domain],
                              NULL,
                   (void ***)&graph_comm_cell_opp);

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      PDM_log_trace_array_long(cell_border_ln_to_gn[i_part], n_face_part_bound[i_domain][i_part], "cell_border_ln_to_gn = ");
      PDM_log_trace_array_long(cell_ln_to_gn_opp[i_part], n_face_part_bound[i_domain][i_part], "cell_ln_to_gn_opp = ");
      PDM_log_trace_array_int(graph_comm_cell_opp[i_part], 3 * n_face_part_bound[i_domain][i_part], "graph_comm_cell_opp = ");
      free(cell_border_ln_to_gn[i_part]);
      free(cell_ln_to_gn_opp[i_part]);
      free(graph_comm_cell_opp[i_part]);
    }
    free(cell_border_ln_to_gn);
    free(cell_ln_to_gn_opp);
    free(graph_comm_cell_opp);

    PDM_distant_neighbor_free(dn);

  }

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(neighbor_idx [i_domain][i_part]);
      free(neighbor_desc[i_domain][i_part]);
      free(graph_comm_cell[i_domain][i_part]);
      free(cell_list[i_domain][i_part]);
    }
    free(neighbor_idx[i_domain]);
    free(neighbor_desc[i_domain]);
    free(n_face_part_bound[i_domain]);
    free(graph_comm_cell[i_domain]);
    free(cell_list[i_domain]);
  }
  free(neighbor_idx);
  free(neighbor_desc);
  free(n_face_part_bound);
  free(graph_comm_cell);
  free(cell_list);

  printf(" PDM_part_extension_compute end \n");
}


/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
)
{
  part_ext->n_part = NULL;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(part_ext->parts[i_domain]);
  }
  free(part_ext->parts);

  free(part_ext);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
