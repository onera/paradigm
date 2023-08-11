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
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_order.h"
#include "pdm_error.h"
#include "pdm_part_extension.h"
#include "pdm_part_to_part.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_part_extension_algorithm.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_domain_interface.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_vtk.h"

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
void
_part_extension_3d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);
}

static
void
_part_extension_2d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);

}

static
void
_part_extension_1d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);

}



/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*
 * Translate and post-treated link by interface with on entity to another
 *  Exemple of use :
 *     - vtx to face
 *     - vtx to edge
 *     - face to cell
 */

void
PDM_part_extension_interface_by_entity1_to_interface_by_entity2
(
  PDM_part_domain_interface_t  *pdi,
  PDM_bound_type_t              entity1_bound,
  int                           n_domain,
  int                          *n_part,
  int                         **pn_entity1_in,
  PDM_g_num_t                ***pentity1_ln_to_gn_in,
  int                        ***pentity1_hint_in,
  int                         **pn_entity2_in,
  PDM_g_num_t                ***pentity2_ln_to_gn_in,
  int                        ***pentity2_entity1_idx_in,
  int                        ***pentity2_entity1_in,
  PDM_MPI_Comm                  comm
)
{
  int n_part_tot = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    n_part_tot += n_part[i_domain];
  }

  int          *pn_entity1            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity1_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity1_hint         = NULL;
  int          *pn_entity2            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_entity1_idx  = malloc(n_part_tot * sizeof(int         *));
  int         **pentity2_entity1      = malloc(n_part_tot * sizeof(int         *));

  if(pentity1_hint_in != NULL) {
    pentity1_hint = malloc(n_part_tot * sizeof(int         *));
  }

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pn_entity1          [ln_part_tot] = pn_entity1_in          [i_dom][i_part];
      pentity1_ln_to_gn   [ln_part_tot] = pentity1_ln_to_gn_in   [i_dom][i_part];
      pn_entity2          [ln_part_tot] = pn_entity2_in          [i_dom][i_part];
      pentity2_ln_to_gn   [ln_part_tot] = pentity2_ln_to_gn_in   [i_dom][i_part];
      pentity2_entity1_idx[ln_part_tot] = pentity2_entity1_idx_in[i_dom][i_part];
      pentity2_entity1    [ln_part_tot] = pentity2_entity1_in    [i_dom][i_part];

      if(pentity1_hint_in != NULL && pentity1_hint_in[i_dom] != NULL) {
        pentity1_hint[ln_part_tot] = pentity1_hint_in[i_dom][i_part];
      }

      ln_part_tot += 1;
    }
  }


  int          *pn_entity1_num             = NULL;
  int         **pentity1_num               = NULL;
  int         **pentity1_opp_location_idx  = NULL;
  int         **pentity1_opp_location      = NULL;
  int         **pentity1_opp_interface_idx = NULL;
  int         **pentity1_opp_interface     = NULL;
  int         **pentity1_opp_sens          = NULL;
  PDM_g_num_t **pentity1_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         entity1_bound,
                                         pn_entity1,
                                         pentity1_ln_to_gn,
                                         &pn_entity1_num,
                                         &pentity1_num,
                                         &pentity1_opp_location_idx,
                                         &pentity1_opp_location,
                                         &pentity1_opp_interface_idx,
                                         &pentity1_opp_interface,
                                         &pentity1_opp_sens,
                                         &pentity1_opp_gnum);

  /*
   * Prepare creation of part_to_part
   */
  int         **part1_to_part2_idx              = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_triplet_idx      = NULL; //malloc(ln_part_tot * sizeof(int *));
  int         **part1_to_part2_triplet          = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_interface        = malloc(ln_part_tot * sizeof(int         *));

  int         **part1_to_part2_edge_n           = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_edge_vtx_n       = malloc(ln_part_tot * sizeof(int         *));
  PDM_g_num_t **part1_to_part2_edge_vtx_gnum    = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **part1_to_part2_edge_gnum        = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **part1_to_part2_edge_vtx_triplet = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_edge_triplet     = malloc(ln_part_tot * sizeof(int         *));

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);


  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(comm, n_part_tot);

  int **ppart_entity1_proc_idx = NULL;
  int **ppart_entity1_part_idx = NULL;
  int **ppart_entity1          = NULL; // (i_vtx, i_proc, i_part, i_vtx_opp)
  PDM_part_generate_entity_graph_comm(comm,
                                      part_distribution,
                                      NULL,
                                      pn_entity1,
                                      pentity1_ln_to_gn,
                                      pentity1_hint,
                                      n_part_tot,
                                      ppart_entity1_proc_idx,
                                      ppart_entity1_part_idx,
                                      ppart_entity1,
                                      NULL);

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    int *n_part_shift = (int *) malloc( (n_rank) * sizeof(int));
    PDM_MPI_Allgather(&li_part,
                      1,
                      PDM_MPI_INT,
                      n_part_shift,
                      1,
                      PDM_MPI_INT,
                      comm);

    if(0 == 1) {
      PDM_log_trace_array_int(n_part_shift, n_rank+1, "n_part_shift ::");
    }

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      li_part += 1;
    }

    free(n_part_shift);
  }



  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(ppart_entity1_proc_idx[i_part]);
    free(ppart_entity1_part_idx[i_part]);
    free(ppart_entity1         [i_part]);
  }
  free(ppart_entity1_proc_idx);
  free(ppart_entity1_part_idx);
  free(ppart_entity1         );



  free(part1_to_part2_idx             );
  // free(part1_to_part2_triplet_idx     );
  free(part1_to_part2_triplet         );
  free(part1_to_part2_interface       );
  free(part1_to_part2_edge_n          );
  free(part1_to_part2_edge_vtx_n      );
  free(part1_to_part2_edge_vtx_gnum   );
  free(part1_to_part2_edge_gnum       );
  free(part1_to_part2_edge_vtx_triplet);
  free(part1_to_part2_edge_triplet    );



  free(pn_entity1          );
  free(pentity1_ln_to_gn   );
  free(pn_entity2          );
  free(pentity2_ln_to_gn   );
  free(pentity2_entity1_idx);
  free(pentity2_entity1    );

  if(pentity1_hint_in != NULL) {
    free(pentity1_hint);
  }



  // On pourra retourné le pdi_entity2 également (pour la recursion par exemple)
}




void
PDM_part_extension_pconnectivity_to_pconnectivity
(
  PDM_part_domain_interface_t  *pdi_entity1,
  PDM_part_domain_interface_t  *pdi_entity2,
  int                           n_domain,
  int                          *n_part,
  int                         **pn_entity1,
  PDM_g_num_t                 **pentity1_ln_to_gn,
  int                         **pn_entity2,
  PDM_g_num_t                 **pentity2_ln_to_gn

)
{




}

void
PDM_part_extension_compute2
(
        PDM_part_extension_t *part_ext,
  const int                   dim
)
{
  // TODO : mv dim in create but break API

  /* Manage shift */
  // _offset_parts_by_domain(part_ext, 1);

  /* Manage dim */
  if(dim == 3) {
    _part_extension_3d(part_ext);
  } else if(dim == 2) {
    _part_extension_2d(part_ext);
  } else if(dim == 1) {
    _part_extension_1d(part_ext);
  } else  {
    PDM_error(__FILE__, __LINE__, 0, "Wrong dim size in PDM_part_extension_compute2 : %d ( Should be >=1 )\n", (int) dim);
  }


  /* Manage loop for depth AND multiple transformation */



  /* Manage unshift */
  // _offset_parts_by_domain(part_ext, -1);
  // _offset_results_by_domain(part_ext);

}



#ifdef __cplusplus
}
#endif /* __cplusplus */

