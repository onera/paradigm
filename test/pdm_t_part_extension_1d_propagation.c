#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_multipart.h"
#include "pdm_array.h"
#include "pdm_part_to_part.h"
#include "pdm_extract_part.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_domain_utils.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_unique.h"
#include "pdm_order.h"
#include "pdm_part_extension_algorithm.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_g_pts,
 int           *n_dom_i,
 int           *periodic_i,
 int           *n_depth
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_g_pts = atol(argv[i]);
        *n_g_pts = (PDM_g_num_t) _n_g_pts;
      }
    }
    else if (strcmp(argv[i], "-ni") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_i = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_depth = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pi") == 0) {
      *periodic_i = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
void
_part_extension_init
(
  PDM_MPI_Comm                  comm,
  int                           n_depth,
  PDM_part_domain_interface_t  *pdi,
  int                           n_domain,
  PDM_g_num_t                  *shift_by_domain_edge,
  PDM_g_num_t                  *shift_by_domain_vtx,
  int                          *n_part,
  int                         **pn_vtx,
  PDM_g_num_t                ***pvtx_ln_to_gn,
  double                     ***pvtx_coord,
  int                         **pn_edge,
  PDM_g_num_t                ***pedge_ln_to_gn,
  int                        ***pedge_vtx_idx,
  int                        ***pedge_vtx,
  int                         **out_pn_edge_extented,
  int                        ***out_pedge_extented_to_pedge_idx,
  int                        ***out_pedge_extented_to_pedge_triplet,
  PDM_g_num_t                ***out_pedge_extented_ln_to_gn,
  int                        ***out_pedge_extented_to_pedge_interface,
  int                         **out_pn_vtx_extented,
  PDM_g_num_t                ***out_pvtx_extented_ln_to_gn,
  int                        ***out_pextented_edge_vtx_idx,
  int                        ***out_pextented_edge_vtx,
  int                        ***out_pvtx_extented_to_pvtx_idx,
  int                        ***out_pvtx_extented_to_pvtx_triplet,
  int                        ***out_pvtx_extented_to_pvtx_interface,
  int                         **out_pn_concat_vtx,
  PDM_g_num_t                ***out_pconcat_vtx_ln_to_gn,
  int                         **out_pn_concat_edge,
  PDM_g_num_t                ***out_pconcat_edge_ln_to_gn,
  int                        ***out_pconcat_edge_vtx,
  int                        ***out_pconcat_edge_vtx_idx,
  double                     ***out_pconcat_vtx_coords,
  int                        ***out_pconcat_pvtx_extented_to_pvtx_idx,
  int                        ***out_pconcat_pvtx_extented_to_pvtx_triplet,
  int                        ***out_pconcat_pvtx_extented_to_pvtx_interface,
  int                        ***out_pconcat_pedge_extented_to_pedge_idx,
  int                        ***out_pconcat_pedge_extented_to_pedge_triplet,
  int                        ***out_pconcat_pedge_extented_to_pedge_interface
)
{
  PDM_UNUSED(n_depth);

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      ln_part_tot += 1;
    }
  }

  int          *pn_edge_extented                = NULL;
  int         **pedge_extented_to_pedge_idx     = NULL;
  int         **pedge_extented_to_pedge_triplet = NULL;
  PDM_g_num_t **pedge_extented_ln_to_gn         = NULL;
  int         **pedge_extented_to_pedge_interface = NULL;

  PDM_part_extension_interface_by_entity1_to_interface_by_entity2(pdi,
                                                                  PDM_BOUND_TYPE_VTX,
                                                                  n_domain,
                                                                  shift_by_domain_edge,
                                                                  n_part,
                                                                  pn_vtx,
                                                                  pvtx_ln_to_gn,
                                                                  NULL, // pvtx_hint,
                                                                  pn_edge,
                                                                  pedge_ln_to_gn,
                                                                  pedge_vtx_idx,
                                                                  pedge_vtx,
                                                                  &pn_edge_extented,
                                                                  &pedge_extented_ln_to_gn,
                                                                  &pedge_extented_to_pedge_idx,
                                                                  &pedge_extented_to_pedge_triplet,
                                                                  &pedge_extented_to_pedge_interface,
                                                                  comm);


  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      int n_triplet = pedge_extented_to_pedge_idx[i_part][ pn_edge_extented[i_part]];
      PDM_log_trace_array_long(pedge_extented_ln_to_gn          [i_part], pn_edge_extented[i_part]  , "pedge_extented_ln_to_gn : ");
      PDM_log_trace_array_int (pedge_extented_to_pedge_idx      [i_part], pn_edge_extented[i_part]+1, "pedge_extented_to_pedge_idx      ::");
      PDM_log_trace_array_int (pedge_extented_to_pedge_interface[i_part], n_triplet/3               , "pedge_extented_to_pedge_interedge       ::");
      PDM_log_trace_array_int (pedge_extented_to_pedge_triplet  [i_part], n_triplet                 , "pedge_extented_to_pedge_triplet ::");
    }
  }

  int          *pn_vtx_extented                 = NULL;
  PDM_g_num_t **pvtx_extented_ln_to_gn          = NULL;
  int         **pextented_edge_vtx_idx          = NULL;
  int         **pextented_edge_vtx              = NULL;
  int         **pvtx_extented_to_pvtx_idx       = NULL;
  int         **pvtx_extented_to_pvtx_triplet   = NULL;
  int         **pvtx_extented_to_pvtx_interface = NULL;
  // Rebuild edge_vtx direchly (after we will do face_edge then edge_vtx)

  PDM_part_extension_pconnectivity_to_extented_pconnectivity(pdi,
                                                             PDM_BOUND_TYPE_VTX,
                                                             n_domain,
                                                             shift_by_domain_vtx,
                                                             n_part,
                                                             pn_edge,
                                                             pedge_ln_to_gn,
                                                             pn_vtx,
                                                             pvtx_ln_to_gn,
                                                             pedge_vtx_idx,
                                                             pedge_vtx,
                                                             pn_edge_extented,
                                                             pedge_extented_ln_to_gn,
                                                             pedge_extented_to_pedge_idx,
                                                             pedge_extented_to_pedge_triplet,
                                                             pedge_extented_to_pedge_interface,
                                                             &pn_vtx_extented,
                                                             &pvtx_extented_ln_to_gn,
                                                             &pextented_edge_vtx_idx,
                                                             &pextented_edge_vtx,
                                                             &pvtx_extented_to_pvtx_idx,
                                                             &pvtx_extented_to_pvtx_triplet,
                                                             &pvtx_extented_to_pvtx_interface,
                                                             comm);

  /*
   * Hook coordinates
   */
  int           *pflat_n_vtx      = NULL;
  double       **pflat_vtx_coords = NULL;
  PDM_malloc(pflat_n_vtx     , ln_part_tot, int     );
  PDM_malloc(pflat_vtx_coords, ln_part_tot, double *);

  int shift_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pflat_n_vtx[shift_part+i_part] = pn_vtx[i_dom][i_part];
      pflat_vtx_coords[shift_part+i_part] = pvtx_coord[i_dom][i_part];
    }
    shift_part += n_part[i_dom];
  }

  PDM_part_to_part_t* ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pvtx_extented_ln_to_gn,
                                                                          (const int          *) pn_vtx_extented,
                                                                          ln_part_tot,
                                                                          (const int          *) pflat_n_vtx,
                                                                          ln_part_tot,
                                                                          (const int         **) pvtx_extented_to_pvtx_idx,
                                                                          (const int         **) NULL,
                                                                          (const int         **) pvtx_extented_to_pvtx_triplet,
                                                                          comm);

  // log_trace(" +++++++++++++++++++++ DEBUG +++++++++++++++++++++ \n");
  // for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

  //   for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]; ++i_vtx) {
  //     log_trace("pvtx_extented_ln_to_gn[%i] = %i \n", i_vtx, pvtx_extented_ln_to_gn[i_part][i_vtx]);
  //     for(int idx_vtx = pvtx_extented_to_pvtx_idx[i_part][i_vtx]/3; idx_vtx < pvtx_extented_to_pvtx_idx[i_part][i_vtx+1]/3; ++idx_vtx) {
  //       log_trace("\t --> (%i, %i , %i) - (%i) \n", pvtx_extented_to_pvtx_triplet  [i_part][3*idx_vtx  ],
  //                                                   pvtx_extented_to_pvtx_triplet  [i_part][3*idx_vtx+1],
  //                                                   pvtx_extented_to_pvtx_triplet  [i_part][3*idx_vtx+2],
  //                                                   pvtx_extented_to_pvtx_interface[i_part][idx_vtx]);
  //     }
  //   }
  // }


  /*
   * Il faut d'abord concatener (ou reallouer ) les connectivités edge_vtx et coordonnées
   * On doit faire pour l'étape d'après 1 part_to_part avec pn_vtx_concat
   * On fait un echange en gnum_come_from
   * On fait un post-tratiement
   * Il faut garder l'info du gnum concaténé et la provenance en interface !!!!
   *
   *
   */

  /*
   *
   */
  int exch_request = -1;
  double      **pextract_vtx_coords           = NULL;
  PDM_part_to_part_reverse_iexch(ptp_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 3 * sizeof(double),
                                 NULL,
                (const void **)  pflat_vtx_coords,
                                 NULL,
                    (void ***)   &pextract_vtx_coords,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp_vtx, exch_request);

  PDM_part_to_part_free(ptp_vtx);

  /*
   * Apply transformation if any
   */
  int n_interface = 0;
  if(pdi != NULL) {
    n_interface = PDM_part_domain_interface_n_interface_get(pdi);
  }
  double  **translation_vector = NULL;
  double ***rotation_matrix    = NULL;
  double  **rotation_direction = NULL;
  double  **rotation_center    = NULL;
  double   *rotation_angle     = NULL;
  PDM_malloc(translation_vector, n_interface, double  *);
  PDM_malloc(rotation_matrix   , n_interface, double **);
  PDM_malloc(rotation_direction, n_interface, double  *);
  PDM_malloc(rotation_center   , n_interface, double  *);
  PDM_malloc(rotation_angle    , n_interface, double   );
  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    translation_vector[i_interf] = NULL;
    PDM_part_domain_interface_translation_get(pdi, i_interf, &translation_vector[i_interf]);

    rotation_matrix[i_interf] = NULL;
    PDM_part_domain_interface_rotation_get   (pdi,
                                              i_interf,
                                              &rotation_direction[i_interf],
                                              &rotation_center   [i_interf],
                                              &rotation_angle    [i_interf]);

    if(rotation_center    [i_interf] != NULL) {
      PDM_malloc(rotation_matrix[i_interf], 3, double *);
      for(int k = 0; k < 3; ++k) {
        PDM_malloc(rotation_matrix[i_interf][k], 3, double);
      }
    }
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]; ++i_vtx) {
      int i_interface   = PDM_ABS (pvtx_extented_to_pvtx_interface[i_part][i_vtx]);
      int sgn_interface = PDM_SIGN(pvtx_extented_to_pvtx_interface[i_part][i_vtx]);
      if(i_interface != 0 && translation_vector[PDM_ABS(i_interface)-1] != NULL) {
        for(int k = 0; k < 3; ++k) {
          pextract_vtx_coords[i_part][3*i_vtx+k] += sgn_interface * translation_vector[PDM_ABS(i_interface)-1][k];
        }
      }
    }
  }


  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    if(translation_vector[i_interf] != NULL) {
      PDM_free(translation_vector[i_interf]);
    }
    if(rotation_center    [i_interf] != NULL) {
      for(int k = 0; k < 3; ++k) {
        PDM_free(rotation_matrix[i_interf][k]);
      }
      PDM_free(rotation_matrix[i_interf]);
    }
  }
  PDM_free(translation_vector);
  PDM_free(rotation_matrix);
  PDM_free(rotation_direction);
  PDM_free(rotation_center);
  PDM_free(rotation_angle);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Concatenation :
   *   - size
   *   - global numbering
   *   - connectivity
   *   - link for part_to_part
   */
  int *pn_concat_edge = NULL;
  int *pn_concat_vtx  = NULL;
  PDM_malloc(pn_concat_edge, ln_part_tot, int);
  PDM_malloc(pn_concat_vtx , ln_part_tot, int);

  PDM_g_num_t **concat_edge_ln_to_gn = NULL;
  PDM_g_num_t **concat_vtx_ln_to_gn  = NULL;
  PDM_malloc(concat_edge_ln_to_gn, ln_part_tot, PDM_g_num_t *);
  PDM_malloc(concat_vtx_ln_to_gn , ln_part_tot, PDM_g_num_t *);

  int    **concat_edge_vtx     = NULL;
  int    **concat_edge_vtx_idx = NULL;
  double **concat_vtx_coords   = NULL;
  PDM_malloc(concat_edge_vtx    , ln_part_tot, int    *);
  PDM_malloc(concat_edge_vtx_idx, ln_part_tot, int    *);
  PDM_malloc(concat_vtx_coords  , ln_part_tot, double *);

  int **edge_kind = NULL;
  PDM_malloc(edge_kind, ln_part_tot, int *);

  int **concat_pvtx_extented_to_pvtx_idx       = NULL;;
  int **concat_pvtx_extented_to_pvtx_triplet   = NULL;;
  int **concat_pvtx_extented_to_pvtx_interface = NULL;;
  PDM_malloc(concat_pvtx_extented_to_pvtx_idx      , ln_part_tot, int *);
  PDM_malloc(concat_pvtx_extented_to_pvtx_triplet  , ln_part_tot, int *);
  PDM_malloc(concat_pvtx_extented_to_pvtx_interface, ln_part_tot, int *);

  int **concat_pedge_extented_to_pedge_idx       = NULL;
  int **concat_pedge_extented_to_pedge_triplet   = NULL;
  int **concat_pedge_extented_to_pedge_interface = NULL;
  PDM_malloc(concat_pedge_extented_to_pedge_idx      , ln_part_tot, int *);
  PDM_malloc(concat_pedge_extented_to_pedge_triplet  , ln_part_tot, int *);
  PDM_malloc(concat_pedge_extented_to_pedge_interface, ln_part_tot, int *);

  shift_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_concat_vtx [shift_part+i_part] = pn_vtx [i_dom][i_part] + pn_vtx_extented [shift_part+i_part];
      pn_concat_edge[shift_part+i_part] = pn_edge[i_dom][i_part] + pn_edge_extented[shift_part+i_part];

      int ln_vtx_extented  = pn_vtx_extented [shift_part+i_part];
      int ln_edge_extented = pn_edge_extented[shift_part+i_part];

      PDM_malloc(edge_kind           [shift_part+i_part],     pn_concat_edge[shift_part+i_part]  , int);
      PDM_malloc(concat_edge_vtx     [shift_part+i_part], 2 * pn_concat_edge[shift_part+i_part]  , int);
      PDM_malloc(concat_edge_vtx_idx [shift_part+i_part],     pn_concat_edge[shift_part+i_part]+1, int);

      for(int i_edge = 0; i_edge < pn_concat_edge[shift_part+i_part]+1; ++i_edge) {
        concat_edge_vtx_idx [shift_part+i_part][i_edge] = 2 * i_edge;
      }

      // PDM_log_trace_array_int(concat_edge_vtx_idx [shift_part+i_part],pn_concat_edge[shift_part+i_part]+1, "concat_edge_vtx_idx ::");

      PDM_malloc(concat_edge_ln_to_gn[shift_part+i_part],     pn_concat_edge[shift_part+i_part], PDM_g_num_t);
      PDM_malloc(concat_vtx_ln_to_gn [shift_part+i_part],     pn_concat_vtx [shift_part+i_part], PDM_g_num_t);
      PDM_malloc(concat_vtx_coords   [shift_part+i_part], 3 * pn_concat_vtx [shift_part+i_part], double     );

      for(int i_edge = 0; i_edge < pn_edge[i_dom][i_part]; ++i_edge) {
        concat_edge_vtx     [shift_part+i_part][2*i_edge  ] = pedge_vtx[i_dom][i_part][2*i_edge  ];
        concat_edge_vtx     [shift_part+i_part][2*i_edge+1] = pedge_vtx[i_dom][i_part][2*i_edge+1];
        concat_edge_ln_to_gn[shift_part+i_part][i_edge] = pedge_ln_to_gn[i_dom][i_part][i_edge];
        edge_kind[shift_part+i_part][i_edge] = 0;
      }

      for(int i_vtx = 0; i_vtx < pflat_n_vtx[shift_part+i_part]; ++i_vtx) {
        concat_vtx_coords[shift_part+i_part][3*i_vtx  ] = pflat_vtx_coords[i_part][3*i_vtx  ];
        concat_vtx_coords[shift_part+i_part][3*i_vtx+1] = pflat_vtx_coords[i_part][3*i_vtx+1];
        concat_vtx_coords[shift_part+i_part][3*i_vtx+2] = pflat_vtx_coords[i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[shift_part+i_part][i_vtx] = pvtx_ln_to_gn[i_dom][i_part][i_vtx];
      }

      for(int i_edge = 0; i_edge < pn_edge_extented[shift_part+i_part]; ++i_edge) {
        int idx_write = pn_edge[i_dom][i_part]+i_edge;
        concat_edge_vtx     [shift_part+i_part][2*idx_write  ] = pextented_edge_vtx     [shift_part+i_part][2*i_edge  ];
        concat_edge_vtx     [shift_part+i_part][2*idx_write+1] = pextented_edge_vtx     [shift_part+i_part][2*i_edge+1];
        concat_edge_ln_to_gn[shift_part+i_part][  idx_write  ] = pedge_extented_ln_to_gn[shift_part+i_part][i_edge];
        edge_kind           [shift_part+i_part][idx_write] = 1;
      }

      for(int i_vtx = 0; i_vtx < pn_vtx_extented[shift_part+i_part]; ++i_vtx) {
        int idx_write = pn_vtx [i_dom][i_part]+i_vtx;
        concat_vtx_coords  [shift_part+i_part][3*idx_write  ] = pextract_vtx_coords   [shift_part+i_part][3*i_vtx  ];
        concat_vtx_coords  [shift_part+i_part][3*idx_write+1] = pextract_vtx_coords   [shift_part+i_part][3*i_vtx+1];
        concat_vtx_coords  [shift_part+i_part][3*idx_write+2] = pextract_vtx_coords   [shift_part+i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[shift_part+i_part][  idx_write  ] = pvtx_extented_ln_to_gn[shift_part+i_part][  i_vtx  ];
      }

      /*
       * Graph between concatenate partition and the extended ones
       */
      int *_pvtx_extented_to_pvtx_idx   = pvtx_extented_to_pvtx_idx  [shift_part+i_part];
      int *_pedge_extented_to_pedge_idx = pedge_extented_to_pedge_idx[shift_part+i_part];

      int n_connect_vtx  = _pvtx_extented_to_pvtx_idx  [ln_vtx_extented ]/3;
      int n_connect_edge = _pedge_extented_to_pedge_idx[ln_edge_extented]/3;

      PDM_malloc(concat_pvtx_extented_to_pvtx_idx  [shift_part+i_part], pn_concat_vtx [shift_part+i_part]+1, int);
      PDM_malloc(concat_pedge_extented_to_pedge_idx[shift_part+i_part], pn_concat_edge[shift_part+i_part]+1, int);

      int *_concat_pvtx_extented_to_pvtx_idx   = concat_pvtx_extented_to_pvtx_idx  [shift_part+i_part];
      int *_concat_pedge_extented_to_pedge_idx = concat_pedge_extented_to_pedge_idx[shift_part+i_part];

      concat_pvtx_extented_to_pvtx_idx[shift_part+i_part][0] = 0;
      for(int i_vtx = 0; i_vtx < pflat_n_vtx[shift_part+i_part]; ++i_vtx) {
        _concat_pvtx_extented_to_pvtx_idx[i_vtx+1] = _concat_pvtx_extented_to_pvtx_idx[i_vtx];
      }

      for(int i_vtx = 0; i_vtx < pn_vtx_extented[shift_part+i_part]; ++i_vtx) {
        int idx_write = pn_vtx [i_dom][i_part]+i_vtx;
        int ln_connect = _pvtx_extented_to_pvtx_idx[i_vtx+1] - _pvtx_extented_to_pvtx_idx[i_vtx];
        _concat_pvtx_extented_to_pvtx_idx[idx_write+1] = _concat_pvtx_extented_to_pvtx_idx[idx_write] + ln_connect;
      }

      concat_pedge_extented_to_pedge_idx        [shift_part+i_part][0] = 0;
      for(int i_edge = 0; i_edge < pn_edge[i_dom][i_part]; ++i_edge) {
        _concat_pedge_extented_to_pedge_idx[i_edge+1] = _concat_pedge_extented_to_pedge_idx[i_edge];
      }

      for(int i_edge = 0; i_edge < pn_edge_extented[shift_part+i_part]; ++i_edge) {
        int idx_write = pn_edge [i_dom][i_part]+i_edge;
        int ln_connect = _pedge_extented_to_pedge_idx[i_edge+1] - _pedge_extented_to_pedge_idx[i_edge];
        _concat_pedge_extented_to_pedge_idx[idx_write+1] = _concat_pedge_extented_to_pedge_idx[idx_write] + ln_connect;
      }

      PDM_malloc(concat_pvtx_extented_to_pvtx_triplet    [shift_part+i_part], 3 * n_connect_vtx , int);
      PDM_malloc(concat_pvtx_extented_to_pvtx_interface  [shift_part+i_part],     n_connect_vtx , int);
      PDM_malloc(concat_pedge_extented_to_pedge_triplet  [shift_part+i_part], 3 * n_connect_edge, int);
      PDM_malloc(concat_pedge_extented_to_pedge_interface[shift_part+i_part],     n_connect_edge, int);

      for(int i = 0; i < n_connect_vtx; ++i) {
        concat_pvtx_extented_to_pvtx_triplet  [shift_part+i_part][3*i  ] = pvtx_extented_to_pvtx_triplet  [shift_part+i_part][3*i  ];
        concat_pvtx_extented_to_pvtx_triplet  [shift_part+i_part][3*i+1] = pvtx_extented_to_pvtx_triplet  [shift_part+i_part][3*i+1];
        concat_pvtx_extented_to_pvtx_triplet  [shift_part+i_part][3*i+2] = pvtx_extented_to_pvtx_triplet  [shift_part+i_part][3*i+2];
        concat_pvtx_extented_to_pvtx_interface[shift_part+i_part][  i  ] = pvtx_extented_to_pvtx_interface[shift_part+i_part][  i  ];
      }

      for(int i = 0; i < n_connect_edge; ++i) {
        concat_pedge_extented_to_pedge_triplet  [shift_part+i_part][3*i  ] = pedge_extented_to_pedge_triplet  [shift_part+i_part][3*i  ];
        concat_pedge_extented_to_pedge_triplet  [shift_part+i_part][3*i+1] = pedge_extented_to_pedge_triplet  [shift_part+i_part][3*i+1];
        concat_pedge_extented_to_pedge_triplet  [shift_part+i_part][3*i+2] = pedge_extented_to_pedge_triplet  [shift_part+i_part][3*i+2];
        concat_pedge_extented_to_pedge_interface[shift_part+i_part][  i  ] = pedge_extented_to_pedge_interface[shift_part+i_part][  i  ];
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_concat_pvtx_extented_to_pvtx_idx  , pn_concat_vtx [shift_part+i_part]+1, "_concat_pvtx_extented_to_pvtx_idx   ::");
        PDM_log_trace_array_int(_concat_pedge_extented_to_pedge_idx, pn_concat_edge[shift_part+i_part]+1, "_concat_pedge_extented_to_pedge_idx ::");
      }


    }
    shift_part += n_part[i_dom];
  }



  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      char filename[999];
      sprintf(filename, "out_part_concate_vtx_i_part=%i_%i.vtk", i_part, i_rank);
      const char* field_name[] = {"edge_kind", 0 };
      int *field[1] = {edge_kind[i_part]};
      PDM_vtk_write_std_elements (filename,
                                  pn_concat_vtx[i_part],
                                  concat_vtx_coords[i_part],
                                  concat_vtx_ln_to_gn[i_part],
                                  PDM_MESH_NODAL_BAR2,
                                  pn_concat_edge[i_part],
                                  concat_edge_vtx[i_part],
                                  concat_edge_ln_to_gn[i_part],
                                  1,
                                  field_name,
                   (const int **) field);
    }
  }


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pextract_vtx_coords[i_part]);
  }
  PDM_free(pflat_n_vtx);
  PDM_free(pflat_vtx_coords);
  PDM_free(pextract_vtx_coords);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(edge_kind           [i_part]);
  }
  PDM_free(edge_kind           );

  *out_pn_edge_extented                  = pn_edge_extented;
  *out_pedge_extented_to_pedge_idx       = pedge_extented_to_pedge_idx;
  *out_pedge_extented_to_pedge_triplet   = pedge_extented_to_pedge_triplet;
  *out_pedge_extented_ln_to_gn           = pedge_extented_ln_to_gn;
  *out_pedge_extented_to_pedge_interface = pedge_extented_to_pedge_interface;
  *out_pn_vtx_extented                   = pn_vtx_extented;
  *out_pvtx_extented_ln_to_gn            = pvtx_extented_ln_to_gn;
  *out_pextented_edge_vtx_idx            = pextented_edge_vtx_idx;
  *out_pextented_edge_vtx                = pextented_edge_vtx;
  *out_pvtx_extented_to_pvtx_idx         = pvtx_extented_to_pvtx_idx;
  *out_pvtx_extented_to_pvtx_triplet     = pvtx_extented_to_pvtx_triplet;
  *out_pvtx_extented_to_pvtx_interface   = pvtx_extented_to_pvtx_interface;

  *out_pn_concat_vtx                             = pn_concat_vtx;
  *out_pconcat_vtx_ln_to_gn                      = concat_vtx_ln_to_gn;
  *out_pn_concat_edge                            = pn_concat_edge;
  *out_pconcat_edge_ln_to_gn                     = concat_edge_ln_to_gn;
  *out_pconcat_edge_vtx                          = concat_edge_vtx;
  *out_pconcat_edge_vtx_idx                      = concat_edge_vtx_idx;
  *out_pconcat_vtx_coords                        = concat_vtx_coords;
  *out_pconcat_pvtx_extented_to_pvtx_idx         = concat_pvtx_extented_to_pvtx_idx;
  *out_pconcat_pvtx_extented_to_pvtx_triplet     = concat_pvtx_extented_to_pvtx_triplet;
  *out_pconcat_pvtx_extented_to_pvtx_interface   = concat_pvtx_extented_to_pvtx_interface;
  *out_pconcat_pedge_extented_to_pedge_idx       = concat_pedge_extented_to_pedge_idx;
  *out_pconcat_pedge_extented_to_pedge_triplet   = concat_pedge_extented_to_pedge_triplet;
  *out_pconcat_pedge_extented_to_pedge_interface = concat_pedge_extented_to_pedge_interface;



}

static
void
_part_extension_one_depth
(
  PDM_MPI_Comm                  comm,
  PDM_part_domain_interface_t  *pdi,
  int                           n_domain,
  PDM_g_num_t                  *shift_by_domain_edge,
  PDM_g_num_t                  *shift_by_domain_vtx,
  int                          *n_part,
  int                         **out_pn_vtx,
  PDM_g_num_t                ***out_pvtx_ln_to_gn,
  int                         **out_pn_edge,
  PDM_g_num_t                ***out_pedge_ln_to_gn,
  int                        ***out_pedge_vtx,
  int                        ***out_pedge_vtx_idx,
  double                     ***out_pvtx_coords,
  int                        ***out_pvtx_extented_to_pvtx_idx,
  int                        ***out_pvtx_extented_to_pvtx_triplet,
  int                        ***out_pvtx_extented_to_pvtx_interface,
  int                        ***out_pedge_extented_to_pedge_idx,
  int                        ***out_pedge_extented_to_pedge_triplet,
  int                        ***out_pedge_extented_to_pedge_interface
)
{

  PDM_UNUSED(pdi);
  PDM_UNUSED(shift_by_domain_edge);
  PDM_UNUSED(shift_by_domain_vtx);

  PDM_UNUSED(out_pvtx_ln_to_gn);
  PDM_UNUSED(out_pedge_ln_to_gn);
  PDM_UNUSED(out_pvtx_coords);
  PDM_UNUSED(out_pedge_extented_to_pedge_idx);
  PDM_UNUSED(out_pedge_extented_to_pedge_triplet);
  PDM_UNUSED(out_pedge_extented_to_pedge_interface);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int          *pn_vtx                            = *out_pn_vtx;
  PDM_g_num_t **pvtx_ln_to_gn                     = *out_pvtx_ln_to_gn;
  int          *pn_edge                           = *out_pn_edge;
  PDM_g_num_t **pedge_ln_to_gn                    = *out_pedge_ln_to_gn;
  int         **pedge_vtx                         = *out_pedge_vtx;
  int         **pedge_vtx_idx                     = *out_pedge_vtx_idx;
  // double      **pvtx_coords                       = *out_pvtx_coords;
  int         **pvtx_extented_to_pvtx_idx         = *out_pvtx_extented_to_pvtx_idx;
  int         **pvtx_extented_to_pvtx_triplet     = *out_pvtx_extented_to_pvtx_triplet;
  int         **pvtx_extented_to_pvtx_interface   = *out_pvtx_extented_to_pvtx_interface;
  int         **pedge_extented_to_pedge_idx       = *out_pedge_extented_to_pedge_idx;
  int         **pedge_extented_to_pedge_triplet   = *out_pedge_extented_to_pedge_triplet;
  int         **pedge_extented_to_pedge_interface = *out_pedge_extented_to_pedge_interface;


  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      ln_part_tot += 1;
    }
  }

  /*
   * Compute offset gnum
   */
  PDM_g_num_t _lshift_edge_gnum = -1;
  PDM_g_num_t _lshift_vtx_gnum  = -1;

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    for(int i_edge = 0; i_edge < pn_edge[i_part]; ++i_edge) {
      _lshift_edge_gnum = PDM_MAX(_lshift_edge_gnum, pedge_ln_to_gn[i_part][i_edge]);
    }
    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {
      _lshift_vtx_gnum = PDM_MAX(_lshift_vtx_gnum, pvtx_ln_to_gn[i_part][i_vtx]);
    }
  }

  PDM_g_num_t shift_edge_gnum = -1;
  PDM_g_num_t shift_vtx_gnum  = -1;

  PDM_MPI_Allreduce(&_lshift_edge_gnum, &shift_edge_gnum, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_MPI_Allreduce(&_lshift_vtx_gnum , &shift_vtx_gnum , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);


  int         **pvtx_edge_idx  = NULL;
  int         **pvtx_edge      = NULL;
  PDM_part_connectivity_transpose(ln_part_tot,
                                  pn_edge,
                                  pn_vtx,
                                  pedge_vtx_idx,
                                  pedge_vtx,
                                  &pvtx_edge_idx,
                                  &pvtx_edge);


  PDM_part_to_part_t* ptp_edge = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pedge_ln_to_gn,
                                                                           (const int          *) pn_edge,
                                                                           ln_part_tot,
                                                                           (const int          *) pn_edge,
                                                                           ln_part_tot,
                                                                           (const int         **) pedge_extented_to_pedge_idx,
                                                                           (const int         **) NULL,
                                                                           (const int         **) pedge_extented_to_pedge_triplet,
                                                                           comm);

  /*
   * Exchange original edge gnum
   */
  int exch_request_edge_gnum = -1;
  PDM_g_num_t **pextented_edge_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp_edge,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                                 NULL,
                (const void **)  pedge_ln_to_gn,
                                 NULL,
                    (void ***)   &pextented_edge_gnum,
                                 &exch_request_edge_gnum);
  PDM_part_to_part_reverse_iexch_wait(ptp_edge, exch_request_edge_gnum);

  PDM_part_to_part_free(ptp_edge);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int n_part1_to_part2 = pedge_extented_to_pedge_idx[i_part][pn_edge[i_part]]/3;
    // PDM_log_trace_array_long(pextented_edge_gnum[i_part], n_part1_to_part2, "pextented_edge_gnum ::");
  }

  /*
   * Prepare previous gnum that alreay computed
   *  Caution : We need to preshift the gnum to have the "true" gnum
   */
  PDM_g_num_t **prev_edge_ln_to_gn_and_interface = NULL;
  PDM_malloc(prev_edge_ln_to_gn_and_interface, ln_part_tot, PDM_g_num_t *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int n_part1_to_part2 = pedge_extented_to_pedge_idx[i_part][pn_edge[i_part]]/3;
    PDM_malloc(prev_edge_ln_to_gn_and_interface[i_part], 2 * n_part1_to_part2, PDM_g_num_t);

    for(int i = 0; i < n_part1_to_part2; ++i) {
      prev_edge_ln_to_gn_and_interface[i_part][2*i  ] = pextented_edge_gnum              [i_part][i];
      prev_edge_ln_to_gn_and_interface[i_part][2*i+1] = pedge_extented_to_pedge_interface[i_part][i];
    }

    int *order = NULL;
    PDM_malloc(order, n_part1_to_part2, int);

    PDM_order_gnum_s(prev_edge_ln_to_gn_and_interface[i_part], 2, order, n_part1_to_part2);

    PDM_order_array(n_part1_to_part2, 2 * sizeof(PDM_g_num_t), order, prev_edge_ln_to_gn_and_interface[i_part]);


    PDM_free(order);
    PDM_free(pextented_edge_gnum[i_part]);


    // PDM_log_trace_array_long(prev_edge_ln_to_gn_and_interface[i_part], 2 * n_part1_to_part2, "prev_edge_ln_to_gn_and_interface ::");


  }
  PDM_free(pextented_edge_gnum);




  PDM_part_to_part_t* ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pvtx_ln_to_gn,
                                                                          (const int          *) pn_vtx,
                                                                          ln_part_tot,
                                                                          (const int          *) pn_vtx,
                                                                          ln_part_tot,
                                                                          (const int         **) pvtx_extented_to_pvtx_idx,
                                                                          (const int         **) NULL,
                                                                          (const int         **) pvtx_extented_to_pvtx_triplet,
                                                                          comm);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_vtx, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp_vtx, &gnum1_come_from_idx, &gnum1_come_from);

  /*
   * Prepare buffer
   */
  int         **gnum1_com_from_triplet_n    = NULL;
  int         **gnum1_com_from_triplet_send = NULL;
  PDM_g_num_t **gnum1_com_from_gnum_send    = NULL;
  PDM_malloc(gnum1_com_from_triplet_n   , ln_part_tot, int         *);
  PDM_malloc(gnum1_com_from_triplet_send, ln_part_tot, int         *);
  PDM_malloc(gnum1_com_from_gnum_send   , ln_part_tot, PDM_g_num_t *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pvtx_edge_idx   = pvtx_edge_idx      [i_part];

    if(0 == 1) {
      printf("n_gnum1_come_from = %i \n", n_gnum1_come_from);
      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 ::");
      PDM_log_trace_array_int(gnum1_come_from_idx[i_part], n_ref_lnum2[i_part]+1, "gnum1_come_from_idx ::");
    }

    // PDM_log_trace_connectivity_int(pedge_vtx_idx[i_part], pedge_vtx[i_part], pn_edge[i_part], "pedge_vtx dbg ::");

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        n_send_part2 += _pvtx_edge_idx[i_entity1+1] - _pvtx_edge_idx[i_entity1];
      }
    }

    printf("n_send_part2 = %i \n", n_send_part2);
    /* Allocate */
    PDM_malloc(gnum1_com_from_triplet_n   [i_part],     n_gnum1_come_from, int        );
    PDM_malloc(gnum1_com_from_triplet_send[i_part], 3 * n_send_part2     , int        );
    PDM_malloc(gnum1_com_from_gnum_send   [i_part],     n_send_part2     , PDM_g_num_t);

    int         *_gnum1_com_from_triplet_n    = gnum1_com_from_triplet_n   [i_part];
    int         *_gnum1_com_from_triplet_send = gnum1_com_from_triplet_send[i_part];
    PDM_g_num_t *_gnum1_com_from_gnum_send    = gnum1_com_from_gnum_send   [i_part];

    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        _gnum1_com_from_triplet_n[k] = _pvtx_edge_idx[i_entity1+1] - _pvtx_edge_idx[i_entity1];
        for(int idx_edge = _pvtx_edge_idx[i_entity1]; idx_edge < _pvtx_edge_idx[i_entity1+1]; ++idx_edge) {

          int i_edge = PDM_ABS(pvtx_edge[i_part][idx_edge])-1;
          _gnum1_com_from_gnum_send   [  n_send_part2  ] = pedge_ln_to_gn[i_part][i_edge];
          _gnum1_com_from_triplet_send[3*n_send_part2  ] = i_rank;
          _gnum1_com_from_triplet_send[3*n_send_part2+1] = i_part;
          _gnum1_com_from_triplet_send[3*n_send_part2+2] = i_edge;
          n_send_part2++;
        }
      }
    }

    // PDM_log_trace_array_long(_gnum1_com_from_gnum_send, n_send_part2, "_gnum1_com_from_gnum_send ::");


  }

  int exch_request = -1;
  int         **pextract_edge_n    = NULL;
  PDM_g_num_t **pextract_edge_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_gnum_send,
                                 &pextract_edge_n,
                    (void ***)   &pextract_edge_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp_vtx, exch_request);

  int **pextract_edge_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int  **)  gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_triplet_send,
                                 &pextract_edge_n,
                    (void ***)   &pextract_edge_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp_vtx, exch_request);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(gnum1_com_from_triplet_n   [i_part]);
    PDM_free(gnum1_com_from_triplet_send[i_part]);
    PDM_free(gnum1_com_from_gnum_send   [i_part]);
  }
  PDM_free(gnum1_com_from_triplet_n   );
  PDM_free(gnum1_com_from_triplet_send);
  PDM_free(gnum1_com_from_gnum_send   );


  if(0 == 1) { // Usefull to know how many data is transfer
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      int n_triplet = pvtx_extented_to_pvtx_idx[i_part][pn_vtx[i_part]];
      PDM_log_trace_array_int(pvtx_extented_to_pvtx_idx       [i_part], pn_vtx[i_part]+1, "pvtx_extented_to_pvtx_idx ::");
      PDM_log_trace_array_int(pvtx_extented_to_pvtx_triplet   [i_part], n_triplet           , "pvtx_extented_to_pvtx_triplet ::");
      PDM_log_trace_array_int(pvtx_extented_to_pvtx_interface [i_part], n_triplet/3         , "pvtx_extented_to_pvtx_interface ::");
      PDM_log_trace_graph_nuplet_int(pvtx_extented_to_pvtx_idx[i_part],
                                     pvtx_extented_to_pvtx_triplet[i_part], 1,
                                     pn_vtx[i_part], "pvtx_extented_to_pvtx_triplet ::");


      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");

      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "extract_lnum2 :");
      PDM_log_trace_connectivity_long(gnum1_come_from_idx[i_part],
                                      gnum1_come_from    [i_part],
                                      n_ref_lnum2  [i_part], "gnum1_come_from ::");
    }
  }

  int **pextract_edge_idx = NULL;
  PDM_malloc(pextract_edge_idx, ln_part_tot, int *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int n_part1_to_part2 = pvtx_extented_to_pvtx_idx[i_part][pn_vtx[i_part]]/3;
    pextract_edge_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_edge_n[i_part], n_part1_to_part2);

    // PDM_log_trace_array_int (pextract_edge_idx[i_part], n_part1_to_part2+1, "pextract_edge_idx ::");
    // PDM_log_trace_array_long(pextract_edge_gnum[i_part], pextract_edge_idx[i_part][n_part1_to_part2], "pextract_edge_gnum ::");
  }


  PDM_part_to_part_free(ptp_vtx);

  /*
   * Post-treatment
   *   - Remove duplicate (same gnum or same equivalent transformation )
   *   - Keep link between edge (extented and current part)
   *   - Create new global numbering for all new entities
   */
  PDM_gen_gnum_t* gen_gnum_edge = PDM_gnum_create(3,
                                                  ln_part_tot,
                                                  PDM_TRUE,
                                                  1.e-6,
                                                  comm,
                                                  PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_edge, 2);

  int          *pn_edge_only_by_interface        = NULL;
  int         **pedge_interface                  = NULL;
  PDM_g_num_t **pedge_ln_to_gn_only_by_interface = NULL;
  PDM_malloc(pn_edge_only_by_interface       , ln_part_tot, int          );
  PDM_malloc(pedge_interface                 , ln_part_tot, int         *);
  PDM_malloc(pedge_ln_to_gn_only_by_interface, ln_part_tot, PDM_g_num_t *);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_edge_gnum    = pextract_edge_gnum   [i_part];
    int         *_pextract_edge_idx     = pextract_edge_idx    [i_part];

    int n_part1_to_part2          = pvtx_extented_to_pvtx_idx[i_part][pn_vtx[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_edge_idx[n_part1_to_part2];

    int n_part1_to_part2_edge = pedge_extented_to_pedge_idx[i_part][pn_edge[i_part]]/3;
    pn_edge_only_by_interface[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(pvtx_extented_to_pvtx_interface[i_part][i] != 0) {

        for(int j = _pextract_edge_idx[i]; j < _pextract_edge_idx[i+1]; ++j) {
          // Remove already gnum generate by interface
          PDM_g_num_t gnum_and_itrf[2] = {_pextract_edge_gnum[j], pvtx_extented_to_pvtx_interface[i_part][i]};
          int pos = PDM_order_binary_search_long(gnum_and_itrf,
                                                 prev_edge_ln_to_gn_and_interface[i_part],
                                                 2,
                                                 n_part1_to_part2_edge);
          log_trace(" ooooooo : %i \n", pos);
          if(pos == -1) {
            pn_edge_only_by_interface[i_part] += 1; // pextract_edge_n[i_part][i];
          }
        }
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_edge_gnum, n_part1_to_part2_recv_tot, "_pextract_edge_gnum ::");
    }

    PDM_malloc(pedge_ln_to_gn_only_by_interface[i_part], 2 * pn_edge_only_by_interface[i_part], PDM_g_num_t);
    PDM_malloc(pedge_interface                 [i_part],     n_part1_to_part2_recv_tot        , int        );
    PDM_g_num_t *_pedge_ln_to_gn_only_by_interface = pedge_ln_to_gn_only_by_interface[i_part];

    pn_edge_only_by_interface[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(pvtx_extented_to_pvtx_interface[i_part][i] != 0) {
        for(int j = _pextract_edge_idx[i]; j < _pextract_edge_idx[i+1]; ++j) {

          // Remove already gnum generate by interface
          PDM_g_num_t gnum_and_itrf[2] = {_pextract_edge_gnum[j], pvtx_extented_to_pvtx_interface[i_part][i]};
          int pos = PDM_order_binary_search_long(gnum_and_itrf,
                                                 prev_edge_ln_to_gn_and_interface[i_part],
                                                 2,
                                                 n_part1_to_part2_edge);
          if(pos == -1) {
            int idx_write = pn_edge_only_by_interface[i_part]++;
            _pedge_ln_to_gn_only_by_interface[2*idx_write  ] = _pextract_edge_gnum[j];
            _pedge_ln_to_gn_only_by_interface[2*idx_write+1] = pvtx_extented_to_pvtx_interface[i_part][i];
          }
        }
      }

      for(int j = _pextract_edge_idx[i]; j < _pextract_edge_idx[i+1]; ++j) {
        pedge_interface[i_part][j] = pvtx_extented_to_pvtx_interface[i_part][i];
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(pedge_ln_to_gn_only_by_interface[i_part], 2 * pn_edge_only_by_interface[i_part], "pedge_ln_to_gn_only_by_interface ::");
    }


    PDM_gnum_set_from_parents(gen_gnum_edge,
                              i_part,
                              pn_edge_only_by_interface[i_part],
                              pedge_ln_to_gn_only_by_interface[i_part]);

  }


  PDM_gnum_compute(gen_gnum_edge);

  /* Update */
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t* extented_edge_ln_to_gn = PDM_gnum_get(gen_gnum_edge, i_part);

    PDM_g_num_t *_pextract_edge_gnum    = pextract_edge_gnum   [i_part];
    int         *_pextract_edge_idx     = pextract_edge_idx    [i_part];

    pn_edge_only_by_interface[i_part] = 0;

    int n_part1_to_part2 = pvtx_extented_to_pvtx_idx[i_part][pn_vtx[i_part]]/3;
    int n_part1_to_part2_edge = pedge_extented_to_pedge_idx[i_part][pn_edge[i_part]]/3;

    int idx_read = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(pvtx_extented_to_pvtx_interface[i_part][i] != 0) {
        for(int j = _pextract_edge_idx[i]; j < _pextract_edge_idx[i+1]; ++j) {
          // Remove already gnum generate by interface
          PDM_g_num_t gnum_and_itrf[2] = {_pextract_edge_gnum[j], pvtx_extented_to_pvtx_interface[i_part][i]};
          int pos = PDM_order_binary_search_long(gnum_and_itrf,
                                                 prev_edge_ln_to_gn_and_interface[i_part],
                                                 2,
                                                 n_part1_to_part2_edge);
          if(pos == -1) {
            _pextract_edge_gnum[j] = extented_edge_ln_to_gn[idx_read++] + shift_edge_gnum;
          }
        }
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_edge_gnum, _pextract_edge_idx[n_part1_to_part2], "_pextract_edge_gnum (Update) : ");
    }

  }


  PDM_gnum_free(gen_gnum_edge);

  /*
   * Local post-treatment
   *   - Remove duplicate
   *   - Remove same entity but wih different path
   */
  int          *pn_edge_extented                       = NULL;
  int         **pnext_edge_extented_to_pedge_idx       = NULL;
  int         **pnextedge_extented_to_pedge_triplet    = NULL;
  PDM_g_num_t **pedge_extented_ln_to_gn                = NULL;
  PDM_g_num_t **extented_edge_orig_gnum                = NULL;
  int         **pnext_edge_extented_to_pedge_interface = NULL;
  PDM_malloc(pn_edge_extented                      , ln_part_tot, int          );
  PDM_malloc(pnext_edge_extented_to_pedge_idx      , ln_part_tot, int         *);
  PDM_malloc(pnextedge_extented_to_pedge_triplet   , ln_part_tot, int         *);
  PDM_malloc(pedge_extented_ln_to_gn               , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(extented_edge_orig_gnum               , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pnext_edge_extented_to_pedge_interface, ln_part_tot, int         *);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_edge_gnum = pextract_edge_gnum[i_part];
    int         *_pextract_edge_idx  = pextract_edge_idx [i_part];

    int n_part1_to_part2 = pvtx_extented_to_pvtx_idx[i_part][pn_vtx[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_edge_idx[n_part1_to_part2];
    int n_part1_to_part2_edge = pedge_extented_to_pedge_idx[i_part][pn_edge[i_part]]/3;


    // All gnum has be unified / shift w/r of interface, only sort along gnum
    int *order = NULL;
    PDM_malloc(order, n_part1_to_part2_recv_tot, int);
    for(int i_edge = 0; i_edge < n_part1_to_part2_recv_tot; ++i_edge) {
      order[i_edge] = i_edge;
    }

    // Maybe faire inplace sort and revert with order after ?
    PDM_g_num_t* sorted_pedge_ln_to_gn = NULL;
    PDM_malloc(sorted_pedge_ln_to_gn, pn_edge[i_part], PDM_g_num_t);
    for(int i = 0; i < pn_edge[i_part]; ++i) {
      sorted_pedge_ln_to_gn[i] = pedge_ln_to_gn[i_part][i];
    }
    PDM_sort_long(sorted_pedge_ln_to_gn, NULL, pn_edge[i_part]);

    // PDM_log_trace_array_long(sorted_pedge_ln_to_gn, pn_edge[i_part], "sorted_pedge_ln_to_gn ::");

    int n_unique = PDM_inplace_unique_long_and_order(_pextract_edge_gnum,
                                                     order,
                                                     0,
                                                     n_part1_to_part2_recv_tot-1);

    /*
     * The new edge can be :
     *  - In current partition
     *  - In interface array
     */

    /*
     * Solution 1 : On mix les numero globaux et les doublets
     * Solution 2 : On se fait un tableau trié complet avec les doublet
     *              --> On genere le gnum a la toute fin ? --> On peut pas car on a besoin de gnum pour le part_to_part
     */


    PDM_malloc(pedge_extented_ln_to_gn[i_part], n_unique, PDM_g_num_t);
    int n_unique2 = 0;
    for(int i = 0; i < n_unique; ++i) {
      int pos = PDM_binary_search_long(_pextract_edge_gnum[i], sorted_pedge_ln_to_gn, pn_edge[i_part]);

      int old_pos = order[i];
      PDM_g_num_t gnum_and_itrf[2] = {_pextract_edge_gnum[i], pedge_interface    [i_part][  old_pos  ]};
      int pos2 = PDM_order_binary_search_long(gnum_and_itrf,
                                              prev_edge_ln_to_gn_and_interface[i_part],
                                              2,
                                              n_part1_to_part2_edge);

      if(pos == -1 && pos2 == -1) {
        pedge_extented_ln_to_gn[i_part][n_unique2] = _pextract_edge_gnum[i];
        n_unique2++;
      }
    }
    pn_edge_extented[i_part] = n_unique2;

    // Extract all data
    PDM_malloc(pnext_edge_extented_to_pedge_idx      [i_part],     n_unique2 + 1, int        );
    PDM_malloc(pnextedge_extented_to_pedge_triplet   [i_part], 3 * n_unique2    , int        );
    PDM_malloc(extented_edge_orig_gnum               [i_part],     n_unique2    , PDM_g_num_t);
    PDM_malloc(pnext_edge_extented_to_pedge_interface[i_part],     n_unique2    , int        );

    /* Count */
    pnext_edge_extented_to_pedge_idx    [i_part][0] = 0;
    int idx_write = 0;
    for(int i_edge = 0; i_edge < n_unique; ++i_edge) {
      int pos = PDM_binary_search_long(_pextract_edge_gnum[i_edge], sorted_pedge_ln_to_gn, pn_edge[i_part]);

      int old_pos = order[i_edge];
      PDM_g_num_t gnum_and_itrf[2] = {_pextract_edge_gnum[i_edge], pedge_interface    [i_part][  old_pos  ]};
      int pos2 = PDM_order_binary_search_long(gnum_and_itrf,
                                              prev_edge_ln_to_gn_and_interface[i_part],
                                              2,
                                              n_part1_to_part2_edge);

      if(pos == -1 && pos2 == -1) {
        pnext_edge_extented_to_pedge_idx[i_part][idx_write+1] = pnext_edge_extented_to_pedge_idx[i_part][idx_write] + 3;

        pnextedge_extented_to_pedge_triplet[i_part][3*idx_write  ] = pextract_edge_triplet[i_part][3*old_pos  ];
        pnextedge_extented_to_pedge_triplet[i_part][3*idx_write+1] = pextract_edge_triplet[i_part][3*old_pos+1];
        pnextedge_extented_to_pedge_triplet[i_part][3*idx_write+2] = pextract_edge_triplet[i_part][3*old_pos+2];

        // Save the link (gnum, interface)
        pnext_edge_extented_to_pedge_interface[i_part][idx_write] = pedge_interface    [i_part][  old_pos  ];
        extented_edge_orig_gnum               [i_part][idx_write] = _pextract_edge_gnum        [  i_edge   ]; // Because it's sorted already
        idx_write++;
      }
    }
    assert(idx_write == n_unique2);
    n_unique = n_unique2;

    int n_triplet = pnext_edge_extented_to_pedge_idx[i_part][n_unique];

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_edge_gnum                           , n_unique   , "_pextract_edge_gnum (UNIQUE)        ::");
      PDM_log_trace_array_int (pnext_edge_extented_to_pedge_idx      [i_part], n_unique+1 , "pnext_edge_extented_to_pedge_idx      ::");
      PDM_log_trace_array_long(extented_edge_orig_gnum               [i_part], n_triplet/3, "extented_edge_orig_gnum             ::");
      PDM_log_trace_array_int (pnext_edge_extented_to_pedge_interface[i_part], n_triplet/3, "pnext_edge_extented_to_pedge_interface::");
      PDM_log_trace_array_int (pnextedge_extented_to_pedge_triplet   [i_part], n_triplet  , "pnextedge_extented_to_pedge_triplet  ::");
    }

    PDM_free(order);
    PDM_free(sorted_pedge_ln_to_gn);

  }


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(prev_edge_ln_to_gn_and_interface[i_part]);
  }
  PDM_free(prev_edge_ln_to_gn_and_interface);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pedge_interface                 [i_part]);
    PDM_free(pedge_ln_to_gn_only_by_interface[i_part]);
  }
  PDM_free(pn_edge_only_by_interface);
  PDM_free(pedge_interface);
  PDM_free(pedge_ln_to_gn_only_by_interface);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pnext_edge_extented_to_pedge_idx      [i_part]);
    PDM_free(pnextedge_extented_to_pedge_triplet   [i_part]);
    PDM_free(extented_edge_orig_gnum               [i_part]);
    PDM_free(pnext_edge_extented_to_pedge_interface[i_part]);
    PDM_free(pedge_extented_ln_to_gn               [i_part]);
  }
  PDM_free(pnext_edge_extented_to_pedge_idx      );
  PDM_free(pnextedge_extented_to_pedge_triplet   );
  PDM_free(extented_edge_orig_gnum               );
  PDM_free(pedge_extented_ln_to_gn               );
  PDM_free(pnext_edge_extented_to_pedge_interface);
  PDM_free(pn_edge_extented);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pextract_edge_idx[i_part]);
  }
  PDM_free(pextract_edge_idx);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pextract_edge_n      [i_part]);
    PDM_free(pextract_edge_gnum   [i_part]);
    PDM_free(pextract_edge_triplet[i_part]);
  }
  PDM_free(pextract_edge_n      );
  PDM_free(pextract_edge_gnum   );
  PDM_free(pextract_edge_triplet);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pvtx_edge_idx[i_part]);
    PDM_free(pvtx_edge    [i_part]);
  }


  PDM_free(pvtx_edge_idx);
  PDM_free(pvtx_edge);



}


static
void
_part_extension
(
  int                          n_depth,
  int                          n_domain,
  int*                         n_part,
  PDM_MPI_Comm                 comm,
  PDM_multipart_t             *mpart,
  PDM_part_domain_interface_t *pdi
)
{
  PDM_UNUSED(n_depth);
  PDM_UNUSED(mpart);
  PDM_UNUSED(pdi);

  // printf("n_domain = %i \n", n_domain);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int* n_part_g = NULL;
  PDM_malloc(n_part_g, n_domain, int);
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      ln_part_tot += 1;
    }
  }



  int          **pn_vtx              = NULL;
  int          **pn_edge             = NULL;
  PDM_g_num_t ***pvtx_ln_to_gn       = NULL;
  PDM_g_num_t ***pedge_ln_to_gn      = NULL;
  double      ***pvtx_coords         = NULL;
  PDM_malloc(pn_vtx        , n_domain, int          *);
  PDM_malloc(pn_edge       , n_domain, int          *);
  PDM_malloc(pvtx_ln_to_gn , n_domain, PDM_g_num_t **);
  PDM_malloc(pedge_ln_to_gn, n_domain, PDM_g_num_t **);
  PDM_malloc(pvtx_coords   , n_domain, double      **);

  int           *pflat_n_vtx         = NULL;
  int           *pflat_n_edge        = NULL;
  int         ***pedge_vtx_idx       = NULL;
  int         ***pedge_vtx           = NULL;
  int          **pflat_edge_vtx      = NULL;
  PDM_g_num_t  **pflat_vtx_ln_to_gn  = NULL;
  PDM_g_num_t  **pflat_edge_ln_to_gn = NULL;
  double       **pflat_vtx_coords    = NULL;
  PDM_malloc(pflat_n_vtx        , n_domain   , int           );
  PDM_malloc(pflat_n_edge       , n_domain   , int           );
  PDM_malloc(pedge_vtx_idx      , n_domain   , int         **);
  PDM_malloc(pedge_vtx          , n_domain   , int         **);
  PDM_malloc(pflat_edge_vtx     , ln_part_tot, int          *);
  PDM_malloc(pflat_vtx_ln_to_gn , ln_part_tot, PDM_g_num_t  *);
  PDM_malloc(pflat_edge_ln_to_gn, ln_part_tot, PDM_g_num_t  *);
  PDM_malloc(pflat_vtx_coords   , ln_part_tot, double       *);
  ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    PDM_malloc(pn_vtx        [i_dom], n_part[i_dom], int          );
    PDM_malloc(pvtx_ln_to_gn [i_dom], n_part[i_dom], PDM_g_num_t *);
    PDM_malloc(pn_edge       [i_dom], n_part[i_dom], int          );
    PDM_malloc(pedge_ln_to_gn[i_dom], n_part[i_dom], PDM_g_num_t *);
    PDM_malloc(pvtx_coords   [i_dom], n_part[i_dom], double      *);

    PDM_malloc(pedge_vtx_idx [i_dom], n_part[i_dom], int         *);
    PDM_malloc(pedge_vtx     [i_dom], n_part[i_dom], int         *);

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VTX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);

      pflat_n_vtx       [ln_part_tot+i_part] = pn_vtx       [i_dom][i_part];
      pflat_vtx_ln_to_gn[ln_part_tot+i_part] = pvtx_ln_to_gn[i_dom][i_part];

      pflat_n_edge      [ln_part_tot+i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                                               i_dom,
                                                                               i_part,
                                                                               PDM_MESH_ENTITY_EDGE,
                                                                               &pflat_edge_ln_to_gn[ln_part_tot+i_part],
                                                                               PDM_OWNERSHIP_KEEP);

      pn_edge       [i_dom][i_part] = pflat_n_edge       [ln_part_tot+i_part];
      pedge_ln_to_gn[i_dom][i_part] = pflat_edge_ln_to_gn[ln_part_tot+i_part];

      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &pedge_vtx_idx[i_dom][i_part],
                                          &pedge_vtx    [i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);

      assert(pedge_vtx_idx[i_dom][i_part] == NULL);
      PDM_malloc(pedge_vtx_idx[i_dom][i_part], pn_edge[i_dom][i_part] + 1, int);
      for(int i_edge = 0; i_edge < pn_edge[i_dom][i_part]+1; ++i_edge) {
        pedge_vtx_idx[i_dom][i_part][i_edge] = 2 * i_edge;
      }
      pflat_edge_vtx[ln_part_tot+i_part] = pedge_vtx    [i_dom][i_part];

      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pflat_vtx_coords[ln_part_tot+i_part],
                                       PDM_OWNERSHIP_KEEP);

      pvtx_coords[i_dom][i_part] = pflat_vtx_coords[ln_part_tot+i_part];

    }
    ln_part_tot += n_part[i_dom];
  }


  PDM_g_num_t* shift_by_domain_vtx = PDM_compute_offset_ln_to_gn_by_domain(n_domain,
                                                                           n_part,
                                                                           pn_vtx,
                                                                           pvtx_ln_to_gn,
                                                                           comm);

  PDM_g_num_t* shift_by_domain_edge = PDM_compute_offset_ln_to_gn_by_domain(n_domain,
                                                                            n_part,
                                                                            pn_edge,
                                                                            pedge_ln_to_gn,
                                                                            comm);

  /* Shift ln_to_gn */
  PDM_offset_ln_to_gn_by_domain(n_domain,
                                n_part,
                                pn_vtx,
                                pvtx_ln_to_gn,
                                shift_by_domain_vtx,
                                1);

  PDM_offset_ln_to_gn_by_domain(n_domain,
                                n_part,
                                pn_edge,
                                pedge_ln_to_gn,
                                shift_by_domain_edge,
                                1);


  /*
   * 1er etape :
   *   - Construire a partir des infos user, une première partition étendu qui va initié une recursion
   *   - On concatene les conectivités de depart et la première couche
   *   - On conserve le lien entre le bout des nouvelles partition étendu et le maillage de base pour preparer le step d'après
   *   - Il faut l'info pour chaque entité
   */
  int          *pn_edge_extented                  = NULL;
  int         **pedge_extented_to_pedge_idx       = NULL;
  int         **pedge_extented_to_pedge_triplet   = NULL;
  PDM_g_num_t **pedge_extented_ln_to_gn           = NULL;
  int         **pedge_extented_to_pedge_interface = NULL;

  int          *pn_vtx_extented                 = NULL;
  PDM_g_num_t **pvtx_extented_ln_to_gn          = NULL;
  int         **pextented_edge_vtx_idx          = NULL;
  int         **pextented_edge_vtx              = NULL;
  int         **pvtx_extented_to_pvtx_idx       = NULL;
  int         **pvtx_extented_to_pvtx_triplet   = NULL;
  int         **pvtx_extented_to_pvtx_interface = NULL;


  /*
   * Concatenate partition
   */
  int *pn_concat_edge = NULL;
  int *pn_concat_vtx  = NULL;

  PDM_g_num_t **pconcat_edge_ln_to_gn = NULL;
  PDM_g_num_t **pconcat_vtx_ln_to_gn  = NULL;
  int         **pconcat_edge_vtx      = NULL;
  int         **pconcat_edge_vtx_idx  = NULL;
  double      **pconcat_vtx_coords    = NULL;

  // int **edge_kind = NULL;
  int **pconcat_pvtx_extented_to_pvtx_idx         = NULL;
  int **pconcat_pvtx_extented_to_pvtx_triplet     = NULL;
  int **pconcat_pvtx_extented_to_pvtx_interface   = NULL;
  int **pconcat_pedge_extented_to_pedge_idx       = NULL;
  int **pconcat_pedge_extented_to_pedge_triplet   = NULL;
  int **pconcat_pedge_extented_to_pedge_interface = NULL;


  /*
   * Init all interface information in block frame to hook easily all data
   * And for all meshes entities
   * One block per interface and by mesh_entities
   * Il faut une version qui update : PDM_domain_interface_make_flat_view
   * Micro-structure =
   *    - n_interface
   *    - part_to_block [PDM_MESH_ENTITY_MAX][n_interface]
   *    - block_data_opp[PDM_MESH_ENTITY_MAX][n_interface]
   *    Une methode de check
   *    Une methode init from domaine_interface / part_domain_interface
   *    Une methode append qui prends l'existant
   *       - Existant = 1 partition
   *       - Nouveau  = 1 nouvelle partition
   *    --> Recreation du ptb + block_data associé
   */




  _part_extension_init(comm,
                       n_depth,
                       pdi,
                       n_domain,
                       shift_by_domain_edge,
                       shift_by_domain_vtx,
                       n_part,
                       pn_vtx,
                       pvtx_ln_to_gn,
                       pvtx_coords,
                       pn_edge,
                       pedge_ln_to_gn,
                       pedge_vtx_idx,
                       pedge_vtx,
                       &pn_edge_extented,
                       &pedge_extented_to_pedge_idx,
                       &pedge_extented_to_pedge_triplet,
                       &pedge_extented_ln_to_gn,
                       &pedge_extented_to_pedge_interface,
                       &pn_vtx_extented,
                       &pvtx_extented_ln_to_gn,
                       &pextented_edge_vtx_idx,
                       &pextented_edge_vtx,
                       &pvtx_extented_to_pvtx_idx,
                       &pvtx_extented_to_pvtx_triplet,
                       &pvtx_extented_to_pvtx_interface,
                       &pn_concat_vtx,
                       &pconcat_vtx_ln_to_gn,
                       &pn_concat_edge,
                       &pconcat_edge_ln_to_gn,
                       &pconcat_edge_vtx,
                       &pconcat_edge_vtx_idx,
                       &pconcat_vtx_coords,
                       &pconcat_pvtx_extented_to_pvtx_idx,
                       &pconcat_pvtx_extented_to_pvtx_triplet,
                       &pconcat_pvtx_extented_to_pvtx_interface,
                       &pconcat_pedge_extented_to_pedge_idx,
                       &pconcat_pedge_extented_to_pedge_triplet,
                       &pconcat_pedge_extented_to_pedge_interface);

  /*
   * At this stage we have a first partition extension but we want to recurse (Realloc inside)
   */
  _part_extension_one_depth(comm,
                            pdi,
                            n_domain,
                            shift_by_domain_edge,
                            shift_by_domain_vtx,
                            n_part,
                            &pn_concat_vtx,
                            &pconcat_vtx_ln_to_gn,
                            &pn_concat_edge,
                            &pconcat_edge_ln_to_gn,
                            &pconcat_edge_vtx,
                            &pconcat_edge_vtx_idx,
                            &pconcat_vtx_coords,
                            &pconcat_pvtx_extented_to_pvtx_idx,
                            &pconcat_pvtx_extented_to_pvtx_triplet,
                            &pconcat_pvtx_extented_to_pvtx_interface,
                            &pconcat_pedge_extented_to_pedge_idx,
                            &pconcat_pedge_extented_to_pedge_triplet,
                            &pconcat_pedge_extented_to_pedge_interface);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pconcat_edge_ln_to_gn[i_part]);
    PDM_free(pconcat_vtx_ln_to_gn [i_part]);
    PDM_free(pconcat_edge_vtx     [i_part]);
    PDM_free(pconcat_edge_vtx_idx [i_part]);
    PDM_free(pconcat_vtx_coords   [i_part]);

    PDM_free(pconcat_pvtx_extented_to_pvtx_idx        [i_part]);
    PDM_free(pconcat_pvtx_extented_to_pvtx_triplet    [i_part]);
    PDM_free(pconcat_pvtx_extented_to_pvtx_interface  [i_part]);
    PDM_free(pconcat_pedge_extented_to_pedge_idx      [i_part]);
    PDM_free(pconcat_pedge_extented_to_pedge_triplet  [i_part]);
    PDM_free(pconcat_pedge_extented_to_pedge_interface[i_part]);
  }
  PDM_free(pn_concat_edge);
  PDM_free(pn_concat_vtx );

  PDM_free(pconcat_edge_ln_to_gn);
  PDM_free(pconcat_vtx_ln_to_gn );
  PDM_free(pconcat_edge_vtx     );
  PDM_free(pconcat_edge_vtx_idx );
  PDM_free(pconcat_vtx_coords   );

  PDM_free(pconcat_pvtx_extented_to_pvtx_idx        );
  PDM_free(pconcat_pvtx_extented_to_pvtx_triplet    );
  PDM_free(pconcat_pvtx_extented_to_pvtx_interface  );
  PDM_free(pconcat_pedge_extented_to_pedge_idx      );
  PDM_free(pconcat_pedge_extented_to_pedge_triplet  );
  PDM_free(pconcat_pedge_extented_to_pedge_interface);


  /* Unshift ln_to_gn */
  PDM_offset_ln_to_gn_by_domain(n_domain,
                                n_part,
                                pn_vtx,
                                pvtx_ln_to_gn,
                                shift_by_domain_vtx,
                                -1);

  PDM_offset_ln_to_gn_by_domain(n_domain,
                                n_part,
                                pn_edge,
                                pedge_ln_to_gn,
                                shift_by_domain_edge,
                                -1);


  /* Free */

  // for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
  //   PDM_free(pextract_vtx_coords[i_part]);
  // }
  // PDM_free(pextract_vtx_coords);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pedge_extented_to_pedge_idx      [i_part]);
    PDM_free(pedge_extented_to_pedge_triplet  [i_part]);
    PDM_free(pedge_extented_ln_to_gn          [i_part]);
    PDM_free(pedge_extented_to_pedge_interface[i_part]);
    PDM_free(pvtx_extented_ln_to_gn           [i_part]);
    PDM_free(pextented_edge_vtx_idx           [i_part]);
    PDM_free(pextented_edge_vtx               [i_part]);
    PDM_free(pvtx_extented_to_pvtx_idx        [i_part]);
    PDM_free(pvtx_extented_to_pvtx_triplet    [i_part]);
    PDM_free(pvtx_extented_to_pvtx_interface  [i_part]);
  }

  PDM_free(pedge_extented_to_pedge_idx      );
  PDM_free(pedge_extented_to_pedge_triplet  );
  PDM_free(pedge_extented_ln_to_gn          );
  PDM_free(pedge_extented_to_pedge_interface);
  PDM_free(pvtx_extented_ln_to_gn           );
  PDM_free(pextented_edge_vtx_idx           );
  PDM_free(pextented_edge_vtx               );
  PDM_free(pvtx_extented_to_pvtx_idx        );
  PDM_free(pvtx_extented_to_pvtx_triplet    );
  PDM_free(pvtx_extented_to_pvtx_interface  );

  PDM_free(pn_edge_extented);
  PDM_free(pn_vtx_extented);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      PDM_free(pedge_vtx_idx[i_dom][i_part]);
    }
    PDM_free(pn_vtx        [i_dom]);
    PDM_free(pvtx_ln_to_gn [i_dom]);
    PDM_free(pn_edge       [i_dom]);
    PDM_free(pedge_ln_to_gn[i_dom]);
    PDM_free(pvtx_coords   [i_dom]);

    PDM_free(pedge_vtx_idx [i_dom]);
    PDM_free(pedge_vtx     [i_dom]);
  }

  PDM_free(pn_vtx             );
  PDM_free(pn_edge            );
  PDM_free(pvtx_ln_to_gn      );
  PDM_free(pedge_ln_to_gn     );
  PDM_free(pflat_n_vtx        );
  PDM_free(pflat_n_edge       );
  PDM_free(pedge_vtx_idx      );
  PDM_free(pedge_vtx          );
  PDM_free(pvtx_coords        );
  PDM_free(pflat_edge_ln_to_gn);
  PDM_free(pflat_vtx_ln_to_gn );
  PDM_free(pflat_vtx_coords   );
  PDM_free(pflat_edge_vtx     );

  PDM_free(shift_by_domain_vtx);
  PDM_free(shift_by_domain_edge);

  PDM_free(n_part_g);


}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */
// mpirun -np 1 ./test/pdm_t_part_extension_1d_propagation -n 5 -depth 1 -pi
int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t n_g_pts   = 10;
  int         n_dom_i    = 1;
  int         periodic_i = 0;
  int         n_depth    = 1;
  _read_args(argc,
             argv,
             &n_g_pts,
             &n_dom_i,
             &periodic_i,
             &n_depth);

  double      **dvtx_coord   = NULL;
  PDM_g_num_t **distrib_edge = NULL;
  PDM_g_num_t **distrib_vtx  = NULL;
  PDM_g_num_t **dedge_vtx    = NULL;

  PDM_domain_interface_t *dom_itrf = NULL;
  PDM_generate_cart_topo_lines(comm,
                               n_dom_i,
                               periodic_i,
                               0.,
                               0.,
                               0.,
                               1.,
                               n_g_pts,
                               &distrib_edge,
                               &distrib_vtx,
                               &dedge_vtx,
                               &dvtx_coord,
                               &dom_itrf);

  PDM_dmesh_t **dm = NULL;
  PDM_malloc(dm, n_dom_i, PDM_dmesh_t *);

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {

    int dn_edge = distrib_edge[i_dom][i_rank+1] - distrib_edge[i_dom][i_rank];
    int dn_vtx  = distrib_vtx [i_dom][i_rank+1] - distrib_vtx [i_dom][i_rank];

    /*
     * Create dmesh
     */
    dm[i_dom] = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                 0,
                                 0,
                                 dn_edge,
                                 dn_vtx,
                                 comm);

    PDM_dmesh_connectivity_set(dm[i_dom],
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               dedge_vtx[i_dom],
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_vtx_coord_set(dm[i_dom],
                            dvtx_coord[i_dom],
                            PDM_OWNERSHIP_USER);
  }


  /*
   * Mulitpart
   */
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_IMPLICIT;
  // PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
  int* n_part = NULL;
  PDM_malloc(n_part, n_dom_i, int);
  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {
    n_part[i_dom] = 1;
  }
  PDM_multipart_t* mpart = PDM_multipart_create(n_dom_i,
                                                n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {
    PDM_multipart_dmesh_set(mpart, i_dom, dm[i_dom]);
  }

  PDM_multipart_compute(mpart);

  int n_domain = n_dom_i;

  int          **pn_vtx         = NULL;
  PDM_g_num_t ***pvtx_ln_to_gn  = NULL;
  PDM_malloc(pn_vtx       , n_domain, int          *);
  PDM_malloc(pvtx_ln_to_gn, n_domain, PDM_g_num_t **);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    PDM_malloc(pn_vtx       [i_dom], n_part[i_dom], int          );
    PDM_malloc(pvtx_ln_to_gn[i_dom], n_part[i_dom], PDM_g_num_t *);

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VTX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);



    }
  }

  if(0 == 1) {
    for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {
      for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

        PDM_g_num_t* lpvtx_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &lpvtx_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t* pedge_ln_to_gn = NULL;
        int pn_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      i_dom,
                                                      i_part,
                                                      PDM_MESH_ENTITY_EDGE,
                                                      &pedge_ln_to_gn,
                                                      PDM_OWNERSHIP_KEEP);

        int *pedge_vtx     = NULL;
        int *pedge_vtx_idx = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            i_dom,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &pedge_vtx_idx,
                                            &pedge_vtx,
                                            PDM_OWNERSHIP_KEEP);

        double *pvtx_coord = NULL;
        int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     i_dom,
                                                     i_part,
                                                     &pvtx_coord,
                                                     PDM_OWNERSHIP_KEEP);
        char filename[999];
        sprintf(filename, "out_part_vtx_i_dom=%i_%i.vtk", i_dom, i_rank);
        PDM_vtk_write_std_elements (filename,
                                    n_vtx,
                                    pvtx_coord,
                                    lpvtx_ln_to_gn,
                                    PDM_MESH_NODAL_BAR2,
                                    pn_edge,
                                    pedge_vtx,
                                    pedge_ln_to_gn,
                                    0,
                                    NULL,
                                    NULL);

      }
    }
  }

  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_itrf,
                                                                                   n_part,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pn_vtx,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pvtx_ln_to_gn);

  PDM_domain_interface_free(dom_itrf);


  /*
   * Extention de partition
   *  - Step 1 : n_dom = 1, n_proc = 2 + rank1
   *  - Step 2 : n_dom = 1, n_proc = 2 + rank2
   *  - Step 3 : n_dom = 1, n_proc = 1 + rank1 + perio
   *  - Step 4 : n_dom = 1, n_proc = 1 + rank2 + perio
   *  - Step 5 : n_dom = 1, n_proc = 1 + rank1 + perio mais qui se recouvre plusieurs fois (genre 2 cellules)
   *  - Step 6 : n_dom = 2, n_proc = 1 + rank1
   *  - Step 7 : n_dom = 2, n_proc = 2 + rank1
   *  - Step 8 : n_dom = 2, n_proc = 2 + rank1 + perio
   *  - Step 9 : n_dom = 2, n_proc = 2 + rank1 + perio mais qui se recouvre plusieurs fois (genre 2 cellules)
   *
   */
  _part_extension(n_depth,
                  n_dom_i,
                  n_part,
                  comm,
                  mpart,
                  pdi);

  PDM_multipart_free(mpart);
  PDM_free(n_part);

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {

    PDM_free (dvtx_coord  [i_dom]);
    PDM_free (distrib_vtx [i_dom]);
    PDM_free (distrib_edge[i_dom]);
    PDM_free (dedge_vtx   [i_dom]);
    PDM_dmesh_free(dm[i_dom]);
  }
  PDM_free (dvtx_coord);
  PDM_free (distrib_vtx);
  PDM_free (distrib_edge);
  PDM_free (dedge_vtx);
  PDM_free (dm);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    PDM_free(pn_vtx       [i_dom]);
    PDM_free(pvtx_ln_to_gn[i_dom]);
  }
  PDM_free(pn_vtx);
  PDM_free(pvtx_ln_to_gn);

  PDM_part_domain_interface_free(pdi);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
