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

  int* n_part_g = malloc(n_domain * sizeof(int));
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      ln_part_tot += 1;
    }
  }



  int          **pn_vtx              = (int          **) malloc( n_domain    * sizeof(int          *));
  int          **pn_edge             = (int          **) malloc( n_domain    * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn       = (PDM_g_num_t ***) malloc( n_domain    * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn      = (PDM_g_num_t ***) malloc( n_domain    * sizeof(PDM_g_num_t **));
  int           *pflat_n_vtx         = (int           *) malloc( n_domain    * sizeof(int           ));
  int           *pflat_n_edge        = (int           *) malloc( n_domain    * sizeof(int           ));
  int         ***pedge_vtx_idx       = (int         ***) malloc( n_domain    * sizeof(int         **));
  int         ***pedge_vtx           = (int         ***) malloc( n_domain    * sizeof(int         **));
  int          **pflat_edge_vtx      = (int          **) malloc( ln_part_tot * sizeof(int          *));
  PDM_g_num_t  **pflat_vtx_ln_to_gn  = (PDM_g_num_t  **) malloc( ln_part_tot * sizeof(PDM_g_num_t  *));
  PDM_g_num_t  **pflat_edge_ln_to_gn = (PDM_g_num_t  **) malloc( ln_part_tot * sizeof(PDM_g_num_t  *));
  double       **pflat_vtx_coords    = (double       **) malloc( ln_part_tot * sizeof(double       *));

  ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx        [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn [i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));
    pn_edge       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pedge_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

    pedge_vtx_idx [i_dom]  = (int         **) malloc( n_part[i_dom] * sizeof(int         *));
    pedge_vtx     [i_dom]  = (int         **) malloc( n_part[i_dom] * sizeof(int         *));

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
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
                                          &pedge_vtx    [i_dom][i_part],
                                          &pedge_vtx_idx[i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);

      assert(pedge_vtx_idx[i_dom][i_part] == NULL);
      pedge_vtx_idx[i_dom][i_part] = malloc((pn_edge[i_dom][i_part] + 1) * sizeof(int));
      for(int i_edge = 0; i_edge < pn_edge[i_dom][i_part]+1; ++i_edge) {
        pedge_vtx_idx[i_dom][i_part][i_edge] = 2 * i_edge;
      }
      pflat_edge_vtx[ln_part_tot+i_part] = pedge_vtx    [i_dom][i_part];

      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pflat_vtx_coords[ln_part_tot+i_part],
                                       PDM_OWNERSHIP_KEEP);

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


  if(1 == 1) {
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
  PDM_part_to_part_t* ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pvtx_extented_ln_to_gn,
                                                                          (const int          *) pn_vtx_extented,
                                                                          ln_part_tot,
                                                                          (const int          *) pflat_n_vtx,
                                                                          ln_part_tot,
                                                                          (const int         **) pvtx_extented_to_pvtx_idx,
                                                                          (const int         **) NULL,
                                                                          (const int         **) pvtx_extented_to_pvtx_triplet,
                                                                          comm);

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
  double  **translation_vector = malloc(n_interface * sizeof(double *  ));
  double ***rotation_matrix    = malloc(n_interface * sizeof(double ** ));
  double  **rotation_direction = malloc(n_interface * sizeof(double *  ));
  double  **rotation_center    = malloc(n_interface * sizeof(double *  ));
  double   *rotation_angle     = malloc(n_interface * sizeof(double    ));
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
      rotation_matrix[i_interf] = malloc(3 * sizeof(double *));
      for(int k = 0; k < 3; ++k) {
        rotation_matrix[i_interf][k] = malloc(3 * sizeof(double));
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
      free(translation_vector[i_interf]);
    }
    if(rotation_center    [i_interf] != NULL) {
      for(int k = 0; k < 3; ++k) {
        free(rotation_matrix[i_interf][k]);
      }
      free(rotation_matrix[i_interf]);
    }
  }
  free(translation_vector);
  free(rotation_matrix);
  free(rotation_direction);
  free(rotation_center);
  free(rotation_angle);


  /*
   * Export vtk
   */
  if(1 == 1) {

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      int pn_concat_vtx  = pflat_n_vtx [i_part] + pn_vtx_extented [i_part];
      int pn_concat_edge = pflat_n_edge[i_part] + pn_edge_extented[i_part];

      int *edge_kind = malloc(pn_concat_edge * sizeof(int));

      int         *concat_edge_vtx      = malloc(2 * pn_concat_edge * sizeof(int         ));
      PDM_g_num_t *concat_edge_ln_to_gn = malloc(    pn_concat_edge * sizeof(PDM_g_num_t ));
      double      *concat_vtx_coord     = malloc(3 * pn_concat_vtx  * sizeof(double      ));
      PDM_g_num_t *concat_vtx_ln_to_gn  = malloc(    pn_concat_vtx  * sizeof(PDM_g_num_t ));

      for(int i_edge = 0; i_edge < pflat_n_edge[i_part]; ++i_edge) {
        concat_edge_vtx[2*i_edge  ] = pflat_edge_vtx[i_part][2*i_edge  ];
        concat_edge_vtx[2*i_edge+1] = pflat_edge_vtx[i_part][2*i_edge+1];
        concat_edge_ln_to_gn[i_edge] = pflat_edge_ln_to_gn[i_part][i_edge];
        edge_kind[i_edge] = 0;
      }

      for(int i_vtx = 0; i_vtx < pflat_n_vtx[i_part]; ++i_vtx) {
        concat_vtx_coord[3*i_vtx  ] = pflat_vtx_coords[i_part][3*i_vtx  ];
        concat_vtx_coord[3*i_vtx+1] = pflat_vtx_coords[i_part][3*i_vtx+1];
        concat_vtx_coord[3*i_vtx+2] = pflat_vtx_coords[i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[i_vtx] = pflat_vtx_ln_to_gn[i_part][i_vtx];
      }


      for(int i_edge = 0; i_edge < pn_edge_extented[i_part]; ++i_edge) {
        int idx_write = pflat_n_edge[i_part]+i_edge;
        concat_edge_vtx     [2*idx_write  ] = pextented_edge_vtx[i_part][2*i_edge  ];
        concat_edge_vtx     [2*idx_write+1] = pextented_edge_vtx[i_part][2*i_edge+1];
        concat_edge_ln_to_gn[  idx_write  ] = pedge_extented_ln_to_gn[i_part][i_edge];
        edge_kind           [idx_write] = 1;
      }

      for(int i_vtx = 0; i_vtx < pn_vtx_extented [i_part]; ++i_vtx) {
        int idx_write = pflat_n_vtx[i_part]+i_vtx;
        concat_vtx_coord   [3*idx_write  ] = pextract_vtx_coords   [i_part][3*i_vtx  ];
        concat_vtx_coord   [3*idx_write+1] = pextract_vtx_coords   [i_part][3*i_vtx+1];
        concat_vtx_coord   [3*idx_write+2] = pextract_vtx_coords   [i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[  idx_write  ] = pvtx_extented_ln_to_gn[i_part][  i_vtx  ];
      }

      char filename[999];
      sprintf(filename, "out_part_concate_vtx_i_part=%i_%i.vtk", i_part, i_rank);
      const char* field_name[] = {"edge_kind", 0 };
      int *field[1] = {edge_kind};
      PDM_vtk_write_std_elements (filename,
                                  pn_concat_vtx,
                                  concat_vtx_coord,
                                  concat_vtx_ln_to_gn,
                                  PDM_MESH_NODAL_BAR2,
                                  pn_concat_edge,
                                  concat_edge_vtx,
                                  concat_edge_ln_to_gn,
                                  1,
                                  field_name,
                                  (const int **) field);

      free(edge_kind);

      free(concat_edge_vtx);
      free(concat_edge_ln_to_gn);
      free(concat_vtx_coord);
      free(concat_vtx_ln_to_gn);

    }
  }


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

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_vtx_coords[i_part]);
  }
  free(pextract_vtx_coords);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pedge_extented_to_pedge_idx      [i_part]);
    free(pedge_extented_to_pedge_triplet  [i_part]);
    free(pedge_extented_ln_to_gn          [i_part]);
    free(pedge_extented_to_pedge_interface[i_part]);
    free(pvtx_extented_ln_to_gn           [i_part]);
    free(pextented_edge_vtx_idx           [i_part]);
    free(pextented_edge_vtx               [i_part]);
    free(pvtx_extented_to_pvtx_idx        [i_part]);
    free(pvtx_extented_to_pvtx_triplet    [i_part]);
    free(pvtx_extented_to_pvtx_interface  [i_part]);
  }

  free(pedge_extented_to_pedge_idx      );
  free(pedge_extented_to_pedge_triplet  );
  free(pedge_extented_ln_to_gn          );
  free(pedge_extented_to_pedge_interface);
  free(pvtx_extented_ln_to_gn           );
  free(pextented_edge_vtx_idx           );
  free(pextented_edge_vtx               );
  free(pvtx_extented_to_pvtx_idx        );
  free(pvtx_extented_to_pvtx_triplet    );
  free(pvtx_extented_to_pvtx_interface  );

  free(pn_edge_extented);
  free(pn_vtx_extented);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      free(pedge_vtx_idx[i_dom][i_part]);
    }
    free(pn_vtx        [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
    free(pn_edge       [i_dom]);
    free(pedge_ln_to_gn[i_dom]);

    free(pedge_vtx_idx [i_dom]);
    free(pedge_vtx     [i_dom]);
  }

  free(pn_vtx             );
  free(pn_edge            );
  free(pvtx_ln_to_gn      );
  free(pedge_ln_to_gn     );
  free(pflat_n_vtx        );
  free(pflat_n_edge       );
  free(pedge_vtx_idx      );
  free(pedge_vtx          );
  free(pflat_edge_ln_to_gn);
  free(pflat_vtx_ln_to_gn );
  free(pflat_vtx_coords   );
  free(pflat_edge_vtx     );

  free(shift_by_domain_vtx);
  free(shift_by_domain_edge);

  free(n_part_g);


}

static
void
_part_extension_old
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

  int* n_part_g = malloc(n_domain * sizeof(int));
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int          **pn_vtx              = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge             = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn       = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn      = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int           *pflat_n_vtx         = (int           *) malloc( n_domain * sizeof(int           ));
  int           *pflat_n_edge        = (int           *) malloc( n_domain * sizeof(int           ));
  int          **pflat_edge_vtx_idx  = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pflat_edge_vtx      = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t  **pflat_vtx_ln_to_gn  = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  PDM_g_num_t  **pflat_edge_ln_to_gn = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  double       **pflat_vtx_coords    = (double       **) malloc( n_domain * sizeof(double       *));

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx        [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn [i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));
    pn_edge       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pedge_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
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
                                          &pflat_edge_vtx    [ln_part_tot+i_part],
                                          &pflat_edge_vtx_idx[ln_part_tot+i_part],
                                          PDM_OWNERSHIP_KEEP);

      assert(pflat_edge_vtx_idx[ln_part_tot+i_part] == NULL);
      pflat_edge_vtx_idx[ln_part_tot+i_part] = malloc((pflat_n_edge[ln_part_tot+i_part] + 1) * sizeof(int));
      for(int i_edge = 0; i_edge < pflat_n_edge[ln_part_tot+i_part]+1; ++i_edge) {
        pflat_edge_vtx_idx[ln_part_tot+i_part][i_edge] = 2 * i_edge;
      }

      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pflat_vtx_coords[ln_part_tot+i_part],
                                       PDM_OWNERSHIP_KEEP);

    }
    ln_part_tot += n_part[i_dom];
  }

  int         **pflat_vtx_edge_idx = NULL;
  int         **pflat_vtx_edge     = NULL;
  PDM_part_connectivity_transpose(ln_part_tot,
                                  pflat_n_edge,
                                  pflat_n_vtx,
                                  pflat_edge_vtx_idx,
                                  pflat_edge_vtx,
                                  &pflat_vtx_edge_idx,
                                  &pflat_vtx_edge);

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
   * Changer les noms ?
   */
  int          *pn_vtx_num             = NULL;
  int         **pvtx_num               = NULL;
  int         **pvtx_opp_location_idx  = NULL;
  int         **pvtx_opp_location      = NULL;
  int         **pvtx_opp_interface_idx = NULL;
  int         **pvtx_opp_interface     = NULL;
  int         **pvtx_opp_sens          = NULL;
  PDM_g_num_t **pvtx_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         PDM_BOUND_TYPE_VTX,
                                         pflat_n_vtx,
                                         pflat_vtx_ln_to_gn,
                                         &pn_vtx_num,
                                         &pvtx_num,
                                         &pvtx_opp_location_idx,
                                         &pvtx_opp_location,
                                         &pvtx_opp_interface_idx,
                                         &pvtx_opp_interface,
                                         &pvtx_opp_sens,
                                         &pvtx_opp_gnum);

  int          *pn_vtx_opp_gnum_and_itrf = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pvtx_opp_gnum_and_itrf   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **vtx_opp_position         = malloc(ln_part_tot * sizeof(int         *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    // int *_pvtx_opp_location = pvtx_opp_location_idx[i_part];
    int *_pvtx_opp_location_idx = pvtx_opp_location_idx[i_part];
    int n_connect_tot = _pvtx_opp_location_idx[pn_vtx_num[i_part]];

    // On garde le lien courant (indirect sort of pvtx_num) + également keep le sens_opp
    pvtx_opp_gnum_and_itrf[i_part] = malloc(2 * n_connect_tot * sizeof(PDM_g_num_t));
    PDM_g_num_t* _pvtx_opp_gnum_and_itrf = pvtx_opp_gnum_and_itrf[i_part];

    for(int idx_entity = 0; idx_entity < pn_vtx_num[i_part]; ++idx_entity) {
      // int i_entity = pvtx_num[i_part][idx_entity];
      for(int idx_opp = _pvtx_opp_location_idx[idx_entity]; idx_opp < _pvtx_opp_location_idx[idx_entity+1]; ++idx_opp) {

        _pvtx_opp_gnum_and_itrf[2*idx_opp  ] = pvtx_opp_gnum     [i_part][idx_opp];
        _pvtx_opp_gnum_and_itrf[2*idx_opp+1] = pvtx_opp_interface[i_part][idx_opp];
      }
    }

    if(1 == 0) {
      PDM_log_trace_array_long(pvtx_opp_gnum_and_itrf[i_part], 2 * n_connect_tot, "pvtx_opp_gnum_and_itrf ::");
    }

    vtx_opp_position[i_part] = malloc(n_connect_tot * sizeof(int));
    pn_vtx_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_connect_tot, 2, pvtx_opp_gnum_and_itrf[i_part], vtx_opp_position[i_part]);


    if(1 == 0) {
      PDM_log_trace_array_int(vtx_opp_position[i_part], pn_vtx_opp_gnum_and_itrf[i_part], "vtx_opp_position = ");
    }


  }

  /*
   * On rajoute l'échange pour le gnum courant
   * Et avec le part_to_part on recupère le gnum opposé
   * Pour chaque entité d'interface
   * Puis on utilise le ptp du extract_part pour tout transferer le gnum et les triplet et le sens
   */
  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pvtx_num[i_part], pn_vtx_num[i_part], "pvtx_num ::");
      PDM_log_trace_graph_nuplet_int(pvtx_opp_location_idx[i_part],
                                     pvtx_opp_location    [i_part],
                                     3,
                                     pn_vtx_num[i_part], "pvtx_opp_location ::");
      PDM_log_trace_array_int (pvtx_opp_interface[i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_interface ::");
      PDM_log_trace_array_int (pvtx_opp_sens     [i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_sens ::");
      PDM_log_trace_array_long(pvtx_opp_gnum     [i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_gnum ::");
    }
  }

  //
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
  /*
   * Pour la recursion, il faut idéalement calculer le graphe avec les partition shifter avec
   * PDM_part_generate_entity_graph_comm avec pentity_hint pour accelerer ou choisir
   * Typiquement pour la deuxième passe on prends que tous les sommets étendu au step d'avant non unifié
   */

  /*
   * Dans tout les cas il ne faut pas oublié de faire le graph de comm proprement
   * avec les connexions entres partitions et entre domaine
   * avec PDM_part_domain_interface_as_graph (à reimplementer)
   */

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


      /*
       * On pourrai le recalculer avec PDM_part_generate_entity_graph_comm avec les données falts au moins pas besoin de gerer des shift ...
       */

      int *ppart_vtx_proc_idx = NULL;
      int *ppart_vtx_part_idx = NULL;
      int *ppart_vtx          = NULL; // (i_vtx, i_proc, i_part, i_vtx_opp)
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_BOUND_TYPE_VTX,
                                        &ppart_vtx_proc_idx,
                                        &ppart_vtx_part_idx,
                                        &ppart_vtx,
                                        PDM_OWNERSHIP_KEEP);

      // PDM_log_trace_array_int(ppart_vtx, 4 * ppart_vtx_part_idx[n_part_g[i_dom]], "ppart_vtx");

      part1_to_part2_idx[li_part] = malloc((pflat_n_vtx[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pflat_n_vtx[li_part]);

      int n_entity_bound = ppart_vtx_part_idx[n_part_g[i_dom]];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_vtx[4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_vtx_num[li_part]; ++idx_entity) {
        int i_entity = pvtx_num[li_part][idx_entity];
        int n_opp = pvtx_opp_location_idx[li_part][idx_entity+1] - pvtx_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }


      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pflat_n_vtx[li_part]];
      part1_to_part2_triplet  [li_part] = malloc(n_connect_tot   * sizeof(int));
      part1_to_part2_interface[li_part] = malloc(n_connect_tot/3 * sizeof(int));

      // printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_vtx[4*idx_entity]-1;
        int i_proc_opp   = ppart_vtx[4*idx_entity+1];
        int i_part_opp   = ppart_vtx[4*idx_entity+2]-1; // A gere -> Le shift
        int i_entity_opp = ppart_vtx[4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp + n_part_shift[i_proc_opp];
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_vtx_num[li_part]; ++idx_entity) {
        int i_entity = pvtx_num[li_part][idx_entity];
        for(int idx_opp = pvtx_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pvtx_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pvtx_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pvtx_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pvtx_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pvtx_opp_interface[li_part][  idx_opp  ];

        }
      }

      if(1 == 0) {
        PDM_log_trace_array_int(part1_to_part2_idx      [li_part], pflat_n_vtx[li_part]+1, "(debug) part1_to_part2_idx       ::");
        PDM_log_trace_array_int(part1_to_part2_triplet  [li_part], n_connect_tot         , "(debug) part1_to_part2_triplet   ::");
        PDM_log_trace_array_int(part1_to_part2_interface[li_part], n_connect_tot/3       , "(debug) part1_to_part2_interface ::");
      }

      /*
       * Creation des buffers d'envoi des connectivités
       */
      int n_send = part1_to_part2_idx[li_part][pflat_n_vtx[li_part]]/3;
      part1_to_part2_edge_vtx_n[li_part] = malloc(n_send * sizeof(int));
      part1_to_part2_edge_n    [li_part] = malloc(n_send * sizeof(int));
      int n_send_edge     = 0;
      int n_send_edge_vtx = 0;

      int *_pflat_vtx_edge_idx        = pflat_vtx_edge_idx       [li_part];
      int *_part1_to_part2_edge_vtx_n = part1_to_part2_edge_vtx_n[li_part];
      int *_part1_to_part2_edge_n     = part1_to_part2_edge_n    [li_part];

      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {

          _part1_to_part2_edge_vtx_n[idx] = 0;
          _part1_to_part2_edge_n    [idx] = 0;

          /* Copy */
          for(int idx_edge = _pflat_vtx_edge_idx[i_entity]; idx_edge < _pflat_vtx_edge_idx[i_entity+1]; ++idx_edge ) {
            int i_edge = PDM_ABS(pflat_vtx_edge[li_part][idx_edge])-1;
            _part1_to_part2_edge_n    [idx] += 1;
            _part1_to_part2_edge_vtx_n[idx] += pflat_edge_vtx_idx[li_part][i_edge+1] - pflat_edge_vtx_idx[li_part][i_edge];
          }
          n_send_edge     += _part1_to_part2_edge_n    [idx];
          n_send_edge_vtx += _part1_to_part2_edge_vtx_n[idx];
        }
      }

      // printf("n_send_edge     = %i \n", n_send_edge);
      // printf("n_send_edge_vtx = %i \n", n_send_edge_vtx);

      part1_to_part2_edge_gnum       [li_part] = malloc(    n_send_edge     * sizeof(PDM_g_num_t));
      part1_to_part2_edge_triplet    [li_part] = malloc(3 * n_send_edge     * sizeof(int        ));
      part1_to_part2_edge_vtx_gnum   [li_part] = malloc(    n_send_edge_vtx * sizeof(PDM_g_num_t));
      part1_to_part2_edge_vtx_triplet[li_part] = malloc(3 * n_send_edge_vtx * sizeof(int        ));
      PDM_g_num_t *_part1_to_part2_edge_gnum        = part1_to_part2_edge_gnum       [li_part];
      PDM_g_num_t *_part1_to_part2_edge_vtx_gnum    = part1_to_part2_edge_vtx_gnum   [li_part];
      int         *_part1_to_part2_edge_triplet     = part1_to_part2_edge_triplet    [li_part];
      int         *_part1_to_part2_edge_vtx_triplet = part1_to_part2_edge_vtx_triplet[li_part];

      /* Fill loop */
      n_send_edge     = 0;
      n_send_edge_vtx = 0;
      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {

          /* Copy */
          for(int idx_edge = _pflat_vtx_edge_idx[i_entity]; idx_edge < _pflat_vtx_edge_idx[i_entity+1]; ++idx_edge ) {
            int i_edge = PDM_ABS(pflat_vtx_edge[li_part][idx_edge])-1;

            _part1_to_part2_edge_gnum    [n_send_edge    ] = pflat_edge_ln_to_gn[li_part][i_edge];
            _part1_to_part2_edge_triplet [3*n_send_edge  ] = i_rank;
            _part1_to_part2_edge_triplet [3*n_send_edge+1] = li_part;
            _part1_to_part2_edge_triplet [3*n_send_edge+2] = i_edge;
            n_send_edge++;

            for(int idx_vtx = pflat_edge_vtx_idx[li_part][i_edge]; idx_vtx < pflat_edge_vtx_idx[li_part][i_edge+1]; ++idx_vtx) {
              int i_vtx = PDM_ABS (pflat_edge_vtx[li_part][idx_vtx])-1;
              int sgn   = PDM_SIGN(pflat_edge_vtx[li_part][idx_vtx]);
              _part1_to_part2_edge_vtx_triplet[3*n_send_edge_vtx  ] = i_rank;
              _part1_to_part2_edge_vtx_triplet[3*n_send_edge_vtx+1] = li_part;
              _part1_to_part2_edge_vtx_triplet[3*n_send_edge_vtx+2] = i_vtx;
              _part1_to_part2_edge_vtx_gnum[n_send_edge_vtx] = sgn * pflat_vtx_ln_to_gn[li_part][i_vtx];
              n_send_edge_vtx++;
            }
          }
        }
      }


      free(part1_to_part2_n);
      li_part += 1;
    }

    free(n_part_shift);
  }


  /*
   * Create part_to_part to exchange all data in opposit part
   */
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pflat_vtx_ln_to_gn,
                                                                      (const int          *) pflat_n_vtx,
                                                                      ln_part_tot,
                                                                      (const int          *) pflat_n_vtx,
                                                                      ln_part_tot,
                                                                      (const int         **) part1_to_part2_idx,
                                                                      (const int         **) part1_to_part2_triplet_idx,
                                                                      (const int         **) part1_to_part2_triplet,
                                                                      comm);


  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  // Il faudrait idéalement refaire un lien explicite face_face / edge_edge puis faire une cascade, pour l'instant on a tout mélanger


  /*
   * Exchange edge gnum
   */
  int exch_request = -1;
  int         **pextract_edge_n    = NULL;
  PDM_g_num_t **pextract_edge_gnum = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(PDM_g_num_t),
        (const int **)   part1_to_part2_edge_n,
        (const void **)  part1_to_part2_edge_gnum,
                         &pextract_edge_n,
            (void ***)   &pextract_edge_gnum,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);

  int **pextract_edge_triplet = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         3 * sizeof(int),
        (const int **)   part1_to_part2_edge_n,
        (const void **)  part1_to_part2_edge_triplet,
                         &pextract_edge_n,
            (void ***)   &pextract_edge_triplet,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);


  /*
   * Exchange connectivity
   */
  int         **pextract_edge_vtx_n    = NULL;
  PDM_g_num_t **pextract_edge_vtx_gnum = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(PDM_g_num_t),
        (const int **)   part1_to_part2_edge_vtx_n,
        (const void **)  part1_to_part2_edge_vtx_gnum,
                         &pextract_edge_vtx_n,
            (void ***)   &pextract_edge_vtx_gnum,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);

  int **pextract_edge_vtx_triplet = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         3 * sizeof(int),
        (const int **)   part1_to_part2_edge_vtx_n,
        (const void **)  part1_to_part2_edge_vtx_triplet,
                         &pextract_edge_vtx_n,
            (void ***)   &pextract_edge_vtx_triplet,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(part1_to_part2_edge_n          [i_part]);
    free(part1_to_part2_edge_gnum       [i_part]);
    free(part1_to_part2_edge_triplet    [i_part]);
    free(part1_to_part2_edge_vtx_n      [i_part]);
    free(part1_to_part2_edge_vtx_gnum   [i_part]);
    free(part1_to_part2_edge_vtx_triplet[i_part]);
  }
  free(part1_to_part2_edge_vtx_n      );
  free(part1_to_part2_edge_vtx_gnum   );
  free(part1_to_part2_edge_vtx_triplet);
  free(part1_to_part2_edge_n          );
  free(part1_to_part2_edge_gnum       );
  free(part1_to_part2_edge_triplet    );

  int **pextract_edge_interface = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(int),
                         NULL,
        (const void **)  part1_to_part2_interface,
                         NULL,
            (void ***)   &pextract_edge_interface,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);

  if(1 == 0) { // Usefull to know how many data is transfer
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");
    }
  }

  int **pextract_edge_vtx_idx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    pextract_edge_vtx_idx[i_part] = malloc((n_ref_lnum2[i_part] + 1) * sizeof(int));

    pextract_edge_vtx_idx[i_part][0] = 0;
    for(int i_edge = 0; i_edge < n_ref_lnum2[i_part]; ++i_edge) {
      pextract_edge_vtx_idx[i_part][i_edge+1] = pextract_edge_vtx_idx[i_part][i_edge] + pextract_edge_vtx_n[i_part][i_edge];
    }

  }

  /* Verbose all recv data */
  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(pextract_edge_interface[i_part], n_ref_lnum2[i_part], "pextract_edge_interface : ");
      PDM_log_trace_array_int(pextract_edge_vtx_n    [i_part], n_ref_lnum2[i_part], "pextract_edge_vtx_n     : ");
      PDM_log_trace_array_int(pextract_edge_n        [i_part], n_ref_lnum2[i_part], "pextract_edge_n         : ");

      int n_recv_edge     = 0;
      int n_recv_edge_vtx = 0;
      for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
        n_recv_edge     += pextract_edge_n    [i_part][i_ref];
        n_recv_edge_vtx += pextract_edge_vtx_n[i_part][i_ref];
      }
      PDM_log_trace_array_long(pextract_edge_vtx_gnum    [i_part], n_recv_edge_vtx    , "pextract_edge_vtx_gnum    : ");
      PDM_log_trace_array_long(pextract_edge_gnum        [i_part], n_recv_edge        , "pextract_edge_gnum        : ");
      PDM_log_trace_array_int (pextract_edge_vtx_triplet [i_part], 3 * n_recv_edge_vtx, "pextract_edge_vtx_triplet : ");
      PDM_log_trace_array_int (pextract_edge_triplet     [i_part], 3 * n_recv_edge    , "pextract_edge_triplet     : ");
    }
  }

  /*
   * Post-treatment
   *   - Create new global numbering with new entity
   *   - Merge connectivity
   *   - Create new global entity for descending connectivity
   */
  PDM_gen_gnum_t* gen_gnum_edge = PDM_gnum_create(3,
                                                  ln_part_tot,
                                                  PDM_TRUE,
                                                  1.e-6,
                                                  comm,
                                                  PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_edge, 2);


  int          *pn_edge_only_by_interface        = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pedge_ln_to_gn_only_by_interface = malloc(ln_part_tot * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    pn_edge_only_by_interface[i_part] = 0;

    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {
        pn_edge_only_by_interface[i_part]++;
      }
    }

    pedge_ln_to_gn_only_by_interface[i_part] = malloc(2 * pn_edge_only_by_interface[i_part] * sizeof(PDM_g_num_t));
    PDM_g_num_t *_pedge_ln_to_gn_only_by_interface = pedge_ln_to_gn_only_by_interface[i_part];
    pn_edge_only_by_interface[i_part] = 0;

    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {
        int idx_write = pn_edge_only_by_interface[i_part]++;
        _pedge_ln_to_gn_only_by_interface[2*idx_write  ] = _pextract_edge_gnum     [i_ref];
        _pedge_ln_to_gn_only_by_interface[2*idx_write+1] = _pextract_edge_interface[i_ref];
      }
    }

    PDM_gnum_set_from_parents(gen_gnum_edge,
                              i_part,
                              pn_edge_only_by_interface[i_part],
                              pedge_ln_to_gn_only_by_interface[i_part]);

  }

  PDM_gnum_compute(gen_gnum_edge);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_g_num_t* extented_edge_ln_to_gn = PDM_gnum_get(gen_gnum_edge, i_part);
    // PDM_log_trace_array_long(extented_edge_ln_to_gn, pn_edge_only_by_interface[i_part], "extented_edge_ln_to_gn ::");

    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    /* Update */
    int idx_read = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {
        _pextract_edge_gnum[i_ref] = extented_edge_ln_to_gn[idx_read++] + shift_by_domain_edge[n_domain];
      }
    }

    // PDM_log_trace_array_long(_pextract_edge_gnum, n_ref_lnum2[i_part], "_pextract_edge_gnum (Update) : ");

  }


  PDM_gnum_free(gen_gnum_edge);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pedge_ln_to_gn_only_by_interface[i_part]);
  }

  free(pn_edge_only_by_interface);
  free(pedge_ln_to_gn_only_by_interface);


  /*
   * ATTENTION : On peut avoir 2 fois le même edge ou face (surtout au deuxieme passage )
   *   Faire un unique puis une indirection (on garde les buffers tel quel c'est plus simple je pense )
   *   Pour chaque sommet connecté par interrface on envoie le gnum pour faire la pair d'arrivé
   *   On doit reconstruire les infos d'interface pour le step d'après (exemple edge --> C'est ultra MANDATORY pour le périodique infini )
   *   Si tout les adjacences sont bonnes on s'arrête sinon on reprendre avec la nouvelle numérotation absolu et les nouvelles interfaces
   *   Attention au multidomaine !!!
   */

  /*
   * Creation numero absolu vtx
   */
  PDM_gen_gnum_t* gen_gnum_vtx = PDM_gnum_create(3,
                                                 ln_part_tot,
                                                 PDM_TRUE,
                                                 1.e-6,
                                                 comm,
                                                 PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_vtx, 2);


  /*
   * Prepare for unify
   */
  int         **vtx_order            = malloc( ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pvtx_ln_to_gn_sorted = malloc( ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    /*
     * Sort current part vtx_ln_to_gn
     */
    vtx_order           [i_part] = malloc( pflat_n_vtx[i_part] * sizeof(int        ));
    pvtx_ln_to_gn_sorted[i_part] = malloc( pflat_n_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < pflat_n_vtx[i_part]; ++i_vtx) {
      pvtx_ln_to_gn_sorted[i_part][i_vtx] = pflat_vtx_ln_to_gn[i_part][i_vtx];
      vtx_order           [i_part][i_vtx] = i_vtx;
    }
    PDM_sort_long(pvtx_ln_to_gn_sorted[i_part], vtx_order[i_part], pflat_n_vtx[i_part]);
  }

  /*
   * We have two kind of extented :
   *   - From partition
   *   - From interfaces (request new gnum génération ...)
   */
  int          *pn_vtx_extented_by_interface         = malloc(ln_part_tot * sizeof(int          ));
  int          *pn_vtx_extented_by_partition         = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pvtx_extented_ln_to_gn_by_interface  = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_extented_ln_to_gn_by_partition  = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pvtx_extented_triplet_by_interface   = malloc(ln_part_tot * sizeof(int         *));
  int         **pvtx_extented_triplet_by_partition   = malloc(ln_part_tot * sizeof(int         *));
  int         **pvtx_extented_interface_by_interface = malloc(ln_part_tot * sizeof(int         *));
  int         **pvtx_extented_interface_by_partition = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t* _pvtx_opp_gnum_and_itrf = pvtx_opp_gnum_and_itrf[i_part];

    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    // PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];
    // int         *_pextract_edge_triplet   = pextract_edge_triplet  [i_part];

    int         *_pextract_edge_vtx_idx     = pextract_edge_vtx_idx    [i_part];
    PDM_g_num_t *_pextract_edge_vtx_gnum    = pextract_edge_vtx_gnum   [i_part];
    int         *_pextract_edge_vtx_triplet = pextract_edge_vtx_triplet[i_part];

    PDM_g_num_t *_pvtx_ln_to_gn_sorted   = pvtx_ln_to_gn_sorted[i_part];

    pn_vtx_extented_by_interface[i_part] = 0;
    pn_vtx_extented_by_partition[i_part] = 0;

    /* Count */
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          PDM_g_num_t gnum_to_find[2] = {vtx_g_num, -_pextract_edge_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pvtx_opp_gnum_and_itrf, 2, pn_vtx_opp_gnum_and_itrf[i_part]);
          if(pos == -1) {
            pn_vtx_extented_by_interface[i_part]++;
          }
          // printf("Search (%i/%i) -> pos = %i \n", vtx_g_num, _pextract_edge_interface[i_ref], pos);
        }
      } else { // Second case : Interior edge

        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          int pos_int = PDM_binary_search_long(vtx_g_num, _pvtx_ln_to_gn_sorted, pflat_n_vtx[i_part]);

          if(pos_int == -1) { // Not in current partition
            pn_vtx_extented_by_partition[i_part]++;
          }
        }
      }
    }

    /* Fill */
    pvtx_extented_ln_to_gn_by_interface [i_part] = malloc(2 * pn_vtx_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    pvtx_extented_ln_to_gn_by_partition [i_part] = malloc(    pn_vtx_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    pvtx_extented_triplet_by_interface  [i_part] = malloc(3 * pn_vtx_extented_by_interface[i_part] * sizeof(int        ));
    pvtx_extented_triplet_by_partition  [i_part] = malloc(3 * pn_vtx_extented_by_partition[i_part] * sizeof(int        ));
    pvtx_extented_interface_by_interface[i_part] = malloc(    pn_vtx_extented_by_interface[i_part] * sizeof(int        ));
    pvtx_extented_interface_by_partition[i_part] = malloc(    pn_vtx_extented_by_partition[i_part] * sizeof(int        ));

    pn_vtx_extented_by_interface[i_part] = 0;
    pn_vtx_extented_by_partition[i_part] = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          PDM_g_num_t gnum_to_find[2] = {vtx_g_num, -_pextract_edge_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pvtx_opp_gnum_and_itrf, 2, pn_vtx_opp_gnum_and_itrf[i_part]);
          if(pos == -1) {
            int idx_write = pn_vtx_extented_by_interface[i_part]++;
            pvtx_extented_ln_to_gn_by_interface[i_part][2*idx_write  ] = vtx_g_num;
            pvtx_extented_ln_to_gn_by_interface[i_part][2*idx_write+1] = -_pextract_edge_interface[i_ref];

            pvtx_extented_triplet_by_interface  [i_part][3*idx_write  ] = _pextract_edge_vtx_triplet[3*idx_edge  ];
            pvtx_extented_triplet_by_interface  [i_part][3*idx_write+1] = _pextract_edge_vtx_triplet[3*idx_edge+1];
            pvtx_extented_triplet_by_interface  [i_part][3*idx_write+2] = _pextract_edge_vtx_triplet[3*idx_edge+2];
            pvtx_extented_interface_by_interface[i_part][  idx_write  ] = -_pextract_edge_interface[i_ref];


          }
        }
      } else { // Second case : Interior edge

        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          int pos_int = PDM_binary_search_long(vtx_g_num, _pvtx_ln_to_gn_sorted, pflat_n_vtx[i_part]);

          if(pos_int == -1) { // Not in current partition
            int idx_write = pn_vtx_extented_by_partition[i_part]++;
            pvtx_extented_ln_to_gn_by_partition [i_part][idx_write] = vtx_g_num;
            pvtx_extented_triplet_by_partition  [i_part][3*idx_write  ] = _pextract_edge_vtx_triplet[3*idx_edge  ];
            pvtx_extented_triplet_by_partition  [i_part][3*idx_write+1] = _pextract_edge_vtx_triplet[3*idx_edge+1];
            pvtx_extented_triplet_by_partition  [i_part][3*idx_write+2] = _pextract_edge_vtx_triplet[3*idx_edge+2];
            pvtx_extented_interface_by_partition[i_part][  idx_write  ] = 0; // Because interior
          }
        }
      }
    }

    PDM_gnum_set_from_parents(gen_gnum_vtx,
                              i_part,
                              pn_vtx_extented_by_interface[i_part],
                              pvtx_extented_ln_to_gn_by_interface[i_part]);

  }
  PDM_gnum_compute(gen_gnum_vtx);

  /*
   * Unify edge_vtx
   */
  int **pextract_edge_vtx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pvtx_opp_gnum_and_itrf  = pvtx_opp_gnum_and_itrf [i_part];
    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    // PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    int         *_pextract_edge_vtx_idx  = pextract_edge_vtx_idx [i_part];
    PDM_g_num_t *_pextract_edge_vtx_gnum = pextract_edge_vtx_gnum[i_part];

    int           n_vtx_opp_position     = pn_vtx_opp_gnum_and_itrf[i_part];
    // int         *_vtx_opp_position       = vtx_opp_position        [i_part];

    PDM_g_num_t *_pvtx_ln_to_gn_sorted   = pvtx_ln_to_gn_sorted[i_part];
    int         *_vtx_order              = vtx_order           [i_part];

    pextract_edge_vtx[i_part] = malloc( _pextract_edge_vtx_idx[n_ref_lnum2[i_part]] * sizeof(int));
    int *_pextract_edge_vtx  = pextract_edge_vtx [i_part];

    PDM_g_num_t* extented_from_itrf_vtx_ln_to_gn = PDM_gnum_get(gen_gnum_vtx, i_part);

    // Erase and realloc :
    pvtx_extented_ln_to_gn_by_interface[i_part] = realloc(pvtx_extented_ln_to_gn_by_interface[i_part], pn_vtx_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < pn_vtx_extented_by_interface[i_part]; ++i_vtx) {
      pvtx_extented_ln_to_gn_by_interface[i_part][i_vtx] = extented_from_itrf_vtx_ln_to_gn[i_vtx] + shift_by_domain_vtx[n_domain];
    }

    // PDM_log_trace_array_long(extented_from_itrf_vtx_ln_to_gn, pn_vtx_extented_by_interface[i_part], "extented_from_itrf_vtx_ln_to_gn ::");

    /* Sort unique the new gnum to unify */
    int         *extented_from_itrf_vtx_order           = malloc( pn_vtx_extented_by_interface[i_part] * sizeof(int        ));
    PDM_g_num_t *extented_from_itrf_vtx_ln_to_gn_sorted = malloc( pn_vtx_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < pn_vtx_extented_by_interface[i_part]; ++i_vtx) {
      extented_from_itrf_vtx_ln_to_gn_sorted[i_vtx] = extented_from_itrf_vtx_ln_to_gn[i_vtx];
      extented_from_itrf_vtx_order          [i_vtx] = i_vtx;
    }
    PDM_sort_long(extented_from_itrf_vtx_ln_to_gn_sorted, extented_from_itrf_vtx_order, pn_vtx_extented_by_interface[i_part]);

    /* Sort unique the new gnum to unify */
    int         *extented_from_part_vtx_order           = malloc( pn_vtx_extented_by_partition[i_part] * sizeof(int        ));
    PDM_g_num_t *extented_from_part_vtx_ln_to_gn_sorted = malloc( pn_vtx_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < pn_vtx_extented_by_partition[i_part]; ++i_vtx) {
      extented_from_part_vtx_ln_to_gn_sorted[i_vtx] = pvtx_extented_ln_to_gn_by_partition[i_part][i_vtx];
      extented_from_part_vtx_order          [i_vtx] = i_vtx;
    }
    PDM_sort_long(extented_from_part_vtx_ln_to_gn_sorted, extented_from_part_vtx_order, pn_vtx_extented_by_partition[i_part]);

    /* */
    int pn_vtx_extented_by_interface2 = 0; // To read in extented_from_itrf_vtx_ln_to_gn
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {

      /*
       * First case :
       *   - Edge is move by interface
       */
      if(_pextract_edge_interface[i_ref] != 0) {
        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          PDM_g_num_t gnum_to_find[2] = {vtx_g_num, -_pextract_edge_interface[i_ref]};

          int pos = PDM_order_binary_search_long(gnum_to_find, _pvtx_opp_gnum_and_itrf, 2, n_vtx_opp_position);
          if(pos == -1) {
            /*
             * Subcase :
             *   - Vtx is not in table of interface : it's a new vtx
             */
            PDM_g_num_t vtx_extented_g_num = extented_from_itrf_vtx_ln_to_gn[pn_vtx_extented_by_interface2++];
            int pos_new = PDM_binary_search_long(vtx_extented_g_num, extented_from_itrf_vtx_ln_to_gn_sorted, pn_vtx_extented_by_interface[i_part]);

            assert(pos_new != -1);
            int shift = pflat_n_vtx[i_part] + pn_vtx_extented_by_partition[i_part];
            _pextract_edge_vtx[idx_edge] = ( shift + extented_from_itrf_vtx_order[pos_new] + 1); // ATTENTION SIGN

          } else {
            /*
             * Subcase :
             *   - Vtx is in table of interface
             */
            int pos2  = vtx_opp_position[i_part][pos ];
            int i_vtx = pvtx_num        [i_part][pos2];

            _pextract_edge_vtx[idx_edge] = ( i_vtx + 1);

          }
        }
      } else { // Second case : Interior edge

        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          int pos_int = PDM_binary_search_long(vtx_g_num, _pvtx_ln_to_gn_sorted, pflat_n_vtx[i_part]);

          if(pos_int != -1) {                   // vtx is already in current partition
            assert(pos_int != -1);
            int orig_pos = _vtx_order[pos_int];
            _pextract_edge_vtx[idx_edge] = ( orig_pos + 1);
          } else {                              // Not in current partition
            int pos_ext = PDM_binary_search_long(vtx_g_num, extented_from_part_vtx_ln_to_gn_sorted, pn_vtx_extented_by_partition[i_part]);
            assert(pos_ext != -1);
            _pextract_edge_vtx[idx_edge] = ( pflat_n_vtx[i_part] + extented_from_part_vtx_order[pos_ext] + 1);

          }
        }
      }
    }

    assert(pn_vtx_extented_by_interface2 == pn_vtx_extented_by_interface[i_part]);

    PDM_log_trace_connectivity_int(_pextract_edge_vtx_idx, _pextract_edge_vtx, n_ref_lnum2[i_part], "pextract_edge_vtx ::");

    free(extented_from_itrf_vtx_order          );
    free(extented_from_itrf_vtx_ln_to_gn_sorted);
    free(extented_from_part_vtx_order          );
    free(extented_from_part_vtx_ln_to_gn_sorted);

  }


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(vtx_order           [i_part]);
    free(pvtx_ln_to_gn_sorted[i_part]);
  }
  free(vtx_order);
  free(pvtx_ln_to_gn_sorted);

  PDM_gnum_free(gen_gnum_vtx);


  /*
   * Reconstruction extented partition
   *      - edge_ln_to_gn   -> pextract_edge_gnum
   *      - edge_interface (il faut toute les infos qu' besoin le part_domain_interface )
   *           -> pextract_edge_interface
   *           -> pextract_edge_opp_gnum    (ne pas updater du coup)
   *           -> pextract_edge_opp_triplet (ne pas updater du coup)
   *      - vtx_ln_to_gn / vtx_triplet_opp
   *           -> pvtx_extented_ln_to_gn puis pvtx_extented_ln_to_gn_by_interface
   *      - pentity_hint (pour le rang 2 par exemple )
   *      - vtx_coords
   *
   *    Refaire un part_domain_interface avec tout -> PDM_part_domain_interface_set
   *      Au premier passage on fait explicitement celui des edges
   *      Pour les autres entités on doit concatener et refaire la structure
   *
   */
  int          *pn_edge_extented         = malloc(ln_part_tot * sizeof(int          ));
  int          *pn_vtx_extented          = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pvtx_extented_ln_to_gn   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pvtx_extented_triplet    = malloc(ln_part_tot * sizeof(int         *));
  int         **pvtx_extented_interface  = malloc(ln_part_tot * sizeof(int         *));
  int         **pvtx_extented_to_vtx_idx = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {


    pn_edge_extented[i_part] = n_ref_lnum2[i_part];

    pn_vtx_extented[i_part] = pn_vtx_extented_by_partition[i_part] + pn_vtx_extented_by_interface[i_part];
    pvtx_extented_ln_to_gn  [i_part] = malloc(    pn_vtx_extented[i_part]      * sizeof(PDM_g_num_t));
    pvtx_extented_triplet   [i_part] = malloc(3 * pn_vtx_extented[i_part]      * sizeof(int        ));
    pvtx_extented_interface [i_part] = malloc(    pn_vtx_extented[i_part]      * sizeof(int        ));
    pvtx_extented_to_vtx_idx[i_part] = malloc( (  pn_vtx_extented[i_part] + 1) * sizeof(int        ));

    for(int i_vtx = 0; i_vtx < pn_vtx_extented_by_partition[i_part]; ++i_vtx) {
      pvtx_extented_ln_to_gn [i_part][  i_vtx  ] = pvtx_extented_ln_to_gn_by_partition [i_part][  i_vtx  ];
      pvtx_extented_triplet  [i_part][3*i_vtx  ] = pvtx_extented_triplet_by_partition  [i_part][3*i_vtx  ];
      pvtx_extented_triplet  [i_part][3*i_vtx+1] = pvtx_extented_triplet_by_partition  [i_part][3*i_vtx+1];
      pvtx_extented_triplet  [i_part][3*i_vtx+2] = pvtx_extented_triplet_by_partition  [i_part][3*i_vtx+2];
      pvtx_extented_interface[i_part][  i_vtx  ] = pvtx_extented_interface_by_partition[i_part][  i_vtx  ];
    }

    for(int i_vtx = 0; i_vtx < pn_vtx_extented_by_interface[i_part]; ++i_vtx) {
      int idx_write = pn_vtx_extented_by_partition[i_part]+i_vtx;
      pvtx_extented_ln_to_gn [i_part][  idx_write  ] = pvtx_extented_ln_to_gn_by_interface [i_part][  i_vtx  ];
      pvtx_extented_triplet  [i_part][3*idx_write  ] = pvtx_extented_triplet_by_interface  [i_part][3*i_vtx  ];
      pvtx_extented_triplet  [i_part][3*idx_write+1] = pvtx_extented_triplet_by_interface  [i_part][3*i_vtx+1];
      pvtx_extented_triplet  [i_part][3*idx_write+2] = pvtx_extented_triplet_by_interface  [i_part][3*i_vtx+2];
      pvtx_extented_interface[i_part][  idx_write  ] = pvtx_extented_interface_by_interface[i_part][  i_vtx  ];
    }

    for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]+1; ++i_vtx) {
      pvtx_extented_to_vtx_idx[i_part][i_vtx] = 3 * i_vtx;
    }

    if(0 == 1) {
      PDM_log_trace_array_int (pvtx_extented_to_vtx_idx[i_part],     pn_vtx_extented[i_part]+1, "pvtx_extented_to_vtx_idx ::" );
      PDM_log_trace_array_int (pvtx_extented_triplet   [i_part], 3 * pn_vtx_extented[i_part]  , "pvtx_extented_triplet    ::" );
      PDM_log_trace_array_int (pvtx_extented_interface [i_part],     pn_vtx_extented[i_part]  , "pvtx_extented_interface  ::" );
      PDM_log_trace_array_long(pvtx_extented_ln_to_gn  [i_part],     pn_vtx_extented[i_part]  , "pvtx_extented_ln_to_gn   ::" );
    }

  }

  /*
   * Create part_to_part to hook coordinates
   */
  PDM_part_to_part_t* ptp_extract_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pvtx_extented_ln_to_gn,
                                                                                  (const int          *) pn_vtx_extented,
                                                                                                         ln_part_tot,
                                                                                  (const int          *) pflat_n_vtx,
                                                                                                         ln_part_tot,
                                                                                  (const int         **) pvtx_extented_to_vtx_idx,
                                                                                  (const int         **) NULL,
                                                                                  (const int         **) pvtx_extented_triplet,
                                                                                                         comm);

  int  *n_extract_lnum2 = NULL;
  int **extract_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_extract_vtx, &n_extract_lnum2, &extract_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp_extract_vtx, &gnum1_come_from_idx, &gnum1_come_from);

  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(extract_lnum2[i_part], n_extract_lnum2[i_part], "extract_lnum2 :");
      PDM_log_trace_connectivity_long(gnum1_come_from_idx[i_part],
                                      gnum1_come_from    [i_part],
                                      n_extract_lnum2  [i_part], "gnum1_come_from ::");
    }
  }

  exch_request = -1;
  double **pextract_vtx_coord = NULL;
  PDM_part_to_part_reverse_iexch(ptp_extract_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 3,
                                 sizeof(double),
                                 NULL,
                (const void **)  pflat_vtx_coords,
                                 NULL,
                    (void ***)   &pextract_vtx_coord,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp_extract_vtx, exch_request);

  /*
   * Apply transformation if any
   */
  int n_interface = 0;
  if(pdi != NULL) {
    n_interface = PDM_part_domain_interface_n_interface_get(pdi);
  }
  double  **translation_vector = malloc(n_interface * sizeof(double *  ));
  double ***rotation_matrix    = malloc(n_interface * sizeof(double ** ));
  double  **rotation_direction = malloc(n_interface * sizeof(double *  ));
  double  **rotation_center    = malloc(n_interface * sizeof(double *  ));
  double   *rotation_angle     = malloc(n_interface * sizeof(double    ));
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
      rotation_matrix[i_interf] = malloc(3 * sizeof(double *));
      for(int k = 0; k < 3; ++k) {
        rotation_matrix[i_interf][k] = malloc(3 * sizeof(double));
      }
    }
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]; ++i_vtx) {
      int i_interface   = PDM_ABS (pvtx_extented_interface[i_part][i_vtx]);
      int sgn_interface = PDM_SIGN(pvtx_extented_interface[i_part][i_vtx]);
      if(i_interface != 0 && translation_vector[PDM_ABS(i_interface)-1] != NULL) {
        for(int k = 0; k < 3; ++k) {
          pextract_vtx_coord[i_part][3*i_vtx+k] += sgn_interface * translation_vector[PDM_ABS(i_interface)-1][k];
        }
      }
    }
  }



  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    if(translation_vector[i_interf] != NULL) {
      free(translation_vector[i_interf]);
    }
    if(rotation_center    [i_interf] != NULL) {
      for(int k = 0; k < 3; ++k) {
        free(rotation_matrix[i_interf][k]);
      }
      free(rotation_matrix[i_interf]);
    }
  }
  free(translation_vector);
  free(rotation_matrix);
  free(rotation_direction);
  free(rotation_center);
  free(rotation_angle);


  PDM_part_to_part_free(ptp_extract_vtx);

  /*
   * Export vtk
   */
  if(1 == 1) {

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      int pn_concat_vtx  = pflat_n_vtx [i_part] + pn_vtx_extented [i_part];
      int pn_concat_edge = pflat_n_edge[i_part] + pn_edge_extented[i_part];

      int *edge_kind = malloc(pn_concat_edge * sizeof(int));

      int         *concat_edge_vtx      = malloc(2 * pn_concat_edge * sizeof(int         ));
      PDM_g_num_t *concat_edge_ln_to_gn = malloc(    pn_concat_edge * sizeof(PDM_g_num_t ));
      double      *concat_vtx_coord     = malloc(3 * pn_concat_vtx  * sizeof(double      ));
      PDM_g_num_t *concat_vtx_ln_to_gn  = malloc(    pn_concat_vtx  * sizeof(PDM_g_num_t ));

      for(int i_edge = 0; i_edge < pflat_n_edge[i_part]; ++i_edge) {
        concat_edge_vtx[2*i_edge  ] = pflat_edge_vtx[i_part][2*i_edge  ];
        concat_edge_vtx[2*i_edge+1] = pflat_edge_vtx[i_part][2*i_edge+1];
        concat_edge_ln_to_gn[i_edge] = pflat_edge_ln_to_gn[i_part][i_edge];
        edge_kind[i_edge] = 0;
      }

      for(int i_vtx = 0; i_vtx < pflat_n_vtx[i_part]; ++i_vtx) {
        concat_vtx_coord[3*i_vtx  ] = pflat_vtx_coords[i_part][3*i_vtx  ];
        concat_vtx_coord[3*i_vtx+1] = pflat_vtx_coords[i_part][3*i_vtx+1];
        concat_vtx_coord[3*i_vtx+2] = pflat_vtx_coords[i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[i_vtx] = pflat_vtx_ln_to_gn[i_part][i_vtx];
      }


      for(int i_edge = 0; i_edge < pn_edge_extented[i_part]; ++i_edge) {
        int idx_write = pflat_n_edge[i_part]+i_edge;
        concat_edge_vtx     [2*idx_write  ] = pextract_edge_vtx [i_part][2*i_edge  ];
        concat_edge_vtx     [2*idx_write+1] = pextract_edge_vtx [i_part][2*i_edge+1];
        concat_edge_ln_to_gn[  idx_write  ] = pextract_edge_gnum[i_part][i_edge];
        edge_kind           [idx_write] = 1;
      }

      for(int i_vtx = 0; i_vtx < pn_vtx_extented [i_part]; ++i_vtx) {
        int idx_write = pflat_n_vtx[i_part]+i_vtx;
        concat_vtx_coord   [3*idx_write  ] = pextract_vtx_coord    [i_part][3*i_vtx  ];
        concat_vtx_coord   [3*idx_write+1] = pextract_vtx_coord    [i_part][3*i_vtx+1];
        concat_vtx_coord   [3*idx_write+2] = pextract_vtx_coord    [i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[  idx_write  ] = pvtx_extented_ln_to_gn[i_part][  i_vtx  ];
      }

      char filename[999];
      sprintf(filename, "out_part_concate_vtx_i_part=%i_%i.vtk", i_part, i_rank);
      const char* field_name[] = {"edge_kind", 0 };
      int *field[1] = {edge_kind};
      PDM_vtk_write_std_elements (filename,
                                  pn_concat_vtx,
                                  concat_vtx_coord,
                                  concat_vtx_ln_to_gn,
                                  PDM_MESH_NODAL_BAR2,
                                  pn_concat_edge,
                                  concat_edge_vtx,
                                  concat_edge_ln_to_gn,
                                  1,
                                  field_name,
                                  (const int **) field);

      free(edge_kind);

      free(concat_edge_vtx);
      free(concat_edge_ln_to_gn);
      free(concat_vtx_coord);
      free(concat_vtx_ln_to_gn);

    }
  }


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    // PDM_log_trace_array_double(pextract_vtx_coord[i_part], 3 * pn_vtx_extented[i_part], "pextract_vtx_coord ::" );
    free(pextract_vtx_coord[i_part]);
  }
  free(pextract_vtx_coord);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_edge_vtx[i_part]);
  }
  free(pextract_edge_vtx);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_extented_ln_to_gn  [i_part]);
    free(pvtx_extented_triplet   [i_part]);
    free(pvtx_extented_interface [i_part]);
    free(pvtx_extented_to_vtx_idx[i_part]);
  }
  free(pn_edge_extented);
  free(pn_vtx_extented);
  free(pvtx_extented_ln_to_gn);
  free(pvtx_extented_triplet);
  free(pvtx_extented_interface);
  free(pvtx_extented_to_vtx_idx);


  /*
   * Pour faire la cascade IL FAUT garder les opp_location !!!!
   */

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_opp_gnum_and_itrf[i_part]);
    free(vtx_opp_position      [i_part]);
  }
  free(pvtx_opp_gnum_and_itrf);
  free(pn_vtx_opp_gnum_and_itrf);
  free(vtx_opp_position);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_extented_ln_to_gn_by_partition [i_part]);
    free(pvtx_extented_ln_to_gn_by_interface [i_part]);
    free(pvtx_extented_triplet_by_partition  [i_part]);
    free(pvtx_extented_triplet_by_interface  [i_part]);
    free(pvtx_extented_interface_by_partition[i_part]);
    free(pvtx_extented_interface_by_interface[i_part]);
  }
  free(pn_vtx_extented_by_interface);
  free(pn_vtx_extented_by_partition);
  free(pvtx_extented_ln_to_gn_by_interface);
  free(pvtx_extented_ln_to_gn_by_partition);
  free(pvtx_extented_triplet_by_interface);
  free(pvtx_extented_triplet_by_partition);
  free(pvtx_extented_interface_by_interface);
  free(pvtx_extented_interface_by_partition);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_edge_interface  [i_part]);
    free(pextract_edge_vtx_n      [i_part]);
    free(pextract_edge_vtx_idx    [i_part]);
    free(pextract_edge_vtx_gnum   [i_part]);
    free(pextract_edge_vtx_triplet[i_part]);
    free(pextract_edge_n          [i_part]);
    free(pextract_edge_gnum       [i_part]);
    free(pextract_edge_triplet    [i_part]);
  }
  free(pextract_edge_interface  );
  free(pextract_edge_vtx_n      );
  free(pextract_edge_vtx_idx    );
  free(pextract_edge_vtx_gnum   );
  free(pextract_edge_vtx_triplet);
  free(pextract_edge_n          );
  free(pextract_edge_gnum       );
  free(pextract_edge_triplet    );



  PDM_part_to_part_free(ptp);

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

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(part1_to_part2_idx        [i_part]);
    free(part1_to_part2_triplet    [i_part]);
    free(part1_to_part2_interface  [i_part]);

    free(pflat_vtx_edge_idx[i_part] );
    free(pflat_vtx_edge    [i_part] );
    free(pflat_edge_vtx_idx[i_part] );
  }
  free(part1_to_part2_idx        );
  free(part1_to_part2_triplet    );
  free(part1_to_part2_interface  );

  free(pflat_vtx_edge_idx );
  free(pflat_vtx_edge     );

  free(shift_by_domain_vtx);
  free(shift_by_domain_edge);
  free(n_part_g);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    free(pn_edge       [i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pn_vtx        [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
  }

  free(pn_edge             );
  free(pedge_ln_to_gn      );
  free(pn_vtx             );
  free(pvtx_ln_to_gn      );
  free(pflat_n_vtx        );
  free(pflat_n_edge       );
  free(pflat_edge_vtx_idx );
  free(pflat_edge_vtx     );
  free(pflat_vtx_ln_to_gn );
  free(pflat_edge_ln_to_gn);
  free(pflat_vtx_coords   );
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_num             [i_part]);
    free(pvtx_opp_location_idx[i_part]);
    free(pvtx_opp_location    [i_part]);
    free(pvtx_opp_interface   [i_part]);
    free(pvtx_opp_sens        [i_part]);
    free(pvtx_opp_gnum        [i_part]);
  }
  free(pn_vtx_num            );
  free(pvtx_num              );
  free(pvtx_opp_location_idx );
  free(pvtx_opp_location     );
  free(pvtx_opp_interface_idx);
  free(pvtx_opp_interface    );
  free(pvtx_opp_sens         );
  free(pvtx_opp_gnum         );



}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

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

  PDM_dmesh_t **dm = malloc(n_dom_i * sizeof(PDM_dmesh_t *));

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
  int* n_part = malloc(n_dom_i * sizeof(int));
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
    PDM_multipart_register_block(mpart, i_dom, dm[i_dom]);
  }

  PDM_multipart_run_ppart(mpart);

  int n_domain = n_dom_i;

  int          **pn_vtx         = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
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
                                        PDM_MESH_ENTITY_VERTEX,
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
                                            &pedge_vtx,
                                            &pedge_vtx_idx,
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
  free(n_part);

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {

    free (dvtx_coord  [i_dom]);
    free (distrib_vtx [i_dom]);
    free (distrib_edge[i_dom]);
    free (dedge_vtx   [i_dom]);
    PDM_dmesh_free(dm[i_dom]);
  }
  free (dvtx_coord);
  free (distrib_vtx);
  free (distrib_edge);
  free (dedge_vtx);
  free (dm);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    free(pn_vtx       [i_dom]);
    free(pvtx_ln_to_gn[i_dom]);
  }
  free(pn_vtx);
  free(pvtx_ln_to_gn);

  PDM_part_domain_interface_free(pdi);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
