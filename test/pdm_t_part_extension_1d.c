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

  int *n_part_g;
  PDM_malloc(n_part_g,n_domain ,int);
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      ln_part_tot += 1;
    }
  }



  int **pn_vtx;
  PDM_malloc(pn_vtx, n_domain    ,int          *);
  int **pn_edge;
  PDM_malloc(pn_edge, n_domain    ,int          *);
  PDM_g_num_t ***pvtx_ln_to_gn;
  PDM_malloc(pvtx_ln_to_gn, n_domain    ,PDM_g_num_t **);
  PDM_g_num_t ***pedge_ln_to_gn;
  PDM_malloc(pedge_ln_to_gn, n_domain    ,PDM_g_num_t **);
  int *pflat_n_vtx;
  PDM_malloc(pflat_n_vtx, n_domain    ,int           );
  int *pflat_n_edge;
  PDM_malloc(pflat_n_edge, n_domain    ,int           );
  int ***pedge_vtx_idx;
  PDM_malloc(pedge_vtx_idx, n_domain    ,int         **);
  int ***pedge_vtx;
  PDM_malloc(pedge_vtx, n_domain    ,int         **);
  int **pflat_edge_vtx;
  PDM_malloc(pflat_edge_vtx, ln_part_tot ,int          *);
  PDM_g_num_t **pflat_vtx_ln_to_gn;
  PDM_malloc(pflat_vtx_ln_to_gn, ln_part_tot ,PDM_g_num_t  *);
  PDM_g_num_t **pflat_edge_ln_to_gn;
  PDM_malloc(pflat_edge_ln_to_gn, ln_part_tot ,PDM_g_num_t  *);
  double **pflat_vtx_coords;
  PDM_malloc(pflat_vtx_coords, ln_part_tot ,double       *);

  ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx        PDM_malloc([i_dom], n_part[i_dom] ,int          );
    pvtx_ln_to_gn PDM_malloc([i_dom], n_part[i_dom] ,PDM_g_num_t *);
    pn_edge       PDM_malloc([i_dom], n_part[i_dom] ,int          );
    PDM_malloc(pedge_ln_to_gn[i_dom], n_part[i_dom] ,PDM_g_num_t *);

    pedge_vtx_idx PDM_malloc([i_dom], n_part[i_dom] ,int         *);
    pedge_vtx     PDM_malloc([i_dom], n_part[i_dom] ,int         *);

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
      PDM_malloc(pedge_vtx_idx[i_dom][i_part],(pn_edge[i_dom][i_part] + 1) ,int);
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
  double **translation_vector;
  PDM_malloc(translation_vector,n_interface ,double *  );
  double ***rotation_matrix;
  PDM_malloc(rotation_matrix,n_interface ,double ** );
  double **rotation_direction;
  PDM_malloc(rotation_direction,n_interface ,double *  );
  double **rotation_center;
  PDM_malloc(rotation_center,n_interface ,double *  );
  double *rotation_angle;
  PDM_malloc(rotation_angle,n_interface ,double    );
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
      PDM_malloc(rotation_matrix[i_interf],3 ,double *);
      for(int k = 0; k < 3; ++k) {
        PDM_malloc(rotation_matrix[i_interf][k],3 ,double);
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


  /*
   * Export vtk
   */
  if(1 == 1) {

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      int pn_concat_vtx  = pflat_n_vtx [i_part] + pn_vtx_extented [i_part];
      int pn_concat_edge = pflat_n_edge[i_part] + pn_edge_extented[i_part];

      int *edge_kind;
      PDM_malloc(edge_kind,pn_concat_edge ,int);

      int *concat_edge_vtx;
      PDM_malloc(concat_edge_vtx,2 * pn_concat_edge ,int         );
      PDM_g_num_t *concat_edge_ln_to_gn;
      PDM_malloc(concat_edge_ln_to_gn,    pn_concat_edge ,PDM_g_num_t );
      double *concat_vtx_coord;
      PDM_malloc(concat_vtx_coord,3 * pn_concat_vtx  ,double      );
      PDM_g_num_t *concat_vtx_ln_to_gn;
      PDM_malloc(concat_vtx_ln_to_gn,    pn_concat_vtx  ,PDM_g_num_t );

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

     PDM_free(edge_kind);

     PDM_free(concat_edge_vtx);
     PDM_free(concat_edge_ln_to_gn);
     PDM_free(concat_vtx_coord);
     PDM_free(concat_vtx_ln_to_gn);

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
   PDM_free(pextract_vtx_coords[i_part]);
  }
 PDM_free(pextract_vtx_coords);


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

  PDM_dmesh_t **dm;
  PDM_malloc(dm,n_dom_i ,PDM_dmesh_t *);

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
  int *n_part;
  PDM_malloc(n_part,n_dom_i ,int);
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

  int **pn_vtx;
  PDM_malloc(pn_vtx, n_domain ,int          *);
  PDM_g_num_t ***pvtx_ln_to_gn;
  PDM_malloc(pvtx_ln_to_gn, n_domain ,PDM_g_num_t **);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx       PDM_malloc([i_dom], n_part[i_dom] ,int          );
    PDM_malloc(pvtx_ln_to_gn[i_dom], n_part[i_dom] ,PDM_g_num_t *);

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

   PDM_free(dvtx_coord  [i_dom]);
   PDM_free(distrib_vtx [i_dom]);
   PDM_free(distrib_edge[i_dom]);
   PDM_free(dedge_vtx   [i_dom]);
    PDM_dmesh_free(dm[i_dom]);
  }
 PDM_free(dvtx_coord);
 PDM_free(distrib_vtx);
 PDM_free(distrib_edge);
 PDM_free(dedge_vtx);
 PDM_free(dm);

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
