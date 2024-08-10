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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_writer.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_extension.h"
#include "pdm_domain_utils.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_vtk.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_logging.h"
#include "pdm_order.h"
#include "pdm_binary_search.h"
#include "pdm_part_connectivity_transform.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_domain_interface.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_part_extension_algorithm.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/
static void
_hexa_ngon_to_nodal
(
 int   n_cell,
 int  *cell_face_idx,
 int  *cell_face,
 int  *face_vtx_idx,
 int  *face_vtx,
 int **cell_vtx
 )
{
  PDM_malloc(*cell_vtx, n_cell * 8, int);

  for (int i_cell = 0; i_cell < n_cell; i_cell++) {

    int *cv = (*cell_vtx) + 8*i_cell;
    int *cf = cell_face + cell_face_idx[i_cell];
    assert(cell_face_idx[i_cell+1] - cell_face_idx[i_cell] == 6);

    int i_face = PDM_ABS(cf[0]) - 1;
    assert(face_vtx_idx[i_face+1] - face_vtx_idx[i_face] == 4);
    if (cf[0] < 0) {
      for (int i = 0; i < 4; i++) {
        cv[  i] = face_vtx[face_vtx_idx[i_face] + i];
      }
    }
    else {
      for (int i = 0; i < 4; i++) {
        cv[3-i] = face_vtx[face_vtx_idx[i_face] + i];
      }
    }

    int found[2] = {0, 0};
    for (int i = 1; i < 6; i++) {
      i_face = PDM_ABS(cf[i]) - 1;
      assert(face_vtx_idx[i_face+1] - face_vtx_idx[i_face] == 4);
      int *fv = face_vtx + face_vtx_idx[i_face];

      int skip = 0;
      for (int j = 0; j < 4; j++) {
        int i_vtx0, i_vtx1;
        if (cf[i] < 0) {
          i_vtx0 = fv[(j+1)%4];
          i_vtx1 = fv[j      ];
        }
        else {
          i_vtx0 = fv[j      ];
          i_vtx1 = fv[(j+1)%4];
        }

        if (!found[0]) {
          if (i_vtx0 == cv[0] && i_vtx1 == cv[1]) {
            if (cf[i] < 0) {
              cv[4] = fv[(j+2)%4];
              cv[5] = fv[(j+3)%4];
            }
            else {
              cv[4] = fv[(j+3)%4];
              cv[5] = fv[(j+2)%4];
            }
            found[0] = 1;
            skip = 1;
            break;
          }
        }
      }
      if (skip) continue;

      if (!found[1]) {
        for (int j = 0; j < 4; j++) {
          int i_vtx0, i_vtx1;
          if (cf[i] < 0) {
            i_vtx0 = fv[j      ];
            i_vtx1 = fv[(j+1)%4];
          }
          else {
            i_vtx0 = fv[(j+1)%4];
            i_vtx1 = fv[j      ];
          }

          if (i_vtx0 == cv[3] && i_vtx1 == cv[2]) {
            if (cf[i] < 0) {
              cv[6] = fv[(j+2)%4];
              cv[7] = fv[(j+3)%4];
            }
            else {
              cv[6] = fv[(j+3)%4];
              cv[7] = fv[(j+2)%4];
            }
            found[1] = 1;
            break;
          }
        }
      }

    }

    assert(found[0] && found[1]);

  }
}



/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
  ("\n"
   "  Usage: \n\n"
   "  -n      <level>  Number of vertices on the cube side.\n\n"
   "  -l      <level>  Cube length.\n\n"
   "  -n_part <level>  Number of partitions par process.\n\n"
   "  -post            Ensight outputs (only if n_part == 1). \n\n"
   "  -parmetis        Call ParMETIS.\n\n"
   "  -pt-scocth       Call PT-Scotch.\n\n"
   "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nx,
           PDM_g_num_t   *ny,
           PDM_g_num_t   *nz,
           int           *n_dom_i,
           int           *n_dom_j,
           int           *n_dom_k,
           int           *periodic_i,
           int           *periodic_j,
           int           *periodic_k,
           int           *n_depth,
           int           *t_elt,
           double        *length,
           int           *n_part,
           int           *with_edges,
           int           *post,
           int           *method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
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
    else if (strcmp(argv[i], "-nj") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_j = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nk") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_k = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pi") == 0) {
      *periodic_i = 1;
    }
    else if (strcmp(argv[i], "-pj") == 0) {
      *periodic_j = 1;
    }
    else if (strcmp(argv[i], "-pk") == 0) {
      *periodic_k = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
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
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-with_edges") == 0) {
      *with_edges = 1;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -t 8
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -ni 2 -pi 1
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi -pj
// Repro LS89 : mpirun -np 2 ./test/pdm_t_join_domains -nx 5 -ny 4 -nz 2 -pt-scotch -pj
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -pt-scotch -n_part 2
int main
(
 int   argc,
 char *argv[]
)
{
  /*
   *  Set default values
   */
  PDM_g_num_t          nx         = 10;
  PDM_g_num_t          ny         = 10;
  PDM_g_num_t          nz         = 10;
  int                  n_dom_i    = 1;
  int                  n_dom_j    = 1;
  int                  n_dom_k    = 1;
  int                  periodic_i = 0;
  int                  periodic_j = 0;
  int                  periodic_k = 0;
  double               length     = 1.;
  int                  n_part     = 1;
  int                  with_edges = 0;
  int                  post       = 0;
  int                  n_depth    = 1;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_HEXA8;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_QUAD4;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TETRA4;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_PRISM6;
  // 2 -> tria
  // 3 -> quad
  // 5 -> tetra
  // 6 -> pyramid
  // 7 -> prism
  // 8 -> hexa

  // PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_IMPLICIT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &n_dom_i,
             &n_dom_j,
             &n_dom_k,
             &periodic_i,
             &periodic_j,
             &periodic_k,
             &n_depth,
     (int *) &t_elt,
             &length,
             &n_part,
             &with_edges,
             &post,
     (int *) &method);


  /**
   * 
   * TOFIX:
   * 
   */


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Initialize structs
   */
  int n_domain = n_dom_i * n_dom_j * n_dom_k;

  // PDM_dcube_nodal_t **dcube = (PDM_dcube_nodal_t **) malloc(sizeof(PDM_dcube_nodal_t *) * n_domain);
  PDM_dmesh_nodal_t **dmn   = NULL;
  PDM_malloc(dmn, n_domain, PDM_dmesh_nodal_t *);
  PDM_dcube_nodal_t **dcube = NULL;
  PDM_domain_interface_t *dom_intrf = NULL;

  const int order = 1;

  PDM_dcube_nodal_cart_topo(comm,
                            n_dom_i,
                            n_dom_j,
                            n_dom_k,
                            periodic_i,
                            periodic_j,
                            periodic_k,
                            nx,
                            ny,
                            nz,
                            length,
                            0.,
                            0.,
                            0.,
                            t_elt,
                            order,
                            &dcube,
                            &dom_intrf,
                            PDM_OWNERSHIP_KEEP);

  /*
   * Partitionnement
   */
  int* n_part_by_domain = NULL;
  PDM_malloc(n_part_by_domain, n_domain, int);
  for(int i = 0; i < n_domain; ++i) {
    n_part_by_domain[i] = n_part;
  }
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                n_part_by_domain,
                                                PDM_FALSE,
                                                method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < n_domain; i++) {
    dmn[i] = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i]);
    PDM_dmesh_nodal_generate_distribution(dmn[i]);
    PDM_multipart_dmesh_nodal_set(mpart, i, dmn[i]);
  }

  // const int renum_properties_cell[6] = {1024, 0, 1, 64, 3, 1};
  // const int renum_properties_cell[6] = {12, 0, 1, 64, 3, 1};
  // const int renum_properties_cell[6] = {256, 0, 1, 64, 1, 1};
  // const int renum_properties_cell[6] = {16, 0, 1, 64, 1, 1};
  // PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_HPC",
  //                                           (void * ) renum_properties_cell,
  //                                                    "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_compute(mpart);


  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part      = NULL;
  int          **pn_cell        = NULL;
  int          **pn_face        = NULL;
  int          **pn_edge        = NULL;
  int          **pn_vtx         = NULL;
  PDM_g_num_t ***pcell_ln_to_gn = NULL;
  PDM_g_num_t ***pface_ln_to_gn = NULL;
  PDM_g_num_t ***pedge_ln_to_gn = NULL;
  PDM_g_num_t ***pvtx_ln_to_gn  = NULL;
  int         ***pcell_face_idx = NULL;
  int         ***pcell_face     = NULL;
  int         ***pcell_vtx_idx  = NULL;
  int         ***pcell_vtx      = NULL;
  int         ***pface_edge_idx = NULL;
  int         ***pface_edge     = NULL;
  int         ***pface_vtx_idx  = NULL;
  int         ***pface_vtx      = NULL;
  int         ***pedge_vtx_idx  = NULL;
  int         ***pedge_vtx      = NULL;
  double      ***pvtx_coords    = NULL;
  PDM_malloc(pn_n_part     , n_domain, int           );
  PDM_malloc(pn_cell       , n_domain, int          *);
  PDM_malloc(pn_face       , n_domain, int          *);
  PDM_malloc(pn_edge       , n_domain, int          *);
  PDM_malloc(pn_vtx        , n_domain, int          *);
  PDM_malloc(pcell_ln_to_gn, n_domain, PDM_g_num_t **);
  PDM_malloc(pface_ln_to_gn, n_domain, PDM_g_num_t **);
  PDM_malloc(pedge_ln_to_gn, n_domain, PDM_g_num_t **);
  PDM_malloc(pvtx_ln_to_gn , n_domain, PDM_g_num_t **);
  PDM_malloc(pcell_face_idx, n_domain, int         **);
  PDM_malloc(pcell_face    , n_domain, int         **);
  PDM_malloc(pcell_vtx_idx , n_domain, int         **);
  PDM_malloc(pcell_vtx     , n_domain, int         **);
  PDM_malloc(pface_edge_idx, n_domain, int         **);
  PDM_malloc(pface_edge    , n_domain, int         **);
  PDM_malloc(pface_vtx_idx , n_domain, int         **);
  PDM_malloc(pface_vtx     , n_domain, int         **);
  PDM_malloc(pedge_vtx_idx , n_domain, int         **);
  PDM_malloc(pedge_vtx     , n_domain, int         **);
  PDM_malloc(pvtx_coords   , n_domain, double      **);

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {

    /**
     * Allocate
     */
    pn_n_part       [i_dom] = n_part;
    PDM_malloc(pn_cell         [i_dom], n_part, int          );
    PDM_malloc(pn_face         [i_dom], n_part, int          );
    PDM_malloc(pn_vtx          [i_dom], n_part, int          );
    PDM_malloc(pcell_ln_to_gn  [i_dom], n_part, PDM_g_num_t *);
    PDM_malloc(pface_ln_to_gn  [i_dom], n_part, PDM_g_num_t *);
    PDM_malloc(pvtx_ln_to_gn   [i_dom], n_part, PDM_g_num_t *);
    PDM_malloc(pcell_face_idx  [i_dom], n_part, int         *);
    PDM_malloc(pcell_face      [i_dom], n_part, int         *);
    PDM_malloc(pcell_vtx_idx   [i_dom], n_part, int         *);
    PDM_malloc(pcell_vtx       [i_dom], n_part, int         *);
    PDM_malloc(pface_vtx_idx   [i_dom], n_part, int         *);
    PDM_malloc(pface_vtx       [i_dom], n_part, int         *);
    if (with_edges==1) {
      PDM_malloc(pn_edge       [i_dom], n_part, int          );
      PDM_malloc(pedge_ln_to_gn[i_dom], n_part, PDM_g_num_t *);
      PDM_malloc(pface_edge_idx[i_dom], n_part, int         *);
      PDM_malloc(pface_edge    [i_dom], n_part, int         *);
      PDM_malloc(pedge_vtx_idx [i_dom], n_part, int         *);
      PDM_malloc(pedge_vtx     [i_dom], n_part, int         *);
    }
    PDM_malloc(pvtx_coords     [i_dom], n_part, double      *);


    /**
     * Get connectivities
     */
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {
      pn_cell[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_CELL,
                                                              &pcell_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_face[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_FACE,
                                                              &pface_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      if (with_edges==1) {
        pn_edge[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                                 i_dom,
                                                                 i_part,
                                                                 PDM_MESH_ENTITY_EDGE,
                                                                &pedge_ln_to_gn[i_dom][i_part],
                                                                 PDM_OWNERSHIP_KEEP);
      }
      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VTX,
                                                             &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                         &pcell_face_idx[i_dom][i_part],
                                         &pcell_face    [i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         &pface_vtx_idx[i_dom][i_part],
                                         &pface_vtx    [i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);

      if (with_edges==1) {
        PDM_multipart_part_connectivity_get(mpart,
                                            i_dom,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                           &pface_edge_idx[i_dom][i_part],
                                           &pface_edge    [i_dom][i_part],
                                            PDM_OWNERSHIP_KEEP);

        PDM_multipart_part_connectivity_get(mpart,
                                            i_dom,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                           &pedge_vtx_idx[i_dom][i_part],
                                           &pedge_vtx    [i_dom][i_part],
                                            PDM_OWNERSHIP_KEEP);
        assert(pedge_vtx_idx[i_dom][i_part]==NULL);
      }


      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pvtx_coords[i_dom][i_part],
                                       PDM_OWNERSHIP_KEEP);


      if (post == 1) {
        char filename[999];

        if (t_elt==PDM_MESH_NODAL_HEXA8) {
          _hexa_ngon_to_nodal(pn_cell       [i_dom][i_part],
                              pcell_face_idx[i_dom][i_part],
                              pcell_face    [i_dom][i_part],
                              pface_vtx_idx [i_dom][i_part],
                              pface_vtx     [i_dom][i_part],
                             &pcell_vtx     [i_dom][i_part]);
          // pcell_vtx_idx[i_dom][i_part] = PDM_array_new_range_with_step_int(pn_cell[i_dom][i_part]+1, 8);
          // PDM_log_trace_array_int(pcell_vtx_idx[i_dom][i_part], pn_cell[i_dom][i_part]+1, "todo::");
        }
        if (pcell_vtx_idx[i_dom][i_part]!=NULL) {
          sprintf(filename, "init_cell_i_part=%i_%i_%i.vtk", i_part, i_dom, i_rank);
          PDM_vtk_write_std_elements(filename,
                                     pn_vtx        [i_dom][i_part],
                                     pvtx_coords   [i_dom][i_part],
                                     pvtx_ln_to_gn [i_dom][i_part],
                                     PDM_MESH_NODAL_HEXA8,
                                     pn_cell       [i_dom][i_part],
                                     pcell_vtx     [i_dom][i_part],
                                     pcell_ln_to_gn[i_dom][i_part],
                                     0,
                                     NULL,
                                     NULL);
        }

        sprintf(filename, "init_face_i_part=%i_%i_%i.vtk", i_part, i_dom, i_rank);
        PDM_vtk_write_polydata(filename,
                               pn_vtx        [i_dom][i_part],
                               pvtx_coords   [i_dom][i_part],
                               pvtx_ln_to_gn [i_dom][i_part],
                               pn_face       [i_dom][i_part],
                               pface_vtx_idx [i_dom][i_part],
                               pface_vtx     [i_dom][i_part],
                               pface_ln_to_gn[i_dom][i_part],
                               NULL);

        if (with_edges==1) {
          sprintf(filename, "init_edge_i_part=%i_%i_%i.vtk", i_part, i_dom, i_rank);
          PDM_vtk_write_std_elements(filename,
                                     pn_vtx        [i_dom][i_part],
                                     pvtx_coords   [i_dom][i_part],
                                     pvtx_ln_to_gn [i_dom][i_part],
                                     PDM_MESH_NODAL_BAR2,
                                     pn_edge       [i_dom][i_part],
                                     pedge_vtx     [i_dom][i_part],
                                     pedge_ln_to_gn[i_dom][i_part],
                                     0,
                                     NULL,
                                     NULL);
        }

        PDM_free(pcell_vtx[i_dom][i_part]);
      }
    } // end loop on partitions
  } // end loop on domains

  /*
   * Deduction en partition du graphe entre domaine
   */
  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   pn_face,
                                                                                   pn_edge,
                                                                                   pn_vtx,
                                                                                   pface_ln_to_gn,
                                                                                   pedge_ln_to_gn,
                                                                                   pvtx_ln_to_gn);



  /**
   * Set partition extension
   */
  PDM_extend_type_t  extend_type = PDM_EXTEND_FROM_VTX;
  PDM_part_extension_t *part_ext = PDM_part_extension_create(n_domain,
                                                             pn_n_part,
                                                             extend_type,
                                                             n_depth,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);



  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++){

      // > Get groups
      int          pn_face_group       = 0;
      int*         group_face_idx      = NULL;
      int*         group_face          = NULL;
      PDM_g_num_t* group_face_ln_to_gn = NULL;
      PDM_multipart_group_get(mpart,
                              i_dom,
                              i_part,
                              PDM_MESH_ENTITY_FACE,
                              &pn_face_group,
                              &group_face_idx,
                              &group_face,
                              &group_face_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);


      // > Get comm graph
      int *face_part_bound_proc_idx = NULL;
      int *face_part_bound_part_idx = NULL;
      int *face_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_part_bound_proc_idx,
                                        &face_part_bound_part_idx,
                                        &face_part_bound,
                                        PDM_OWNERSHIP_KEEP);

      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &vtx_part_bound_proc_idx,
                                        &vtx_part_bound_part_idx,
                                        &vtx_part_bound,
                                        PDM_OWNERSHIP_KEEP);


    // > Set connectivities in part_extension
    PDM_part_extension_connectivity_set(part_ext,
                                        i_dom,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        pcell_face_idx[i_dom][i_part],
                                        pcell_face    [i_dom][i_part]);

    if (with_edges==0) {
      PDM_part_extension_connectivity_set(part_ext,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          pface_vtx_idx[i_dom][i_part],
                                          pface_vtx    [i_dom][i_part]);
    }
    else {
      PDM_part_extension_connectivity_set(part_ext,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          pface_edge_idx[i_dom][i_part],
                                          pface_edge    [i_dom][i_part]);

      PDM_part_extension_connectivity_set(part_ext,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          NULL,
                                          pedge_vtx[i_dom][i_part]);
    }


    // > Set gnum in part_extension
    PDM_part_extension_ln_to_gn_set(part_ext,
                                    i_dom,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    pn_cell       [i_dom][i_part],
                                    pcell_ln_to_gn[i_dom][i_part]);
    PDM_part_extension_ln_to_gn_set(part_ext,
                                    i_dom,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    pn_face       [i_dom][i_part],
                                    pface_ln_to_gn[i_dom][i_part]);
    if (with_edges==1) {
      PDM_part_extension_ln_to_gn_set(part_ext,
                                      i_dom,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      pn_edge       [i_dom][i_part],
                                      pedge_ln_to_gn[i_dom][i_part]);
    }

    PDM_part_extension_ln_to_gn_set(part_ext,
                                    i_dom,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    pn_vtx       [i_dom][i_part],
                                    pvtx_ln_to_gn[i_dom][i_part]);

    PDM_part_extension_vtx_coord_set(part_ext,
                                     i_dom,
                                     i_part,
                                     pvtx_coords[i_dom][i_part]);

    }
  }

  /* Set interface */
  PDM_part_extension_part_domain_interface_shared_set(part_ext, pdi);


  PDM_part_extension_compute2(part_ext, 3);


  /**
   * Get some result
   */
  int l_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++){
      if (0) {
        log_trace("\ni_dom = %d; i_part = %d\n", i_dom, i_part);
      }

      /**
       * Get coords, connectitvity and gnum
       */
      int          extended_n_cell        = 0;
      PDM_g_num_t *extended_cell_gnum     = NULL;
      int         *extended_cell_face_idx = NULL;
      int         *extended_cell_face     = NULL;

      int          extended_n_face        = 0;
      PDM_g_num_t *extended_face_gnum     = NULL;
      int         *extended_face_vtx_idx  = NULL;
      int         *extended_face_vtx      = NULL;
      int         *extended_face_edge_idx = NULL;
      int         *extended_face_edge     = NULL;

      int          extended_n_edge        = 0;
      PDM_g_num_t *extended_edge_gnum     = NULL;
      int         *extended_edge_vtx_idx  = NULL;
      int         *extended_edge_vtx      = NULL;

      int          extended_n_vtx     = 0;
      PDM_g_num_t *extended_vtx_gnum  = NULL;
      double      *extended_vtx_coord = NULL;

      extended_n_cell = PDM_part_extension_ln_to_gn_get2(part_ext, i_dom, i_part,
                                                         PDM_MESH_ENTITY_CELL,
                                                        &extended_cell_gnum);
      PDM_part_extension_connectivity_get2(part_ext, i_dom, i_part,
                                           PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &extended_cell_face_idx,
                                          &extended_cell_face);

      if (with_edges==1) {

        extended_n_face = PDM_part_extension_ln_to_gn_get2(part_ext, i_dom, i_part,
                                                           PDM_MESH_ENTITY_FACE,
                                                          &extended_face_gnum);
        PDM_part_extension_connectivity_get2(part_ext, i_dom, i_part,
                                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &extended_face_edge_idx,
                                            &extended_face_edge);

        extended_n_edge = PDM_part_extension_ln_to_gn_get2(part_ext, i_dom, i_part,
                                                           PDM_MESH_ENTITY_EDGE,
                                                          &extended_edge_gnum);
        PDM_part_extension_connectivity_get2(part_ext, i_dom, i_part,
                                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &extended_edge_vtx_idx,
                                            &extended_edge_vtx);
        if (1 == 0) {
          log_trace("\n");
          log_trace("FACE::\n");
          log_trace("extended_n_face = %d\n", extended_n_face);
          int size_face_conn = extended_face_edge_idx[extended_n_face];
          PDM_log_trace_array_long(extended_face_gnum    , extended_n_face  , "\t extended_face_gnum     ::");
          PDM_log_trace_array_int (extended_face_edge_idx, extended_n_face+1, "\t extended_face_edge_idx ::");
          PDM_log_trace_array_int (extended_face_edge    , size_face_conn   , "\t extended_face_edge     ::");

          log_trace("\n");
          log_trace("EDGE::\n");
          log_trace("extended_n_edge = %d\n", extended_n_edge);
          int size_edge_conn = extended_edge_vtx_idx[extended_n_edge];
          PDM_log_trace_array_long(extended_edge_gnum   , extended_n_edge  , "\t extended_edge_gnum    ::");
          PDM_log_trace_array_int (extended_edge_vtx_idx, extended_n_edge+1, "\t extended_edge_vtx_idx ::");
          PDM_log_trace_array_int (extended_edge_vtx    , size_edge_conn   , "\t extended_edge_vtx     ::");
        }
      }
      else {
        extended_n_face = PDM_part_extension_ln_to_gn_get2(part_ext, i_dom, i_part,
                                                           PDM_MESH_ENTITY_FACE,
                                                          &extended_face_gnum);
        PDM_part_extension_connectivity_get2(part_ext, i_dom, i_part,
                                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &extended_face_vtx_idx,
                                            &extended_face_vtx);

        if (1 == 0) {
          log_trace("\n");
          log_trace("FACE::\n");
          log_trace("extended_n_face = %d\n", extended_n_face);
          int size_face_conn = extended_face_vtx_idx[extended_n_face];
          PDM_log_trace_array_long(extended_face_gnum   , extended_n_face  , "\t extended_face_gnum    ::");
          PDM_log_trace_array_int (extended_face_vtx_idx, extended_n_face+1, "\t extended_face_vtx_idx ::");
          PDM_log_trace_array_int (extended_face_vtx    , size_face_conn   , "\t extended_face_vtx     ::");
        }
      }

      extended_n_vtx = PDM_part_extension_ln_to_gn_get2(part_ext, i_dom, i_part,
                                                        PDM_MESH_ENTITY_VTX,
                                                       &extended_vtx_gnum);
      PDM_part_extension_vtx_coord_get2(part_ext, i_dom, i_part,
                                       &extended_vtx_coord);
      if (1 == 0) {
        log_trace("\n");
        log_trace("VTX::\n");
        log_trace("extended_n_vtx = %d\n", extended_n_vtx);
        PDM_log_trace_array_long  (extended_vtx_gnum ,   extended_n_vtx, "\t extended_vtx_gnum  ::");
        PDM_log_trace_array_double(extended_vtx_coord, 3*extended_n_vtx, "\t extended_vtx_coord ::");
      }


      /**
       * Get graph to init partition
       */
      // int *cell_cell_extended_idx  = NULL;
      int *cell_cell_extended      = NULL;
      int *cell_cell_path_itrf_idx = NULL;
      int *cell_cell_path_itrf     = NULL;
      PDM_part_extension_graph_get(part_ext, i_dom, i_part,
                                   PDM_MESH_ENTITY_CELL,
                                  // &cell_cell_extended_idx,
                                  &cell_cell_extended);
      PDM_part_extension_path_interface_get(part_ext, i_dom, i_part,
                                            PDM_MESH_ENTITY_CELL,
                                           &cell_cell_path_itrf_idx,
                                           &cell_cell_path_itrf);
      if (1 == 0) {
        log_trace("CELL::\n");
        // int n_extended = cell_cell_extended_idx[pn_cell[i_dom][i_part]];
        int n_extended = extended_n_cell;
        int n_itrf     = cell_cell_path_itrf_idx[n_extended];
        // PDM_log_trace_array_int(cell_cell_extended_idx , pn_cell[i_dom][i_part]+1, "\t cell_cell_extended_idx  ::");
        PDM_log_trace_array_int(cell_cell_extended     ,3*n_extended             , "\t cell_cell_extended      ::");
        PDM_log_trace_array_int(cell_cell_path_itrf_idx, n_extended              , "\t cell_cell_path_itrf_idx ::");
        PDM_log_trace_array_int(cell_cell_path_itrf    , n_itrf                  , "\t cell_cell_path_itrf     ::");
      }

      // int *face_face_extended_idx  = NULL;
      int *face_face_extended      = NULL;
      int *face_face_path_itrf_idx = NULL;
      int *face_face_path_itrf     = NULL;
      PDM_part_extension_graph_get(part_ext, i_dom, i_part,
                                   PDM_MESH_ENTITY_FACE,
                                  // &face_face_extended_idx,
                                  &face_face_extended);
      PDM_part_extension_path_interface_get(part_ext, i_dom, i_part,
                                            PDM_MESH_ENTITY_FACE,
                                           &face_face_path_itrf_idx,
                                           &face_face_path_itrf);
      if (1 == 0) {
        log_trace("FACE::\n");
        // int n_extended = face_face_extended_idx[pn_face[i_dom][i_part]];
        int n_extended = extended_n_face;
        int n_itrf     = face_face_path_itrf_idx[n_extended];
        // PDM_log_trace_array_int(face_face_extended_idx , pn_face[i_dom][i_part]+1, "\t face_face_extended_idx  ::");
        PDM_log_trace_array_int(face_face_extended     ,3*n_extended             , "\t face_face_extended      ::");
        PDM_log_trace_array_int(face_face_path_itrf_idx, n_extended              , "\t face_face_path_itrf_idx ::");
        PDM_log_trace_array_int(face_face_path_itrf    , n_itrf                  , "\t face_face_path_itrf     ::");
      }

      // int *edge_edge_extended_idx  = NULL;
      int *edge_edge_extended      = NULL;
      int *edge_edge_path_itrf_idx = NULL;
      int *edge_edge_path_itrf     = NULL;
      if (with_edges==1) {
        PDM_part_extension_graph_get(part_ext, i_dom, i_part,
                                     PDM_MESH_ENTITY_EDGE,
                                    // &edge_edge_extended_idx,
                                    &edge_edge_extended);
        PDM_part_extension_path_interface_get(part_ext, i_dom, i_part,
                                              PDM_MESH_ENTITY_EDGE,
                                             &edge_edge_path_itrf_idx,
                                             &edge_edge_path_itrf);
        if (1 == 0) {
          log_trace("EDGE::\n");
          // int n_extended = edge_edge_extended_idx[pn_vtx[i_dom][i_part]];
          int n_extended = extended_n_edge;
          int n_itrf     = edge_edge_path_itrf_idx[n_extended];
          // PDM_log_trace_array_int(edge_edge_extended_idx , pn_edge[i_dom][i_part]+1, "\t edge_edge_extended_idx  ::");
          PDM_log_trace_array_int(edge_edge_extended     ,3*n_extended             , "\t edge_edge_extended      ::");
          PDM_log_trace_array_int(edge_edge_path_itrf_idx, n_extended              , "\t edge_edge_path_itrf_idx ::");
          PDM_log_trace_array_int(edge_edge_path_itrf    , n_itrf                  , "\t edge_edge_path_itrf     ::");
        }
      }

      // int *vtx_vtx_extended_idx  = NULL;
      int *vtx_vtx_extended      = NULL;
      int *vtx_vtx_path_itrf_idx = NULL;
      int *vtx_vtx_path_itrf     = NULL;
      PDM_part_extension_graph_get(part_ext, i_dom, i_part,
                                   PDM_MESH_ENTITY_VTX,
                                  // &vtx_vtx_extended_idx,
                                  &vtx_vtx_extended);
      PDM_part_extension_path_interface_get(part_ext, i_dom, i_part,
                                            PDM_MESH_ENTITY_VTX,
                                           &vtx_vtx_path_itrf_idx,
                                           &vtx_vtx_path_itrf);
      if (1 == 0) {
        log_trace("VTX::\n");
        // int n_extended = vtx_vtx_extended_idx[pn_vtx[i_dom][i_part]];
        int n_extended = extended_n_vtx;
        int n_itrf     = vtx_vtx_path_itrf_idx[n_extended];
        // PDM_log_trace_array_int(vtx_vtx_extended_idx , pn_vtx[i_dom][i_part]+1, "\t vtx_vtx_extended_idx  ::");
        PDM_log_trace_array_int(vtx_vtx_extended     ,3*n_extended            , "\t vtx_vtx_extended      ::");
        PDM_log_trace_array_int(vtx_vtx_path_itrf_idx, n_extended             , "\t vtx_vtx_path_itrf_idx ::");
        PDM_log_trace_array_int(vtx_vtx_path_itrf    , n_itrf                 , "\t vtx_vtx_path_itrf     ::");
      }

      l_part++;
    }
  }


  PDM_part_extension_free(part_ext);
  /*
   *  Pour le debug : extration des faces avec le extract part + vtk
   */
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_free(pn_vtx        [i_dom]);
    PDM_free(pvtx_ln_to_gn [i_dom]);
    PDM_free(pvtx_coords   [i_dom]);

    PDM_free(pn_face       [i_dom]);
    PDM_free(pface_ln_to_gn[i_dom]);
    PDM_free(pface_vtx_idx [i_dom]);
    PDM_free(pface_vtx     [i_dom]);

    PDM_free(pn_cell       [i_dom]);
    PDM_free(pcell_ln_to_gn[i_dom]);
    PDM_free(pcell_face_idx[i_dom]);
    PDM_free(pcell_face    [i_dom]);
    PDM_free(pcell_vtx_idx [i_dom]);
    PDM_free(pcell_vtx     [i_dom]);


    if (with_edges==1) {
      PDM_free(pn_edge       [i_dom]);
      PDM_free(pedge_ln_to_gn[i_dom]);
      PDM_free(pface_edge_idx[i_dom]);
      PDM_free(pface_edge    [i_dom]);
      PDM_free(pedge_vtx_idx [i_dom]);
      PDM_free(pedge_vtx     [i_dom]);
    }
  }

  PDM_free(pn_vtx);
  PDM_free(pvtx_ln_to_gn);
  PDM_free(pn_edge);
  PDM_free(pedge_ln_to_gn);
  PDM_free(pn_face);
  PDM_free(pface_ln_to_gn);
  PDM_free(pn_cell);
  PDM_free(pcell_ln_to_gn);

  PDM_free(pcell_face_idx);
  PDM_free(pcell_face);
  PDM_free(pcell_vtx_idx);
  PDM_free(pcell_vtx);
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge);
  PDM_free(pface_vtx_idx);
  PDM_free(pface_vtx);
  PDM_free(pedge_vtx_idx);
  PDM_free(pedge_vtx);
  PDM_free(pvtx_coords);
  PDM_free(pn_n_part);

  /*
   *  Free memory
   */
  // PDM_free(i_period);

  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }
  PDM_multipart_free(mpart);

  if(pdi != NULL) {
    PDM_part_domain_interface_free(pdi);
  }
  PDM_UNUSED(pdi);

  PDM_domain_interface_free(dom_intrf);
  PDM_free(n_part_by_domain);

  PDM_free(dcube);
  PDM_free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
