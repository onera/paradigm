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
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_extension.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_vtk.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_logging.h"
#include "pdm_order.h"
#include "pdm_binary_search.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_domain_interface.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
        *t_elt = atoi(argv[i]);
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

static void
_compute_face_vtx
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
 )
{
  int dbg = 0;

  *face_vtx = malloc (sizeof(int) * face_edge_idx[n_face]);

  for(int i_face = 0; i_face <  face_edge_idx[n_face]; ++i_face) {
    (*face_vtx)[i_face] = 10000;
  }

  int n_edge = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_edge = PDM_MAX(n_edge, PDM_ABS(face_edge[i]));
  }


  int *edge_tag = PDM_array_zeros_int(n_edge);

  for (int iface = 0; iface < n_face; iface++) {
    int *_face_vtx  = *face_vtx  + face_edge_idx[iface];
    int *_face_edge =  face_edge + face_edge_idx[iface];

    if (dbg) {
      log_trace("\nFace %d\n", iface);
      for (int idx_edge = face_edge_idx[iface]; idx_edge < face_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  face_edge[idx_edge],
                  edge_vtx[2*iedge], edge_vtx[2*iedge+1]);
      }
    }

    int _n_edge = face_edge_idx[iface+1] - face_edge_idx[iface];
    // first edge
    int iedge = PDM_ABS(_face_edge[0]) - 1;
    edge_tag[iedge] = 1;
    _face_vtx[0] = PDM_ABS(edge_vtx[2*iedge  ]);
    _face_vtx[1] = PDM_ABS(edge_vtx[2*iedge+1]);

    for (int i = 2; i < _n_edge; i++) {

      for (int j = 1; j < _n_edge; j++) {
        iedge = PDM_ABS(_face_edge[j]) - 1;

        if (edge_tag[iedge]) {
          continue;
        }

        if (edge_vtx[2*iedge] == _face_vtx[i-1]) {
          _face_vtx[i] = PDM_ABS(edge_vtx[2*iedge+1]);
          edge_tag[iedge] = 1;
          break;
        }
        else if (edge_vtx[2*iedge+1] == _face_vtx[i-1]) {
          _face_vtx[i] = PDM_ABS(edge_vtx[2*iedge]);
          edge_tag[iedge] = 1;
          break;
        }
      }
    }

    if (dbg) {
      log_trace("  face_vtx = ");
      for (int ivtx = 0; ivtx < face_edge_idx[iface+1] - face_edge_idx[iface]; ivtx++) {
        log_trace("%d ", _face_vtx[ivtx]);
      }
      log_trace("\n");
    }

    // reset tags
    for (int i = 0; i < _n_edge; i++) {
      iedge = PDM_ABS(_face_edge[i]) - 1;
      edge_tag[iedge] = 0;
    }
  }
  free(edge_tag);
}


static
PDM_g_num_t*
_compute_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_MPI_Comm     comm
)
{

  PDM_g_num_t *shift_by_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *shift_by_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        shift_by_domain_loc[i_domain] = PDM_MAX(shift_by_domain_loc[i_domain], _pentity_ln_to_gn[i]);
      }
    }
  }

  shift_by_domain[0] = 0;
  PDM_MPI_Allreduce(shift_by_domain_loc, &shift_by_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(shift_by_domain, n_domain+1);

  free(shift_by_domain_loc);

  return shift_by_domain;
}

// static
// void
// _shift_ln_to_gn
// (
//   int          n_entity,
//   PDM_g_num_t *entity_ln_to_gn,
//   PDM_g_num_t  shift,
//   int          sens
// )
// {
//   for(int i = 0; i < n_entity; ++i) {
//     entity_ln_to_gn[i] = entity_ln_to_gn[i] + sens * shift;
//   }
// }

static
void
_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_g_num_t     *shift_by_domain,
  int              sens
)
{
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        _pentity_ln_to_gn[i] = _pentity_ln_to_gn[i] + sens * shift_by_domain[i_domain];
      }
    }
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

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int* n_part_g = malloc(n_domain * sizeof(int));
  printf("n_domain = %i \n", n_domain);
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int          **pn_vtx              = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge             = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_face             = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn       = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn      = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pface_ln_to_gn      = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));

  int           *pflat_n_vtx         = (int           *) malloc( n_domain * sizeof(int           ));
  int           *pflat_n_edge        = (int           *) malloc( n_domain * sizeof(int           ));
  int           *pflat_n_face        = (int           *) malloc( n_domain * sizeof(int           ));

  PDM_g_num_t  **pflat_vtx_ln_to_gn  = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  PDM_g_num_t  **pflat_edge_ln_to_gn = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  PDM_g_num_t  **pflat_face_ln_to_gn = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  double       **pflat_vtx_coords    = (double       **) malloc( n_domain * sizeof(double       *));


  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx        [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pn_edge       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pn_face       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn [i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));
    pedge_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));
    pface_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

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

      pflat_n_face      [ln_part_tot+i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                                               i_dom,
                                                                               i_part,
                                                                               PDM_MESH_ENTITY_FACE,
                                                                               &pflat_face_ln_to_gn[ln_part_tot+i_part],
                                                                               PDM_OWNERSHIP_KEEP);

      pn_face       [i_dom][i_part] = pflat_n_face       [ln_part_tot+i_part];
      pface_ln_to_gn[i_dom][i_part] = pflat_face_ln_to_gn[ln_part_tot+i_part];

      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pflat_vtx_coords[ln_part_tot+i_part],
                                       PDM_OWNERSHIP_KEEP);

    }
    ln_part_tot += n_part[i_dom];
  }


  // int         **pflat_vtx_edge_idx = NULL;
  // int         **pflat_vtx_edge     = NULL;
  // PDM_part_connectivity_transpose(ln_part_tot,
  //                                 pflat_n_edge,
  //                                 pflat_n_vtx,
  //                                 pflat_edge_vtx_idx,
  //                                 pflat_edge_vtx,
  //                                 &pflat_vtx_edge_idx,
  //                                 &pflat_vtx_edge);

  PDM_g_num_t* shift_by_domain_vtx = _compute_offset_ln_to_gn_by_domain(n_domain,
                                                                        n_part,
                                                                        pn_vtx,
                                                                        pvtx_ln_to_gn,
                                                                        comm);

  PDM_g_num_t* shift_by_domain_edge = _compute_offset_ln_to_gn_by_domain(n_domain,
                                                                         n_part,
                                                                         pn_edge,
                                                                         pedge_ln_to_gn,
                                                                         comm);

  PDM_g_num_t* shift_by_domain_face = _compute_offset_ln_to_gn_by_domain(n_domain,
                                                                         n_part,
                                                                         pn_face,
                                                                         pface_ln_to_gn,
                                                                         comm);

  /* Shift ln_to_gn */
  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                             shift_by_domain_vtx,
                             1);

  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_edge,
                             pedge_ln_to_gn,
                             shift_by_domain_edge,
                             1);

  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_face,
                             pface_ln_to_gn,
                             shift_by_domain_face,
                             1);

  int  *pn_vtx_num             = NULL;
  int **pvtx_num               = NULL;
  int **pvtx_opp_location_idx  = NULL;
  int **pvtx_opp_location      = NULL;
  int **pvtx_opp_interface_idx = NULL;
  int **pvtx_opp_interface     = NULL;
  int **pvtx_opp_sens          = NULL;

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
                                         &pvtx_opp_sens);


  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pvtx_num[i_part], pn_vtx_num[i_part], "pvtx_num ::");
      PDM_log_trace_graph_nuplet_int(pvtx_opp_location_idx[i_part],
                                     pvtx_opp_location    [i_part],
                                     3,
                                     pn_vtx_num[i_part], "pvtx_opp_location ::");
      PDM_log_trace_array_int(pvtx_opp_interface[i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_interface ::");
      PDM_log_trace_array_int(pvtx_opp_sens     [i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_sens ::");
    }
  }

  int  *pn_edge_num             = NULL;
  int **pedge_num               = NULL;
  int **pedge_opp_location_idx  = NULL;
  int **pedge_opp_location      = NULL;
  int **pedge_opp_interface_idx = NULL;
  int **pedge_opp_interface     = NULL;
  int **pedge_opp_sens          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         PDM_BOUND_TYPE_EDGE,
                                         pflat_n_edge,
                                         pflat_edge_ln_to_gn,
                                         &pn_edge_num,
                                         &pedge_num,
                                         &pedge_opp_location_idx,
                                         &pedge_opp_location,
                                         &pedge_opp_interface_idx,
                                         &pedge_opp_interface,
                                         &pedge_opp_sens);


  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pedge_num[i_part], pn_edge_num[i_part], "pedge_num ::");
      PDM_log_trace_graph_nuplet_int(pedge_opp_location_idx[i_part],
                                     pedge_opp_location    [i_part],
                                     3,
                                     pn_edge_num[i_part], "pedge_opp_location ::");
      PDM_log_trace_array_int(pedge_opp_interface[i_part], pedge_opp_location_idx[i_part][pn_edge_num[i_part]], "pedge_opp_interface ::");
      PDM_log_trace_array_int(pedge_opp_sens     [i_part], pedge_opp_location_idx[i_part][pn_edge_num[i_part]], "pedge_opp_sens ::");
    }
  }


  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    free(pn_vtx        [i_dom]);
    free(pn_edge       [i_dom]);
    free(pn_face       [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pface_ln_to_gn[i_dom]);
  }
  free(pn_vtx        );
  free(pn_edge       );
  free(pn_face       );
  free(pvtx_ln_to_gn );
  free(pedge_ln_to_gn);
  free(pface_ln_to_gn);

  free(pflat_vtx_ln_to_gn );
  free(pflat_edge_ln_to_gn);
  free(pflat_face_ln_to_gn);
  free(pflat_vtx_coords   );

  free(pflat_n_vtx );
  free(pflat_n_edge);
  free(pflat_n_face);

  free(shift_by_domain_vtx);
  free(shift_by_domain_edge);
  free(shift_by_domain_face);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pedge_num             [i_part]);
    free(pedge_opp_location_idx[i_part]);
    free(pedge_opp_location    [i_part]);
    free(pedge_opp_interface   [i_part]);
    free(pedge_opp_sens        [i_part]);
  }
  free(pn_edge_num            );
  free(pedge_num              );
  free(pedge_opp_location_idx );
  free(pedge_opp_location     );
  free(pedge_opp_interface_idx);
  free(pedge_opp_interface    );
  free(pedge_opp_sens         );


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_num             [i_part]);
    free(pvtx_opp_location_idx[i_part]);
    free(pvtx_opp_location    [i_part]);
    free(pvtx_opp_interface   [i_part]);
    free(pvtx_opp_sens        [i_part]);
  }
  free(pn_vtx_num            );
  free(pvtx_num              );
  free(pvtx_opp_location_idx );
  free(pvtx_opp_location     );
  free(pvtx_opp_interface_idx);
  free(pvtx_opp_interface    );
  free(pvtx_opp_sens         );

  free(n_part_g);
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
  int                  post       = 0;
  int                  n_depth    = 1;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_HEXA8;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_QUAD4;
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
             &post,
     (int *) &method);

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
  int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);

  if (dim == 2) {
    n_dom_k    = 1;
    nz         = 1;
    periodic_k = 0;
  }

  int n_domain = n_dom_i * n_dom_j * n_dom_k;

  // PDM_dcube_nodal_t **dcube = (PDM_dcube_nodal_t **) malloc(sizeof(PDM_dcube_nodal_t *) * n_domain);
  PDM_dmesh_nodal_t **dmn   = (PDM_dmesh_nodal_t **) malloc(sizeof(PDM_dmesh_nodal_t *) * n_domain);
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
  int* n_part_by_domain = (int *) malloc(n_domain * sizeof(int));
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
    PDM_multipart_register_dmesh_nodal(mpart, i, dmn[i]);
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
  PDM_multipart_run_ppart(mpart);


  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part      = (int           *) malloc( n_domain * sizeof(int          *));
  int          **pn_face        = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge        = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int          **pn_vtx         = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int         ***pface_vtx      = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pedge_vtx      = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pedge_vtx_idx  = (int         ***) malloc( n_domain * sizeof(int         **));
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    pn_n_part     [i_dom] = n_part;
    pn_face       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pn_edge       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pface_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pedge_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pn_vtx        [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pvtx_ln_to_gn [i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pface_vtx     [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pedge_vtx     [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pedge_vtx_idx [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {
      pn_face[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_FACE,
                                                               &pface_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_edge[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_EDGE,
                                                               &pedge_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &pedge_vtx     [i_dom][i_part],
                                          &pedge_vtx_idx [i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);

      assert(pedge_vtx_idx [i_dom][i_part] == NULL);

      pedge_vtx_idx [i_dom][i_part] = malloc((pn_edge[i_dom][i_part]+1) * sizeof(int));
      for(int i_edge = 0; i_edge < pn_edge[i_dom][i_part]+1; ++i_edge) {
        pedge_vtx_idx [i_dom][i_part][i_edge] = 2*i_edge;
      }


    }
  }

  /*
   * Deduction en partition du graphe entre domaine
   */

  // PDM_domain_interface_translate_vtx2edge(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_edge,
  //                                         dedge_vtx_idx,
  //                                         dedge_vtx);


  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   pn_face,
                                                                                   pn_edge,
                                                                                   pn_vtx,
                                                                                   pface_ln_to_gn,
                                                                                   pedge_ln_to_gn,
                                                                                   pvtx_ln_to_gn);

  PDM_part_domain_interface_add(pdi,
                                PDM_BOUND_TYPE_VTX,
                                PDM_BOUND_TYPE_EDGE,
                                pn_n_part,
                                pn_vtx,
                                pvtx_ln_to_gn,
                                pn_edge,
                                pedge_ln_to_gn,
                                pedge_vtx_idx,
                                pedge_vtx,
                                1); // Connectivity_is_signed

  _part_extension(n_depth,
                  n_domain,
                  pn_n_part,
                  comm,
                  mpart,
                  pdi);

  /*
   * Extension
   */
  int shift_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++){

      PDM_g_num_t* cell_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_dom,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      int *cell_face     = NULL;
      int *cell_face_idx = NULL;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       i_dom,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face,
                                                       &cell_face_idx,
                                                       PDM_OWNERSHIP_KEEP);

      int *face_vtx     = NULL;
      int *face_vtx_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx,
                                          &face_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* face_ln_to_gn = NULL;
      int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                   i_dom,
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE,
                                                   &face_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  i_dom,
                                                  i_part,
                                                  PDM_MESH_ENTITY_VERTEX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);
      PDM_UNUSED(n_cell);
      PDM_UNUSED(n_vtx);
      PDM_UNUSED(n_face);

      double *vtx = NULL;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &vtx,
                                       PDM_OWNERSHIP_KEEP);


      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_BOUND_TYPE_VTX,
                                        &vtx_part_bound_proc_idx,
                                        &vtx_part_bound_part_idx,
                                        &vtx_part_bound,
                                        PDM_OWNERSHIP_KEEP);

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face2 = PDM_multipart_part_connectivity_get(mpart,
                                                        i_dom,
                                                        i_part,
                                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                        &face_edge,
                                                        &face_edge_idx,
                                                        PDM_OWNERSHIP_KEEP);
      // assert(n_face == n_face2);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                       i_dom,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);
      PDM_UNUSED(n_edge);
      assert(edge_vtx_idx == NULL);

      int *face_cell     = NULL;
      int *face_cell_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                          &face_cell,
                                          &face_cell_idx,
                                          PDM_OWNERSHIP_KEEP);
      assert(face_cell_idx == NULL);

      PDM_g_num_t* edge_ln_to_gn = NULL;
      int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    i_dom,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &edge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);

      PDM_UNUSED(n_face2);
      PDM_UNUSED(n_edge2);

      int          pn_face_group       = 0;
      int*         group_face_idx      = NULL;
      int*         group_face          = NULL;
      PDM_g_num_t* group_face_ln_to_gn = NULL;
      PDM_multipart_bound_get(mpart,
                              i_dom,
                              i_part,
                              PDM_BOUND_TYPE_FACE,
                              &pn_face_group,
                              &group_face_idx,
                              &group_face,
                              &group_face_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);

      // Concatenate face_vtx
      pface_vtx[i_dom][i_part] = NULL;
      _compute_face_vtx(n_face2,
                        face_edge_idx,
                        face_edge,
                        edge_vtx,
                        &pface_vtx[i_dom][i_part]);
      free(pface_vtx[i_dom][i_part]);
      // _compute_face_vtx2(n_face2,
      //                   face_edge_idx,
      //                   face_edge,
      //                   edge_vtx,
      //                   &pface_vtx[i_dom][i_part]);

      int *face_part_bound_proc_idx = NULL;
      int *face_part_bound_part_idx = NULL;
      int *face_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_BOUND_TYPE_FACE,
                                        &face_part_bound_proc_idx,
                                        &face_part_bound_part_idx,
                                        &face_part_bound,
                                        PDM_OWNERSHIP_KEEP);

    }
    shift_part += pn_n_part[i_dom];
  }
  PDM_UNUSED(shift_part);

  /*
   *  Pour le debug : extration des faces avec le extract part + vtk
   */
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {

    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++){
      free(pedge_vtx_idx [i_dom][i_part]);
    }
    free(pn_face       [i_dom]);
    free(pn_edge       [i_dom]);
    free(pface_ln_to_gn[i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pn_vtx        [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
    free(pface_vtx     [i_dom]);
    free(pedge_vtx     [i_dom]);
    free(pedge_vtx_idx [i_dom]);
  }
  free(pn_face       );
  free(pn_edge       );
  free(pface_ln_to_gn);
  free(pedge_ln_to_gn);
  free(pn_vtx       );
  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(pface_vtx);
  free(pn_n_part);

  /*
   *  Free memory
   */
  // free(i_period);

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
  free(n_part_by_domain);

  free(dcube);
  free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
