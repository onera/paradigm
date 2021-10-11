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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_nodal_gen2.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

#define _MIN(a,b) ((a) < (b) ? (a) : (b))
#define _MAX(a,b) ((a) > (b) ? (a) : (b))


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
           int           *order,
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
    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = atoi(argv[i]);
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
_dmesh_nodal_dump_vtk
(
       PDM_dmesh_nodal_t   *dmn,
       int                  order,
       PDM_geometry_kind_t  geom_kind,
 const char                *filename_patter
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int* sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int n_section    = PDM_DMesh_nodal_n_section_get(dmn, geom_kind);


  const char *field_name = "group";
  int n_field = 0;
  _pdm_dmesh_nodal_elts_t *dmne = NULL;
  int *delt_group_idx = NULL;
  int *delt_group     = NULL;
  double **field = NULL;
  if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
    dmne = dmn->ridge;
  } else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC && dmn->mesh_dimension == 3) {
    dmne = dmn->surfacic;
  }

  PDM_g_num_t *distrib_elt = NULL;
  if (dmne != NULL) {
    distrib_elt = PDM_compute_uniform_entity_distribution(dmn->comm,
                                                          dmne->n_g_elmts);
    PDM_log_trace_array_long(distrib_elt, n_rank+1, "distrib_elt : ");

    PDM_dgroup_entity_transpose(dmne->n_group_elmt,
                                dmne->dgroup_elmt_idx,
                                dmne->dgroup_elmt,
                (PDM_g_num_t *) distrib_elt,
                                &delt_group_idx,
                                &delt_group,
                                dmn->comm);
  }

  PDM_g_num_t shift = 0;
  for(int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);

    int         *dconnec_idx    = (int         * ) malloc( (n_elt+1) * sizeof(int        ));
    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));

    int strid = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    dconnec_idx[0] = 0;
    for(int i = 0; i < n_elt; ++i) {
      dconnec_idx[i+1] = dconnec_idx[i] + strid;
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    log_trace("section %d (%d) :\n", i_section, id_section);
    PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn : ");

    PDM_g_num_t *pvtx_ln_to_gn;
    int         *pcell_vtx_idx;
    int         *pcell_vtx;
    int          pn_vtx;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(dmn->comm,
                                                             delmt_distribution,
                                                             dconnec_idx,
                                                             dconnec,
                                                             n_elt,
                                    (const PDM_g_num_t *)    delmt_ln_to_gn,
                                                            &pn_vtx,
                                                            &pvtx_ln_to_gn,
                                                            &pcell_vtx_idx,
                                                            &pcell_vtx);

    /*
     * Coordinates
     */
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    // int          dn_vtx   = PDM_DMesh_nodal_n_vtx_get(dln->dmesh_nodal_in);
    // assert(dn_vtx == (vtx_distrib[i_rank+1]-vtx_distrib[i_rank]));
    double** tmp_pvtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(dmn->comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_vtx,
                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                          &tmp_pvtx_coord);

    double* pvtx_coord_out = tmp_pvtx_coord[0];

    /*
     * Groups
     */
    if (dmne != NULL) {
      for(int i = 0; i < n_elt; ++i) {
        delmt_ln_to_gn[i] += shift;
      }

      int **tmp_elt_group_idx = NULL;
      int **tmp_elt_group     = NULL;
      PDM_part_dentity_group_to_pentity_group(dmn->comm,
                                              1,
                                              distrib_elt,
                                              delt_group_idx,
                                              delt_group,
                                              &n_elt,
                      (const PDM_g_num_t **)  &delmt_ln_to_gn,
                                              &tmp_elt_group_idx,
                                              &tmp_elt_group);
      int *pelt_group_idx = tmp_elt_group_idx[0];
      int *pelt_group     = tmp_elt_group    [0];
      PDM_log_trace_connectivity_int(pelt_group_idx, pelt_group, n_elt, "pelt_group : ");
      free (tmp_elt_group_idx);
      free (tmp_elt_group);
      PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn (shifted) : ");

      n_field = 1;
      field = malloc (sizeof(double *) * n_field);
      field[0] = malloc (sizeof(double) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        assert (pelt_group_idx[i+1] == pelt_group_idx[i] + 1);
        field[0][i] = (double) pelt_group[i];
      }
      free (pelt_group);
      free (pelt_group_idx);
    }

    /*
     *  Dump
     */
    char filename[999];
    sprintf(filename, "%s_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    PDM_vtk_write_std_elements_ho(filename,
                                  order,
                                  pn_vtx,
                                  pvtx_coord_out,
                                  pvtx_ln_to_gn,
                                  t_elt,
                                  n_elt,
                                  pcell_vtx,
                                  delmt_ln_to_gn,
                                  n_field,
                (const char   **) &field_name,
                (const double **) field);

    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];

    if (dmne != NULL) {
      free (field[0]);
      free (field);
    }
  }

  if (dmne != NULL) {
    free (delt_group_idx);
    free (delt_group);
    free (distrib_elt);
  }
}


static void
_lagrange_to_bezier_tria
(
 const int     order,
 const double *lag,
 double       *bez
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIA3, order);
  //double *bez = malloc (sizeof(double) * n_nodes * 3);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {
    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = -0.5*lag[j] + 2*lag[9+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -0.5*lag[6+j] + 2*lag[12+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = lag[15+j];
    }
  }

  else if (order == 3) {
    double f5_6 = 5. / 6.;
    double f1_3 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f5_6*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f1_3*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f1_3*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f5_6*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -f5_6*lag[j] + 3*lag[12+j] - 1.5*lag[21+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = f1_3*lag[j] - 0.75*lag[3+j] - 0.75*lag[6+j] + f1_3*lag[9+j] - 0.75*lag[12+j] + 4.5*lag[15+j] - 0.75*lag[18+j] - 0.75*lag[21+j] - 0.75*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = -f5_6*lag[9+j] + 3*lag[18+j] - 1.5*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = f1_3*lag[j] - 1.5*lag[12+j] + 3*lag[21+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = f1_3*lag[9+j] - 1.5*lag[18+j] + 3*lag[24+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[27+j] = lag[27+j];
    }
  }
  /*int *M = malloc (sizeof(double) * n_nodes * n_nodes);

  int idx = 0;
  if (order == 1) {
    M[idx++] = 1; M[idx++] = 0; M[idx++] = 0;
    M[idx++] = 0; M[idx++] = 1; M[idx++] = 0;
    M[idx++] = 0; M[idx++] = 0; M[idx++] = 1;
  }
  else if (order == 2) {
    M[idx++] =  1;   M[idx++] = 0; M[idx++] =  0;   M[idx++] = 0; M[idx++] = 0; M[idx++] =  0;
    M[idx++] = -0.5; M[idx++] = 2; M[idx++] = -0.5; M[idx++] = 0; M[idx++] = 0; M[idx++] =  0;
    M[idx++] =  0.;  M[idx++] = 0; M[idx++] =  1;   M[idx++] = 0; M[idx++] = 0; M[idx++] =  0;
    M[idx++] = -0.5; M[idx++] = 0; M[idx++] =  0;   M[idx++] = 2; M[idx++] = 0; M[idx++] = -0.5;
    M[idx++] =  0;   M[idx++] = 0; M[idx++] = -0.5; M[idx++] = 0; M[idx++] = 2; M[idx++] = -0.5;
    M[idx++] =  1;   M[idx++] = 0; M[idx++] =  0;   M[idx++] = 0; M[idx++] = 0; M[idx++] =  1;
  }
  else if (order == 3) {
    M[idx++] =  1;   M[idx++] = 0; M[idx++] =  0;   M[idx++] = 0; M[idx++] = 0; M[idx++] =  0;
  }*/
}


static void
_bezier_bounding_boxes
(
 PDM_dmesh_nodal_t   *dmn,
 int                  order,
 PDM_geometry_kind_t  geom_kind
 )
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, geom_kind);

  PDM_g_num_t shift = 0;
  for (int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);

    if (t_elt != PDM_MESH_NODAL_TRIA3 || order > 3) continue;

    int         *dconnec_idx    = (int         * ) malloc( (n_elt+1) * sizeof(int        ));
    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));

    int strid = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    dconnec_idx[0] = 0;
    for(int i = 0; i < n_elt; ++i) {
      dconnec_idx[i+1] = dconnec_idx[i] + strid;
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    int *ijk_to_vtk = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                      t_elt,
                                                      order);


    PDM_g_num_t *pvtx_ln_to_gn;
    int         *pcell_vtx_idx;
    int         *pcell_vtx;
    int          pn_vtx;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(dmn->comm,
                                                             delmt_distribution,
                                                             dconnec_idx,
                                                             dconnec,
                                                             n_elt,
                                    (const PDM_g_num_t *)    delmt_ln_to_gn,
                                                            &pn_vtx,
                                                            &pvtx_ln_to_gn,
                                                            &pcell_vtx_idx,
                                                            &pcell_vtx);

    /*
     * Coordinates
     */
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    // int          dn_vtx   = PDM_DMesh_nodal_n_vtx_get(dln->dmesh_nodal_in);
    // assert(dn_vtx == (vtx_distrib[i_rank+1]-vtx_distrib[i_rank]));
    double** tmp_pvtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(dmn->comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_vtx,
                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                          &tmp_pvtx_coord);

    double* pvtx_coord_out = tmp_pvtx_coord[0];

    //double *bezier_coord = _lagrange_to_bezier_tria (order, pvtx_coord_out);
    int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    double *lagrange_coord = malloc (sizeof(double) * n_nodes * 3);
    double *bezier_coord   = malloc (sizeof(double) * n_nodes * 3);
    double *elt_coord      = malloc (sizeof(double) * n_elt * n_nodes * 3);
    int    *elt_vtx        = malloc (sizeof(int)    * n_elt * n_nodes);

    double *extents = malloc (sizeof(double) * n_elt * 6);
    int idx2 = 0;
    for (int i = 0; i < n_elt; i++) {
      double *_min = extents + 6*i;
      double *_max = _min + 3;

      for (int j = 0; j < 3; j++) {
        _min[j] =  1e30;
        _max[j] = -1e30;
      }

      int idx = 0;
      for (int k = pcell_vtx_idx[i]; k < pcell_vtx_idx[i+1]; k++) {
        int ivtx = pcell_vtx[k] - 1;
        for (int j = 0; j < 3; j++) {
          lagrange_coord[idx++] = pvtx_coord_out[3*ivtx + j];
        }
      }

      _lagrange_to_bezier_tria (order, lagrange_coord, bezier_coord);

      for (int k = 0; k < n_nodes; k++) {
        for (int j = 0; j < 3; j++) {
          elt_coord[3*idx2 + j] = bezier_coord[3*k + j];//lagrange_coord[3*k + j];//
          _min[j] = _MIN(_min[j], bezier_coord[3*k + j]);
          _max[j] = _MAX(_max[j], bezier_coord[3*k + j]);
        }
        elt_vtx[n_nodes*i + ijk_to_vtk[k]] = ++idx2;
      }
    }

    /*
     *  Dump
     */
    char filename[999];
    sprintf(filename, "bezier_section_%2.2d_%2.2d.vtk", i_section, i_rank);
    PDM_vtk_write_std_elements_ho(filename,
                                  order,
                                  n_elt * n_nodes,
                                  elt_coord,
                                  NULL,
                                  t_elt,
                                  n_elt,
                                  elt_vtx,
                                  delmt_ln_to_gn,
                                  0,
                                  NULL,
                                  NULL);
    free(elt_vtx);
    free(elt_coord);

    sprintf(filename, "boxes_section_%2.2d_%2.2d.vtk", i_section, i_rank);
    PDM_vtk_write_boxes(filename,
                        n_elt,
                        extents,
                        delmt_ln_to_gn);
    free (extents);
    free (bezier_coord);
    free (lagrange_coord);

    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];
  }
}



/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */
  PDM_g_num_t          nx     = 10;
  PDM_g_num_t          ny     = 10;
  PDM_g_num_t          nz     = 10;
  int                  order  = 1;
  double               length = 1.;
  int                  n_part = 1;
  int                  post   = 0;
  PDM_Mesh_nodal_elt_t t_elt  = PDM_MESH_NODAL_TRIA3;
  // 2 -> tria
  // 3 -> quad
  // 5 -> tetra
  // 6 -> pyramid
  // 7 -> prism
  // 8 -> hexa

  #ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif

/*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &order,
             (int *) &t_elt,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */
  struct timeval t_elaps_debut;

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dim = 3;
  if (t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_QUAD4 ) {
    dim = 2;
  }

  if (order > 3) {
    int *ijk = NULL;

    for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_BAR2;
         type <= PDM_MESH_NODAL_HEXA8;
         type++) {

      if (type == PDM_MESH_NODAL_POLY_2D ||
          type == PDM_MESH_NODAL_PYRAMID5) continue;

      ijk = PDM_vtk_lagrange_to_ijk(type, order);
      PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                       type,
                                       order,
                                       PDM_Mesh_nodal_n_vtx_elt_get(type, order),
                                       ijk);
      free (ijk);
    }
  }

  /*
   *  Create distributed cube
   */
  PDM_dcube_nodal2_t *dcube = PDM_dcube_nodal_gen2_create(comm,
                                                          nx,
                                                          ny,
                                                          nz,
                                                          length,
                                                          0.,
                                                          0.,
                                                          0.,
                                                          t_elt,
                                                          order,
                                                          PDM_OWNERSHIP_KEEP);

  /*PDM_dcube_nodal_gen2_ordering_set (dcube,
    "PDM_HO_ORDERING_VTK");*/

  PDM_dcube_nodal_gen2_build (dcube);

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen2_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);


  /* Deform */
  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  double amplitude = 0.2;
  double frequence = 3.;

  for (int i = 0; i < dn_vtx; i++) {
    double x = (dvtx_coord[3*i    ] - 0.5) / length;
    double y = (dvtx_coord[3*i + 1] - 0.5) / length;
    double z = (dvtx_coord[3*i + 2] - 0.5) / length;

    //double scale = length * pow(2, order-1);

    if (dim == 2) {
      //dvtx_coord[3*i + 2] = scale * (pow(x, order) + pow(y, order));
      dvtx_coord[3*i + 2] = length * (x*x + y*y);
    } else {
      dvtx_coord[3*i    ] += amplitude*length*cos(frequence*y);
      dvtx_coord[3*i + 1] += amplitude*length*cos(frequence*z);
      dvtx_coord[3*i + 2] += amplitude*length*cos(frequence*x);
    }
  }

  /* Bounding boxes */
  //_bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_SURFACIC);


  /* Reorder */
  PDM_dmesh_nodal_reorder (dmn,
                           "PDM_HO_ORDERING_VTK",
                           order);

  if (dim == 3) {
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic");
  }
  //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
  //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_RIDGE,    "out_ridge");
  _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_CORNER,   "out_corner");


  //PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  gettimeofday(&t_elaps_debut, NULL);
  PDM_dcube_nodal_gen2_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
