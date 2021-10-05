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
                                  0,
                                  NULL,
                                  NULL);

    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);
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

  /*
   *  Create distributed cube
   */
  PDM_dcube_nodal2_t *dcube = PDM_dcube_nodal_gen2_init(comm,
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

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen2_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);


  /* Deform */
  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);

  for (int i = 0; i < dn_vtx; i++) {
    double x = (dvtx_coord[3*i    ] - 0.5) / length;
    double y = (dvtx_coord[3*i + 1] - 0.5) / length;
    double z = (dvtx_coord[3*i + 2] - 0.5) / length;

    if (dim == 2) {
      dvtx_coord[3*i + 2] = length * (x*x + y*y);
    } else {
      dvtx_coord[3*i    ] += 0.1*length*cos(3*y);
      dvtx_coord[3*i + 1] += 0.1*length*cos(3*z);
      dvtx_coord[3*i + 2] += 0.1*length*cos(3*x);
    }
  }

  if (dim == 3) {
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic");
  } else {
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }
  //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
  //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_CORNER, "out_corner");


  //PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  gettimeofday(&t_elaps_debut, NULL);
  //PDM_dcube_nodal_gen2_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
