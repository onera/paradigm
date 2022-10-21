#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_block_to_block.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh_nodal.h"

#include "pdm_vtk.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

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
           PDM_g_num_t   *n_vtx_seg,
           char         **filename,
           PDM_g_num_t   *n_g_pts_clouds,
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
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_g_pts_clouds = atol(argv[i]);
        *n_g_pts_clouds = (PDM_g_num_t) _n_g_pts_clouds;
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
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static void
_read_surface_mesh
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 double       **dvtx_coord,
 PDM_g_num_t  **distrib_vtx,
 PDM_g_num_t  **distrib_face
)
{

  PDM_UNUSED(comm);
  PDM_UNUSED(filename);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(dvtx_coord);
  PDM_UNUSED(distrib_vtx);
  PDM_UNUSED(distrib_face);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // assert(n_rank == 1);
  int n_vtx  = 0;
  int n_face = 0;
  double      *face_normal    = NULL;
  double      *vtx_coord      = NULL;
  int         *face_vtx_idx   = NULL;
  int         *face_vtx_n     = NULL;
  // PDM_g_num_t *face_vtx       = NULL;
  int *face_vtx       = NULL;
  double      *face_vtx_coord = NULL;

  if (i_rank == 0) {
    FILE *f = fopen(filename, "r");

    if (f == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
    }

    char header[999];
    fscanf(f, "%[^\n]", header);

    printf("header = %s \n", header);

    int end = 0;
    while(end == 0) {

      int stat = fscanf(f, "%s", header);

      // printf("header = %s \n", header);

      if(strcmp(header, "vertex") == 0) {
        n_vtx += 1;
      }
      if(strcmp(header, "facet") == 0) {
        n_face += 1;
      }

      // if(strcmp(header, "OpenSCAD_Model") == 0) {
      //   end = 1;
      // }
      if(stat == EOF) {
        end = 1;
      }
    }

    printf("n_face = %i \n", n_face);
    printf("n_vtx  = %i \n", n_vtx );


    face_normal    = (double * ) malloc(3 * n_face    * sizeof(double));
    vtx_coord      = (double * ) malloc(3 * n_vtx     * sizeof(double));
    face_vtx_idx   = (int    * ) malloc( (n_face + 1) * sizeof(int   ));
    face_vtx_n     = (int    * ) malloc( (n_face + 1) * sizeof(int   ));
    face_vtx       = (int    * ) malloc(3 * n_face    * sizeof(int   ));
    face_vtx_coord = (double * ) malloc(9 * n_face    * sizeof(double));


    fseek(f, 0, SEEK_SET);
    fscanf(f, "%[^\n]", header);

    int i_vtx  = 0;
    int i_face = 0;
    end = 0;
    face_vtx_idx[0] = 0;
    while(end == 0) {

      int stat = fscanf(f, "%s", header);

      // printf("header = %s \n", header);

      if(strcmp(header, "vertex") == 0) {
        fscanf(f, "%lf %lf %lf \n", &vtx_coord[3*i_vtx  ],
                                          &vtx_coord[3*i_vtx+1],
                                          &vtx_coord[3*i_vtx+2]);
        // printf("%f %f %f \n", vtx_coord[3*i_vtx  ],
        //                                   vtx_coord[3*i_vtx+1],
        //                                   vtx_coord[3*i_vtx+2]);

        face_vtx[face_vtx_idx[i_face]] = i_vtx+1;//(PDM_g_num_t) i_vtx+1;

        face_vtx_coord[3*face_vtx_idx[i_face]  ] = vtx_coord[3*i_vtx  ];
        face_vtx_coord[3*face_vtx_idx[i_face]+1] = vtx_coord[3*i_vtx+1];
        face_vtx_coord[3*face_vtx_idx[i_face]+2] = vtx_coord[3*i_vtx+2];

        face_vtx_idx[i_face]++;
        i_vtx += 1;

        // if(i_vtx == 12) {
        //   return;
        // }
      }
      if(strcmp(header, "facet") == 0) {
        face_vtx_idx[i_face+1] = face_vtx_idx[i_face];
        i_face += 1;
      }

      // if(strcmp(header, "OpenSCAD_Model") == 0) {
      // if(strcmp(header, "endsolid") == 0) {
      //   end = 1;
      // }
      if(stat == EOF) {
        end = 1;
      }

    }


    // const char outfilename[999] = "debug_stl.vtk";
    // PDM_vtk_write_point_cloud(outfilename,
    //                           n_vtx,
    //                           vtx_coord,
    //                           NULL,
    //                           NULL);

    // const char outfilenamep[999] = "debug_stl_face_vtx.vtk";
    // PDM_vtk_write_polydata(outfilenamep,
    //                        n_vtx,
    //                        vtx_coord,
    //                        NULL,
    //                        n_face,
    //                        face_vtx_idx,
    //                        face_vtx,
    //                        NULL,
    //                        NULL);

    for(int i = 0; i < n_face; ++i) {
      face_vtx_n[i] = face_vtx_idx[i+1] - face_vtx_idx[i];
    }

    fclose(f);
  }

  // Re-repart surface among all procs
  PDM_g_num_t* init_distrib_vtx  = PDM_compute_entity_distribution(comm, n_vtx );
  PDM_g_num_t* init_distrib_face = PDM_compute_entity_distribution(comm, n_face);

  PDM_g_num_t _n_g_face = n_face;
  PDM_MPI_Bcast(&_n_g_face, 1, PDM__PDM_MPI_G_NUM, 0, comm);

  // PDM_g_num_t* _distrib_vtx  = PDM_compute_uniform_entity_distribution(comm, n_vtx );
  PDM_g_num_t* _distrib_face = PDM_compute_uniform_entity_distribution(comm, _n_g_face);

  PDM_block_to_block_t *btb_face = PDM_block_to_block_create(init_distrib_face,
                                                             _distrib_face,
                                                             comm);

  int dn_face_end = _distrib_face[i_rank+1] - _distrib_face[i_rank];
  int *tmp_dface_vtx_n = (int * ) malloc( dn_face_end * sizeof(int));
  PDM_g_num_t *tmp_dface_vtx   = NULL;

  PDM_block_to_block_exch(btb_face,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          -1,
                          face_vtx_n,
                (void *)  face_vtx,
                          tmp_dface_vtx_n,
                (void **) &tmp_dface_vtx);

  // PDM_log_trace_array_int   (face_vtx_n    , n_face, "face_vtx_n : ");
  // PDM_log_trace_array_double(face_vtx_coord, 3 * 3 * n_face, "face_vtx_coord : ");

  double *dface_vtx_coord = NULL;
  PDM_block_to_block_exch(btb_face,
                          3 * sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          -1,
                          face_vtx_n,
                (void *)  face_vtx_coord,
                          tmp_dface_vtx_n,
                (void **) &dface_vtx_coord);


  PDM_block_to_block_free(btb_face);

  // *distrib_vtx  = _distrib_vtx;
  *distrib_face = _distrib_face;

  int *tmp_dface_vtx_idx = (int * ) malloc( (dn_face_end + 1) * sizeof(int));
  tmp_dface_vtx_idx[0] = 0;
  for(int i = 0; i < dn_face_end; ++i) {
    tmp_dface_vtx_idx[i+1] = tmp_dface_vtx_idx[i] + tmp_dface_vtx_n[i];
  }

  // PDM_log_trace_array_long(tmp_dface_vtx_idx, dn_face_end+1, "tmp_dface_vtx_idx ::");

  /*
   * Create gnum
   */
  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3, 1, PDM_TRUE, 1.e-6, comm, PDM_OWNERSHIP_USER);

  double *char_length = malloc( tmp_dface_vtx_idx[dn_face_end] * sizeof(double));

  for (int i = 0; i < tmp_dface_vtx_idx[dn_face_end]; ++i) {
    char_length[i] = HUGE_VAL;//1.e-6;
  }

  double tol = 1e-6;
  const double eps_base = 1e-12;
  for(int i = 0; i < dn_face_end; ++i) {
    // double tmp_char_lenght = 1e30;     /* Big value */
    int n_vtx_on_face = tmp_dface_vtx_idx[i+1] - tmp_dface_vtx_idx[i];
    double *fvc = dface_vtx_coord + tmp_dface_vtx_idx[i]*3;
    double *chl = char_length     + tmp_dface_vtx_idx[i];
    for (int j = 0; j < n_vtx_on_face; j++) {
      int ivtx1 = j;
      int ivtx2 = (j+1) % n_vtx_on_face;

      double length2 = 0.;
      for (int k = 0; k < 3; k++) {
        double delta = fvc[3*ivtx1 + k] - fvc[3*ivtx2 + k];
        length2 += delta * delta;
      }

      chl[ivtx1] = PDM_MIN(chl[ivtx1], length2);
      chl[ivtx2] = PDM_MIN(chl[ivtx2], length2);
    }
    // for(int idx_vtx = tmp_dface_vtx_idx[i]; idx_vtx < tmp_dface_vtx_idx[i+1]; ++idx_vtx) {
    //   int idx_vtx1 = dface_vtx[idx_vtx] -1;
    //   double x1 = dface_vtx_coord[3*idx_vtx1  ];
    //   double y1 = dface_vtx_coord[3*idx_vtx1+1];
    //   double z1 = dface_vtx_coord[3*idx_vtx1+2];
    //   int idx_vtx2 = dface_vtx[(idx_vtx +1) % n_vtx_on_face] -1;
    //   double x2 = dface_vtx_coord[3*idx_vtx2  ];
    //   double y2 = dface_vtx_coord[3*idx_vtx2+1];
    //   double z2 = dface_vtx_coord[3*idx_vtx2+2];
    // }
  }

  for (int i = 0; i < tmp_dface_vtx_idx[dn_face_end]; ++i) {
    char_length[i] = PDM_MAX(eps_base, tol*sqrt(char_length[i]));
  }

  PDM_gnum_set_from_coords(gen_gnum,
                           0,
                           tmp_dface_vtx_idx[dn_face_end],
                           dface_vtx_coord,
                           char_length);

  PDM_gnum_compute(gen_gnum);
  free(char_length);

  PDM_g_num_t* _dface_vtx = PDM_gnum_get(gen_gnum, 0);
  PDM_gnum_free(gen_gnum);

  // PDM_log_trace_array_long(_dface_vtx, tmp_dface_vtx_idx[dn_face_end], "_dface_vtx");

  /*
   * Unified vtx
   */
  PDM_part_to_block_t* ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          &_dface_vtx,
                                                          NULL,
                                                          &tmp_dface_vtx_idx[dn_face_end],
                                                          1,
                                                          comm);

  PDM_g_num_t* _distrib_vtx = PDM_part_to_block_distrib_index_get(ptb_vtx);

  // PDM_log_trace_array_long(init_distrib_vtx, n_rank+1, "init_distrib_vtx : ");
  // PDM_log_trace_array_long(_distrib_vtx    , n_rank+1, "_distrib_vtx : ");

  // int    *dvtx_coord_n  = NULL;
  double *_dvtx_coord   = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
                (void**) &dface_vtx_coord,
                         NULL,
                (void**) &_dvtx_coord);

  free(dface_vtx_coord);

  char outfilename[999];
  sprintf(outfilename, "debug_stl_dvtx_coord_%i.vtk", i_rank);
  int dn_vtx = _distrib_vtx[i_rank+1] - _distrib_vtx[i_rank];
  // PDM_vtk_write_point_cloud(outfilename,
  //                           dn_vtx,
  //                           _dvtx_coord,
  //                           NULL,
  //                           NULL);

  PDM_g_num_t *_tmp_distrib_vtx = malloc( (n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    _tmp_distrib_vtx[i] = _distrib_vtx[i];
  }

  PDM_part_to_block_free(ptb_vtx);

  *dvtx_coord   = _dvtx_coord;
  *dface_vtx    = _dface_vtx;

  *distrib_vtx  = _tmp_distrib_vtx;
  *distrib_face = _distrib_face;


  // free(_dvtx_coord);
  // free(_dface_vtx);

  free(tmp_dface_vtx_idx);
  free(tmp_dface_vtx);

  free(tmp_dface_vtx_n  );
  // free(tmp_dvtx_coord  );

  free(init_distrib_vtx );
  free(init_distrib_face);

  if(i_rank == 0) {
    free(face_normal  );
    free(vtx_coord    );
    free(face_vtx_idx );
    free(face_vtx     );
    free(face_vtx_n   );
    free(face_vtx_coord);
  }
}

static
PDM_dmesh_nodal_t*
_create_dmesh_nodal
(
 PDM_MPI_Comm  comm,
 PDM_g_num_t  *distrib_vtx,
 PDM_g_num_t  *distrib_face,
 double       *dvtx_coord,
 PDM_g_num_t  *dface_vtx
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Create dmesh nodal
   */
  PDM_g_num_t gn_vtx  = distrib_vtx [n_rank];
  PDM_g_num_t gn_face = distrib_face[n_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  // dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_section_add(dmn,
                                               PDM_GEOMETRY_KIND_SURFACIC,
                                               PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_SURFACIC,
                                  id_section,
                                  dn_face,
                                  dface_vtx,
                                  PDM_OWNERSHIP_KEEP);


  PDM_dmesh_nodal_generate_distribution(dmn);

  return dmn;
}

static
void
_split_surface_mesh
(
 const PDM_MPI_Comm        comm,
 const PDM_split_dual_t    part_method,
 const int                 n_part,
       PDM_dmesh_nodal_t  *dmn,
       PDM_multipart_t   **_mpart
)
{
  int n_zone = 1;
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);
  PDM_multipart_run_ppart(mpart);

  free(n_part_zones);

  *_mpart = mpart;
}



int main(int argc, char *argv[])
{
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
   *  Set default values
   */

  PDM_g_num_t           n_vtx_seg      = 100;
  PDM_g_num_t           n_g_pts_clouds = 10;
  double                length         = 1.;
  int                   n_part         = 1;
  int                   post           = 1;
  char                 *filename  = NULL;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &filename,
             &n_g_pts_clouds,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  if(filename == NULL) {
    PDM_MPI_Finalize ();
    return 0;
  }
  double radius         = length;//2*length;

  /*
   *  Generate cloud
   */
  if (i_rank == 0) {
    printf("-- Generate volume mesh\n");
    fflush(stdout);
  }
  int          n_pts_clouds;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_random (comm,
                              0, // seed
                              0, // geometric_g_num
                              n_g_pts_clouds,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_pts_clouds,
                              &pts_coord,
                              &pts_g_num);


  /*
   *  Read surface mesh
   */
  if (i_rank == 0) {
    printf("-- Read surface mesh : %s \n", filename);
    fflush(stdout);
  }

  int         *dsurf_face_vtx_idx = NULL;
  PDM_g_num_t *dsurf_face_vtx     = NULL;
  double      *dsurf_vtx_coord    = NULL;
  PDM_g_num_t *surf_distrib_vtx   = NULL;
  PDM_g_num_t *surf_distrib_face  = NULL;

  _read_surface_mesh(comm,
                     filename,
                     &dsurf_face_vtx_idx,
                     &dsurf_face_vtx,
                     &dsurf_vtx_coord,
                     &surf_distrib_vtx,
                     &surf_distrib_face);

  // free(dsurf_face_vtx_idx); // Unused
  PDM_g_num_t  nsurf_g_vtx  = surf_distrib_vtx [n_rank];
  // PDM_g_num_t  nsurf_g_face = surf_distrib_face[n_rank];

  if (i_rank == 0) {
    printf("-- Done \n");
    fflush(stdout);
  }

  log_trace("nsurf_g_vtx = %i \n", nsurf_g_vtx);

  PDM_dmesh_nodal_t* dmn = _create_dmesh_nodal(comm,
                                               surf_distrib_vtx,
                                               surf_distrib_face,
                                               dsurf_vtx_coord,
                                               dsurf_face_vtx);

  PDM_multipart_t       *mpart_surf = NULL;
  _split_surface_mesh(comm,
                      part_method,
                      n_part,
                      dmn,
                      &mpart_surf);

  int          *surf_pn_vtx          = (int          *) malloc(sizeof(int          ) * n_part);
  int          *surf_pn_face         = (int          *) malloc(sizeof(int          ) * n_part);
  int          *surf_pn_edge         = (int          *) malloc(sizeof(int          ) * n_part);
  int         **surf_pface_edge_idx  = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pface_edge      = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pedge_vtx       = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pface_vtx       = (int         **) malloc(sizeof(int         *) * n_part);
  double      **surf_pvtx_coord      = (double      **) malloc(sizeof(double      *) * n_part);
  PDM_g_num_t **surf_pvtx_ln_to_gn   = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surf_pface_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    surf_pn_vtx[i_part] = PDM_multipart_part_vtx_coord_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           &surf_pvtx_coord[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    &surf_pvtx_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           PDM_MESH_ENTITY_FACE,
                                                           &surf_pface_ln_to_gn[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_connectivity_get(mpart_surf,
                                                               0,
                                                               i_part,
                                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                               &surf_pface_edge[i_part],
                                                               &surf_pface_edge_idx[i_part],
                                                               PDM_OWNERSHIP_KEEP);

    int* tmp_pedge_vtx_idx = NULL;
    surf_pn_edge[i_part] = PDM_multipart_part_connectivity_get(mpart_surf,
                                                               0,
                                                               i_part,
                                                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                               &surf_pedge_vtx[i_part],
                                                               &tmp_pedge_vtx_idx,
                                                               PDM_OWNERSHIP_KEEP);
    assert(tmp_pedge_vtx_idx == NULL);


    /*
     * Compute face_vtx
     */
    PDM_compute_face_vtx_from_face_and_edge(surf_pn_face       [i_part],
                                            surf_pface_edge_idx[i_part],
                                            surf_pface_edge    [i_part],
                                            surf_pedge_vtx     [i_part],
                                            &surf_pface_vtx    [i_part]);

    if(0 == 1) {
      PDM_log_trace_array_long(surf_pface_ln_to_gn[i_part], surf_pn_face[i_part]                            , "surf_pface_ln_to_gn :: ");
      PDM_log_trace_array_int (surf_pface_edge_idx[i_part], surf_pn_face[i_part]+1                          , "surf_pface_edge_idx :: ");
      // PDM_log_trace_array_int(surf_pface_vtx[i_part]     , surf_pface_edge_idx[i_part][surf_pn_face[i_part]], "surf_pface_vtx     :: ");
    }
  }


  /*
   * Identity for each point in cloud if it's inside or outside surf
   */
  PDM_dist_cloud_surf_t* ics = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED,
                                                          1,
                                                          comm,
                                                          PDM_OWNERSHIP_KEEP);

  int n_part_cloud  = 1;
  int i_point_cloud = 0;
  int i_part_cloud  = 0;
  PDM_dist_cloud_surf_n_part_cloud_set(ics, i_part_cloud, n_part_cloud);

  PDM_dist_cloud_surf_cloud_set(ics,
                                  i_point_cloud,
                                  i_part_cloud,
                                  n_pts_clouds,
                                  pts_coord,
                                  pts_g_num);

  PDM_dist_cloud_surf_surf_mesh_global_data_set(ics, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_dist_cloud_surf_surf_mesh_part_set(ics,
                                             i_part,
                                             surf_pn_face       [i_part],
                                             surf_pface_edge_idx[i_part],
                                             surf_pface_vtx     [i_part],
                                             surf_pface_ln_to_gn[i_part],
                                             surf_pn_vtx        [i_part],
                                             surf_pvtx_coord    [i_part],
                                             surf_pvtx_ln_to_gn [i_part]);
  }

  PDM_dist_cloud_surf_compute_optim(ics);

  PDM_dist_cloud_surf_dump_times(ics);

  double      *distance = NULL;
  double      *projected = NULL;
  PDM_g_num_t *closest_elt_gnum = NULL;
  PDM_dist_cloud_surf_get(ics,
                          i_point_cloud,
                          i_part_cloud,
                          &distance,
                          &projected,
                          &closest_elt_gnum);


  // if (1 == 1) {
  //   for (int i_part = 0; i_part < n_part_cloud; ++i_part)
  //   {
  //     char filename_out[999];

  //     sprintf(filename_out, "raytracing_distance_%2.2d.vtk", i_rank);
  //     PDM_vtk_write_point_cloud(filename_out,
  //                               n_pts_clouds,
  //                               pts_coord,
  //                               NULL,
  //                               distance);
  //                               // &distance[0],
  //                               // NULL);
  //   }
  // }

  PDM_dist_cloud_surf_free(ics);

  free(pts_coord);
  free(pts_g_num);
  free(surf_pn_edge       );
  free(surf_pn_vtx        );
  free(surf_pn_face       );
  free(surf_pface_edge    );
  free(surf_pface_edge_idx);
  for(int i_part = 0; i_part < n_part; i_part++) {
    free(surf_pface_vtx[i_part]);
  }
  free(surf_pface_vtx     );
  free(surf_pedge_vtx     );
  free(surf_pvtx_coord    );
  free(surf_pvtx_ln_to_gn );
  free(surf_pface_ln_to_gn);




  // free(dsurf_face_vtx);
  // free(dsurf_vtx_coord);
  free(surf_distrib_vtx);
  free(surf_distrib_face);

  PDM_DMesh_nodal_free(dmn);

  PDM_multipart_free(mpart_surf);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
