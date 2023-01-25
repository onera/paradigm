/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dist_cloud_surf_priv.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_surf_mesh.h"
#include "pdm_octree.h"
#include "pdm_para_octree.h"
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"
#include "pdm_hash_tab.h"
#include "pdm_sort.h"
#include "pdm_logging.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_unique.h"
#include "pdm_array.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 8


/*============================================================================
 * Global variable
 *============================================================================*/

static int idebug = 0;

/*=============================================================================
 * Private type
 *============================================================================*/

//--->>>
typedef enum {
  PDM_OCTREE_SERIAL,
  PDM_OCTREE_PARALLEL,
} _octree_type_t;
//<<<---

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Compute distance
 *
 * \param [in]   dist  Pointer to \ref PDM_dist_cloud_surf object
 *
 */

static void
_dist_cloud_surf_compute
(
 PDM_dist_cloud_surf_t *dist
)
{
  const int n_point_cloud      = dist->n_point_cloud;
  PDM_Mesh_nodal_t *mesh_nodal = NULL;//dist->mesh_nodal;
  PDM_surf_mesh_t  *surf_mesh  = dist->_surf_mesh;
  PDM_MPI_Comm comm            = dist->comm;

  int rank;
  PDM_MPI_Comm_rank (comm, &rank);
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  //--->>>
  _octree_type_t octree_type = PDM_OCTREE_PARALLEL;
  char *env_octree_type = getenv ("PDM_OCTREE_TYPE");
  if (env_octree_type != NULL) {
    if (atoi(env_octree_type) == 0) {
      octree_type = PDM_OCTREE_SERIAL;
    }
    else if (atoi(env_octree_type) == 1) {
      octree_type = PDM_OCTREE_PARALLEL;
    }
  }
  if (idebug && rank == 0) printf("octree_type = %d\n", octree_type);
  //<<<---

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  //PDM_timer_hang_on(dist->timer);
  dist->times_elapsed[BEGIN] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[BEGIN]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);



  const double tolerance = 1e-4;
  // const int depth_max = 35;
  const int depth_max = 31;
  const int points_in_leaf_max = 4;

  int n_part_mesh = 0;
  if (mesh_nodal != NULL) {
    n_part_mesh = PDM_Mesh_nodal_n_part_get (mesh_nodal);
  }
  else if (surf_mesh != NULL) {
    n_part_mesh = PDM_surf_mesh_n_part_get (surf_mesh);
  }
  else {
    PDM_error(__FILE__, __LINE__, 0,
        "PDM_dist_cloud_surf error : The surface mesh is not defined. "
        "To do that : \n"
        "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
        "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
        " PDM_dist_cloud_surf_surf_mesh_part_set\n");
  }

  PDM_octree_t      *octree      = NULL;
  PDM_para_octree_t *para_octree = NULL;
  if (octree_type == PDM_OCTREE_SERIAL) {
    octree = PDM_octree_create (n_part_mesh,
        depth_max,
        points_in_leaf_max,
        tolerance,
        comm);
  } else {
    para_octree = PDM_para_octree_create (n_part_mesh,
        depth_max,
        points_in_leaf_max,
        0,
        comm);
  }

  for (int i_part = 0; i_part < n_part_mesh; i_part++) {

    int n_vertices = 0;
    const double *vertices_coords = NULL;
    const PDM_g_num_t *vertices_gnum = NULL;

    if (mesh_nodal != NULL) {
      n_vertices      = PDM_Mesh_nodal_n_vertices_get (mesh_nodal, i_part);
      vertices_coords = PDM_Mesh_nodal_vertices_get   (mesh_nodal, i_part);
      vertices_gnum   = PDM_Mesh_nodal_vertices_g_num_get (mesh_nodal, i_part);
    } else if (surf_mesh != NULL) {
      n_vertices      = PDM_surf_mesh_part_n_vtx_get(surf_mesh, i_part);
      vertices_coords = PDM_surf_mesh_part_vtx_get  (surf_mesh, i_part);
      vertices_gnum   = PDM_surf_mesh_part_vtx_g_num_get (surf_mesh, i_part);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0,
          "PDM_dist_cloud_surf error : The surface mesh is not defined. "
          "To do that : \n"
          "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
          "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
          " PDM_dist_cloud_surf_surf_mesh_part_set\n");
    }

    if (octree_type == PDM_OCTREE_SERIAL) {
      PDM_octree_point_cloud_set (octree, i_part, n_vertices,
          vertices_coords, vertices_gnum);
    } else {
      PDM_para_octree_point_cloud_set (para_octree, i_part, n_vertices,
          vertices_coords, vertices_gnum);
    }
  }

  /*
   * Build octree
   */
  int                *part_n_elt       = malloc (sizeof(int          ) * n_part_mesh);
  const double      **part_elt_extents = malloc (sizeof(double      *) * n_part_mesh);
  const PDM_g_num_t **part_elt_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);

  if (mesh_nodal != NULL) {
    //...
  }
  else if (surf_mesh != NULL) {
    PDM_surf_mesh_compute_faceExtentsMesh (surf_mesh, 1e-8);
    for (int i_part = 0; i_part < n_part_mesh; i_part++) {
      part_n_elt[i_part] = PDM_surf_mesh_part_n_face_get (surf_mesh, i_part);

      part_elt_g_num[i_part] = PDM_surf_mesh_part_face_g_num_get (surf_mesh, i_part);

      part_elt_extents[i_part] = PDM_surf_mesh_part_extents_get (surf_mesh, i_part);

    }
  }

  /* Compute local extents */
  double my_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    for (int i = 0; i < part_n_elt[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        my_extents[j]   = PDM_MIN (my_extents[j],   part_elt_extents[ipart][6*i + j]);
        my_extents[j+3] = PDM_MAX (my_extents[j+3], part_elt_extents[ipart][6*i + 3 + j]);
      }
    }
  }

  /* Compute global extents */
  double global_extents[6];
  PDM_MPI_Allreduce (my_extents,   global_extents,   3,
      PDM_MPI_DOUBLE, PDM_MPI_MIN, dist->comm);
  PDM_MPI_Allreduce (my_extents+3, global_extents+3, 3,
      PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i]   -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }

  if (octree_type == PDM_OCTREE_SERIAL) {
    PDM_octree_build (octree);
  } else {
    PDM_para_octree_build (para_octree, NULL);//global_extents);
    if (0) {
      PDM_para_octree_dump_times (para_octree);
    }
  }


  /***************************************************************************
   *
   * Compute bounding box structure to find candidates closest
   *     to the upper bound distance
   *
   **************************************************************************/

  PDM_timer_hang_on(dist->timer);
  b_t_elapsed = PDM_timer_elapsed(dist->timer);
  b_t_cpu     = PDM_timer_cpu(dist->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create (dist->comm, 3, global_extents);

  PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set (dbbt,
      n_part_mesh,
      part_n_elt,
      part_elt_extents,
      part_elt_g_num);

  if (idebug) {
    printf ("surf_mesh_boxes->n_boxes : %d\n", PDM_box_set_get_size (surf_mesh_boxes));
    for (int i_part = 0; i_part < n_part_mesh; i_part++) {
      printf (" PDM_dbbtree_boxes_set n_elmts %d : %d\n", i_part, part_n_elt[i_part]);
      for (int i = 0; i < part_n_elt[i_part]; i++) {
        printf ("%d : extents %12.5e %12.5e %12.5e / %12.5e %12.5e %12.5e gnum "PDM_FMT_G_NUM"\n",
            i,
            part_elt_extents[i_part][6*i  ],
            part_elt_extents[i_part][6*i+1],
            part_elt_extents[i_part][6*i+2],
            part_elt_extents[i_part][6*i+3],
            part_elt_extents[i_part][6*i+4],
            part_elt_extents[i_part][6*i+5],
            part_elt_g_num[i_part][i]);
      }
    }
  }

  PDM_timer_hang_on(dist->timer);
  e_t_elapsed = PDM_timer_elapsed(dist->timer);
  e_t_cpu     = PDM_timer_cpu(dist->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

  dist->times_elapsed[BBTREE_CREATE] += e_t_elapsed - b_t_elapsed;
  dist->times_cpu[BBTREE_CREATE]     += e_t_cpu - b_t_cpu;
  dist->times_cpu_u[BBTREE_CREATE]   += e_t_cpu_u - b_t_cpu_u;
  dist->times_cpu_s[BBTREE_CREATE]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(dist->timer);

  PDM_timer_hang_on(dist->timer);
  b_t_elapsed = PDM_timer_elapsed(dist->timer);
  b_t_cpu     = PDM_timer_cpu(dist->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);


  /*
   * For each cloud
   */
  for (int i_point_cloud = 0; i_point_cloud < n_point_cloud; i_point_cloud++) {

    _points_cloud_t *pt_cloud = &(dist->points_cloud[i_point_cloud]);
    const int n_part = pt_cloud->n_part;

    /***************************************************************************
     *
     * Compute the upper bound distance. It is the distance from the closest
     * vertex
     *      - Store mesh vetices in a octree
     *      - Compute the closest distance of points to vertices stored in
     *        the octree
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);


    /*
     * Concatenation of the partitions
     */

    int n_pts_rank = 0;

    for (int i_part = 0; i_part < n_part; i_part++) {
      n_pts_rank += pt_cloud->n_points[i_part];
    }

    double      *pts_rank       = malloc (sizeof(double)      * n_pts_rank * 3);
    PDM_g_num_t *pts_g_num_rank = malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    n_pts_rank = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < pt_cloud->n_points[i_part]; i++) {
        for (int k = 0; k < 3; k++) {
          pts_rank[3*n_pts_rank + k] = pt_cloud->coords[i_part][3*i + k];
        }
        pts_g_num_rank[n_pts_rank++] = pt_cloud->gnum[i_part][i];
      }
    }

    /*
     * Look for closest surface mesh vertices
     */
    PDM_g_num_t *closest_vertices_gnum =
      malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    double *closest_vertices_dist2 =  malloc (sizeof(double) * n_pts_rank);
    // log_trace("n_pts_rank:: %d\n", n_pts_rank);

    if (octree_type == PDM_OCTREE_SERIAL) {
      // log_trace("PDM_OCTREE_SERIAL \n");
      PDM_octree_closest_point (octree,
                                n_pts_rank,
                                pts_rank,
                                pts_g_num_rank,
                                closest_vertices_gnum,
                                closest_vertices_dist2);
    } else {
      // log_trace("PDM_OCTREE_PARALLEL \n");
      PDM_para_octree_single_closest_point (para_octree,
                                            n_pts_rank,
                                            pts_rank,
                                            pts_g_num_rank,
                                            closest_vertices_gnum,
                                            closest_vertices_dist2);
    }
    // PDM_log_trace_array_long(closest_vertices_gnum, n_pts_rank, "closest_vertices_gnum::");
    // PDM_log_trace_array_double(closest_vertices_dist2, n_pts_rank, "closest_vertices_dist2::");
    free (closest_vertices_gnum);

    if (i_point_cloud == n_point_cloud -1) { //Octree is not needed anymore
      if (octree_type == PDM_OCTREE_SERIAL) {
        PDM_octree_free (octree);
      } else {
        PDM_para_octree_free (para_octree);
      }
    }


    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[UPPER_BOUND_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[UPPER_BOUND_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[UPPER_BOUND_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[UPPER_BOUND_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);




    /*
     * Find elements closer than closest_vertices_dist2 distance
     */

    int         *part_pts_elt_idx;
    PDM_g_num_t *part_pts_elt_g_num;
    PDM_dbbtree_closest_upper_bound_dist_boxes_get (dbbt,
                                                    n_pts_rank,
                                                    pts_rank,
                                                    pts_g_num_rank,
                                                    closest_vertices_dist2,
                                                    &part_pts_elt_idx,
                                                    &part_pts_elt_g_num);
    if (0) {
      int nmax = 0;
      int imax = 0;
      for (int i = 0; i < n_pts_rank; i++) {
        int n = part_pts_elt_idx[i+1] - part_pts_elt_idx[i];
        if (n > nmax) {
          nmax = n;
          imax = i;
        }
      }

      printf("[%3d] pt %6d ("PDM_FMT_G_NUM") : %5d candidates (dist = %f)\n",
             rank, imax, pts_g_num_rank[imax], nmax, sqrt(closest_vertices_dist2[imax]));
    }
    if (idebug) {
      printf (" PDM_dbbtree_closest_upper_bound_dist_boxes_get n_pts_rank : %d\n", n_pts_rank);
      for (int i = 0; i < n_pts_rank; i++) {
        printf (PDM_FMT_G_NUM" : (%12.5e %12.5e %12.5e) %12.5e\n", pts_g_num_rank[i],
                pts_rank[3*i], pts_rank[3*i+1], pts_rank[3*i+2],
                closest_vertices_dist2[i]);
        printf ("  boxes %d :" , part_pts_elt_idx[i+1] - part_pts_elt_idx[i]);
        for (int j = part_pts_elt_idx[i]; j < part_pts_elt_idx[i+1]; j++) {
          printf (" "PDM_FMT_G_NUM, part_pts_elt_g_num[j]);
        }
        printf ("\n");
      }
    }

    free (closest_vertices_dist2);

    if (i_point_cloud == n_point_cloud -1 ) { //Now useless
      PDM_dbbtree_free (dbbt);
      PDM_box_set_destroy (&surf_mesh_boxes);
    }


    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[CANDIDATE_SELECTION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[CANDIDATE_SELECTION]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[CANDIDATE_SELECTION]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[CANDIDATE_SELECTION]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);


    /*******************************************************************
     *  Adopt SOURCE point-of-view
     *******************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    double *elt_weight  = malloc (sizeof(double) * part_pts_elt_idx[n_pts_rank]);
    int    *part_stride = malloc (sizeof(int)    * part_pts_elt_idx[n_pts_rank]);

    for (int i = 0; i < part_pts_elt_idx[n_pts_rank]; i++) {
      elt_weight[i] = 1.;
      part_stride[i] = 1;
    }

    PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_MERGE,
                                                         1.,
                                                         (PDM_g_num_t **) &part_pts_elt_g_num,
                                                         &elt_weight,
                                                         &(part_pts_elt_idx[n_pts_rank]),
                                                         1,
                                                         comm);
    free (elt_weight);

    /*******************************************************************
     *  Send pts coords & gnum from TARGET parts to SOURCE blocks
     *******************************************************************/
    PDM_g_num_t *part_pts_g_num = malloc (sizeof(PDM_g_num_t) * part_pts_elt_idx[n_pts_rank]);
    double      *part_pts_coord = malloc (sizeof(double) *      part_pts_elt_idx[n_pts_rank] * 3);
    for (int i = 0; i < n_pts_rank; i++) {
      for (int j = part_pts_elt_idx[i]; j < part_pts_elt_idx[i+1]; j++) {
        part_pts_g_num[j] = pts_g_num_rank[i];
        for (int k = 0; k < 3; k++) {
          part_pts_coord[3*j + k] = pts_rank[3*i + k];
        }
      }
    }
    free (pts_rank);

    int *block_elt_pts_n = NULL;
    PDM_g_num_t *block_elt_pts_g_num = NULL;
    PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &part_pts_g_num,
                            &block_elt_pts_n,
                            (void **) &block_elt_pts_g_num);
    free (part_pts_g_num);


    /*for (int i = 0; i < part_pts_elt_idx[n_pts_rank]; i++) {
      part_stride[i] = 3;
      }*/

    int *block_elt_pts_n3 = NULL;
    double *block_elt_pts_coord = NULL;
    PDM_part_to_block_exch (ptb,
                            3*sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &part_pts_coord,
                            &block_elt_pts_n3,
                            (void **) &block_elt_pts_coord);
    free (part_pts_coord);
    free (block_elt_pts_n3);
    free (part_stride);
    free (part_pts_elt_g_num);

    /*******************************************************************
     *  Transfer element coords from parts to blocks
     *******************************************************************/
    PDM_g_num_t *block_elt_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);
    PDM_g_num_t *_block_elt_distrib_idx = block_elt_distrib_idx;

    /* Fix incomplete distribution */
    /*PDM_g_num_t l_max_elt_g_num = 0;
    for (int ipart = 0; ipart < n_part_mesh; ipart++) {
      for (int i = 0; i < part_n_elt[ipart]; i++) {
        l_max_elt_g_num = PDM_MAX (l_max_elt_g_num, part_elt_g_num[ipart][i]);
      }
    }
    PDM_g_num_t g_max_elt_g_num;
    PDM_MPI_Allreduce (&l_max_elt_g_num, &g_max_elt_g_num, 1,
    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);*/
    PDM_g_num_t g_max_elt_g_num = PDM_surf_mesh_n_g_face_get (surf_mesh);

    if (block_elt_distrib_idx[n_rank] < g_max_elt_g_num) {
      _block_elt_distrib_idx = malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
      for (int i = 0; i < n_rank; i++) {
        _block_elt_distrib_idx[i] = block_elt_distrib_idx[i];
      }
      _block_elt_distrib_idx[n_rank] = g_max_elt_g_num;
    }

    PDM_part_to_block_t *ptb_elt = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_MERGE,
                                                              1.,
                                                              (PDM_g_num_t **) part_elt_g_num,
                                                              _block_elt_distrib_idx,
                                                              part_n_elt,
                                                              n_part_mesh,
                                                              comm);

    int **part_elt_vtx_n = malloc (sizeof(int *) * n_part_mesh);
    double **part_elt_vtx_coord = malloc (sizeof(double *) * n_part_mesh);
    for (int ipart = 0; ipart < n_part_mesh; ipart++) {

      const int *part_face_vtx = PDM_surf_mesh_part_face_vtx_get (surf_mesh, ipart);
      const int *part_face_vtx_idx = PDM_surf_mesh_part_face_vtx_idx_get (surf_mesh, ipart);
      const double *part_vtx = PDM_surf_mesh_part_vtx_get (surf_mesh, ipart);

      part_elt_vtx_n[ipart] = malloc (sizeof(int) * part_n_elt[ipart]);
      part_elt_vtx_coord[ipart] =
        malloc (sizeof(double) * part_face_vtx_idx[part_n_elt[ipart]] * 3);

      for (int i = 0; i < part_n_elt[ipart]; i++) {
        int face_vtx_n = part_face_vtx_idx[i+1] - part_face_vtx_idx[i];
        part_elt_vtx_n[ipart][i] = face_vtx_n;//3 * face_vtx_n;
        for (int j = part_face_vtx_idx[i]; j < part_face_vtx_idx[i+1]; j++) {
          int ivtx = part_face_vtx[j] - 1;
          for (int k = 0; k < 3; k++) {
            part_elt_vtx_coord[ipart][3*j + k] = part_vtx[3*ivtx + k];
          }
        }
      }
    }

    int *block_elt_vtx_n = NULL;
    double *block_elt_vtx_coord = NULL;
    PDM_part_to_block_exch (ptb_elt,
                            3*sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            part_elt_vtx_n,
                            (void **) part_elt_vtx_coord,
                            &block_elt_vtx_n,
                            (void **) &block_elt_vtx_coord);
    for (int ipart = 0; ipart < n_part_mesh; ipart++) {
      free (part_elt_vtx_n[ipart]);
      free (part_elt_vtx_coord[ipart]);
    }
    free (part_elt_vtx_n);
    free (part_elt_vtx_coord);


    /*******************************************************************
     *  Compute element-point distances from SOURCE point-of-view
     *******************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);
    PDM_g_num_t *block_elt_g_num = PDM_part_to_block_block_gnum_get (ptb);
    PDM_g_num_t *block_elt_g_num_full = PDM_part_to_block_block_gnum_get (ptb_elt);

    int l_block_elt_pts = 0;
    for (int ielt = 0; ielt < n_elt_block; ielt++) {
      l_block_elt_pts += block_elt_pts_n[ielt];
    }

    double *block_elt_pts_dist2 = malloc (sizeof(double) * l_block_elt_pts);
    double *block_elt_pts_proj  = malloc (sizeof(double) * l_block_elt_pts * 3);

    double *_vtx_coord = block_elt_vtx_coord;
    double *_pts_coord = block_elt_pts_coord;
    double *_dist2     = block_elt_pts_dist2;
    double *_proj      = block_elt_pts_proj;

    int i1 = 0;
    for (int ielt = 0; ielt < n_elt_block; ielt++) {

      while (block_elt_g_num_full[i1] < block_elt_g_num[ielt]) {
        _vtx_coord += 3*block_elt_vtx_n[i1];//block_elt_vtx_n[i1];
        i1++;
      }

      int elt_vtx_n = block_elt_vtx_n[i1];// / 3;

      /* Line */
      if (elt_vtx_n == 2) {
        for (int i = 0; i < block_elt_pts_n[ielt]; i++) {
          double t;
          _dist2[i] = PDM_line_distance (_pts_coord + 3*i,
                                         _vtx_coord,
                                         _vtx_coord + 3,
                                         &t,
                                         _proj + 3*i);
        } // End of loop on points
      }

      /* Triangle */
      else if (elt_vtx_n == 3) {
        for (int i = 0; i < block_elt_pts_n[ielt]; i++) {
          PDM_triangle_status_t status =
          PDM_triangle_evaluate_position(_pts_coord + 3*i,
                                         _vtx_coord,
                                         _proj + 3*i,
                                         _dist2 + i,
                                         NULL);

          if (status == PDM_TRIANGLE_DEGENERATED) {
            for (int j = 0; j < block_elt_pts_n[ielt]; j++) {
              _dist2[j] = HUGE_VAL;
            }
            break;
          }
        } // End of loop on points
      }

      /* Polygon */
      else {
        for (int i = 0; i < block_elt_pts_n[ielt]; i++) {
          PDM_polygon_status_t status =
            PDM_polygon_evaluate_position (_pts_coord + 3*i,
                                           elt_vtx_n,
                                           _vtx_coord,
                                           _proj + 3*i,
                                           _dist2 + i);

          if (status == PDM_POLYGON_DEGENERATED) {
            for (int j = 0; j < block_elt_pts_n[ielt]; j++) {
              _dist2[j] = HUGE_VAL;
            }
            break;
          }
        } // End of loop on points
      }

      _pts_coord += block_elt_pts_n[ielt] * 3;
      _dist2     += block_elt_pts_n[ielt];
      _proj      += block_elt_pts_n[ielt] * 3;

    } // End of loop on elements

    free (block_elt_vtx_n);
    free (block_elt_vtx_coord);
    free (block_elt_pts_coord);
    ptb_elt = PDM_part_to_block_free (ptb_elt);


    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[COMPUTE_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[COMPUTE_ELEM_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[COMPUTE_ELEM_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[COMPUTE_ELEM_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);



    /*******************************************************************
     *  Back to TARGET point-of-view
     *******************************************************************/
    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);

    double *part_pts_weight = malloc (sizeof(double) * n_pts_rank);
    for (int i = 0; i < n_pts_rank; i++) {
      part_pts_weight[i] = (double) (part_pts_elt_idx[i+1] - part_pts_elt_idx[i]);
    }
    free (part_pts_elt_idx);

    /* Build block distribution for points, weighted by number of candidate elements */
    PDM_part_to_block_t *ptb_pts = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                             PDM_PART_TO_BLOCK_POST_MERGE,
                                                             1.,
                                                             &pts_g_num_rank,
                                                             &part_pts_weight,
                                                             &n_pts_rank,
                                                             1,
                                                             comm);
    free (part_pts_weight);


    PDM_g_num_t *block_pts_distrib_idx = PDM_part_to_block_distrib_index_get (ptb_pts);

    /* Exchange results from elements to points */
    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           (PDM_g_num_t **) &block_elt_pts_g_num,
                                                           block_pts_distrib_idx,
                                                           &l_block_elt_pts,
                                                           1,
                                                           comm);

    PDM_g_num_t *part_block_elt_g_num = malloc (sizeof(PDM_g_num_t) * l_block_elt_pts);
    int idx = 0;
    for (int ielt = 0; ielt < n_elt_block; ielt++) {
      for (int i = 0; i < block_elt_pts_n[ielt]; i++) {
        part_block_elt_g_num[idx++] = block_elt_g_num[ielt];
      }
    }
    free (block_elt_pts_n);


    part_stride = malloc (sizeof(int) * l_block_elt_pts);
    for (int i = 0; i < l_block_elt_pts; i++) {
      part_stride[i] = 1;
    }
    int *block_pts_elt_n = NULL;
    PDM_g_num_t *tmp_block_pts_elt_g_num = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &part_block_elt_g_num,
                            &block_pts_elt_n,
                            (void **) &tmp_block_pts_elt_g_num);
    free (block_pts_elt_n);
    free (part_block_elt_g_num);


    double *tmp_block_pts_elt_dist2 = NULL;
    PDM_part_to_block_exch (ptb2,
                            sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &block_elt_pts_dist2,
                            &block_pts_elt_n,
                            (void **) &tmp_block_pts_elt_dist2);
    free (block_elt_pts_dist2);


    /*for (int i = 0; i < l_block_elt_pts; i++) {
      part_stride[i] = 3;
      }*/
    double *tmp_block_pts_elt_proj = NULL;
    int *block_pts_elt_n3 = NULL;
    PDM_part_to_block_exch (ptb2,
                            3*sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                            (void **) &block_elt_pts_proj,
                            &block_pts_elt_n3,
                            (void **) &tmp_block_pts_elt_proj);
    free (block_elt_pts_proj);
    free (part_stride);
    free (block_pts_elt_n3);


    /* Merge block data (keep closest element for each point) */
    int n_pts_block = PDM_part_to_block_n_elt_block_get (ptb_pts);
    PDM_g_num_t *block_pts_elt_g_num = malloc (sizeof(PDM_g_num_t) * n_pts_block);
    double      *block_pts_elt_dist2 = malloc (sizeof(double)      * n_pts_block);
    double      *block_pts_elt_proj  = malloc (sizeof(double)      * n_pts_block * 3);

    idx = 0;
    int n_max = 0;
    for (int i = 0; i < n_pts_block; i++) {
      n_max = PDM_MAX (n_max, block_pts_elt_n[i]);
    }
    int *order = malloc (sizeof(int) * n_max);

    for (int i = 0; i < n_pts_block; i++) {
      block_pts_elt_dist2[i] = HUGE_VAL;

      PDM_g_num_t *_tmp_g_num = tmp_block_pts_elt_g_num + idx;
      double      *_tmp_dist2 = tmp_block_pts_elt_dist2 + idx;
      double      *_tmp_proj  = tmp_block_pts_elt_proj  + idx * 3;

      for (int j = 0; j < block_pts_elt_n[i]; j++) {
        order[j] = j;
      }
      PDM_sort_long (_tmp_g_num, order, block_pts_elt_n[i]);
      int jmin = 0;
      for (int j = 0; j < block_pts_elt_n[i]; j++) {
        if (_tmp_dist2[order[j]] < block_pts_elt_dist2[i]) {
          block_pts_elt_dist2[i] = _tmp_dist2[order[j]];
          jmin = j;
        }
        idx++;
      }

      block_pts_elt_g_num[i] = _tmp_g_num[jmin];
      for (int k = 0; k < 3; k++) {
        block_pts_elt_proj[3*i + k] = _tmp_proj[3*order[jmin] + k];
      }
    }
    free (order);
    free (tmp_block_pts_elt_g_num);
    free (tmp_block_pts_elt_dist2);
    free (tmp_block_pts_elt_proj);
    free (block_pts_elt_n);

    ptb2 = PDM_part_to_block_free (ptb2);
    free (block_elt_pts_g_num);




    /*******************************************************************
     *  Final Block-to-Part transfer
     *******************************************************************/
    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_pts_distrib_idx,
                                                         (const PDM_g_num_t **) pt_cloud->gnum,
                                                         pt_cloud->n_points,
                                                         n_part,
                                                         comm);

    for (int i = 0; i < n_part; i++) {
      int n_pts = pt_cloud->n_points[i];
      pt_cloud->dist[i] = malloc (sizeof(double) * n_pts);
      pt_cloud->proj[i] = malloc (sizeof(double) * n_pts * 3);
      pt_cloud->closest_elt_gnum[i] = malloc (sizeof(PDM_g_num_t) * n_pts);
    }

    int one = 1;
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            &one,
                            block_pts_elt_dist2,
                            NULL,
                            (void **) pt_cloud->dist);
    free (block_pts_elt_dist2);

    int three = 3;
    PDM_block_to_part_exch_in_place (btp,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            &three,
                            block_pts_elt_proj,
                            NULL,
                            (void **) pt_cloud->proj);
    free (block_pts_elt_proj);

    PDM_block_to_part_exch_in_place (btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            &one,
                            block_pts_elt_g_num,
                            NULL,
                            (void **) pt_cloud->closest_elt_gnum);
    free (block_pts_elt_g_num);


    if (_block_elt_distrib_idx != block_elt_distrib_idx) free (_block_elt_distrib_idx);
    btp = PDM_block_to_part_free (btp);
    ptb_pts = PDM_part_to_block_free (ptb_pts);
    ptb = PDM_part_to_block_free (ptb);

    free (pts_g_num_rank);

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[RESULT_TRANSMISSION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[RESULT_TRANSMISSION]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[RESULT_TRANSMISSION]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[RESULT_TRANSMISSION]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

  } // End of loop on point clouds
  free (part_n_elt);
  free(part_elt_g_num);
  free(part_elt_extents);

  PDM_timer_hang_on(dist->timer);
  dist->times_elapsed[END] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[END]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[END]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[END]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

}


/**
 *
 * \brief Compute distance
 *
 * \param [in]   dist  Pointer to \ref PDM_dist_cloud_surf object
 *
 */

static void
_dist_cloud_surf_compute_optim
(
 PDM_dist_cloud_surf_t *dist
)
{
  int dbg_enabled = 1;

  const int n_point_cloud      = dist->n_point_cloud;
  // PDM_Mesh_nodal_t *mesh_nodal = dist->mesh_nodal;
  PDM_part_mesh_nodal_t *mesh_nodal = dist->mesh_nodal;
  PDM_surf_mesh_t  *surf_mesh  = dist->_surf_mesh;
  PDM_MPI_Comm comm            = dist->comm;

  int rank;
  PDM_MPI_Comm_rank (comm, &rank);
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  //--->>>
  _octree_type_t octree_type = PDM_OCTREE_PARALLEL;
  char *env_octree_type = getenv ("PDM_OCTREE_TYPE");
  if (env_octree_type != NULL) {
    if (atoi(env_octree_type) == 0) {
      octree_type = PDM_OCTREE_SERIAL;
    }
    else if (atoi(env_octree_type) == 1) {
      octree_type = PDM_OCTREE_PARALLEL;
    }
  }
  if (idebug && rank == 0) printf("octree_type = %d\n", octree_type);
  //<<<---

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  //PDM_timer_hang_on(dist->timer);
  dist->times_elapsed[BEGIN] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[BEGIN]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  const double tolerance          = 1e-4;
  const int    depth_max          = 31;
  const int    points_in_leaf_max = 4;

  int n_part_mesh = 0;
  if (mesh_nodal != NULL) {
    n_part_mesh = PDM_part_mesh_nodal_n_part_get (mesh_nodal);
  }
  else if (surf_mesh != NULL) {
    n_part_mesh = PDM_surf_mesh_n_part_get (surf_mesh);
  }
  else {
    PDM_error(__FILE__, __LINE__, 0,
        "PDM_dist_cloud_surf error : The surface mesh is not defined. "
        "To do that : \n"
        "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
        "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
        " PDM_dist_cloud_surf_surf_mesh_part_set\n");
  }

  PDM_octree_t      *octree      = NULL;
  PDM_para_octree_t *para_octree = NULL;
  if (octree_type == PDM_OCTREE_SERIAL) {
    octree = PDM_octree_create (n_part_mesh,
        depth_max,
        points_in_leaf_max,
        tolerance,
        comm);
  } else {
    para_octree = PDM_para_octree_create (n_part_mesh,
        depth_max,
        points_in_leaf_max,
        0,
        comm);
  }

  for (int i_part = 0; i_part < n_part_mesh; i_part++) {

    int n_vertices = 0;
    const double *vertices_coords = NULL;
    const PDM_g_num_t *vertices_gnum = NULL;

    if (mesh_nodal != NULL) {
      n_vertices      = PDM_part_mesh_nodal_n_vtx_get    (mesh_nodal, i_part);
      vertices_coords = PDM_part_mesh_nodal_vtx_coord_get(mesh_nodal, i_part);
      vertices_gnum   = PDM_part_mesh_nodal_vtx_g_num_get(mesh_nodal, i_part);
    } else if (surf_mesh != NULL) {
      n_vertices      = PDM_surf_mesh_part_n_vtx_get(surf_mesh, i_part);
      vertices_coords = PDM_surf_mesh_part_vtx_get  (surf_mesh, i_part);
      vertices_gnum   = PDM_surf_mesh_part_vtx_g_num_get (surf_mesh, i_part);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0,
          "PDM_dist_cloud_surf error : The surface mesh is not defined. "
          "To do that : \n"
          "        Call PDM_dist_cloud_surf_nodal_mesh_set or\n"
          "        Call PDM_dist_cloud_surf_surf_mesh_global_data_set +"
          " PDM_dist_cloud_surf_surf_mesh_part_set\n");
    }

    if (octree_type == PDM_OCTREE_SERIAL) {
      PDM_octree_point_cloud_set (octree, i_part, n_vertices,
                                  vertices_coords, vertices_gnum);
    } else {
      PDM_para_octree_point_cloud_set (para_octree, i_part, n_vertices,
                                       vertices_coords, vertices_gnum);
    }
  }

  /*
   * Build octree
   */
  int                *part_n_elt       = malloc (sizeof(int          ) * n_part_mesh);
  // const double      **part_elt_extents = malloc (sizeof(double      *) * n_part_mesh);
  // const PDM_g_num_t **part_elt_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);
  double      **part_elt_extents = malloc (sizeof(double      *) * n_part_mesh);
  PDM_g_num_t **part_elt_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);

  PDM_geometry_kind_t geom_kind;

  if (mesh_nodal != NULL) {
    /* Infer mesh dimension from mesh_nodal (< 3?) */
    PDM_geometry_kind_t l_geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
    for (l_geom_kind = PDM_GEOMETRY_KIND_SURFACIC; l_geom_kind < PDM_GEOMETRY_KIND_MAX; l_geom_kind++) {
      int n_section = PDM_part_mesh_nodal_n_section_get(mesh_nodal,
                                                        l_geom_kind);

      if (n_section > 0) {
        break;
      }
    }

    PDM_MPI_Allreduce(&l_geom_kind, &geom_kind, 1, PDM_MPI_INT, PDM_MPI_MIN, dist->comm);

    if (dbg_enabled) {
      PDM_part_mesh_nodal_dump_vtk(mesh_nodal,
                                   geom_kind,
                                   "dist_cloud_surf_nodal_");
    }

    for (int i_part = 0; i_part < n_part_mesh; i_part++) {
      part_n_elt[i_part] = 0;

      int n_section = PDM_part_mesh_nodal_n_section_get(mesh_nodal, geom_kind);

      int *sections_id = PDM_part_mesh_nodal_sections_id_get(mesh_nodal, geom_kind);

      int max_n_elt = 0;
      for (int isection = 0; isection < n_section; isection++) {
        int id_section = sections_id[isection];

        int n_elt = PDM_part_mesh_nodal_block_n_elt_get(mesh_nodal,
                                                        geom_kind,
                                                        id_section,
                                                        i_part);

        part_n_elt[i_part] += n_elt;

        max_n_elt = PDM_MAX(max_n_elt, n_elt);
      }

      part_elt_g_num  [i_part] = malloc(sizeof(PDM_g_num_t) * part_n_elt[i_part]);
      part_elt_extents[i_part] = malloc(sizeof(double     ) * part_n_elt[i_part] * 6);

      double *_extents = malloc(sizeof(double) * max_n_elt * 6);

      int idx = 0;
      for (int isection = 0; isection < n_section; isection++) {
        int id_section = sections_id[isection];

        int *parent_num = PDM_part_mesh_nodal_block_parent_num_get(mesh_nodal,
                                                                   geom_kind,
                                                                   id_section,
                                                                   i_part);

        PDM_g_num_t *_elt_g_num = PDM_part_mesh_nodal_block_g_num_get(mesh_nodal,
                                                                      geom_kind,
                                                                      id_section,
                                                                      i_part);

        int n_elt = PDM_part_mesh_nodal_block_n_elt_get(mesh_nodal,
                                                        geom_kind,
                                                        id_section,
                                                        i_part);

        PDM_part_mesh_nodal_block_elt_extents_compute(mesh_nodal,
                                                      geom_kind,
                                                      id_section,
                                                      i_part,
                                                      1e-8,
                                                      _extents);

        for (int i = 0; i < n_elt; i++) {
          idx = i;
          if (parent_num) {
            idx = parent_num[i];
          }
          part_elt_g_num[i_part][idx] = _elt_g_num[i];
          memcpy(part_elt_extents[i_part] + 6*idx, _extents + 6*i, sizeof(double)*6);
        }
      }
      free(_extents);
    }
  }
  else if (surf_mesh != NULL) {
    PDM_surf_mesh_compute_faceExtentsMesh (surf_mesh, 1e-8);
    for (int i_part = 0; i_part < n_part_mesh; i_part++) {
      part_n_elt[i_part] = PDM_surf_mesh_part_n_face_get (surf_mesh, i_part);

      part_elt_g_num[i_part] = (PDM_g_num_t *) PDM_surf_mesh_part_face_g_num_get (surf_mesh, i_part);

      part_elt_extents[i_part] = (double *) PDM_surf_mesh_part_extents_get (surf_mesh, i_part);
    }
  }

  if (dbg_enabled) {
    for (int i_part = 0; i_part < n_part_mesh; i_part++) {
      char filename[999];
      sprintf(filename, "dist_cloud_surf_bbox_%d_%d.vtk", i_part, rank);

      PDM_vtk_write_boxes(filename,
                          part_n_elt[i_part],
                          part_elt_extents[i_part],
                          part_elt_g_num[i_part]);
    }
  }

  /* Compute local extents */
  double my_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    for (int i = 0; i < part_n_elt[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        my_extents[j]   = PDM_MIN (my_extents[j],   part_elt_extents[ipart][6*i + j]);
        my_extents[j+3] = PDM_MAX (my_extents[j+3], part_elt_extents[ipart][6*i + 3 + j]);
      }
    }
  }

  /* Compute global extents */
  double global_extents[6];
  PDM_MPI_Allreduce (my_extents,   global_extents,   3,
      PDM_MPI_DOUBLE, PDM_MPI_MIN, dist->comm);
  PDM_MPI_Allreduce (my_extents+3, global_extents+3, 3,
      PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i]   -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }

  if (octree_type == PDM_OCTREE_SERIAL) {
    PDM_octree_build (octree);
  } else {
    PDM_para_octree_build (para_octree, NULL);//global_extents);
    if (0) {
      PDM_para_octree_dump_times (para_octree);
    }
  }


  /***************************************************************************
   *
   * Compute bounding box structure to find candidates closest
   *     to the upper bound distance
   *
   **************************************************************************/

  PDM_timer_hang_on(dist->timer);
  b_t_elapsed = PDM_timer_elapsed(dist->timer);
  b_t_cpu     = PDM_timer_cpu(dist->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create (dist->comm, 3, global_extents);

  PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set (dbbt,
                                                           n_part_mesh,
                                                           part_n_elt,
                                    (const double      **) part_elt_extents,
                                    (const PDM_g_num_t **) part_elt_g_num);

  if (idebug) {
    printf ("surf_mesh_boxes->n_boxes : %d\n", PDM_box_set_get_size (surf_mesh_boxes));
    for (int i_part = 0; i_part < n_part_mesh; i_part++) {
      printf (" PDM_dbbtree_boxes_set n_elmts %d : %d\n", i_part, part_n_elt[i_part]);
      for (int i = 0; i < part_n_elt[i_part]; i++) {
        printf ("%d : extents %12.5e %12.5e %12.5e / %12.5e %12.5e %12.5e gnum "PDM_FMT_G_NUM"\n",
            i,
            part_elt_extents[i_part][6*i  ],
            part_elt_extents[i_part][6*i+1],
            part_elt_extents[i_part][6*i+2],
            part_elt_extents[i_part][6*i+3],
            part_elt_extents[i_part][6*i+4],
            part_elt_extents[i_part][6*i+5],
            part_elt_g_num[i_part][i]);
      }
    }
  }

  PDM_timer_hang_on(dist->timer);
  e_t_elapsed = PDM_timer_elapsed(dist->timer);
  e_t_cpu     = PDM_timer_cpu(dist->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

  dist->times_elapsed[BBTREE_CREATE] += e_t_elapsed - b_t_elapsed;
  dist->times_cpu[BBTREE_CREATE]     += e_t_cpu - b_t_cpu;
  dist->times_cpu_u[BBTREE_CREATE]   += e_t_cpu_u - b_t_cpu_u;
  dist->times_cpu_s[BBTREE_CREATE]   += e_t_cpu_s - b_t_cpu_s;

  PDM_timer_resume(dist->timer);

  PDM_timer_hang_on(dist->timer);
  b_t_elapsed = PDM_timer_elapsed(dist->timer);
  b_t_cpu     = PDM_timer_cpu(dist->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);


  /*
   * For each cloud
   */
  for (int i_point_cloud = 0; i_point_cloud < n_point_cloud; i_point_cloud++) {

    _points_cloud_t *pt_cloud = &(dist->points_cloud[i_point_cloud]);
    const int n_part = pt_cloud->n_part;

    /***************************************************************************
     *
     * Compute the upper bound distance. It is the distance from the closest
     * vertex
     *      - Store mesh vetices in a octree
     *      - Compute the closest distance of points to vertices stored in
     *        the octree
     *
     **************************************************************************/

    PDM_timer_hang_on(dist->timer);
    b_t_elapsed = PDM_timer_elapsed(dist->timer);
    b_t_cpu     = PDM_timer_cpu(dist->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
    PDM_timer_resume(dist->timer);



    /*
     * Concatenation of the partitions
     */

    int n_pts_rank = 0;

    for (int i_part = 0; i_part < n_part; i_part++) {
      n_pts_rank += pt_cloud->n_points[i_part];
    }

    double      *pts_rank       = malloc (sizeof(double)      * n_pts_rank * 3);
    PDM_g_num_t *pts_g_num_rank = malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    n_pts_rank = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < pt_cloud->n_points[i_part]; i++) {
        for (int k = 0; k < 3; k++) {
          pts_rank[3*n_pts_rank + k] = pt_cloud->coords[i_part][3*i + k];
        }
        pts_g_num_rank[n_pts_rank++] = pt_cloud->gnum[i_part][i];
      }
    }

    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "dist_cloud_surf_pts_%d.vtk", rank);

      PDM_vtk_write_point_cloud(filename,
                                n_pts_rank,
                                pts_rank,
                                pts_g_num_rank,
                                NULL);
    }

    /*
     * Look for closest surface mesh vertices
     */
    PDM_g_num_t *closest_vertices_gnum =
      malloc (sizeof(PDM_g_num_t) * n_pts_rank);

    double *closest_vertices_dist2 =  malloc (sizeof(double) * n_pts_rank);
    // log_trace("n_pts_rank:: %d\n", n_pts_rank);

    if (octree_type == PDM_OCTREE_SERIAL) {
      // log_trace("PDM_OCTREE_SERIAL \n");
      PDM_octree_closest_point (octree,
                                n_pts_rank,
                                pts_rank,
                                pts_g_num_rank,
                                closest_vertices_gnum,
                                closest_vertices_dist2);
    } else {
      // log_trace("PDM_OCTREE_PARALLEL \n");
      PDM_para_octree_single_closest_point (para_octree,
                                            n_pts_rank,
                                            pts_rank,
                                            pts_g_num_rank,
                                            closest_vertices_gnum,
                                            closest_vertices_dist2);
    }
    // PDM_log_trace_array_long(closest_vertices_gnum, n_pts_rank, "closest_vertices_gnum::");
    // PDM_log_trace_array_double(closest_vertices_dist2, n_pts_rank, "closest_vertices_dist2::");
    free (closest_vertices_gnum);

    if (i_point_cloud == n_point_cloud -1) { //Octree is not needed anymore
      if (octree_type == PDM_OCTREE_SERIAL) {
        PDM_octree_free (octree);
      } else {
        PDM_para_octree_free (para_octree);
      }
    }


    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[UPPER_BOUND_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu[UPPER_BOUND_DIST]     += e_t_cpu - b_t_cpu;
    dist->times_cpu_u[UPPER_BOUND_DIST]   += e_t_cpu_u - b_t_cpu_u;
    dist->times_cpu_s[UPPER_BOUND_DIST]   += e_t_cpu_s - b_t_cpu_s;

    PDM_timer_resume(dist->timer);

    /*
     * Find elements closer than closest_vertices_dist2 distance
     */
    int          n_extract_boxes = 0;
    PDM_g_num_t *box_gnum          = NULL;
    int         *box_init_location = NULL;
    int         *dbox_pts_idx      = NULL;
    PDM_g_num_t *dbox_pts_g_num    = NULL;
    double      *dbox_pts_coord    = NULL;

    // printf("PDM_dbbtree_closest_upper_bound_dist_boxes_pts_shared_get ... \n");
    PDM_dbbtree_closest_upper_bound_dist_boxes_pts_shared_get (dbbt,
                                                               n_pts_rank,
                                                               pts_rank,
                                                               pts_g_num_rank,
                                                               closest_vertices_dist2,
                                                               &n_extract_boxes,
                                                               &box_gnum,
                                                               &box_init_location,
                                                               &dbox_pts_idx,
                                                               &dbox_pts_g_num,
                                                               &dbox_pts_coord);
    free(pts_rank);
    free(pts_g_num_rank);
    // printf("PDM_dbbtree_closest_upper_bound_dist_boxes_pts_shared_get END \n");

    if (0 == 1) {
      PDM_log_trace_array_long(box_gnum,
                               n_extract_boxes,
                               "box_gnum :");
      PDM_log_trace_connectivity_long(dbox_pts_idx,
                                      dbox_pts_g_num,
                                      n_extract_boxes,
                                      "dbox_pts_g_num : ");
    }

    free (closest_vertices_dist2);

    if (i_point_cloud == n_point_cloud -1 ) { //Now useless
      PDM_dbbtree_free (dbbt);
      PDM_box_set_destroy (&surf_mesh_boxes);
      // free(part_n_elt);
      // if (mesh_nodal != NULL) {
      //   for (int i = 0; i < n_part_mesh; i++) {
      //     free(part_elt_g_num  [i]);
      //     free(part_elt_extents[i]);
      //   }
      // }
      // free(part_elt_g_num);
      // free(part_elt_extents);
    }

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[CANDIDATE_SELECTION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu    [CANDIDATE_SELECTION] += e_t_cpu     - b_t_cpu;
    dist->times_cpu_u  [CANDIDATE_SELECTION] += e_t_cpu_u   - b_t_cpu_u;
    dist->times_cpu_s  [CANDIDATE_SELECTION] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(dist->timer);

    /*
     * Setup almost partition for elementary operation
     *    - Make pts as partition
     *    - Extract part on boxes
     */
    PDM_g_num_t *pts_ln_to_gn = dbox_pts_g_num;
    int* unique_order = malloc(dbox_pts_idx[n_extract_boxes] * sizeof(int));
    int n_extract_pts = PDM_inplace_unique_long2(pts_ln_to_gn,
                                                 unique_order,
                                                 0,
                                                 dbox_pts_idx[n_extract_boxes]-1);
    pts_ln_to_gn = realloc(pts_ln_to_gn, n_extract_pts * sizeof(PDM_g_num_t));

    // log_trace("Reduce to pts connectivity = %i / %i \n", dbox_pts_idx[n_extract_boxes], n_extract_pts);

    /*
     * Deduce local numbering and coords
     */
    int    *box_pts    = malloc(    dbox_pts_idx[n_extract_boxes] * sizeof(int   ));
    double *pts_coords = malloc(3 * n_extract_pts                 * sizeof(double));
    for(int idx = 0; idx < dbox_pts_idx[n_extract_boxes] ; ++idx) {
      int l_elmt = unique_order[idx];
      box_pts[idx] = (l_elmt + 1);
      pts_coords[3*l_elmt  ] = dbox_pts_coord[3*idx  ];
      pts_coords[3*l_elmt+1] = dbox_pts_coord[3*idx+1];
      pts_coords[3*l_elmt+2] = dbox_pts_coord[3*idx+2];
    }
    free(unique_order);
    free(dbox_pts_coord);


    PDM_extract_part_t *extrp = PDM_extract_part_create(2, // dim 1?
                                                        n_part_mesh,
                                                        1,                                 // n_part_out
                                                        PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                        PDM_SPLIT_DUAL_WITH_PTSCOTCH,      // unused in this case
                                                        PDM_FALSE,                         // compute_child_gnum
                                                        PDM_OWNERSHIP_KEEP,
                                                        dist->comm);
    if (mesh_nodal != NULL) {
      PDM_part_mesh_nodal_elmts_t *pmne = NULL;
      if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
        pmne = mesh_nodal->ridge;
      }
      else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
        pmne = mesh_nodal->surfacic;
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Invalid geom_kind %d\n", (int) geom_kind);
      }

      PDM_extract_part_part_nodal_set(extrp, pmne);

      for (int i_part = 0; i_part < n_part_mesh; i_part++) {

        int n_vtx = PDM_part_mesh_nodal_n_vtx_get(mesh_nodal, i_part);
        const double      *vtx_coord    = PDM_part_mesh_nodal_vtx_coord_get(mesh_nodal, i_part);
        const PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(mesh_nodal, i_part);

        int n_cell = 0;
        int n_face = 0;
        int n_edge = 0;
        PDM_g_num_t *cell_ln_to_gn = NULL;
        PDM_g_num_t *face_ln_to_gn = NULL;
        PDM_g_num_t *edge_ln_to_gn = NULL;
        switch (geom_kind) {
          case PDM_GEOMETRY_KIND_VOLUMIC: {
            n_cell = part_n_elt[i_part];
            cell_ln_to_gn = part_elt_g_num[i_part];
            break;
          }
          case PDM_GEOMETRY_KIND_SURFACIC: {
            n_face = part_n_elt[i_part];
            face_ln_to_gn = part_elt_g_num[i_part];
            break;
          }
          case PDM_GEOMETRY_KIND_RIDGE: {
            n_edge = part_n_elt[i_part];
            edge_ln_to_gn = part_elt_g_num[i_part];
            break;
          }
          default:
            PDM_error(__FILE__, __LINE__, 0, "incorrect geom_kind %d\n", (int) geom_kind);
        }

        PDM_extract_part_part_set(extrp,
                                  i_part,
                                  n_cell,
                                  n_face,
                                  n_edge,
                                  n_vtx,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  cell_ln_to_gn,
                                  face_ln_to_gn,
                                  edge_ln_to_gn,
                  (PDM_g_num_t *) vtx_ln_to_gn,
                  (double      *) vtx_coord);

      }
    }
    else {
      /* Set vtx_coord */
      for (int i_part = 0; i_part < n_part_mesh; i_part++) {


        int n_face = PDM_surf_mesh_part_n_face_get(surf_mesh, i_part);
        int n_vtx  = PDM_surf_mesh_part_n_vtx_get (surf_mesh, i_part);

        const PDM_g_num_t* face_ln_to_gn = PDM_surf_mesh_part_face_g_num_get(surf_mesh, i_part);
        const PDM_g_num_t* vtx_ln_to_gn  = PDM_surf_mesh_part_vtx_g_num_get (surf_mesh, i_part);

        const int *part_face_vtx = PDM_surf_mesh_part_face_vtx_get (surf_mesh, i_part);
        const int *part_face_vtx_idx = PDM_surf_mesh_part_face_vtx_idx_get (surf_mesh, i_part);
        const double *part_vtx_coords = PDM_surf_mesh_part_vtx_get (surf_mesh, i_part);

        // PDM_log_trace_array_long(vtx_ln_to_gn, n_vtx, "vtx_ln_to_gn :: ");

        PDM_extract_part_part_set(extrp,
                                  i_part,
                                  0,
                                  n_face,
                                  0,
                                  n_vtx,
                                NULL, // pcell_face_idx[i_part],
                                NULL, // pcell_face[i_part],
                                NULL, // pface_edge_idx[i_part],
                                NULL, // pface_edge[i_part],
                                NULL, // pedge_vtx[i_part],
                  (int *)       part_face_vtx_idx, // pface_vtx_idx[i_part],
                  (int *)       part_face_vtx, // pface_vtx[i_part],
                  NULL,
                  (PDM_g_num_t *) face_ln_to_gn,
                  NULL,
                  (PDM_g_num_t *) vtx_ln_to_gn,
                  (double *) part_vtx_coords);

      }
    }

    PDM_extract_part_target_set(extrp,
                                0,
                                n_extract_boxes,
                                box_gnum,
                                box_init_location);

    if(0 == 1) {
      PDM_log_trace_array_long(box_gnum, n_extract_boxes, "box_gnum :: ");
      PDM_log_trace_array_int(box_init_location, 3 * n_extract_boxes, "box_init_location :: ");
    }

    PDM_extract_part_compute(extrp);
    free(box_init_location);

    if (i_point_cloud == n_point_cloud -1 ) { //Now useless
      free(part_n_elt);
      if (mesh_nodal != NULL) {
        for (int i = 0; i < n_part_mesh; i++) {
          free(part_elt_g_num  [i]);
          free(part_elt_extents[i]);
        }
      }
      free(part_elt_g_num);
      free(part_elt_extents);
    }

    if (mesh_nodal != NULL) {

      PDM_part_mesh_nodal_elmts_t *extract_pmne = NULL;
      PDM_extract_part_part_mesh_nodal_get(extrp,
                                           &extract_pmne,
                                           PDM_OWNERSHIP_USER);

      int mesh_dim = 0;
      if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
        mesh_dim = 1;
      }
      else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
        mesh_dim = 2;
      }

      PDM_part_mesh_nodal_t *extract_pmn = PDM_part_mesh_nodal_create(mesh_dim,
                                                                      n_part_mesh,
                                                                      dist->comm);

      for (int i_part = 0; i_part < n_part_mesh; i_part++) {
        double *_vtx_coord = NULL;
        PDM_extract_part_vtx_coord_get(extrp,
                                       i_part,
                                       &_vtx_coord,
                                       PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *_vtx_ln_to_gn = NULL;
        int _n_vtx = PDM_extract_part_parent_ln_to_gn_get(extrp,
                                                          i_part,
                                                          PDM_MESH_ENTITY_VERTEX,
                                                          &_vtx_ln_to_gn,
                                                          PDM_OWNERSHIP_KEEP);

        PDM_part_mesh_nodal_coord_set(extract_pmn,
                                      i_part,
                                      _n_vtx,
                                      _vtx_coord,
                                      _vtx_ln_to_gn,
                                      PDM_OWNERSHIP_USER);
      }
      PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(extract_pmn, extract_pmne, PDM_OWNERSHIP_KEEP);

      if (dbg_enabled) {
        PDM_part_mesh_nodal_dump_vtk(extract_pmn,
                                     geom_kind,
                                     "dist_cloud_surf_extrp_nodal_");
      }

      PDM_MPI_Barrier(dist->comm);
      log_trace("OK jusqu'ici :D\n");

      abort();
    }
    else {
      //??
    }

    int pextract_n_vtx = PDM_extract_part_n_entity_get(extrp,
                                                       0,
                                                       PDM_MESH_ENTITY_VERTEX);

    int *pextract_face_vtx     = NULL;
    int *pextract_face_vtx_idx = NULL;
    int pn_extract_face = PDM_extract_part_connectivity_get(extrp,
                                                            0,
                                                            PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                            &pextract_face_vtx,
                                                            &pextract_face_vtx_idx,
                                                            PDM_OWNERSHIP_KEEP);
    PDM_g_num_t *pextract_face_ln_to_gn = NULL;
    PDM_extract_part_ln_to_gn_get(extrp,
                                  0,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

    double *pextract_vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp,
                                   0,
                                   &pextract_vtx_coord,
                                   PDM_OWNERSHIP_KEEP);

    if(dbg_enabled) {
      char filename[999];
      int i_rank;
      PDM_MPI_Comm_rank (comm, &i_rank);
      sprintf(filename, "out_extrac_part_%i.vtk", i_rank);
      PDM_vtk_write_polydata(filename,
                             pextract_n_vtx,
                             pextract_vtx_coord,
                             NULL,
                             pn_extract_face,
                             pextract_face_vtx_idx,
                             pextract_face_vtx,
                             box_gnum,// pextract_face_ln_to_gn,
                             NULL);

      sprintf(filename, "out_extrac_pts_%i.vtk", i_rank);
      PDM_vtk_write_point_cloud(filename,
                                n_extract_pts,
                                pts_coords,
                                pts_ln_to_gn,
                                NULL);


    }

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[LOAD_BALANCING_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu    [LOAD_BALANCING_ELEM_DIST] += e_t_cpu     - b_t_cpu;
    dist->times_cpu_u  [LOAD_BALANCING_ELEM_DIST] += e_t_cpu_u   - b_t_cpu_u;
    dist->times_cpu_s  [LOAD_BALANCING_ELEM_DIST] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(dist->timer);
    // if (1) {
    //   PDM_log_trace_connectivity_int(dbox_pts_idx,
    //                                  box_pts,
    //                                  n_extract_boxes,
    //                                  "box_pts : ");
    // }


    double      *pts_dist2        = malloc(    n_extract_pts * sizeof(double));
    double      *pts_proj         = malloc(3*  n_extract_pts * sizeof(double));
    PDM_g_num_t *pts_closest_face_g_num = malloc(    n_extract_pts * sizeof(PDM_g_num_t));
    // int         *pts_closest_face_init_loc = malloc(3*  n_extract_pts * sizeof(int));

    for(int i_pts = 0; i_pts < n_extract_pts; ++i_pts) {
      pts_dist2[i_pts] = HUGE_VAL;
      pts_closest_face_g_num[i_pts] = -1;
    }

    int max_elmt_vtx = 0;
    for(int i_elmt = 0; i_elmt < n_extract_boxes; ++i_elmt) {
      max_elmt_vtx = PDM_MAX(max_elmt_vtx, pextract_face_vtx_idx[i_elmt+1] - pextract_face_vtx_idx[i_elmt]);
    }
    double *lvtx_coords = malloc(3 * max_elmt_vtx * sizeof(double));

    for(int i_elmt = 0; i_elmt < n_extract_boxes; ++i_elmt) {

      int n_elmt_vtx = pextract_face_vtx_idx[i_elmt+1] - pextract_face_vtx_idx[i_elmt];

      int idx_write = 0;
      for(int idx_vtx = pextract_face_vtx_idx[i_elmt]; idx_vtx < pextract_face_vtx_idx[i_elmt+1]; ++idx_vtx) {
        int i_vtx = pextract_face_vtx[idx_vtx]-1;
        lvtx_coords[3*idx_write  ] = pextract_vtx_coord[3*i_vtx  ];
        lvtx_coords[3*idx_write+1] = pextract_vtx_coord[3*i_vtx+1];
        lvtx_coords[3*idx_write+2] = pextract_vtx_coord[3*i_vtx+2];
        idx_write++;
      }

      /* Line */
      if (n_elmt_vtx == 2) {
        for(int idx_pts = dbox_pts_idx[i_elmt]; idx_pts < dbox_pts_idx[i_elmt+1]; ++idx_pts) {
          int i_pts = box_pts[idx_pts]-1;
          double t;
          double lproj[3];
          double ldist = PDM_line_distance(&pts_coords[3*i_pts],
                                           lvtx_coords,
                                           lvtx_coords + 3,
                                           &t,
                                           lproj);
          if (ldist < pts_dist2[i_pts]) {
            pts_dist2[i_pts]     = ldist;
            pts_proj [3*i_pts  ] = lproj[0];
            pts_proj [3*i_pts+1] = lproj[1];
            pts_proj [3*i_pts+2] = lproj[2];
            pts_closest_face_g_num[i_pts] = box_gnum[i_elmt];
          }
        } // End of loop on points
      }

      /* Triangle */
      else if (n_elmt_vtx == 3) {
        for(int idx_pts = dbox_pts_idx[i_elmt]; idx_pts < dbox_pts_idx[i_elmt+1]; ++idx_pts) {
          int i_pts = box_pts[idx_pts]-1;
          double lproj[3];
          double ldist;
          PDM_triangle_status_t status =
          PDM_triangle_evaluate_position(&pts_coords[3*i_pts],
                                         lvtx_coords,
                                         lproj,
                                         &ldist,
                                         NULL);

          if (status == PDM_TRIANGLE_DEGENERATED) {
            continue;
          }

          if (ldist < pts_dist2[i_pts]) {
            pts_dist2[i_pts]     = ldist;
            pts_proj [3*i_pts  ] = lproj[0];
            pts_proj [3*i_pts+1] = lproj[1];
            pts_proj [3*i_pts+2] = lproj[2];
            pts_closest_face_g_num[i_pts] = box_gnum[i_elmt];
          }
        } // End of loop on points
      }

      /* Polygon */
      else {
        for(int idx_pts = dbox_pts_idx[i_elmt]; idx_pts < dbox_pts_idx[i_elmt+1]; ++idx_pts) {
          int i_pts = box_pts[idx_pts]-1;
          double lproj[3];
          double ldist;
          PDM_polygon_status_t status =
          PDM_polygon_evaluate_position(&pts_coords[3*i_pts],
                                        n_elmt_vtx,
                                        lvtx_coords,
                                        lproj,
                                        &ldist);

          // if (pts_ln_to_gn[i_pts] == 429) {
          //   log_trace("pt "PDM_FMT_G_NUM" (%f %f %f) at dist2 %e from face "PDM_FMT_G_NUM", proj = %f %f %f\n",
          //             pts_ln_to_gn[i_pts],
          //             pts_coords[3*i_pts], pts_coords[3*i_pts+1], pts_coords[3*i_pts+2],
          //             ldist, box_gnum[i_elmt], lproj[0], lproj[1], lproj[2]);
          // }

          if (status == PDM_POLYGON_DEGENERATED) {
            continue;
          }

          if (ldist < pts_dist2[i_pts]) {
            pts_dist2[i_pts]     = ldist;
            pts_proj [3*i_pts  ] = lproj[0];
            pts_proj [3*i_pts+1] = lproj[1];
            pts_proj [3*i_pts+2] = lproj[2];
            pts_closest_face_g_num[i_pts] = box_gnum[i_elmt];
          }
        } // End of loop on points
      }
    }
    free(lvtx_coords);
    free(box_pts);
    free(pts_coords);


    PDM_extract_part_free(extrp);

    free(box_gnum         );
    free(dbox_pts_idx     );

    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[COMPUTE_ELEM_DIST] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu    [COMPUTE_ELEM_DIST] += e_t_cpu     - b_t_cpu;
    dist->times_cpu_u  [COMPUTE_ELEM_DIST] += e_t_cpu_u   - b_t_cpu_u;
    dist->times_cpu_s  [COMPUTE_ELEM_DIST] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(dist->timer);

    /*
     *  Pass in block frame for points to select closest element
     */
    int    *part_stride = PDM_array_const_int(n_extract_pts, 1);
    double *part_weight = malloc(sizeof(double) * n_extract_pts);
    for (int i = 0; i < n_extract_pts; i++) {
      part_weight[i] = 1.;
    }

    // if (1) {
    //   for (int i = 0; i < n_extract_pts; i++) {
    //     log_trace("point %d ("PDM_FMT_G_NUM") -> face "PDM_FMT_G_NUM" at dist2 = %f, proj = %f %f %f\n",
    //               i, pts_ln_to_gn[i], pts_closest_face_g_num[i],
    //               pts_dist2[i],
    //               pts_proj[3*i], pts_proj[3*i+1], pts_proj[3*i+2]);
    //   }
    // }

    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &pts_ln_to_gn,
                                                        &part_weight,
                                                        &n_extract_pts,
                                                        1,
                                                        dist->comm);
    free(part_weight);
    free(pts_ln_to_gn);

    int    *block_pts_elt_n     = NULL;
    double *block_pts_elt_dist2 = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 (void **) &pts_dist2,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_dist2);
    free(pts_dist2);
    free(block_pts_elt_n);

    PDM_g_num_t *block_pts_elt_g_num = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 (void **) &pts_closest_face_g_num,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_g_num);
    free(pts_closest_face_g_num);
    free(block_pts_elt_n);

    double *block_pts_elt_proj = NULL;
    PDM_part_to_block_exch(ptb,
                           3*sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 (void **) &pts_proj,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_proj);
    free(pts_proj);
    free(part_stride);


    /* Pick closest elt for each point in current block */
    int block_n_pts = PDM_part_to_block_n_elt_block_get(ptb);
    PDM_g_num_t *block_pts_g_num = PDM_part_to_block_block_gnum_get(ptb);

    int idx = 0;
    for (int i = 0; i < block_n_pts; i++) {
      double min_dist2 = HUGE_VAL;
      for (int j = 0; j < block_pts_elt_n[i]; j++) {
        if (block_pts_elt_dist2[idx] < min_dist2) {
          min_dist2              = block_pts_elt_dist2[idx];
          block_pts_elt_g_num[i] = block_pts_elt_g_num[idx];
          // memcpy(block_pts_elt_proj + 3*i, block_pts_elt_proj + 3*idx, sizeof(double) * 3);
          for (int k = 0; k < 3; k++) {
            block_pts_elt_proj[3*i+k] = block_pts_elt_proj[3*idx+k];
          }
        }
        idx++;
      }
      block_pts_elt_dist2[i] = min_dist2;
    }
    free(block_pts_elt_n);
    block_pts_elt_dist2 = realloc(block_pts_elt_dist2, sizeof(double     ) * block_n_pts);
    block_pts_elt_g_num = realloc(block_pts_elt_g_num, sizeof(PDM_g_num_t) * block_n_pts);
    block_pts_elt_proj  = realloc(block_pts_elt_proj,  sizeof(double     ) * block_n_pts * 3);

    /* Send back */
    /* Or send directly to user frame ?? */
    PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_pts_g_num,
                                                                          block_n_pts,
                                                   (const PDM_g_num_t **) pt_cloud->gnum,
                                                   (const int          *) pt_cloud->n_points,
                                                                          pt_cloud->n_part,
                                                                          dist->comm);

    PDM_part_to_block_free(ptb);
    for (int i = 0; i < n_part; i++) {
      int n_pts = pt_cloud->n_points[i];
      pt_cloud->dist[i] = malloc (sizeof(double) * n_pts);
      pt_cloud->proj[i] = malloc (sizeof(double) * n_pts * 3);
      pt_cloud->closest_elt_gnum[i] = malloc (sizeof(PDM_g_num_t) * n_pts);
    }

    int one = 1;
    PDM_block_to_part_exch_in_place(btp,
                                    sizeof(double),
                                    PDM_STRIDE_CST_INTERLACED,
                                    &one,
                         (void   *) block_pts_elt_dist2,
                                    NULL,
                          (void **) pt_cloud->dist);
    free(block_pts_elt_dist2);

    PDM_block_to_part_exch_in_place(btp,
                                    sizeof(PDM_g_num_t),
                                    PDM_STRIDE_CST_INTERLACED,
                                    &one,
                         (void   *) block_pts_elt_g_num,
                                    NULL,
                         (void **)  pt_cloud->closest_elt_gnum);
    free(block_pts_elt_g_num);

    PDM_block_to_part_exch_in_place(btp,
                                    3 * sizeof(double),
                                    PDM_STRIDE_CST_INTERLACED,
                                    &one,
                         (void   *) block_pts_elt_proj,
                                    NULL,
                         (void **)  pt_cloud->proj);
    free(block_pts_elt_proj);

    PDM_block_to_part_free(btp);


    PDM_timer_hang_on(dist->timer);
    e_t_elapsed = PDM_timer_elapsed(dist->timer);
    e_t_cpu     = PDM_timer_cpu(dist->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

    dist->times_elapsed[RESULT_TRANSMISSION] += e_t_elapsed - b_t_elapsed;
    dist->times_cpu    [RESULT_TRANSMISSION] += e_t_cpu     - b_t_cpu;
    dist->times_cpu_u  [RESULT_TRANSMISSION] += e_t_cpu_u   - b_t_cpu_u;
    dist->times_cpu_s  [RESULT_TRANSMISSION] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(dist->timer);

  }

  PDM_timer_hang_on(dist->timer);
  dist->times_elapsed[END] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[END]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[END]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[END]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

}



/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Pointer to \ref PDM_dist_cloud_surf object
 */

PDM_dist_cloud_surf_t*
PDM_dist_cloud_surf_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
)
{
  PDM_dist_cloud_surf_t *dist = (PDM_dist_cloud_surf_t *) malloc(sizeof(PDM_dist_cloud_surf_t));

  dist->comm              = comm;
  dist->owner             = owner;
  dist->results_is_getted = PDM_FALSE;

  dist->mesh_nature   = mesh_nature;
  dist->mesh_nodal    = NULL;
  dist->surf_mesh     = NULL;
  dist->_surf_mesh    = NULL;
  dist->n_point_cloud = n_point_cloud;
  dist->points_cloud  = (_points_cloud_t*) malloc (sizeof(_points_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    dist->points_cloud[i].n_part           = -1;
    dist->points_cloud[i].n_points         = NULL;
    dist->points_cloud[i].coords           = NULL;
    dist->points_cloud[i].gnum             = NULL;
    dist->points_cloud[i].dist             = NULL;
    dist->points_cloud[i].proj             = NULL;
    dist->points_cloud[i].closest_elt_gnum = NULL;
  }

  dist->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    dist->times_elapsed[i] = 0.;
    dist->times_cpu    [i] = 0.;
    dist->times_cpu_u  [i] = 0.;
    dist->times_cpu_s  [i] = 0.;
  }

  return dist;
}



/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_dist_cloud_surf_n_part_cloud_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    n_part
)
{

  dist->points_cloud[i_point_cloud].n_part = n_part;
  dist->points_cloud[i_point_cloud].n_points =
    realloc(dist->points_cloud[i_point_cloud].n_points, n_part * sizeof(int));
  dist->points_cloud[i_point_cloud].coords =
    realloc(dist->points_cloud[i_point_cloud].coords,
            n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].gnum =
    realloc(dist->points_cloud[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));
  dist->points_cloud[i_point_cloud].dist =
    realloc(dist->points_cloud[i_point_cloud].dist, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].proj =
    realloc(dist->points_cloud[i_point_cloud].proj, n_part * sizeof(double *));
  dist->points_cloud[i_point_cloud].closest_elt_gnum =
    realloc(dist->points_cloud[i_point_cloud].closest_elt_gnum,
            n_part * sizeof(PDM_g_num_t * ));

  for (int i = 0; i < n_part; i++) {
    dist->points_cloud[i_point_cloud].n_points[i] = -1;
    dist->points_cloud[i_point_cloud].coords[i] = NULL;
    dist->points_cloud[i_point_cloud].gnum[i] = NULL;
    dist->points_cloud[i_point_cloud].dist[i] = NULL;
    dist->points_cloud[i_point_cloud].proj[i] = NULL;
    dist->points_cloud[i_point_cloud].closest_elt_gnum[i] = NULL;
  }
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_dist_cloud_surf_cloud_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    i_part,
 const int                    n_points,
       double                *coords,
       PDM_g_num_t           *gnum
)
{

  dist->points_cloud[i_point_cloud].n_points[i_part] = n_points;
  dist->points_cloud[i_point_cloud].coords[i_part] = coords;
  dist->points_cloud[i_point_cloud].gnum[i_part] = gnum;
}


/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   dist           Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   mesh_nodal_id  Mesh nodal Pointer to \ref PDM_dist_cloud_surf object
 *
 */

void
PDM_dist_cloud_surf_nodal_mesh_set
(
 PDM_dist_cloud_surf_t *dist,
 PDM_part_mesh_nodal_t *mesh_nodal
 // PDM_Mesh_nodal_t      *mesh_nodal
)
{
  dist->mesh_nodal = mesh_nodal;
}




/**
 *
 * \brief Map a surface mesh
 *
 * \param [in]   dist       Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   surf_mesh  Surface mesh pointer
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_map
(
 PDM_dist_cloud_surf_t *dist,
 PDM_surf_mesh_t       *surf_mesh
)
{
  dist->_surf_mesh = surf_mesh;
}


/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   dist           Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_global_data_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    n_part
)
{
  assert (dist->surf_mesh == NULL);

  dist->surf_mesh = PDM_surf_mesh_create (n_part, dist->comm);
  dist->_surf_mesh = dist->surf_mesh;
}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   dist          Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_part        Partition to define
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering
 *                             to global vertex numbering
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_part_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_part,
 const int                    n_face,
 const int                   *face_vtx_idx,
 const int                   *face_vtx,
 const PDM_g_num_t           *face_ln_to_gn,
 const int                    n_vtx,
 const double                *coords,
 const PDM_g_num_t           *vtx_ln_to_gn
)
{
  PDM_surf_mesh_part_input (dist->surf_mesh,
                            i_part,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);
}

/**
 *
 * \brief Compute distance
 *
 * \param [in]   dist  Pointer to \ref PDM_dist_cloud_surf object
 *
 */

void
PDM_dist_cloud_surf_compute
(
 PDM_dist_cloud_surf_t *dist
)
{
  int use_optim = 0;

  char *env_var = NULL;
  env_var = getenv ("PDM_DIST_CLOUD_SURF_OPTIM");


  if (env_var != NULL) {
    use_optim = atoi(env_var);
  }


  if (use_optim) {
    _dist_cloud_surf_compute_optim(dist);
  }
  else {
    _dist_cloud_surf_compute(dist);
  }

}

/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   dist              Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud     Current cloud
 * \param [in]   i_part            Index of partition of the cloud
 * \param [out]  distance          Distance
 * \param [out]  projected         Projected point coordinates
 * \param [out]  closest_elt_g_num Global number of the closest element
 *
 */

void
PDM_dist_cloud_surf_get
(
       PDM_dist_cloud_surf_t  *dist,
 const int                     i_point_cloud,
 const int                     i_part,
       double                **distance,
       double                **projected,
       PDM_g_num_t           **closest_elt_gnum
)
{

  *distance         = dist->points_cloud[i_point_cloud].dist            [i_part];
  *projected        = dist->points_cloud[i_point_cloud].proj            [i_part];
  *closest_elt_gnum = dist->points_cloud[i_point_cloud].closest_elt_gnum[i_part];

  dist->results_is_getted = PDM_TRUE;
}


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  dist     Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_dist_cloud_surf_free
(
 PDM_dist_cloud_surf_t  *dist
)
{

  if(( dist->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dist->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !dist->results_is_getted)){
    for (int i_point_cloud = 0;
         i_point_cloud < dist->n_point_cloud;
         i_point_cloud++) {

      for (int i = 0; i < (dist->points_cloud[i_point_cloud]).n_part; i++) {
        free (dist->points_cloud[i_point_cloud].dist[i]);
        free (dist->points_cloud[i_point_cloud].proj[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_gnum[i]);
      }
    }
  }

  for (int i_point_cloud = 0;
       i_point_cloud < dist->n_point_cloud;
       i_point_cloud++) {

    free (dist->points_cloud[i_point_cloud].n_points);
    free (dist->points_cloud[i_point_cloud].coords);
    free (dist->points_cloud[i_point_cloud].gnum);
    free (dist->points_cloud[i_point_cloud].dist);
    free (dist->points_cloud[i_point_cloud].proj);
    free (dist->points_cloud[i_point_cloud].closest_elt_gnum);
  }

  free (dist->points_cloud);

  PDM_timer_free(dist->timer);

  if (dist->_surf_mesh != NULL) {
    if (dist->surf_mesh != NULL) {
      PDM_surf_mesh_free (dist->surf_mesh);
    }
  }

  free (dist);
}


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  dist     Pointer to \ref PDM_dist_cloud_surf object
 *
 */

void
PDM_dist_cloud_surf_dump_times
(
 PDM_dist_cloud_surf_t  *dist
)
{
  double t1 = dist->times_elapsed[END] - dist->times_elapsed[BEGIN];
  double t2 = dist->times_cpu    [END] - dist->times_cpu    [BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (dist->times_elapsed, t_elaps_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (dist->times_cpu, t_cpu_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  //-->>
  double t_elaps_min[NTIMER];
  PDM_MPI_Allreduce (dist->times_elapsed, t_elaps_min, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MIN, dist->comm);

  double t_cpu_min[NTIMER];
  PDM_MPI_Allreduce (dist->times_cpu, t_cpu_min, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MIN, dist->comm);
  //<<--

  int rank;
  PDM_MPI_Comm_rank (dist->comm, &rank);

  if (rank == 0) {


    PDM_printf( "distance timer : all                  (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "distance timer : Upper bound distance (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[UPPER_BOUND_DIST],
                t_cpu_max[UPPER_BOUND_DIST]);
    PDM_printf( "distance timer : Bbtree building      (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BBTREE_CREATE],
                t_cpu_max[BBTREE_CREATE]);
    PDM_printf( "distance timer : Candidate selection  (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[CANDIDATE_SELECTION],
                t_cpu_max[CANDIDATE_SELECTION]);
    PDM_printf( "distance timer : Load balacing of elementary computations of distance"
                " from the points to the candidates  (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[LOAD_BALANCING_ELEM_DIST],
                t_cpu_max[LOAD_BALANCING_ELEM_DIST]);
    PDM_printf( "distance timer : Computations of the distance"
                " from the points to the candidates    (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[COMPUTE_ELEM_DIST],
                t_cpu_max[COMPUTE_ELEM_DIST]);
    PDM_printf( "distance timer : Results exchange     (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[RESULT_TRANSMISSION],
                t_cpu_max[RESULT_TRANSMISSION]);
    PDM_printf_flush();
  }
}



/**
 *
 * \brief Get the dimension of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  n_points        Number of points
 *
 */

void
PDM_dist_cloud_surf_cloud_dim_get
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    i_part,
       int                   *n_points
)
{
  assert(dist != NULL);

  *n_points = dist->points_cloud[i_point_cloud].n_points[i_part];
}


#ifdef	__cplusplus
}
#endif
