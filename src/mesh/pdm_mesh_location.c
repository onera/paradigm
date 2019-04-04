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
#include "pdm_mesh_dist.h"
#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
#include "pdm_octree.h"
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"
#include "pdm_hash_tab.h"
#include "pdm_mesh_location.h"

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

#define NTIMER 2
  
/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ol_timer_step_t
 *
 */
 
typedef enum {

  BEGIN                         = 0,
  END                           = 1,        

} _ol_timer_step_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 * 
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;
  PDM_g_num_t **location;
  double      **uvw;
  int         **weights_idx;
  double      **weights;
  
} _point_cloud_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 * 
 */

typedef struct {

  int  n_point_cloud; /*!< Number of point clouds */
  PDM_MPI_Comm comm;  /*!< MPI communicator */ 

  PDM_mesh_nature_t mesh_nature;  /*!< Nature of the mesh */
  
  int  shared_nodal;   /*!< 1 if mesh nodal is shared, 0 otherwise */
  int  mesh_nodal_id;  /*!< Volume mesh identifier */

  _point_cloud_t *point_clouds; /*!< Point clouds */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */
  
  double times_cpu[NTIMER];     /*!< CPU time */
  
  double times_cpu_u[NTIMER];  /*!< User CPU time */
  
  double times_cpu_s[NTIMER];  /*!< System CPU time */

  
} _PDM_location_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_locations   = NULL;

static int idebug = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/
  
/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _PDM_location_t *
_get_from_id
(
 int  id
)
{
  _PDM_location_t *location = (_PDM_location_t *) PDM_Handles_get (_locations, id);
    
  if (location == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad identifier\n");
  }

  return location;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute the location of point clouds inta a mesh
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 *
 */

int
PDM_mesh_location_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
)
{
  if (_locations == NULL) {
    _locations = PDM_Handles_create (4);
  }

  _PDM_location_t *location = (_PDM_location_t *) malloc(sizeof(_PDM_location_t));

  int id = PDM_Handles_store (_locations, location);

 location->n_point_cloud = n_point_cloud;
 location->comm = comm;
 location->mesh_nature = mesh_nature;

 location->shared_nodal  = 0;
 location->mesh_nodal_id = -1;

  location->point_clouds =
    (_point_cloud_t*) malloc (sizeof(_point_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    location->point_clouds[i].n_part = -1;
    location->point_clouds[i].n_points = NULL;
    location->point_clouds[i].coords = NULL;
    location->point_clouds[i].gnum = NULL;
    location->point_clouds[i].location = NULL;
    location->point_clouds[i].uvw = NULL;
    location->point_clouds[i].weights = NULL;
    location->point_clouds[i].weights_idx = NULL;
  }

  location->timer = PDM_timer_create ();
  
  for (int i = 0; i < NTIMER; i++) {
    location->times_elapsed[i] = 0.;
    location->times_cpu[i] = 0.;
    location->times_cpu_u[i] = 0.;
    location->times_cpu_s[i] = 0.;
  }
  
  return id;

}

void
PDM_mesh_location_create_cf 
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Fint comm,
 int *id
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_mesh_location_create (mesh_nature, n_point_cloud, _comm);
  
}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_location_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
)
{
  _PDM_location_t *location = _get_from_id (id);

  location->point_clouds[i_point_cloud].n_part = n_part;
  location->point_clouds[i_point_cloud].n_points =
    realloc(location->point_clouds[i_point_cloud].n_points, n_part * sizeof(int));
  location->point_clouds[i_point_cloud].coords =
    realloc(location->point_clouds[i_point_cloud].coords,
            n_part * sizeof(double *));
  location->point_clouds[i_point_cloud].gnum =
    realloc(location->point_clouds[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));
  
  for (int i = 0; i < n_part; i++) {
    location->point_clouds[i_point_cloud].n_points[i] = -1;
    location->point_clouds[i_point_cloud].coords[i] = NULL;
    location->point_clouds[i_point_cloud].gnum[i] = NULL;
  }

}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number 
 *
 */

void
PDM_mesh_location_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
)
{
  _PDM_location_t *location = _get_from_id (id);

  location->point_clouds[i_point_cloud].n_points[i_part] = n_points;
  location->point_clouds[i_point_cloud].coords[i_part] = coords;
  location->point_clouds[i_point_cloud].gnum[i_part] = gnum;

}


/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Identifier
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 *
 */

void
PDM_mesh_location_shared_nodal_mesh_set
(
 const int  id,
 const int  mesh_nodal_id
)
{
  _PDM_location_t *location = _get_from_id (id);
  
  location->mesh_nodal_id = mesh_nodal_id;
  location->shared_nodal = 1;
}


/**
 *
 * \brief Set global data of a volume mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_mesh_location_volume_mesh_global_data_set
(
 const int  id,
 const int  n_part
)
{
  _PDM_location_t *location = _get_from_id (id);

  if ((location->shared_nodal == 0) && (location->mesh_nodal_id != -1)) {
    PDM_Mesh_nodal_free (location->mesh_nodal_id);
  }

    location->mesh_nodal_id = PDM_Mesh_nodal_create (n_part, location->comm);
}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define  
 * \param [in]   n_cell        Number of cells                     
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering 
 * \param [in]   n_face        Number of faces                     
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering 
 * \param [in]   n_vtx         Number of vertices              
 * \param [in]   coords        Coordinates       
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering 
 *
 */

void
PDM_mesh_location_volume_part_set
(
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_face_idx,
 const int         *cell_face,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx, 
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
)
{
  /* PDM_l_num_t      *face_vtx_nb; */
  /* PDM_l_num_t      *cell_face_nb; */

  /* PDM_Mesh_nodal_cell3d_cellface_add (location->mesh_nodal_id, */
  /*                                     i_part, */
  /*                                     n_cell, */
  /*                                     n_face, */

  /*                                     face_vtx_idx, */
  /*                                     face_vtx_nb, */
  /*                                     face_vtx, */
  /*                                     cell_face_idx, */
  /*                                     cell_face_nb, */
  /*                                     cell_face, */
  /*                                     PDM_g_num_t      *numabs);  */
}


/**
 *
 * \brief Compute point location
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_location_compute
(
 const int id
)
{
}


/**
 *
 * \brief Get mesh location
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  location_elt_g_num    The global number of the closest element if the point is located,
 *                                     -1 otherwise
 *
 */

void
PDM_mesh_location_get
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
       PDM_g_num_t **location_elt_gnum
)
{
}


/**
 *
 * \brief Free a locationd mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed. 
 *                       Otherwise, results are kept. 
 *
 */

void
PDM_mesh_location_free
(
 const int id,
 const int partial
 );

  
/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_mesh_location_dump_times
(
 const int id
);

#ifdef	__cplusplus
}
#endif
