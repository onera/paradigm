
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_BLOCK_ID_BLOCK_STD    = 0,     
  PDM_BLOCK_ID_BLOCK_POLY2D = 1000000,
  PDM_BLOCK_ID_BLOCK_POLY3D = 2000000

} PDM_section_id_section_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * \brief Storage of mesh handles
 *
 */

static PDM_Handles_t *mesh_handles = NULL; 

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * 
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 * 
 */

static 
PDM_DMesh_nodal_vtx_t *
_vtx_free
(
 PDM_DMesh_nodal_vtx_t *vtx
)
{
  
  if (vtx != NULL) {

    if (vtx->distrib != NULL) {
      free (vtx->distrib);
      vtx->distrib = NULL;
    }

    free (vtx);
  }
  return NULL;
}


/**
 * 
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void   
_mesh_init
(
PDM_DMesh_nodal_t *mesh,
const PDM_MPI_Comm comm        
)
{

  mesh->n_som_abs                = -1;          
  mesh->n_cell_abs               = -1;          

  mesh->vtx                      = malloc(sizeof(PDM_DMesh_nodal_vtx_t ));
  mesh->vtx->_coords             = NULL;
  mesh->vtx->n_vtx               = 0;

  mesh->sections_std             = NULL;  
  mesh->sections_poly2d          = NULL;            
  mesh->sections_poly3d          = NULL;           

  mesh->pdm_mpi_comm             = comm;
  mesh->sections_id              = NULL;
  mesh->n_sections               = 0;

  mesh->n_dcell                  = -1;
  mesh->dcell_face               = NULL;
  mesh->dcell_face_idx           = NULL;
  mesh->cell_distrib             = NULL;
  
}

/**
 * 
 * \brief Update sections identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static void   
_update_sections_id
(
PDM_DMesh_nodal_t *mesh
)
{
  int n_sections = 0;
  
  if (mesh->sections_std != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_std);
  }
    
  if (mesh->sections_poly2d != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_poly2d);
  }

  if (mesh->sections_poly3d != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_poly3d);
  }
  
  if (mesh->n_sections < n_sections) {
    mesh->sections_id = (int *) realloc(mesh->sections_id, sizeof(int) * n_sections);
  }
  
  int k = 0;
  if (mesh->sections_std != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_std);
    int n = PDM_Handles_n_get (mesh->sections_std);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_STD;
    }
  }
    
  if (mesh->sections_poly2d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_poly2d);
    int n = PDM_Handles_n_get (mesh->sections_poly2d);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }
  
  if (mesh->sections_poly3d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_poly3d);
    int n = PDM_Handles_n_get (mesh->sections_poly3d);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }
  
  mesh->n_sections = n_sections;  

}


/**
 * 
 * \brief Free a standard section
 *
 * \param [inout]  _bloc_std    Standard section
 *
 * \return         Null
 *   
 */

static
PDM_DMesh_nodal_section_std_t *
_section_std_free
(
PDM_DMesh_nodal_section_std_t *_section_std
)
{

  if (_section_std == NULL) {
    return NULL;
  }
  
  if (_section_std->distrib != NULL) {
    free (_section_std->distrib);
    _section_std->distrib = NULL;
  }
 
  free(_section_std);
  return NULL;
}


/**
 * 
 * \brief Free a polygon section
 *
 * \param [inout]  _bloc_poly2d    Polygon section
 *
 * \return         Null
 *   
 */

static
PDM_DMesh_nodal_section_poly2d_t *
_section_poly2d_free
(
PDM_DMesh_nodal_section_poly2d_t *_section_poly2d
)
{
  
  if (_section_poly2d == NULL) {
    return NULL;
  }
  
  if (_section_poly2d->distrib != NULL) {
    free (_section_poly2d->distrib);
    _section_poly2d->distrib = NULL;
  }
 
  
  free(_section_poly2d);

  return NULL;
}


/**
 * 
 * \brief Free a polyhedron section
 *
 * \param [inout]  _section_poly3d    Polyhedron section
 *
 * \return         Null
 *   
 */

static
PDM_DMesh_nodal_section_poly3d_t *
_section_poly3d_free
(
PDM_DMesh_nodal_section_poly3d_t *_section_poly3d
)
{
  if (_section_poly3d == NULL) {
    return NULL;
  }
  
  if (_section_poly3d->distrib != NULL) {
    free (_section_poly3d->distrib);
    _section_poly3d->distrib = NULL;
  }
 
  
  free(_section_poly3d);

  return NULL;

}




/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

int 
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm        
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) malloc (sizeof(PDM_DMesh_nodal_t));
  
  _mesh_init (mesh, comm);
  
  if (mesh_handles == NULL) {
    mesh_handles = PDM_Handles_create (4);
  }
  
  return PDM_Handles_store (mesh_handles, mesh);
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_DMesh_nodal_free
(
const int hdl
)
{
  
  PDM_DMesh_nodal_t * mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh != NULL) {

    if (mesh->sections_id != NULL) {
      free (mesh->sections_id);
    }
    
    mesh->sections_id = NULL;
    
     _vtx_free (mesh->vtx);

    /* free standard sections */

    if (mesh->sections_std != NULL) {
      int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

      while (n_sections_std > 0) {
        PDM_DMesh_nodal_section_std_t *_bloc_std = 
          (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[0]);
        _section_std_free(_bloc_std);
        PDM_Handles_handle_free (mesh->sections_std, list_ind[0], PDM_FALSE);
        n_sections_std = PDM_Handles_n_get (mesh->sections_std);
      }

      mesh->sections_std = PDM_Handles_free (mesh->sections_std); 
    }

    /* Free polygon sections */ 

    if (mesh->sections_poly2d != NULL) {
      int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

      while (n_sections_poly2d > 0) {
        PDM_DMesh_nodal_section_poly2d_t *_bloc_poly2d = 
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[0]);
        _section_poly2d_free(_bloc_poly2d);
        PDM_Handles_handle_free (mesh->sections_poly2d, list_ind[0], PDM_FALSE);
        n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
      }

      mesh->sections_poly2d = PDM_Handles_free (mesh->sections_poly2d); 
    }

    /* Free polyhedron sections */ 

    if (mesh->sections_poly3d != NULL) {
      int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

      while (n_sections_poly3d > 0) {
        PDM_DMesh_nodal_section_poly3d_t *_bloc_poly3d = 
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[0]);
        _section_poly3d_free(_bloc_poly3d);
        PDM_Handles_handle_free (mesh->sections_poly3d, list_ind[0], PDM_FALSE);
        n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
      }

      mesh->sections_poly3d = PDM_Handles_free (mesh->sections_poly3d); 
    }
    
    if (mesh->sections_id != NULL) {
      free (mesh->sections_id);
    }

    if (mesh->dcell_face_idx != NULL) {
      free (mesh->dcell_face_idx);
    }
     
    if (mesh->dcell_face != NULL) {
      free (mesh->dcell_face);
    }
     
    if (mesh->cell_distrib != NULL) {
      free (mesh->cell_distrib);
    }
     
    free(mesh);

    PDM_Handles_handle_free (mesh_handles, hdl, PDM_FALSE);
  
    int n_mesh_array = PDM_Handles_n_get (mesh_handles); 

    if (n_mesh_array == 0) {
      mesh_handles = PDM_Handles_free (mesh_handles);
    }
  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_DMesh_nodal_coord_set
(
 const int          hdl,
 const int          n_vtx,  
 const PDM_real_t  *coords  
)
{
  PDM_DMesh_nodal_t * mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  if (vtx->_coords != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;

}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Number of vertices
 *
 */

int
PDM_DMesh_nodal_n_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->n_vtx;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_DMesh_nodal_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->_coords;
}


/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_sections_get
(
 const int   hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  
  return mesh->n_sections;

}


/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */

int *
PDM_DMesh_nodal_sections_id_get
(
const int   hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  return mesh->sections_id;

}


/**
 * \brief  Return type of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */

PDM_Mesh_nodal_elt_t
PDM_Mesh_nodal_section_type_get
(
const int   hdl,
const int   id_section     
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  PDM_Mesh_nodal_elt_t t_elt;
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    t_elt = PDM_MESH_NODAL_POLY_3D;
    PDM_DMesh_nodal_section_std_t *section = 
            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");  
    }

    t_elt = section->t_elt;
  }
  
  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {
    
    t_elt = PDM_MESH_NODAL_POLY_2D;

  }
  
  else {
    
    t_elt = PDM_MESH_NODAL_POLY_3D;

  }
  
  return t_elt;
    
}


/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory 
 *                             when the section is destroyed
 * \param [in]  id_section       Block identifier
 *
 * \return Block identifier     
 *
 */

int 
PDM_Mesh_nodal_section_add 
(
const int                    hdl,
const PDM_Mesh_nodal_elt_t   t_elt
)
{

  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int id_section;
  
  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :    
  case PDM_MESH_NODAL_BAR2     :   
  case PDM_MESH_NODAL_TRIA3    :    
  case PDM_MESH_NODAL_QUAD4    :    
  case PDM_MESH_NODAL_TETRA4   :     
  case PDM_MESH_NODAL_PYRAMID5 :     
  case PDM_MESH_NODAL_PRISM6   :     
  case PDM_MESH_NODAL_HEXA8    : 
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_std == NULL) {
        mesh->sections_std = PDM_Handles_create (4);
      } 
      
      /* Allocation du bloc */
      
      PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) malloc(sizeof(PDM_DMesh_nodal_section_std_t));

      id_section = PDM_Handles_store (mesh->sections_std, section_std);

      /* Intialisation du bloc */

      section_std->t_elt = t_elt;

      section_std->n_elt      = -1;
      section_std->_connec    = NULL;
      section_std->distrib    = NULL;
      
      id_section += PDM_BLOCK_ID_BLOCK_STD;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard sections must be less than %d\n", 
               PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }
    
    break;

  case PDM_MESH_NODAL_POLY_2D  :    
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_poly2d == NULL) {
        mesh->sections_poly2d = PDM_Handles_create (4);
      }
      
      /* Allocation du bloc */
      
      PDM_DMesh_nodal_section_poly2d_t *section_poly2d =
              (PDM_DMesh_nodal_section_poly2d_t *) malloc(sizeof(PDM_DMesh_nodal_section_poly2d_t));

      id_section = PDM_Handles_store (mesh->sections_poly2d, section_poly2d);
      
      /* Intialisation du bloc */

      section_poly2d->n_elt       = -1;
      section_poly2d->_connec_idx = NULL;
      section_poly2d->_connec     = NULL;
      section_poly2d->distrib     = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_POLY2D;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of polygon sections must be less than %d\n",
               PDM_BLOCK_ID_BLOCK_POLY3D - PDM_BLOCK_ID_BLOCK_POLY2D);
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :     
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_poly3d == NULL) {
        mesh->sections_poly3d = PDM_Handles_create (4);
      }
      
      /* Allocation du bloc */
      
      PDM_DMesh_nodal_section_poly3d_t *section_poly3d =
              (PDM_DMesh_nodal_section_poly3d_t *) malloc(sizeof(PDM_DMesh_nodal_section_poly3d_t));

      id_section = PDM_Handles_store (mesh->sections_poly3d, section_poly3d);

      /* Intialisation du bloc */

      section_poly3d->n_elt        = -1;
      section_poly3d->n_face       = -1;
      section_poly3d->_facvtx_idx  = NULL;
      section_poly3d->_facvtx      = NULL;
      section_poly3d->_cellfac_idx = NULL;
      section_poly3d->_cellfac     = NULL;
      section_poly3d->distrib      = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");

  }

  _update_sections_id (mesh);
  return id_section ;
  
}


/**
 * \brief Define a standard section
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 *
 */

void
PDM_Mesh_nodal_section_std_set 
(
const int          hdl,
const int          id_section,     
const int          n_elt,    
      PDM_g_num_t *connec   
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
  PDM_DMesh_nodal_section_std_t *section = (PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }
   
  /* Mapping */
  
  section->n_elt = n_elt;
  section->_connec = connec;
 
}


/**
 * \brief Return standard section description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *   
PDM_DMesh_nodal_section_std_get 
(   
const int            hdl,
const int            id_section     
) 
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
  const PDM_DMesh_nodal_section_std_t *section = (const PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }
  
  return section->_connec;
}


/**
 * \brief Get number of section elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Number of elements
 *  
 */

//int
//PDM_DMesh_nodal_section_n_elt_get 
//(   
//const int            hdl,
//const int            id_section,     
//const int            id_part 
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//
//  int _id_section;
//  
//  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
//  
//    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
//  
//    const PDM_DMesh_nodal_section_poly3d_t *section = (const PDM_DMesh_nodal_section_poly3d_t *) 
//     PDM_Handles_get (mesh->sections_poly3d, _id_section);
//  
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//    }
//  
//    if (id_part >= section->n_part) {
//      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//    }
//    
//    return section->n_elt[id_part];
//  }
//  
//  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
//  
//    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//  
//    const PDM_DMesh_nodal_section_poly2d_t *section = (const PDM_DMesh_nodal_section_poly2d_t *) 
//     PDM_Handles_get (mesh->sections_poly2d, _id_section);
//  
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//    }
//  
//    if (id_part >= section->n_part) {
//      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//    }
//
//    return section->n_elt[id_part];
//  }
//  
//  else {
//  
//    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
//  
//    const PDM_DMesh_nodal_section_std_t *section = (const PDM_DMesh_nodal_section_std_t *) 
//     PDM_Handles_get (mesh->sections_std, _id_section);
//  
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//    }
//  
//    if (id_part >= section->n_part) {
//      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//    }
//
//    return section->n_elt[id_part];
//  }
//  
//}
//
//
///**
// * \brief Get global numbering of section elements
// *
// * \param [in]  idx            Nodal mesh handle
// * \param [in]  id_section       Block identifier
// * \param [in]  id_part        Partition identifier
// *
// * \return      Return global numbering of section elements
// *  
// */
//
//PDM_g_num_t *
//PDM_DMesh_nodal_section_g_num_get 
//(   
//const int            hdl,
//const int            id_section,     
//const int            id_part 
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//
//  int _id_section;
//  
//  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
//  
//    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
//  
//    const PDM_DMesh_nodal_section_poly3d_t *section = (const PDM_DMesh_nodal_section_poly3d_t *) 
//     PDM_Handles_get (mesh->sections_poly3d, _id_section);
//  
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//    }
//  
//    if (id_part >= section->n_part) {
//      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//    }
//    
//    return section->numabs_int[id_part];
//  }
//  
//  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
//  
//    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//  
//    const PDM_DMesh_nodal_section_poly2d_t *section = (const PDM_DMesh_nodal_section_poly2d_t *) 
//     PDM_Handles_get (mesh->sections_poly2d, _id_section);
//  
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//    }
//  
//    if (id_part >= section->n_part) {
//      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//    }
//
//    return section->numabs_int[id_part];
//  }
//  
//  else {
//  
//    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
//  
//    const PDM_DMesh_nodal_section_std_t *section = (const PDM_DMesh_nodal_section_std_t *) 
//     PDM_Handles_get (mesh->sections_std, _id_section);
//  
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//    }
//  
//    if (id_part >= section->n_part) {
//      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//    }
//
//    return section->numabs_int[id_part];
//  }
//}


/**
 * \brief Define a polygon section
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */
 
//void
//PDM_DMesh_nodal_section_poly2d_set 
//(
//const int            hdl,
//const int            id_section, 
//const int            id_part, 
//const PDM_l_num_t    n_elt,    
//      PDM_l_num_t   *connec_idx,   
//      PDM_l_num_t   *connec,
//      PDM_g_num_t   *numabs,
//      PDM_l_num_t   *parent_num
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//  
//  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//  
//  PDM_DMesh_nodal_section_poly2d_t *section = 
//          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);
//
//  if (section == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//  }
//  
//  if (id_part >= section->n_part) {
//    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//  }
//
//  /* Mapping */
//
//  mesh->n_cell[id_part]      += -section->n_elt[id_part]; 
//  mesh->n_cell[id_part]      += n_elt;
//  section->n_elt[id_part]       = n_elt;
//  section->_connec_idx[id_part] = connec_idx;
//  section->_connec[id_part]     = connec;
//  section->_numabs[id_part]     = numabs;
//
//  for (int i = 0; i < n_elt; i++) {
//    mesh->n_elt_abs = PDM_MAX(mesh->n_elt_abs, numabs[i]);
//  }
//
//  if (parent_num != NULL) {
//    if (section->_parent_num == NULL) {
//      section->_parent_num = malloc (sizeof(PDM_l_num_t *) * section->n_part);
//    }
//    section->_parent_num[id_part] = parent_num;
//  }
//  
//}
//
//
//
///**
// * \brief Return a polygon section description
// *
// * \param [in]  idx            Nodal mesh handle
// * \param [in]  id_section       Block identifier
// * \param [in]  id_part        Partition identifier
// * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
// * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
// *
// */
// 
//void
//PDM_DMesh_nodal_section_poly2d_get 
//(
// const int          hdl,
// const int          id_section, 
// const int          id_part, 
//       PDM_l_num_t  **connec_idx,   
//       PDM_l_num_t  **connec
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//  
//  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//  
//  PDM_DMesh_nodal_section_poly2d_t *section = 
//          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);
//
//  if (section == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//  }
//  
//  if (id_part >= section->n_part) {
//    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//  }
//
//  *connec_idx = section->_connec_idx[id_part];
//  *connec     = section->_connec[id_part];
// 
//}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

//void
//PDM_DMesh_nodal_section_poly3d_set 
//(
//const int            hdl,
//const int            id_section, 
//const int            id_part, 
//const PDM_l_num_t    n_elt,    
//const PDM_l_num_t    n_face,   
//      PDM_l_num_t   *facvtx_idx,   
//      PDM_l_num_t   *facvtx,
//      PDM_l_num_t   *cellfac_idx,   
//      PDM_l_num_t   *cellfac,
//      PDM_g_num_t   *numabs,
//      PDM_l_num_t   *parent_num
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//  
//  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//  
//  PDM_DMesh_nodal_section_poly3d_t *section = 
//          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);
//
//  if (section == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//  }
//  
//  if (id_part >= section->n_part) {
//    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//  }
//
//  mesh->n_cell[id_part]       += -section->n_elt[id_part]; 
//  mesh->n_cell[id_part]       += n_elt;
//  section->n_elt[id_part]        = n_elt;
//  section->n_face[id_part]       = n_face;
//  section->_facvtx_idx[id_part]  = facvtx_idx;
//  section->_facvtx[id_part]      = facvtx;
//  section->_cellfac_idx[id_part] = cellfac_idx;
//  section->_cellfac[id_part]     = cellfac;
//  section->_numabs[id_part]      = numabs;
//
//  for (int i = 0; i < n_elt; i++) {
//    mesh->n_elt_abs = PDM_MAX (mesh->n_elt_abs, numabs[i]);
//  }
//
//  if (parent_num != NULL) {
//    if (section->_parent_num == NULL) {
//      section->_parent_num = malloc (sizeof(PDM_l_num_t *) * section->n_part);
//    }
//    section->_parent_num[id_part] = parent_num;
//  }
//  
//}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

//void
//PDM_DMesh_nodal_section_poly3d_get 
//(
//const int            hdl,
//const int            id_section, 
//const int            id_part, 
//      PDM_l_num_t   *n_face,   
//      PDM_l_num_t  **facvtx_idx,   
//      PDM_l_num_t  **facvtx,
//      PDM_l_num_t  **cellfac_idx,   
//      PDM_l_num_t  **cellfac
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//  
//  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//  
//  PDM_DMesh_nodal_section_poly3d_t *section = 
//          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);
//
//  if (section == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//  }
//  
//  if (id_part >= section->n_part) {
//    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//  }
//
//  *n_face     = section->n_face[id_part];
//  *facvtx_idx  = section->_facvtx_idx[id_part];
//  *facvtx      = section->_facvtx[id_part];
//  *cellfac_idx = section->_cellfac_idx[id_part];
//  *cellfac     = section->_cellfac[id_part];
//
//}


/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Return number elements of a partition
 * 
 */

//int
//PDM_DMesh_nodal_n_cell_get
//(
//const int  hdl,
//const int  id_part 
//)
//{
//  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, idx);
//  
//  if (mesh == NULL) {
//    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
//  }
//
//  if (id_part >= mesh->n_part) {
//    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
//  }
//
//  return mesh->n_cell[id_part];
//
//}

#ifdef __cplusplus
}
#endif /* __cplusplus */

