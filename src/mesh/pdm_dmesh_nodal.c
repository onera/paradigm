
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
#include "pdm_geom_elem.h"

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
  int n_proc = 0;
  PDM_MPI_Comm_size (mesh->pdm_mpi_comm, &n_proc);

  mesh->n_proc  = n_proc;

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

  mesh->n_dface                  = -1;
  mesh->dface_vtx               = NULL;
  mesh->dface_vtx_idx           = NULL;
  mesh->face_distrib             = NULL;
  
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



/**
 * 
 * \brief _section_elt_faces_get
 *
 * \param [in]     mesh               Current mesh
 * \param [in]     id_section         Section identifier
 * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
 * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
 *   
 */

static int
_section_size_elt_faces_get
(
      PDM_DMesh_nodal_t *mesh,
      int               *s_elt_face_vtx_idx,
      int               *s_elt_face_vtx,
      int               *s_elt_face_cell
)
{
  
  int _s_elt_face_vtx_idx = 0;
  int _s_elt_face_vtx = 0;
  
  
  int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);
  
  for (int i = 0; i < n_sections_std; i++) {
    PDM_DMesh_nodal_section_std_t *section = 
      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
    int n_face_elt = 0;
    int n_sum_vtx_face = 0;
    
    switch (section->t_elt) {
    case PDM_MESH_NODAL_TRIA3:
      n_face_elt = 3;
      n_sum_vtx_face = 6; 
      break;
    case PDM_MESH_NODAL_TETRA4:
      n_face_elt = 4;
      n_sum_vtx_face = 12; 
      break;
    case PDM_MESH_NODAL_QUAD4:
      n_face_elt = 4;
      n_sum_vtx_face = 8; 
      break;
    case PDM_MESH_NODAL_HEXA8:
      n_face_elt = 6;
      n_sum_vtx_face = 24; 
      break;
    case PDM_MESH_NODAL_PYRAMID5:
      n_face_elt = 5;
      n_sum_vtx_face = 16; 
      break;
    case PDM_MESH_NODAL_PRISM6:
      n_face_elt = 5;
      n_sum_vtx_face = 18; 
      break;
    }
    
    _s_elt_face_vtx_idx += section->n_elt * n_face_elt;
    _s_elt_face_vtx += section->n_elt * n_sum_vtx_face;
  }
      
  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);
  
  for (int i = 0; i < n_sections_poly2d; i++) {
    PDM_DMesh_nodal_section_poly2d_t *section = 
      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
    _s_elt_face_vtx_idx += section->_connec_idx[section->n_elt];
    _s_elt_face_vtx += 2 * section->_connec_idx[section->n_elt]; 
  }
      
  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);
  
  for (int i = 0; i < n_sections_poly3d; i++) {
    PDM_DMesh_nodal_section_poly3d_t *section = 
      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
    _s_elt_face_vtx_idx += section->n_face;
    _s_elt_face_vtx += section->_facvtx[section->_facvtx_idx[section->n_face]];
  }

  *s_elt_face_cell = 2 * _s_elt_face_vtx_idx;
  *s_elt_face_vtx_idx = _s_elt_face_vtx_idx + 1;
  *s_elt_face_vtx = _s_elt_face_vtx + 1;
  

  return *s_elt_face_vtx - 1;
  
}  


/**
 * 
 * \brief _section_elt_faces_get
 *
 * \param [in]     mesh               Current mesh
 * \param [in]     id_section         Section identifier
 * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
 * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
 *   
 */

static void
_section_elt_faces_get
(
      PDM_DMesh_nodal_t *mesh,
const int                id_section,
      int               *n_elt_current,  
      int               *n_face_current,  
      int               *elt_face_vtx_idx,
      PDM_g_num_t       *elt_face_vtx,        
      PDM_g_num_t       *elt_face_cell        
)
{
  int _n_face_current = *n_face_current;
  int _n_elt_current = *n_elt_current;
  
  int         *_current_elt_face_vtx_idx = elt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elt_face_vtx     = elt_face_vtx + elt_face_vtx_idx[_n_face_current];        
  PDM_g_num_t *_current_elt_face_cell    = elt_face_cell + 2 * _n_face_current;        
  
  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    PDM_DMesh_nodal_section_std_t *section = 
            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");  
    }

    switch (section->t_elt) {

    case PDM_MESH_NODAL_TRIA3:
      {
        const int n_face_elt        = 3;
        const int n_sum_vtx_face    = 6;
        const int n_sum_vtx_elt     = 3;

        for (int ielt = 0; ielt < section->n_elt; ielt++) {
          for (int iface = 0; iface < n_face_elt; iface++) {
            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] = 
              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 3;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface)    ] = *n_elt_current + ielt + 1;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface) + 1] = 0;
          }

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
        }
        
        *n_elt_current += section->n_elt;

      }
      break;

    case PDM_MESH_NODAL_TETRA4: 
      {
        const int n_face_elt        = 4;
        const int n_sum_vtx_face    = 12;
        const int n_sum_vtx_elt     = 4;

        for (int ielt = 0; ielt < section->n_elt; ielt++) {
          for (int iface = 0; iface < n_face_elt; iface++) {
            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] = 
              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 3;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface)    ] = *n_elt_current + ielt + 1;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface) + 1] = 0;
          }

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt    ];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 2];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 1];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt    ];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 3];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];

          *n_elt_current += section->n_elt;

        }
      }
      break;
    case PDM_MESH_NODAL_QUAD4: 
      {
        const int n_face_elt        = 4;
        const int n_sum_vtx_face    = 8;
        const int n_sum_vtx_elt     = 4;

        for (int ielt = 0; ielt < section->n_elt; ielt++) {
          for (int iface = 0; iface < n_face_elt; iface++) {
            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] = 
              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 3;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface)    ] = *n_elt_current + ielt + 1;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface) + 1] = 0;
          }

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
        }
        
        *n_elt_current += section->n_elt;

      }
      break;
    
    case PDM_MESH_NODAL_HEXA8: 
      {
        const int n_face_elt        = 6;
        const int n_sum_vtx_face    = 24;
        const int n_sum_vtx_elt     = 8;

        for (int ielt = 0; ielt < section->n_elt; ielt++) {
          for (int iface = 0; iface < n_face_elt; iface++) {
            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] = 
              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 4;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface)    ] = *n_elt_current + ielt + 1;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface) + 1] = 0;
          }

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 6];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 7];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 4];        
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 5];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 7];        
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt    ];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 7];        
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 6];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 3];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt + 6];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 18] = section->_connec[n_sum_vtx_elt * ielt + 5];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 19] = section->_connec[n_sum_vtx_elt * ielt + 1];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 20] = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 21] = section->_connec[n_sum_vtx_elt * ielt + 5];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 22] = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 23] = section->_connec[n_sum_vtx_elt * ielt + 0];
        }
        
        *n_elt_current += section->n_elt;
        
      }
      break;
  
    case PDM_MESH_NODAL_PYRAMID5: 
      {
        const int n_face_elt        = 5;
        const int n_sum_vtx_face    = 1*4 + 4*3;
        const int n_sum_vtx_elt     = 5;

        elt_face_vtx_idx[0] = 0;

        for (int ielt = 0; ielt < section->n_elt; ielt++) {

          for (int iface = 0; iface < n_face_elt; iface++) {
            _current_elt_face_cell[2*(ielt * n_face_elt + iface)    ] = *n_elt_current + ielt + 1;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface) + 1] = 0;
          }

          _current_elt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 4;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 3;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 3;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 3;

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt    ];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 1];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 3];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt    ];

        }
        *n_elt_current += section->n_elt;
      }
      break;
    
    case PDM_MESH_NODAL_PRISM6: 
      {
        const int n_face_elt        = 5;
        const int n_sum_vtx_face    = 3*4 + 2*3;
        const int n_sum_vtx_elt     = 6;

        for (int ielt = 0; ielt < section->n_elt; ielt++) {

          for (int iface = 0; iface < n_face_elt; iface++) {
            _current_elt_face_cell[2*(ielt * n_face_elt + iface)    ] = *n_elt_current + ielt + 1;
            _current_elt_face_cell[2*(ielt * n_face_elt + iface) + 1] = 0;
          }

          _current_elt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 3;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 4;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 4;
          _current_elt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 4;

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt    ];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt + 4];        
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 5];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 3];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 5];        
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt    ];        
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 1];

          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 5];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
          _current_elt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt    ];

        }
        *n_elt_current += section->n_elt;
    
      }
      break;
    }
  }
  
  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    int _id = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
    PDM_DMesh_nodal_section_poly2d_t *section = 
            (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id);
    
    int idx = 0;
    elt_face_vtx_idx[0] = 0;

    for (int ielt = 0; ielt < section->n_elt; ielt++) {
      int n_face_elt = section->_connec_idx[section->n_elt];
      int idx2 = section->_connec_idx[ielt];
      for (int iface = 0; iface < n_face_elt; iface++) {
        elt_face_vtx_idx[idx + 1] = elt_face_vtx_idx[idx] + 2;
        int inext = (iface + 1) % n_face_elt;
        elt_face_vtx[2 * idx    ] = section->_connec[idx2 + iface];
        elt_face_vtx[2 * idx + 1] = section->_connec[idx2 + inext];
        elt_face_cell[idx] = *n_elt_current + ielt + 1;
      }
    }
    
    *n_elt_current += section->n_elt;
    
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");  
    }
    
  }
  
  else {
    
    int _id = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
    PDM_DMesh_nodal_section_poly3d_t *section = 
            (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");  
    }

  }
  
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
     
    if (mesh->dface_vtx_idx != NULL) {
      free (mesh->dcell_face_idx);
    }
     
    if (mesh->dface_vtx != NULL) {
      free (mesh->dface_vtx);
    }
     
    if (mesh->face_distrib != NULL) {
      free (mesh->face_distrib);
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

  vtx->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  PDM_g_num_t *_distrib = vtx->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_vtx = n_vtx;
  
  PDM_MPI_Scan (&_n_vtx, _distrib, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);


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

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  PDM_g_num_t *_distrib = section->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_elt = n_elt;
  
  PDM_MPI_Scan (&_n_elt, _distrib, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
 
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
  
  const PDM_DMesh_nodal_section_std_t *section = 
  (const PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, _id_section);
  
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


/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *  
 */

int
PDM_DMesh_nodal_section_n_elt_get 
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

  int _id_section;
  
  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_DMesh_nodal_section_poly3d_t *section = 
    (const PDM_DMesh_nodal_section_poly3d_t *) 
     PDM_Handles_get (mesh->sections_poly3d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }
      
    return section->n_elt;
  }
  
  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_DMesh_nodal_section_poly2d_t *section = 
    (const PDM_DMesh_nodal_section_poly2d_t *) 
     PDM_Handles_get (mesh->sections_poly2d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }
  
    return section->n_elt;
  }
  
  else {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_DMesh_nodal_section_std_t *section = 
    (const PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }
  
    return section->n_elt;
  }
  
}


/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */
 
void
PDM_DMesh_nodal_section_poly2d_set 
(
const int            hdl,
const int            id_section, 
const PDM_l_num_t    n_elt,    
      PDM_l_num_t   *connec_idx,   
      PDM_g_num_t   *connec
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly2d_t *section = 
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  /* Mapping */

  section->n_elt       = n_elt;
  section->_connec_idx = connec_idx;
  section->_connec     = connec;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  PDM_g_num_t *_distrib = section->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_elt = n_elt;
  
  PDM_MPI_Scan (&_n_elt, _distrib, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  
}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */
 
void
PDM_DMesh_nodal_section_poly2d_get 
(
 const int          hdl,
 const int          id_section, 
       PDM_l_num_t  **connec_idx,   
       PDM_g_num_t  **connec
) 
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly2d_t *section = 
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }
  
  *connec_idx = section->_connec_idx;
  *connec     = section->_connec;
 
}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_set 
(
const int            hdl,
const int            id_section, 
const PDM_l_num_t    n_elt,    
const PDM_l_num_t    n_face,   
      PDM_l_num_t   *facvtx_idx,   
      PDM_g_num_t   *facvtx,
      PDM_l_num_t   *cellfac_idx,   
      PDM_g_num_t   *cellfac
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly3d_t *section = 
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  section->n_elt        = n_elt;
  section->n_face       = n_face;
  section->_facvtx_idx  = facvtx_idx;
  section->_facvtx      = facvtx;
  section->_cellfac_idx = cellfac_idx;
  section->_cellfac     = cellfac;
  
  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  PDM_g_num_t *_distrib = section->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_elt = n_elt;
  
  PDM_MPI_Scan (&_n_elt, _distrib, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  
  
}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_get 
(
const int            hdl,
const int            id_section, 
      PDM_l_num_t   *n_face,   
      PDM_l_num_t  **facvtx_idx,   
      PDM_g_num_t  **facvtx,
      PDM_l_num_t  **cellfac_idx,   
      PDM_g_num_t  **cellfac
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly3d_t *section = 
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *n_face      = section->n_face;
  *facvtx_idx  = section->_facvtx_idx;
  *facvtx      = section->_facvtx;
  *cellfac_idx = section->_cellfac_idx;
  *cellfac     = section->_cellfac;

}


/**
 * \brief  Return total number of elements of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return number elements of a partition
 * 
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_cell_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  PDM_g_num_t total_n_cell = 0;
  
  int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);
  
  for (int i = 0; i < n_sections_std; i++) {
    PDM_DMesh_nodal_section_std_t *_section_std = 
      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
    total_n_cell += _section_std->distrib[mesh->n_proc];
  }
      
  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);
  
  for (int i = 0; i < n_sections_poly2d; i++) {
    PDM_DMesh_nodal_section_poly2d_t *_section_poly2d = 
      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
    total_n_cell += _section_poly2d->distrib[mesh->n_proc];
  }
      
  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);
  
  for (int i = 0; i < n_sections_poly3d; i++) {
    PDM_DMesh_nodal_section_poly3d_t *_section_poly3d = 
      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
    total_n_cell += _section_poly3d->distrib[mesh->n_proc];
  }
  
  return total_n_cell;

}


/**
 * \brief  Return total number of faces of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of faces
 * 
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_face_get
(
const int  hdl
) 
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  if (mesh->dcell_face == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
  }
  
  return mesh->face_distrib[mesh->n_proc];
  
}


/**
 * \brief  Return total number of vertices of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of vertices
 * 
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_vtx_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->distrib[mesh->n_proc];
}


/**
 * \brief  Compute cell->face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * 
 */

void
PDM_DMesh_nodal_cell_face_compute
(
const int   hdl
)
{
  PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
}


/**
 * \brief  Return cell \rightarrow face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dcell_face_idx   Index of distributed cell->face connectivity
 * \param [out]  dcell_face       Distributed cell->face connectivity
 *
 * \return     Number of cells on the current process
 *  
 */

int
PDM_DMesh_nodal_cell_face_get
(
const int   hdl,
      int   **dcell_face_idx,  
PDM_g_num_t **dcell_face  
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  *dcell_face_idx = mesh->dcell_face_idx;  
  *dcell_face = mesh->dcell_face;

  return mesh->n_dcell;

}


/**
 * \brief  Return face \rightarrow vertex connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dcell_face_idx   Index of distributed cell->face connectivity
 * \param [out]  dcell_face       Distributed cell->face connectivity
 *
 * \return     Number of faces on the current process
 *  
 */

int
PDM_DMesh_nodal_face_vtx_get
(
const int   hdl,
      int   **dface_vtx_idx,  
PDM_g_num_t **dface_vtx  
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  *dface_vtx_idx = mesh->dface_vtx_idx;  
  *dface_vtx     = mesh->dface_vtx;

  return mesh->n_dface;

}

/**
 * \brief  Return vertices distribution
 * 
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_vtx_get
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

  return vtx->distrib;

}


/**
 * \brief  Return section distribution
 * 
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_section_get
(
 const int   hdl,
 const int   id_section     
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_section;
  
  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_DMesh_nodal_section_poly3d_t *section = 
    (const PDM_DMesh_nodal_section_poly3d_t *) 
     PDM_Handles_get (mesh->sections_poly3d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }
      
    return section->distrib;
  }
  
  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_DMesh_nodal_section_poly2d_t *section = 
    (const PDM_DMesh_nodal_section_poly2d_t *) 
     PDM_Handles_get (mesh->sections_poly2d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }
  
    return section->distrib;
  }
  
  else {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_DMesh_nodal_section_std_t *section = 
    (const PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }
  
    return section->distrib;
  }
}


/**
 * \brief  Return cell distribution
 * 
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_cell_get
(
 const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  return mesh->cell_distrib;
  
}


/**
 * \brief  Return face distribution
 * 
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_face_get
(
 const int hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  return mesh->face_distrib;
  
}

        


#ifdef __cplusplus
}
#endif /* __cplusplus */

