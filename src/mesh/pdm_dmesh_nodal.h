#ifndef __PDM_DMESH_NODAL_H__
#define __PDM_DMESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

///*----------------------------------------------------------------------------
// * Geometric element type
// *----------------------------------------------------------------------------*/
//
//typedef enum {
//
//  PDM_MESH_NODAL_POINT,     
//  PDM_MESH_NODAL_BAR2,     
//  PDM_MESH_NODAL_TRIA3,     
//  PDM_MESH_NODAL_QUAD4,     
//  PDM_MESH_NODAL_POLY_2D,     
//  PDM_MESH_NODAL_TETRA4,     
//  PDM_MESH_NODAL_PYRAMID5,     
//  PDM_MESH_NODAL_PRISM6,     
//  PDM_MESH_NODAL_HEXA8,     
//  PDM_MESH_NODAL_POLY_3D     
//
//} PDM_Mesh_nodal_elt_t;


typedef struct _PDM_DMesh_nodal_t PDM_DMesh_nodal_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

int 
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm        
);


void
PDM_DMesh_nodal_free
(
const int hdl
);

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
);


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
);


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
);


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
);


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
);


/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_sections_get
(
const int   hdl
);


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
);


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
);


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


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
); 


/**
 * \brief  Return cell->face connectivity
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
      int   **cell_face_idx,  
PDM_g_num_t **cell_face  
); 


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
);


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
);


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
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */