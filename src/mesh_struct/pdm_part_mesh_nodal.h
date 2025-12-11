/*
 * \file
 */

#ifndef __PDM_PART_MESH_NODAL_H__
#define __PDM_PART_MESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_comm_graph.h"

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

typedef struct _pdm_part_mesh_nodal_t PDM_part_mesh_nodal_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a new \ref PDM_part_mesh_nodal_t instance
 *
 * \param [in]   mesh_dimension   Mesh dimension
 * \param [in]   n_part           Number of partitions on current process
 * \param [in]   comm             MPI communicator
 *
 * \return       Pointer to new \ref PDM_part_mesh_nodal_t instance
 *
 */
PDM_part_mesh_nodal_t*
PDM_part_mesh_nodal_create
(
  const int          mesh_dimension,
  const int          n_part,
  const PDM_MPI_Comm comm
);

/**
 * \brief Define partition vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  owner     Vertices ownership
 *
 */
void
PDM_part_mesh_nodal_coord_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const int                    n_vtx,
  const double                *coords,
        PDM_ownership_t        owner
);


/**
 * \brief Define vertices globals ids
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part   Partition identifier
 * \param [in]  numabs    Global IDs
 * \param [in]  owner     Vertices ownership
 *
 */
void
PDM_part_mesh_nodal_vtx_gnum_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const PDM_g_num_t           *numabs,
        PDM_ownership_t        owner
);

/**
 * \brief Define partition vertices from parents
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part       Partition identifier
 * \param [in]  n_vtx         Number of vertices
 * \param [in]  n_vtx_parent  Number of parent vertices
 * \param [in]  numabs        Global IDs (size = \ref n_vtx)
 * \param [in]  num_parent    Child to parent (local IDs, size = \ref n_vtx)
 * \param [in]  coords_parent Interlaced coordinates (size = 3 * \ref n_vtx_parent)
 * \param [in]  numabs_parent Global IDs of parent vertices (size = \ref n_vtx_parent)
 * \param [in]  owner         Vertices ownership
 *
 */
void
PDM_part_mesh_nodal_coord_from_parent_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const int                    n_vtx,
  const int                    n_vtx_parent,
  const PDM_g_num_t           *numabs,
  const int                   *num_parent,
  const PDM_real_t            *coords_parent,
  const PDM_g_num_t           *numabs_parent,
  const PDM_ownership_t        ownership
);

/**
 * \brief  Return number of partitions
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return  Number of partitions
 *
 */
int
PDM_part_mesh_nodal_n_part_get
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Return the mesh dimension
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return  Mesh dimension
 *
 */

int
PDM_part_mesh_nodal_mesh_dimension_get
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Return number of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */
int
PDM_part_mesh_nodal_n_vtx_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part
);

/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part   Partition identifier
 * \param [in]  ownership Ownership
 *
 * \return  Coordinates of vertices
 *
 */
double*
PDM_part_mesh_nodal_vtx_coord_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
        PDM_ownership_t        ownership
);

/**
 * \brief  Return global ids of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part   Partition identifier
 * \param [in]  ownership Ownership
 *
 * \return  Global ids of vertices
 *
 */
PDM_g_num_t*
PDM_part_mesh_nodal_vtx_g_num_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
        PDM_ownership_t        ownership
);

/**
 * \brief  Return number of sections in a specific geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Number of sections
 *
 */
int
PDM_part_mesh_nodal_n_section_in_geom_kind_get
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 * \brief  Return ids of sections in a specific geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Ids of sections
 *
 */
int *
PDM_part_mesh_nodal_sections_id_in_geom_kind_get
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 * \brief  Return type of section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section    Section identifier
 *
 * \return  Type of section
 *
 */
PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section
);

/**
 * \brief  Return type of section in a specific geometry kind
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section   Section identifier
 *
 * \return  Type of section
 *
 */
PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_in_geom_kind_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_section
);

/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  t_elt        Section type
 *
 * \return Section identifier
 *
 */
int
PDM_part_mesh_nodal_section_add
(
        PDM_part_mesh_nodal_t *pmn,
  const PDM_Mesh_nodal_elt_t   t_elt
);

/**
 * \brief Define a standard section
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  n_elt                   Number of elements
 * \param [in]  connec                  Connectivity
 * \param [in]  numabs                  Global numbering
 * \param [in]  parent_num              Parent numbering or NULL
 * \param [in]  parent_entity_g_num     Parent global numbering or NULL
 * \param [in]  owner                   Ownership
 *
 */
void
PDM_part_mesh_nodal_section_std_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    id_part,
  const int                    n_elt,
  const int                   *connec,
  const PDM_g_num_t           *numabs,
  const int                   *parent_num,
  const PDM_g_num_t           *parent_entity_g_num,
        PDM_ownership_t        owner
);

/**
 * \brief Define a standard high-order section
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  n_elt                   Number of elements
 * \param [in]  connec                  Connectivity
 * \param [in]  numabs                  Global numbering
 * \param [in]  parent_num              Parent numbering or NULL
 * \param [in]  parent_entity_g_num     Parent global numbering or NULL
 * \param [in]  order                   Element order
 * \param [in]  ho_ordering             HO ordering
 * \param [in]  owner                   Ownership
 *
 */
void
PDM_part_mesh_nodal_section_std_ho_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    id_part,
  const int                    n_elt,
  const int                   *connec,
  const PDM_g_num_t           *numabs,
  const int                   *parent_num,
  const PDM_g_num_t           *parent_entity_g_num,
  const int                    order,
  const char                  *ho_ordering,
        PDM_ownership_t        owner
);

/**
 * \brief Get number of elements in section
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section  Section identifier
 * \param [in]  id_part    Partition identifier
 *
 * \return      Number of elements in section
 *
 */
int
PDM_part_mesh_nodal_section_n_elt_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part
);

/**
 * \brief Return standard section description
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [out] connec                  Connectivity
 * \param [out] numabs                  Global ids
 * \param [out] parent_num              Parent local ids or NULL
 * \param [out] parent_entity_g_num     Parent global ids or NULL
 * \param [in]  ownership               Data ownership
 *
 */
void
PDM_part_mesh_nodal_section_std_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        int                   **connec,
        PDM_g_num_t           **numabs,
        int                   **parent_num,
        PDM_g_num_t           **parent_entity_g_num,
        PDM_ownership_t         ownership
);

/**
 * \brief Return standard high-order section description
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [out] connec                  Connectivity
 * \param [out] numabs                  Global ids
 * \param [out] parent_num              Parent local ids or NULL
 * \param [out] parent_entity_g_num     Parent global ids or NULL
 * \param [out] order                   Element order
 * \param [out] ho_ordering             HO ordering
 * \param [in]  ownership               Data ownership
 *
 */
void
PDM_part_mesh_nodal_section_std_ho_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        int                   **connec,
        PDM_g_num_t           **numabs,
        int                   **parent_num,
        PDM_g_num_t           **parent_entity_g_num,
        int                    *order,
  const char                  **ho_ordering,
        PDM_ownership_t         ownership
);

/**
 * \brief Get parent numbering of elements in section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return      Return parent numbering of elements in section
 *
 */
int *
PDM_part_mesh_nodal_section_parent_num_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        PDM_ownership_t         ownership
);

/**
 * \brief Get element global IDs of section elements
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return      Return element global IDs of section elements
 *
 */
PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        PDM_ownership_t         ownership
);

/**
 * \brief Add a \ref PDM_part_mesh_nodal_elmts_t to a \ref PDM_part_mesh_nodal_t
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  pmne         Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 *
 */
void
PDM_part_mesh_nodal_add_part_mesh_nodal_elmts
(
  PDM_part_mesh_nodal_t       *pmn,
  PDM_part_mesh_nodal_elmts_t *pmne
);

/**
 * \brief Free a \ref PDM_part_mesh_nodal_t instance
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */
void
PDM_part_mesh_nodal_free
(
  PDM_part_mesh_nodal_t* pmn
);


/**
 * \brief Export the current nodal mesh in vtk format with scalar fields attached to elements and vertices
 *
 * \param [in]  pmn              Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind        Geometry kind (corner, ridge, surface or volume)
 * \param [in]  filename_patter  Pattern for file naming (the function will append i_rank and i_part to this current pattern)
 * \param [in]  n_elt_field      Number of fields attached to elements
 * \param [in]  elt_field_name   Names of fields attached to elements
 * \param [in]  elt_field        Values of fields attached to elements (for each field, size = n_part)
 * \param [in]  n_vtx_field      Number of fields attached to vertices
 * \param [in]  vtx_field_name   Names of fields attached to vertices
 * \param [in]  vtx_field        Values of fields attached to vertices (for each field, size = n_part)
 *
 */
void
PDM_part_mesh_nodal_dump_vtk_with_fields
(
        PDM_part_mesh_nodal_t  *pmn,
        PDM_geometry_kind_t     geom_kind,
  const char                   *filename_pattern,
  const int                     n_elt_field,
  const char                   *elt_field_name[],
  const double                **elt_field     [],
  const int                     n_vtx_field,
  const char                   *vtx_field_name[],
  const double                **vtx_field     []
);

/**
 * \brief Export the current nodal mesh in vtk format
 *
 * \param [in]  pmn              Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind        Geometry kind (corner, ridge, surface or volume)
 * \param [in]  filename_patter  Pattern for file naming (the function will append i_rank and i_part to this current pattern)
 *
 */
void
PDM_part_mesh_nodal_dump_vtk
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind,
  const char            *filename_pattern
);

/**
 * \brief Compute element extents of a part of a section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  tolerance      Expansion tolerance for bounding boxes
 * \param [out] extents        Extents of mesh elements in current part of current section
 *
 */
void
PDM_part_mesh_nodal_section_elt_extents_compute
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    i_part,
  const double                 tolerance,
        double                *extents
);

/**
 * \brief Compute cell centers of a part of section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  ownership      Ownership
 *
 */
void
PDM_part_mesh_nodal_section_elt_center_compute
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    i_part,
  const PDM_ownership_t        ownership
);

/**
 * \brief  Return cell centers
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  ownership      Data ownership
 *
 * \return  Return cell centers
 *
 */
const double *
PDM_part_mesh_nodal_section_elt_center_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    i_part,
        PDM_ownership_t        ownership
);


/**
 * \brief Reset cell center computation
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 *
 */
void
PDM_part_mesh_nodal_section_elt_center_reset
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    i_part
);


/**
 * \brief Define a polygon section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connec_idx     Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connec         Connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 * \param [in]  owner          Ownership
 *
 */
void
PDM_part_mesh_nodal_section_poly2d_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    id_part,
  const int                    n_elt,
  const int                   *connec_idx,
  const int                   *connec,
  const PDM_g_num_t           *numabs,
  const int                   *parent_num,
        PDM_ownership_t        owner
);

/**
 * \brief Return a polygon section description
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connec_idx     Connectivity index
 * \param [out] connec         Connectivity
 * \param [in]  ownership      Data ownership
 *
 */
void
PDM_part_mesh_nodal_section_poly2d_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        int                   **connec_idx,
        int                   **connec,
        PDM_ownership_t         ownership
);


/**
 * \brief Define a polyhedron section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  n_face         Number of faces
 * \param [in]  facvtx_idx     Face->vertex connectivity index (size = \ref n_face + 1)
 * \param [in]  facvtx         Face->vertex connectivity
 * \param [in]  face_ln_to_gn  Face global ids
 * \param [in]  cellfac_idx    Cell->face connectivity index (size = \ref n_cell + 1)
 * \param [in]  cellfac        Cell->face connectivity
 * \param [in]  numabs         Cell global ids
 * \param [in]  parent_num     Cell parent numbering or NULL
 * \param [in]  owner          Ownership
 *
 */
void
PDM_part_mesh_nodal_section_poly3d_set
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    id_part,
  const int                    n_elt,
  const int                    n_face,
  const int                   *facvtx_idx,
  const int                   *facvtx,
  const PDM_g_num_t           *face_ln_to_gn,
  const int                   *cellfac_idx,
  const int                   *cellfac,
  const PDM_g_num_t           *numabs,
  const int                   *parent_num,
  const PDM_g_num_t           *parent_entity_g_num,
        PDM_ownership_t        owner
);

/**
 * \brief Return a polyhedron section
 *
 * \param [in]  pmn                  Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section            Section identifier
 * \param [in]  id_part              Partition identifier
 * \param [out] n_face               Number of faces
 * \param [out] face_ln_to_gn        Face global ids
 * \param [out] facvtx_idx           Face->vertex connectivity index (size = \ref n_face + 1)
 * \param [out] facvtx               Face->vertex connectivity
 * \param [out] numabs               Cell global ids
 * \param [out] cell_face_idx        Cell->face connectivity index (size = \ref n_cell + 1)
 * \param [out] cell_face            Cell->face connectivity
 * \param [out] parent_num           Cell parent numbering or NULL
 * \param [out] parent_entity_g_num  Cell parent global ids or NULL
 * \param [in]  ownership            Data ownership
 *
 */
void
PDM_part_mesh_nodal_section_poly3d_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        int                    *n_face,
        PDM_g_num_t           **face_ln_to_gn,
        int                   **face_vtx_idx,
        int                   **face_vtx,
        PDM_g_num_t           **numabs,
        int                   **cell_face_idx,
        int                   **cell_face,
        int                   **parent_num,
        PDM_g_num_t           **parent_entity_g_num,
        PDM_ownership_t         ownership
);

/**
 * \brief Get the cell->vertex connectivity of a polyhedron section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [out] cellvtx_idx  Index of cell->vertex connectivity
 * \param [out] cellvtx      Cell->vertex connectivity
 * \param [in]  ownership    Data ownership
 *
 */
void
PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
(
        PDM_part_mesh_nodal_t  *pmn,
  const int                     i_section,
  const int                     id_part,
        int                   **cellvtx_idx,
        int                   **cellvtx,
        PDM_ownership_t         ownership
);

/**
 * \brief Reset a nodal mesh structure
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return      NULL
 *
 */
void
PDM_part_mesh_nodal_reset
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Compute global IDs in a section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section    Section identifier
 * \param [in]  ownership    Ownership
 *
 */
void
PDM_part_mesh_nodal_g_num_in_section_compute
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const PDM_ownership_t        ownership
);

/**
 * \brief  Return number of elements in a partition
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Return number of elements in a partition
 *
 */
int
PDM_part_mesh_nodal_n_elmts_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_part
);

/**
 * \brief Get the element global numbering taking into account parent_num
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return  Global ids of element in current partition
 *
 */
PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get_from_part
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_part,
        PDM_ownership_t        ownership
);

/**
 * \brief Free partially a \ref PDM_part_mesh_nodal_t instance
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return      NULL
 *
 */
void
PDM_part_mesh_nodal_partial_free
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return true if the vertices are defined from parents
 */
int
PDM_part_mesh_nodal_is_set_coord_from_parent
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief Get global element numbering of section elements inside the section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return      Return global numbering of section elements inside the section
 *
 */
PDM_g_num_t *
PDM_part_mesh_nodal_section_g_num_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    id_part,
        PDM_ownership_t        ownership
);

/**
 * \brief  Return parent element number to local number
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent element number to local number
 *
 */
int *
PDM_part_mesh_nodal_num_elmt_parent_to_local_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_part
);


/**
 * \brief Return element to entity indirection for a section
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  ownership               Data ownership
 *
 */
int *
PDM_part_mesh_nodal_section_elmt_to_entity_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
  const int                    id_part,
        PDM_ownership_t        ownership
);

/**
 * \brief  Return parent num of vertices
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent of vertices
 *
 */
const int *
PDM_part_mesh_nodal_vertices_parent_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part
);

/**
 * \brief  Return parent global IDs
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent of vertices
 *
 */
const PDM_g_num_t *
PDM_part_mesh_nodal_vertices_g_num_parent_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part
);

/**
 * \brief  Add some 3D cells from cell face connectivity.
 *
 * For each cell, this function determines the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding section.
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_cell         Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  face_vtx_idx   Index of face vertex connectivity
 * \param [in]  face_vtx       Face vertex connectivity
 * \param [in]  face_ln_to_gn  Face global numbering
 * \param [in]  cell_face_idx  Index of cell face connectivity
 * \param [in]  cell_face      Cell face connectivity
 * \param [in]  cell_ln_to_gn  Global numbering
 * \param [in]  ownership      Ownership
 *
 */
void
PDM_part_mesh_nodal_cell3d_cellface_add
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const int                    n_cell,
  const int                    n_face,
  const int                   *face_vtx_idx,
  const int                   *face_vtx,
  const PDM_g_num_t           *face_ln_to_gn,
  const int                   *cell_face_idx,
  const int                   *cell_face,
  const PDM_g_num_t           *cell_ln_to_gn,
  const PDM_ownership_t        ownership
);

/**
 * \brief  Add some 2D faces from face edge connectivity.
 *
 * For each face, this function determines the type of the face (triangles, quadrangles, ...)
 * and stores it in the corresponding section.
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygons
 * \param [in]  n_edge         Number of edges used to describe polygons
 * \param [in]  edge_vtx       edge vertex connectivity
 * \param [in]  face_edge_idx  Index of face edge connectivity
 * \param [in]  face_edge      face edge connectivity
 * \param [in]  face_ln_to_gn  Global numbering
 * \param [in]  ownership      Ownership
 *
 */
void
PDM_part_mesh_nodal_face2d_faceedge_add
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const int                    n_face,
  const int                    n_edge,
  const int                   *edge_vtx,
  const int                   *face_edge_idx,
  const int                   *face_edge,
  const PDM_g_num_t           *face_ln_to_gn,
  const PDM_ownership_t        ownership
);

/**
 * \brief  Add some standard 3D cells from cell vertex connectivity.
 *
 * For each cell, this function determines the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding section.
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts instance
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_cell         Number of cells
 * \param [in]  cell_vtx_idx   Index of cell vertex connectivity
 * \param [in]  cell_vtx       Cell vertex connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  ownership      Ownership
 *
 */
void
PDM_part_mesh_nodal_cells_cellvtx_add
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const int                    n_cell,
  const int                   *cell_vtx_idx,
  const int                   *cell_vtx,
  const PDM_g_num_t           *numabs,
  const PDM_ownership_t        ownership
);

/**
 * \brief  Add some 2D faces from face vertex connectivity.
 *
 * For each face, this function determines the type of the cell (triangles, quadrangles, ...)
 * and stores it in the corresponding section.
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts instance
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygons
 * \param [in]  face_vtx_idx   Index of edge vertex connectivity
 * \param [in]  face_vtx       Edge vertex connectivity
 * \param [in]  ownership      Ownership
 *
 */
void
PDM_part_mesh_nodal_faces_facevtx_add
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    id_part,
  const int                    n_face,
  const int                   *face_vtx_idx,
  const int                   *face_vtx,
  const PDM_g_num_t           *numabs,
  const PDM_ownership_t        ownership
);


/**
 * \brief  Return geom_kind and identifier (local to this geom_kind) of a section
 *
 * \param [in]  pmn                      Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  i_section                Unique section identifier
 * \param [out] geom_kind                Geometry kind (corner, ridge, surface or volume)
 * \param [out] id_section_in_geom_kind  Section identifier local to the geometry kind
 *
 */
void
PDM_part_mesh_nodal_section_id_and_geom_kind_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section,
        PDM_geometry_kind_t   *geom_kind,
        int                   *id_section_in_geom_kind
);

/**
 * \brief  Return unique identifier of a section
 *
 * \param [in]  pmn                      Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind                Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section_in_geom_kind  Section identifier local to the geometry kind
 *
 * \return   Unique section identifier
 *
 */
int
PDM_part_mesh_nodal_section_id_from_geom_kind_get
(
        PDM_part_mesh_nodal_t *pmn,
  const PDM_geometry_kind_t    geom_kind,
  const int                    id_section_in_geom_kind
);

/**
 * \brief  Return number of sections
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return  Number of sections
 *
 */
int
PDM_part_mesh_nodal_n_section_get
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Return ids of sections
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return  Ids of sections
 *
 */
int *
PDM_part_mesh_nodal_sections_id_get
(
  PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Set number of group for a current geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 * \param [in]  n_group    Number of group in geom_kind
 */
void
PDM_part_mesh_nodal_n_group_set
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    n_group
);

/**
 * \brief  Set partition group
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  i_part         Partition identifier
 * \param [in]  i_group        Group identifier
 * \param [in]  n_group_elmt   Number of element in current group for current part
 * \param [in]  group_elmt     List of entity in group (size = \p n_group_elmt)
 * \param [in]  group_ln_to_gn List of global IDs in group (size = \p n_group_elmt)
 * \param [in]  ownership      Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
void
PDM_part_mesh_nodal_group_set
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    i_part,
  const int                    i_group,
        int                    n_group_elmt,
        int                   *group_elmt,
        PDM_g_num_t           *group_ln_to_gn,
        PDM_ownership_t        ownership
);

/**
 * \brief  Get partition group
 *
 * \param [in]   pmn            Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]   geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]   i_part         Partition identifier
 * \param [in]   i_group        Group identifier
 * \param [out]  n_group_elmt   Number of element in current group for current part
 * \param [out]  group_elmt     List of entity in group (size = \p n_group_elmt)
 * \param [out]  group_ln_to_gn List of global IDs in group (size = \p n_group_elmt)
 * \param [in]   ownership      Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
void
PDM_part_mesh_nodal_group_get
(
        PDM_part_mesh_nodal_t  *pmn,
        PDM_geometry_kind_t     geom_kind,
  const int                     i_part,
  const int                     i_group,
        int                    *n_group_elmt,
        int                   **group_elmt,
        PDM_g_num_t           **group_ln_to_gn,
        PDM_ownership_t         ownership
);

/**
 * \brief Get the section index for a current kind for a partition
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 * \param [in]  i_part     Partition identifier
 *
 * \return Index of sections (size=n_sections)
 */
int*
PDM_part_mesh_nodal_compute_sections_idx
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind,
  const int              id_part
);

/**
 * \brief Get number of group for a current geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Number of group in geom_kind
 */
int
PDM_part_mesh_nodal_n_group_get
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 * \brief Get the substructure \ref PDM_part_mesh_nodal_elmts_t for a current geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  \ref PDM_part_mesh_nodal_elmts_t of geometry kind
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_part_mesh_nodal_elmts_get
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
);

/**
 * \brief Get the substructure \ref PDM_part_mesh_nodal_elmts_t at principal dimension of the current \ref PDM_part_mesh_nodal_t
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return  \ref PDM_part_mesh_nodal_elmts_t of principal dimension
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_part_mesh_nodal_elmts_principal_dim_get
(
  PDM_part_mesh_nodal_t  *pmn
);

/**
 * \brief Return the geometry kind of highest dimension
 * for a given \ref PDM_part_mesh_nodal_t instance
 *
 * \param [in] pmn    Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return Geometry kind of highest dimension
 *
 */
PDM_geometry_kind_t
PDM_part_mesh_nodal_principal_geom_kind_get
(
  PDM_part_mesh_nodal_t *pmn
);


/**
 * \brief Return the cell->vertex connectivity
 * The output pointers are owned by the user.
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  geom_kind     Geometry kind (corner, ridge, surface or volume)
 * \param [in]  i_part        Partition identifier
 * \param [out] cell_vtx_idx  Index for the cell->vertex connectivity
 * \param [out] cell_vtx      Cell->vertex connectivity
 *
 * \return Number of cells in current partition
 *
 */
int
PDM_part_mesh_nodal_cell_vtx_connect_get
(
        PDM_part_mesh_nodal_t  *pmn,
        PDM_geometry_kind_t     geom_kind,
  const int                     i_part,
        int                   **cell_vtx_idx,
        int                   **cell_vtx
);

/**
 *
 * \brief Set part_comm_graph onto part_mesh_nodal struct
 *
 * \param [in]  pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  pcg         Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  geom_kind   Geometry kind (see \ref PDM_geometry_kind_t )
 * \param [in]  ownership   part_mesh_nodal ownership on given part_comm_graph
 *
 */
void
PDM_part_mesh_nodal_part_comm_graph_set
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_part_comm_graph_t *pcg,
  PDM_geometry_kind_t    geom_kind,
  PDM_ownership_t        ownership
);


/**
 *
 * \brief Get part_mesh_nodal's part_comm_graph
 *
 * \param [in]   pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]   geom_kind   Geometry kind (see \ref PDM_geometry_kind_t )
 * \param [out]  pcg         Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]   ownership   part_mesh_nodal ownership on returned part_comm_graph
 */
void
PDM_part_mesh_nodal_part_comm_graph_get
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  PDM_part_comm_graph_t **pcg,
  PDM_ownership_t         ownership
);

/**
 *
 * \brief Set part_comm_graph onto part_mesh_nodal struct for vertices
 *
 * \param [in]  pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [in]  pcg         Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]  ownership   part_mesh_nodal ownership on given part_comm_graph
 *
 */
void
PDM_part_mesh_nodal_part_comm_graph_vtx_set
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_part_comm_graph_t *pcg,
  PDM_ownership_t        ownership
);

/**
 *
 * \brief Get part_mesh_nodal's part_comm_graph for vertices
 *
 * \param [in]   pmn         Pointer to \ref PDM_part_mesh_nodal_t instance
 * \param [out]  pcg         Pointer to \ref PDM_part_comm_graph_t instance
 * \param [in]   ownership   part_mesh_nodal ownership on returned part_comm_graph
 */
void
PDM_part_mesh_nodal_part_comm_graph_vtx_get
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_part_comm_graph_t **pcg,
  PDM_ownership_t         ownership
);

/**
 * \brief Transform group information inside a \ref PDM_part_mesh_nodal_t to tag for all elements
 *
 * \param [in]  pmne            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind       Geometry kind (corner, ridge, surface or volume)
 * \param [in]  allow_multiple  If PDM_TRUE allows that one element can be referenced by more than one group or not referenced at all
 * \param [out] out_tag_idx     Identifier index if allow_multiple is PDM_TRUE, else NULL
 * \param [out] out_tag         Identifier for all elements in current \ref PDM_part_mesh_nodal_t
 *                              that follows the natural order of the elements (size = n_part)
 *                              For each part size is : n_elmts or out_tag_idx[n_elmts]
 *                              Value is between [0, n_group-1]
 */
void
PDM_part_mesh_nodal_group_to_tag
(
  PDM_part_mesh_nodal_t   *pmn,
  PDM_geometry_kind_t      geom_kind,
  PDM_bool_t               allow_multiple,
  int                   ***out_tag_idx,
  int                   ***out_tag
);

/**
 * \brief Transform tag for all elements into group information inside a \ref PDM_part_mesh_nodal_t
 *
 * \param [in]  pmne      Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind Geometry kind (corner, ridge, surface or volume)
 * \param [in]  n_group   Number of groups, if the specify n_group is negative, n_group is automatically compute
 * \param [in]  tag_idx   Identifier index or NULL if no multiplicity in element tag
 * \param [in]  tag       Identifier for all elements in current \ref PDM_part_mesh_nodal_elmts_t
 *                        that follows the natural order of the elements (size = n_part)
 *                        For each part size is : n_elmts or out_tag_idx[n_elmts]
 *                        Value is between [0, n_group-1]
 */
void
PDM_part_mesh_nodal_tag_to_group
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  int                     n_group,
  int                   **tag_idx,
  int                   **tag
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_H__ */
