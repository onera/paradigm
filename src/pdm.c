
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_error.h"
#include "pdm_part.h"
#include "pdm_part_renum.h"
#include "pdm_part_coarse_mesh.h"


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

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_entities_t
PDM_connectivity_type_to_entity_type
(
 PDM_connectivity_type_t   connectivity_type
)
{
  PDM_mesh_entities_t mesh_entity = PDM_MESH_ENTITY_MAX;
  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    mesh_entity = PDM_MESH_ENTITY_CELL;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    mesh_entity = PDM_MESH_ENTITY_FACE;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    mesh_entity = PDM_MESH_ENTITY_EDGE;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    mesh_entity = PDM_MESH_ENTITY_VTX;
  }
  return mesh_entity;
}


PDM_mesh_entities_t
PDM_bound_type_to_entity_type
(
 PDM_bound_type_t   bound_type
)
{
  switch (bound_type) {

    case PDM_BOUND_TYPE_CELL: {
      return PDM_MESH_ENTITY_CELL;
    }
    case PDM_BOUND_TYPE_FACE: {
      return PDM_MESH_ENTITY_FACE;
    }
    case PDM_BOUND_TYPE_EDGE: {
      return PDM_MESH_ENTITY_EDGE;
    }
    case PDM_BOUND_TYPE_VTX: {
      return PDM_MESH_ENTITY_VTX;
    }
    default: {
      PDM_error(__FILE__, __LINE__, 0, "Invalid bound_type %d\n", bound_type);
    }

  }

  return PDM_MESH_ENTITY_MAX;
}


PDM_mesh_entities_t
PDM_geometry_kind_to_entity_type
(
 PDM_geometry_kind_t   geom_kind
)
{
  switch (geom_kind) {

    case PDM_GEOMETRY_KIND_VOLUMIC: {
      return PDM_MESH_ENTITY_CELL;
    }
    case PDM_GEOMETRY_KIND_SURFACIC: {
      return PDM_MESH_ENTITY_FACE;
    }
    case PDM_GEOMETRY_KIND_RIDGE: {
      return PDM_MESH_ENTITY_EDGE;
    }
    case PDM_GEOMETRY_KIND_CORNER: {
      return PDM_MESH_ENTITY_VTX;
    }
    default: {
      PDM_error(__FILE__, __LINE__, 0, "Invalid geom_kind %d\n", geom_kind);
    }

  }

  return PDM_MESH_ENTITY_MAX;
}


PDM_geometry_kind_t
PDM_entity_type_to_geometry_kind
(
 PDM_mesh_entities_t   entity_type
)
{
  switch (entity_type) {

    case PDM_MESH_ENTITY_CELL: {
      return PDM_GEOMETRY_KIND_VOLUMIC;
    }
    case PDM_MESH_ENTITY_FACE: {
      return PDM_GEOMETRY_KIND_SURFACIC;
    }
    case PDM_MESH_ENTITY_EDGE: {
      return PDM_GEOMETRY_KIND_RIDGE;
    }
    case PDM_MESH_ENTITY_VTX: {
      return PDM_GEOMETRY_KIND_CORNER;
    }
    default: {
      PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
    }

  }

  return PDM_GEOMETRY_KIND_MAX;
}


PDM_connectivity_type_t
PDM_entity_pair_to_connectivity_type
(
  PDM_mesh_entities_t entity_type1,
  PDM_mesh_entities_t entity_type2
)
{
  if (entity_type1 == PDM_MESH_ENTITY_CELL) {
    switch (entity_type2) {
      case PDM_MESH_ENTITY_CELL: return PDM_CONNECTIVITY_TYPE_CELL_CELL;
      case PDM_MESH_ENTITY_FACE: return PDM_CONNECTIVITY_TYPE_CELL_FACE;
      case PDM_MESH_ENTITY_EDGE: return PDM_CONNECTIVITY_TYPE_CELL_EDGE;
      case PDM_MESH_ENTITY_VTX : return PDM_CONNECTIVITY_TYPE_CELL_VTX;
      default: PDM_error(__FILE__, __LINE__, 0, "Invalid entity pair %d %d\n", entity_type1, entity_type2);
    }
  }
  else if (entity_type1 == PDM_MESH_ENTITY_FACE) {
    switch (entity_type2) {
      case PDM_MESH_ENTITY_CELL: return PDM_CONNECTIVITY_TYPE_FACE_CELL;
      case PDM_MESH_ENTITY_FACE: return PDM_CONNECTIVITY_TYPE_FACE_FACE;
      case PDM_MESH_ENTITY_EDGE: return PDM_CONNECTIVITY_TYPE_FACE_EDGE;
      case PDM_MESH_ENTITY_VTX : return PDM_CONNECTIVITY_TYPE_FACE_VTX;
      default: PDM_error(__FILE__, __LINE__, 0, "Invalid entity pair %d %d\n", entity_type1, entity_type2);
    }
  }
  else if (entity_type1 == PDM_MESH_ENTITY_EDGE) {
    switch (entity_type2) {
      case PDM_MESH_ENTITY_CELL: return PDM_CONNECTIVITY_TYPE_EDGE_CELL;
      case PDM_MESH_ENTITY_FACE: return PDM_CONNECTIVITY_TYPE_EDGE_FACE;
      case PDM_MESH_ENTITY_EDGE: return PDM_CONNECTIVITY_TYPE_EDGE_EDGE;
      case PDM_MESH_ENTITY_VTX : return PDM_CONNECTIVITY_TYPE_EDGE_VTX;
      default: PDM_error(__FILE__, __LINE__, 0, "Invalid entity pair %d %d\n", entity_type1, entity_type2);
    }
  }
  else if (entity_type1 == PDM_MESH_ENTITY_VTX) {
    switch (entity_type2) {
      case PDM_MESH_ENTITY_CELL: return PDM_CONNECTIVITY_TYPE_VTX_CELL;
      case PDM_MESH_ENTITY_FACE: return PDM_CONNECTIVITY_TYPE_VTX_FACE;
      case PDM_MESH_ENTITY_EDGE: return PDM_CONNECTIVITY_TYPE_VTX_EDGE;
      case PDM_MESH_ENTITY_VTX : return PDM_CONNECTIVITY_TYPE_VTX_VTX;
      default: PDM_error(__FILE__, __LINE__, 0, "Invalid entity pair %d %d\n", entity_type1, entity_type2);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid entity pair %d %d\n", entity_type1, entity_type2);
  }

  return PDM_CONNECTIVITY_TYPE_MAX;
}


int
PDM_connectivity_type_to_entity_pair
(
  PDM_connectivity_type_t  connectivity_type,
  PDM_mesh_entities_t     *entity_type1,
  PDM_mesh_entities_t     *entity_type2
)
{
  switch (connectivity_type) {
    case PDM_CONNECTIVITY_TYPE_CELL_CELL: {
      *entity_type1 = PDM_MESH_ENTITY_CELL;
      *entity_type2 = PDM_MESH_ENTITY_CELL;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_CELL_FACE: {
      *entity_type1 = PDM_MESH_ENTITY_CELL;
      *entity_type2 = PDM_MESH_ENTITY_FACE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_CELL_EDGE: {
      *entity_type1 = PDM_MESH_ENTITY_CELL;
      *entity_type2 = PDM_MESH_ENTITY_EDGE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_CELL_VTX: {
      *entity_type1 = PDM_MESH_ENTITY_CELL;
      *entity_type2 = PDM_MESH_ENTITY_VTX;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_CELL: {
      *entity_type1 = PDM_MESH_ENTITY_FACE;
      *entity_type2 = PDM_MESH_ENTITY_CELL;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_FACE_FACE: {
      *entity_type1 = PDM_MESH_ENTITY_FACE;
      *entity_type2 = PDM_MESH_ENTITY_FACE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_FACE_EDGE: {
      *entity_type1 = PDM_MESH_ENTITY_FACE;
      *entity_type2 = PDM_MESH_ENTITY_EDGE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_FACE_VTX: {
      *entity_type1 = PDM_MESH_ENTITY_FACE;
      *entity_type2 = PDM_MESH_ENTITY_VTX;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_EDGE_CELL: {
      *entity_type1 = PDM_MESH_ENTITY_EDGE;
      *entity_type2 = PDM_MESH_ENTITY_CELL;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_EDGE_FACE: {
      *entity_type1 = PDM_MESH_ENTITY_EDGE;
      *entity_type2 = PDM_MESH_ENTITY_FACE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_EDGE_EDGE: {
      *entity_type1 = PDM_MESH_ENTITY_EDGE;
      *entity_type2 = PDM_MESH_ENTITY_EDGE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_EDGE_VTX: {
      *entity_type1 = PDM_MESH_ENTITY_EDGE;
      *entity_type2 = PDM_MESH_ENTITY_VTX;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_VTX_CELL: {
      *entity_type1 = PDM_MESH_ENTITY_VTX;
      *entity_type2 = PDM_MESH_ENTITY_CELL;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_VTX_FACE: {
      *entity_type1 = PDM_MESH_ENTITY_VTX;
      *entity_type2 = PDM_MESH_ENTITY_FACE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_VTX_EDGE: {
      *entity_type1 = PDM_MESH_ENTITY_VTX;
      *entity_type2 = PDM_MESH_ENTITY_EDGE;
      break;
    }
    case PDM_CONNECTIVITY_TYPE_VTX_VTX: {
      *entity_type1 = PDM_MESH_ENTITY_VTX;
      *entity_type2 = PDM_MESH_ENTITY_VTX;
      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid connectivity type %d\n", connectivity_type);
      // return 0;
    }
  }

  return 1;
}


/**
 * \brief Finalize PDM
 *
 * This function frees all allocated global variables
 *
 */

void
PDM_Finalize
(
void
)
{

 /**
  *  Free global array inside part_renum
  */

 PDM_part_renum_method_purge();
 PDM_coarse_mesh_method_purge();

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
