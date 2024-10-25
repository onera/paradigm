/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_ho_ordering.h"

#include "pdm_writer.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Tables for decomposing standard elments in edges and faces */
// https://cgns.github.io/CGNS_docs_current/sids/conv.html

static const int bar_edge_vtx[] = {
  0, 1
};

static const int tria_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 0
};

static const int quad_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0
};

static const int tetra_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 0,
  0, 3,
  1, 3,
  2, 3
};

static const int pyramid_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  0, 4,
  1, 4,
  2, 4,
  3, 4
};

static const int prism_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 0,
  0, 3,
  1, 4,
  2, 5,
  3, 4,
  4, 5,
  5, 3
};

static const int hexa_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  0, 4,
  1, 5,
  2, 6,
  3, 7,
  4, 5,
  5, 6,
  6, 7,
  7, 4
};


static const int tria_face_vtx_idx[] = {
  0, 3
};

static const int tria_face_vtx[] = {
  0, 1, 2
};

static const int quad_face_vtx_idx[] = {
  0, 4
};

static const int quad_face_vtx[] = {
  0, 1, 2, 3
};

static const int tetra_face_vtx_idx[] = {
  0, 3, 6, 9, 12
};

static const int tetra_face_vtx[] = {
  0, 2, 1,
  0, 1, 3,
  0, 3, 2,
  1, 2, 3
};

static const int pyramid_face_vtx_idx[] = {
  0, 4, 7, 10, 13, 16
};

static const int pyramid_face_vtx[] = {
  3, 2, 1, 0,
  4, 0, 1,
  4, 1, 2,
  4, 2, 3,
  4, 3, 0
};

static const int prism_face_vtx_idx[] = {
  0, 3, 6, 10, 14, 18
};

static const int prism_face_vtx[] = {
  2, 1, 0,
  4, 5, 3,
  5, 4, 1, 2,
  4, 3, 0, 1,
  3, 5, 2, 0
};

static const int hexa_face_vtx_idx[] = {
  0, 4, 8, 12, 16, 20, 24
};

static const int hexa_face_vtx[] = {
  3, 2, 1, 0,
  6, 7, 4, 5,
  4, 7, 3, 0,
  7, 6, 2, 3,
  2, 6, 5, 1,
  1, 5, 4, 0
};

static inline int
ij2idx_tria
(
 const int i,
 const int j,
 const int order
 )
{
    return i + j*(order+1) - (j-1)*j/2;
}

static inline int
ij2idx_quad
(
 const int i,
 const int j,
 const int order
 )
{
    return i + j*(order+1);
}

static inline int
ijk2idx_tetra
(
 const int i,
 const int j,
 const int k,
 const int order
 )
{
    return i + j*(order + 1 - k) - j*(j-1)/2 + (k*(k*(k - 3*order - 6) + 3*order*(order + 4) + 11)) / 6;
}

static inline int
ijk2idx_pyramid
(
 const int i,
 const int j,
 const int k,
 const int order
 )
{
    return i + j*(order+1-k) + (k*(k*(2*k - 6*order - 9) + 6*order*(order + 3) + 13)) / 6;
}

static inline int
ijk2idx_prism
(
 const int i,
 const int j,
 const int k,
 const int order
 )
{
    return i + j*(order+1) - j*(j-1)/2 + k*(order+1)*(order+2)/2;
}

static inline int
ijk2idx_hexa
(
 const int i,
 const int j,
 const int k,
 const int order
 )
{
    return i + (order+1)*(j + (order+1)*k);
}

static void
_principal_to_ijk
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 int                        *principal
 )
{
  if (t_elt == PDM_MESH_NODAL_BARHO ||
      t_elt == PDM_MESH_NODAL_BARHO_BEZIER) {
    principal[0] = 0;
    principal[1] = order;
  }

  else if (t_elt == PDM_MESH_NODAL_TRIAHO ||
           t_elt == PDM_MESH_NODAL_TRIAHO_BEZIER) {
    principal[0] = ij2idx_tria(0,     0,     order);
    principal[1] = ij2idx_tria(order, 0,     order);
    principal[2] = ij2idx_tria(0,     order, order);
  }

  else if (t_elt == PDM_MESH_NODAL_QUADHO) {
    principal[0] = ij2idx_quad(0,     0,     order);
    principal[1] = ij2idx_quad(order, 0,     order);
    principal[2] = ij2idx_quad(order, order, order);
    principal[3] = ij2idx_quad(0,     order, order);
  }

  else if (t_elt == PDM_MESH_NODAL_TETRAHO) {
    principal[0] = ijk2idx_tetra(0,     0,     0,     order);
    principal[1] = ijk2idx_tetra(order, 0,     0,     order);
    principal[2] = ijk2idx_tetra(0,     order, 0,     order);
    principal[3] = ijk2idx_tetra(0,     0,     order, order);
  }

  else if (t_elt == PDM_MESH_NODAL_PYRAMIDHO) {
    principal[0] = ijk2idx_pyramid(0,     0,     0,     order);
    principal[1] = ijk2idx_pyramid(order, 0,     0,     order);
    principal[2] = ijk2idx_pyramid(order, order, 0,     order);
    principal[3] = ijk2idx_pyramid(0,     order, 0,     order);
    principal[4] = ijk2idx_pyramid(0,     0,     order, order);
  }

  else if (t_elt == PDM_MESH_NODAL_PRISMHO) {
    principal[0] = ijk2idx_prism(0,     0,     0,     order);
    principal[1] = ijk2idx_prism(order, 0,     0,     order);
    principal[2] = ijk2idx_prism(0,     order, 0,     order);
    principal[3] = ijk2idx_prism(0,     0,     order, order);
    principal[4] = ijk2idx_prism(order, 0,     order, order);
    principal[5] = ijk2idx_prism(0,     order, order, order);
  }

  else if (t_elt == PDM_MESH_NODAL_HEXAHO) {
    principal[0] = ijk2idx_hexa(0,     0,     0,     order);
    principal[1] = ijk2idx_hexa(order, 0,     0,     order);
    principal[2] = ijk2idx_hexa(order, order, 0,     order);
    principal[3] = ijk2idx_hexa(0,     order, 0,     order);
    principal[4] = ijk2idx_hexa(0,     0,     order, order);
    principal[5] = ijk2idx_hexa(order, 0,     order, order);
    principal[6] = ijk2idx_hexa(order, order, order, order);
    principal[7] = ijk2idx_hexa(0,     order, order, order);
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid element type %d\n", (int) t_elt);
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Get the dimension of an element
 *
 * \param[in]  type    Element type
 *
 * \return     Dimension of the element
 *
 */

int
PDM_Mesh_nodal_elt_dim_get
(
 PDM_Mesh_nodal_elt_t type
 )
{
  int elt_dim = -1;

  switch(type) {
    case PDM_MESH_NODAL_POINT:
      elt_dim = 0;
      break;
    case PDM_MESH_NODAL_BAR2:
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER:
      elt_dim = 1;
      break;
    case PDM_MESH_NODAL_TRIA3:
    case PDM_MESH_NODAL_QUAD4:
    case PDM_MESH_NODAL_POLY_2D:
    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_QUADHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
      elt_dim = 2;
      break;
    case PDM_MESH_NODAL_TETRA4:
    case PDM_MESH_NODAL_PYRAMID5:
    case PDM_MESH_NODAL_PRISM6:
    case PDM_MESH_NODAL_HEXA8:
    case PDM_MESH_NODAL_POLY_3D:
    case PDM_MESH_NODAL_TETRAHO:
    case PDM_MESH_NODAL_PYRAMIDHO:
    case PDM_MESH_NODAL_PRISMHO:
    case PDM_MESH_NODAL_HEXAHO:
      elt_dim = 3;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt %d\n", (int) type);
  }

  return elt_dim;
}


/**
 *
 * \brief Check if an element is two-dimensional
 *
 * \param [in]  type     Element type
 *
 * \return    1 if the element is 2D, 0 else
 *
 */

int
PDM_Mesh_nodal_is_2D_element
(
  PDM_Mesh_nodal_elt_t type
)
{
  return (PDM_Mesh_nodal_elt_dim_get(type) == 2);
}


/**
 *
 * \brief Check if an element is three-dimensional
 *
 * \param [in]  type     Element type
 *
 * \return    1 if the element is 3D, 0 else
 *
 */

int
PDM_Mesh_nodal_is_3D_element
(
  PDM_Mesh_nodal_elt_t type
)
{
  return (PDM_Mesh_nodal_elt_dim_get(type) == 3);
}


/**
 * \brief Get the number of vertices of an element type
 *
 * \param [in]   type     Element type
 * \param [in]   order    Element order
 *
 * \return       Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vtx_elt_get
(
  PDM_Mesh_nodal_elt_t type,
  const int            order
)
{
  switch (type) {
  case PDM_MESH_NODAL_POINT :
    return 1;
    break;
  case PDM_MESH_NODAL_BAR2 :
    return 2;
    break;
  case PDM_MESH_NODAL_TRIA3 :
    return 3;
    break;
  case PDM_MESH_NODAL_QUAD4 :
    return 4;
    break;
  case PDM_MESH_NODAL_TETRA4 :
    return 4;
    break;
  case PDM_MESH_NODAL_PYRAMID5 :
    return 5;
    break;
  case PDM_MESH_NODAL_PRISM6 :
    return 6;
    break;
  case PDM_MESH_NODAL_HEXA8 :
    return 8;
    break;

  case PDM_MESH_NODAL_BARHO :
    return order + 1;
    break;
  case PDM_MESH_NODAL_TRIAHO :
    return (order + 1) * (order + 2) / 2;
    break;
  case PDM_MESH_NODAL_QUADHO :
    return (order + 1) * (order + 1);
    break;
  case PDM_MESH_NODAL_TETRAHO :
    return (order + 1) * (order + 2) * (order + 3) / 6;
    break;
  case PDM_MESH_NODAL_PYRAMIDHO :
    return (order + 1) * (order + 2) * (2*order + 3) / 6;
    break;
  case PDM_MESH_NODAL_PRISMHO :
    return (order + 1) * (order + 1) * (order + 2) / 2;
    break;
  case PDM_MESH_NODAL_HEXAHO :
    return (order + 1) * (order + 1) * (order + 1);
    break;

  case PDM_MESH_NODAL_BARHO_BEZIER:
    return order + 1;
    break;
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
   return (order + 1) * (order + 2) / 2;
    break;
  default :
    PDM_error (__FILE__, __LINE__, 0, "Unknown order for Poly2D and Poly3D (type %d)\n", type);
  }
  return -1;
}

int
PDM_Mesh_nodal_elmt_is_ho
(
  PDM_Mesh_nodal_elt_t type
)
{  switch (type) {
  case PDM_MESH_NODAL_POINT :
  case PDM_MESH_NODAL_BAR2 :
  case PDM_MESH_NODAL_TRIA3 :
  case PDM_MESH_NODAL_QUAD4 :
  case PDM_MESH_NODAL_TETRA4 :
  case PDM_MESH_NODAL_PYRAMID5 :
  case PDM_MESH_NODAL_PRISM6 :
  case PDM_MESH_NODAL_HEXA8 :
  case PDM_MESH_NODAL_POLY_2D :
  case PDM_MESH_NODAL_POLY_3D :
    return 0;
    break;

  case PDM_MESH_NODAL_BARHO :
  case PDM_MESH_NODAL_TRIAHO :
  case PDM_MESH_NODAL_QUADHO :
  case PDM_MESH_NODAL_TETRAHO :
  case PDM_MESH_NODAL_PYRAMIDHO :
  case PDM_MESH_NODAL_PRISMHO :
  case PDM_MESH_NODAL_HEXAHO :
  case PDM_MESH_NODAL_BARHO_BEZIER:
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
    return 1;
    break;
  default :
    PDM_error (__FILE__, __LINE__, 0, "Unknown elt type %d\n", (int) type);
  }
  return -1;
}


void
PDM_Mesh_nodal_ho_parent_node
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const char                 *ho_ordering,
       int                  *parent_node
 )
{
  int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

  int principal_to_ijk[8];
  _principal_to_ijk(t_elt,
                    order,
                    principal_to_ijk);

  if (ho_ordering == NULL) {
    for (int i = 0; i < elt_vtx_n; i++) {
      parent_node[i] = principal_to_ijk[i];
    }
  }
  else {
    int *ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering,
                                                       t_elt,
                                                       order);
    assert(ijk_to_user != NULL);

    for (int i = 0; i < elt_vtx_n; i++) {
      parent_node[i] = ijk_to_user[principal_to_ijk[i]];
    }
  }

}



void
PDM_Mesh_nodal_reorder_elt_vtx
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const char                 *ho_ordering_in,
 const char                 *ho_ordering_out,
 const int                   n_elt,
       int                  *elt_vtx_in,
       int                  *elt_vtx_out
 )
{
  int stride = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                           order);

  int *ijk_to_in = NULL;
  if (ho_ordering_in != NULL) {
    ijk_to_in = PDM_ho_ordering_ijk_to_user_get(ho_ordering_in,
                                                t_elt,
                                                order);
    assert(ijk_to_in != NULL);
  }

  int *ijk_to_out = NULL;
  if (ho_ordering_out != NULL) {
    ijk_to_out = PDM_ho_ordering_ijk_to_user_get(ho_ordering_out,
                                                 t_elt,
                                                 order);
    assert(ijk_to_out != NULL);
  }


  int *tmp = NULL;
  PDM_malloc(tmp, stride, int);
  for (int ielt = 0; ielt < n_elt; ielt++) {

    int *ev_in  = elt_vtx_in  + stride*ielt;
    int *ev_out = elt_vtx_out + stride*ielt;
    memcpy(tmp, ev_in, sizeof(int) * stride);

    /* In --> IJK */
    if (ijk_to_in != NULL) {
      for (int i = 0; i < stride; i++) {
        tmp[i] = ev_in[ijk_to_in[i]];
      }
    }

    /* IJK --> Out */
    if (ijk_to_out != NULL) {
      for (int i = 0; i < stride; i++) {
        ev_out[ijk_to_out[i]] = tmp[i];
      }
    }
    else {
      memcpy(ev_out, tmp, sizeof(int) * stride);
    }

  }

  // int *ijk_to_user = NULL;

  // /* In --> IJK */
  // if (ho_ordering_in != NULL){
  //   ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering_in,
  //                                                 t_elt,
  //                                                 order);
  //   assert(ijk_to_user != NULL);

  //   for (int ielt = 0; ielt < n_elt; ielt++) {

  //     int *ev_in  = elt_vtx_in  + ielt*stride;
  //     int *ev_out = elt_vtx_out + ielt*stride;

  //     for (int i = 0; i < stride; i++) {
  //       ev_out[i] = ev_in[ijk_to_user[i]];
  //     }

  //   }
  // }

  // else {
  //   if (elt_vtx_in != elt_vtx_out) {
  //     memcpy(elt_vtx_out, elt_vtx_in, sizeof(int) * n_elt * stride);
  //   }
  // }


  // /* IJK --> Out */
  // int *ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering_out,
  //                                                    t_elt,
  //                                                    elt_order);
  // assert(ijk_to_user != NULL;)

  // int *tmp;
  // PDM_malloc(tmp,stride,int);

  // for (int ielt = 0; ielt < n_elt; ielt++) {
  //   int *ev = elt_vtx_out + stride*i;
  //   memcpy(tmv, ev, sizeof(int) * stride);

  //   for (int i = 0; i < stride; i++) {
  //     ev[ijk_to_user[i]] = tmp[i];
  //   }
  // }

  PDM_free(tmp);
}


PDM_geometry_kind_t
PDM_Mesh_nodal_geom_kind_from_elt_type
(
 PDM_Mesh_nodal_elt_t t_elt
 )
{
  switch (PDM_Mesh_nodal_elt_dim_get(t_elt)) {
  case 0:
    return PDM_GEOMETRY_KIND_CORNER;
    break;
  case 1:
    return PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 2:
    return PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 3:
    return PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid elt type %d\n", (int) t_elt);
  }

  return PDM_GEOMETRY_KIND_MAX;
}



/**
 * \brief Return for standard elements the number of face that build this element
 *
 */
int
PDM_n_face_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_face_elt = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
   case PDM_MESH_NODAL_TRIAHO_BEZIER:
     n_face_elt = 1;
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     n_face_elt = 1;
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     n_face_elt = 4;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     n_face_elt = 5;
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     n_face_elt = 5;
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     n_face_elt = 6;
     break;
   default:
     n_face_elt = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_face_elt_per_elmt : Element type is supported\n");
  }
  return n_face_elt;
}


int
PDM_face_vtx_per_elmt
(
  PDM_Mesh_nodal_elt_t   t_elt,
  const int            **face_vtx_idx,
  const int            **face_vtx
)
{
  switch (t_elt) {
    case PDM_MESH_NODAL_POINT:
    case PDM_MESH_NODAL_BAR2:
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER: {
      *face_vtx_idx = NULL;
      *face_vtx     = NULL;
      break;
    }
    case PDM_MESH_NODAL_TRIA3:
    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER: {
      *face_vtx_idx = tria_face_vtx_idx;
      *face_vtx     = tria_face_vtx;
      break;
    }
    case PDM_MESH_NODAL_QUAD4:
    case PDM_MESH_NODAL_QUADHO: {
      *face_vtx_idx = quad_face_vtx_idx;
      *face_vtx     = quad_face_vtx;
      break;
    }
    case PDM_MESH_NODAL_TETRA4:
    case PDM_MESH_NODAL_TETRAHO: {
      *face_vtx_idx = tetra_face_vtx_idx;
      *face_vtx     = tetra_face_vtx;
      break;
    }
    case PDM_MESH_NODAL_PYRAMID5:
    case PDM_MESH_NODAL_PYRAMIDHO: {
      *face_vtx_idx = pyramid_face_vtx_idx;
      *face_vtx     = pyramid_face_vtx;
      break;
    }
    case PDM_MESH_NODAL_PRISM6:
    case PDM_MESH_NODAL_PRISMHO: {
      *face_vtx_idx = prism_face_vtx_idx;
      *face_vtx     = prism_face_vtx;
      break;
    }
    case PDM_MESH_NODAL_HEXA8:
    case PDM_MESH_NODAL_HEXAHO: {
      *face_vtx_idx = hexa_face_vtx_idx;
      *face_vtx     = hexa_face_vtx;
      break;
    }
    default : {
      PDM_error(__FILE__, __LINE__, 0, "PDM_face_vtx_per_elmt : Invalid t_elt %d\n", t_elt);
    }
  }

  return PDM_n_face_elt_per_elmt(t_elt);
}


/**
 * \brief Return for standard elements the number of edge that build this element
 *
 */
int
PDM_n_edge_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_nedge_elt = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     n_nedge_elt = 0;
     break;
   case PDM_MESH_NODAL_BAR2:
   case PDM_MESH_NODAL_BARHO:
   case PDM_MESH_NODAL_BARHO_BEZIER:
     n_nedge_elt = 1;
     break;
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
   case PDM_MESH_NODAL_TRIAHO_BEZIER:
     n_nedge_elt = 3;
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     n_nedge_elt = 4;
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     n_nedge_elt = 6;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     n_nedge_elt = 8;
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     n_nedge_elt = 9;
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     n_nedge_elt = 12;
     break;
   default:
     n_nedge_elt = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_edge_elt_per_elmt : Element type is not supported\n");
  }
  return n_nedge_elt;
}

int
PDM_edge_vtx_per_elmt
(
  PDM_Mesh_nodal_elt_t   t_elt,
  const int            **edge_vtx
)
{
  switch (t_elt) {
    case PDM_MESH_NODAL_POINT: {
      *edge_vtx = NULL;
      break;
    }
    case PDM_MESH_NODAL_BAR2:
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER: {
      *edge_vtx = bar_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_TRIA3:
    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER: {
      *edge_vtx = tria_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_QUAD4:
    case PDM_MESH_NODAL_QUADHO: {
      *edge_vtx = quad_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_TETRA4:
    case PDM_MESH_NODAL_TETRAHO: {
      *edge_vtx = tetra_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_PYRAMID5:
    case PDM_MESH_NODAL_PYRAMIDHO: {
      *edge_vtx = pyramid_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_PRISM6:
    case PDM_MESH_NODAL_PRISMHO: {
      *edge_vtx = prism_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_HEXA8:
    case PDM_MESH_NODAL_HEXAHO: {
      *edge_vtx = hexa_edge_vtx;
      break;
    }
    default : {
      PDM_error(__FILE__, __LINE__, 0, "PDM_edge_vtx_per_elmt : Invalid t_elt %d\n", t_elt);
    }
  }

  return PDM_n_edge_elt_per_elmt(t_elt);
}



/**
 * \brief Return for standard elements the total number of face vtx connectivity that build this element
 *
 */
int
PDM_n_sum_vtx_face_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_sum_vtx_face = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
   case PDM_MESH_NODAL_TRIAHO_BEZIER:
     n_sum_vtx_face = 3;
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     n_sum_vtx_face = 4;
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     n_sum_vtx_face = 12;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     n_sum_vtx_face = 16;
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     n_sum_vtx_face = 18;
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     n_sum_vtx_face = 24;
     break;
   default:
     n_sum_vtx_face = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_sum_vtx_face_per_elmt : Element type is not supported\n");
  }
  return n_sum_vtx_face;
}


/**
 * \brief Return for standard elements the total number of edge vtx connectivity that build this element
 *
 */
int
PDM_n_sum_vtx_edge_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_sum_vtx_edge = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     n_sum_vtx_edge = 0;
     break;
   case PDM_MESH_NODAL_BAR2:
   case PDM_MESH_NODAL_BARHO:
   case PDM_MESH_NODAL_BARHO_BEZIER:
     n_sum_vtx_edge = 2;
     break;
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
   case PDM_MESH_NODAL_TRIAHO_BEZIER:
     n_sum_vtx_edge = 6;
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     n_sum_vtx_edge = 8;
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     n_sum_vtx_edge = 12;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     n_sum_vtx_edge = 16;
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     n_sum_vtx_edge = 18;
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     n_sum_vtx_edge = 24;
     break;
   default:
     n_sum_vtx_edge = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_sum_vtx_edge_per_elmt : Element type is supported\n");
  }
  return n_sum_vtx_edge;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
