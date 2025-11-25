/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_lagrange_to_bezier.h"
#include "pdm_mem_tool.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm.h"

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

#define CHECK_PMNE(pmne)                                                            \
  if ((pmne) == NULL) {                                                             \
    PDM_error(__FILE__, __LINE__, 0, "Undefined Part Mesh Nodal Elmts instance\n"); \
  }

#define CHECK_I_PART(pmne, i_part)                                                          \
  if ((i_part) < 0 || (i_part) >= (pmne)->n_part) {                                         \
    PDM_error(__FILE__, __LINE__, 0, "Invalid i_part (%d / %d)", (i_part), (pmne)->n_part); \
  }

#define CHECK_BLOCK(block)                                                \
  if ((block) == NULL) {                                                  \
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n"); \
  }

#define CHECK_GROUP(pmne, i_group)                                                             \
  if ((i_group) < 0 || (i_group) >= (pmne)->n_group) {                                         \
    PDM_error(__FILE__, __LINE__, 0, "Invalid i_group (%d / %d)", (i_group), (pmne)->n_group); \
  }

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_block_std_free_partial
(
  PDM_Mesh_nodal_block_std_t *_block_std
)
{

  if (_block_std == NULL) {
    return;
  }

  if (_block_std->_connec != NULL) {
    if (_block_std->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_connec[i] != NULL)
          PDM_free(_block_std->_connec[i]);
      }
    }
    PDM_free(_block_std->_connec);
  }

  if (_block_std->_numabs != NULL) {

    if (_block_std->numabs_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        PDM_free(_block_std->_numabs[i]);
      }
    }
    PDM_free(_block_std->_numabs);
  }

  if (_block_std->ho_ordering != NULL) {
    PDM_free(_block_std->ho_ordering);
  }

  if (_block_std->_elt_to_entity != NULL) {
    if (_block_std->elt_to_entity_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        PDM_free(_block_std->_elt_to_entity[i]);
      }
    }
    PDM_free(_block_std->_elt_to_entity);
  }

}


static
PDM_Mesh_nodal_block_std_t *
_block_std_free
(
  PDM_Mesh_nodal_block_std_t *_block_std
)
{

  if (_block_std == NULL) {
    return NULL;
  }

  _block_std_free_partial(_block_std);

  PDM_free(_block_std->n_elt);

  if (_block_std->numabs_int != NULL) {
    if (_block_std->numabs_int_owner == PDM_OWNERSHIP_KEEP) {
      for (int j = 0; j < _block_std->n_part; j++) {
        PDM_free(_block_std->numabs_int[j]);
      }
      PDM_free(_block_std->numabs_int);
    }
  }

  if (_block_std->cell_centers != NULL) {
    if (_block_std->cell_centers_owner == PDM_OWNERSHIP_KEEP) {
      for (int j = 0; j < _block_std->n_part; j++) {
        PDM_free(_block_std->cell_centers[j]);
      }
      PDM_free(_block_std->cell_centers);
    }
  }

  if (_block_std->cell_centers_to_compute != NULL) {
    PDM_free(_block_std->cell_centers_to_compute);
  }

  if (_block_std->_parent_num != NULL) {
    if (_block_std->parent_num_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        PDM_free(_block_std->_parent_num[i]);
      }
    }
    PDM_free(_block_std->_parent_num);
  }

  if (_block_std->_parent_entity_g_num != NULL) {
    if (_block_std->parent_num_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        PDM_free(_block_std->_parent_entity_g_num[i]);
      }
    }
    PDM_free(_block_std->_parent_entity_g_num);
  }

  PDM_free(_block_std);
  return NULL;
}

/**
 *
 * \brief Free partially a polygon block
 *
 * \param [inout]  _block_poly2d   polygon block
 *
 */

static
void
_block_poly2d_free_partial
(
  PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
)
{

  if (_block_poly2d->_connec_idx != NULL) {
    if (_block_poly2d->elt_vtx_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        PDM_free(_block_poly2d->_connec_idx[i]);
      }
    }
    PDM_free(_block_poly2d->_connec_idx);
  }

  if (_block_poly2d->_connec != NULL) {
    if (_block_poly2d->elt_vtx_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        PDM_free(_block_poly2d->_connec[i]);
      }
    }
    PDM_free(_block_poly2d->_connec);
  }

  if (_block_poly2d->_numabs != NULL) {
    if (_block_poly2d->numabs_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        PDM_free(_block_poly2d->_numabs[i]);
      }
    }
    PDM_free(_block_poly2d->_numabs);
  }

  if (_block_poly2d->_elt_to_entity != NULL) {
    if (_block_poly2d->elt_to_entity_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        PDM_free(_block_poly2d->_elt_to_entity[i]);
      }
    }
    PDM_free(_block_poly2d->_elt_to_entity);
  }

}


/**
 *
 * \brief Free a polygon block
 *
 * \param [inout]  _bloc_poly2d    Polygon block
 *
 * \return         Null
 *
 */

static
PDM_Mesh_nodal_block_poly2d_t *
_block_poly2d_free
(
  PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
)
{
  _block_poly2d_free_partial(_block_poly2d);

  PDM_free(_block_poly2d->n_elt);

  if (_block_poly2d->numabs_int != NULL) {
    if (_block_poly2d->numabs_int_owner == PDM_OWNERSHIP_KEEP) {

      for (int j = 0; j < _block_poly2d->n_part; j++) {
        if (_block_poly2d->numabs_int[j] != NULL) {
          PDM_free(_block_poly2d->numabs_int[j]);
        }
      }
    }
    PDM_free(_block_poly2d->numabs_int);
  }

  if (_block_poly2d->cell_centers != NULL) {
    if (_block_poly2d->cell_centers_owner == PDM_OWNERSHIP_KEEP) {
      for (int j = 0; j < _block_poly2d->n_part; j++) {
        PDM_free(_block_poly2d->cell_centers[j]);
      }
    }
    PDM_free(_block_poly2d->cell_centers);
  }

  if (_block_poly2d->cell_centers_to_compute != NULL) {
    PDM_free(_block_poly2d->cell_centers_to_compute);
  }

  if (_block_poly2d->_parent_num != NULL) {
    if (_block_poly2d->parent_num_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        PDM_free(_block_poly2d->_parent_num[i]);
      }
    }
    PDM_free(_block_poly2d->_parent_num);
  }

  PDM_free(_block_poly2d);

  return NULL;
}


/**
 *
 * \brief Free partially a polyhedron block
 *
 * \param [inout]  _block_poly3d   polyhedron block
 *
 */

static
void
_block_poly3d_free_partial
(
  PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
)
{

  if (_block_poly3d->_facvtx_idx != NULL) {
    if (_block_poly3d->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_facvtx_idx[i]);
      }
    }
    PDM_free(_block_poly3d->_facvtx_idx);
  }

  if (_block_poly3d->_facvtx != NULL) {
    if (_block_poly3d->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_facvtx[i]);
      }
    }
    PDM_free(_block_poly3d->_facvtx);
  }

  if (_block_poly3d->_cellfac_idx != NULL) {
    if (_block_poly3d->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_cellfac_idx[i]);
      }
    }
    PDM_free(_block_poly3d->_cellfac_idx);
  }

  if (_block_poly3d->_cellfac != NULL) {
    if (_block_poly3d->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_cellfac[i]);
      }
    }
    PDM_free(_block_poly3d->_cellfac);
  }

  if (_block_poly3d->_cellvtx_idx != NULL) {
    if (_block_poly3d->elt_vtx_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_cellvtx_idx[i]);
      }
    }
    PDM_free(_block_poly3d->_cellvtx_idx);
  }

  if (_block_poly3d->_cellvtx != NULL) {
    if (_block_poly3d->elt_vtx_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_cellvtx[i]);
      }
    }
    PDM_free(_block_poly3d->_cellvtx);
  }

  if (_block_poly3d->_numabs != NULL) {
    if (_block_poly3d->numabs_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_numabs[i]);
      }
    }
    PDM_free(_block_poly3d->_numabs);
  }

  if (_block_poly3d->_face_ln_to_gn != NULL) {
    if (_block_poly3d->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_face_ln_to_gn[i]);
      }
    }
    PDM_free(_block_poly3d->_face_ln_to_gn);
  }

  if (_block_poly3d->_elt_to_entity != NULL) {
    if (_block_poly3d->elt_to_entity_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_elt_to_entity[i]);
      }
    }
    PDM_free(_block_poly3d->_elt_to_entity);
  }

}


/**
 *
 * \brief Free a polyhedron block
 *
 * \param [inout]  _block_poly3d    Polyhedron block
 *
 * \return         Null
 *
 */

static
PDM_Mesh_nodal_block_poly3d_t *
_block_poly3d_free
(
  PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
)
{
  _block_poly3d_free_partial(_block_poly3d);

  PDM_free(_block_poly3d->n_elt);

  PDM_free(_block_poly3d->n_face);

  if (_block_poly3d->numabs_int != NULL) {
    if (_block_poly3d->numabs_int_owner == PDM_OWNERSHIP_KEEP) {
      for (int j = 0; j < _block_poly3d->n_part; j++) {
        PDM_free(_block_poly3d->numabs_int[j]);
      }
    }
    PDM_free(_block_poly3d->numabs_int);
  }

  if (_block_poly3d->cell_centers != NULL) {
    if (_block_poly3d->cell_centers_owner == PDM_OWNERSHIP_KEEP) {
      for (int j = 0; j < _block_poly3d->n_part; j++) {
        PDM_free(_block_poly3d->cell_centers[j]);
      }
    }
    PDM_free(_block_poly3d->cell_centers);
  }

  PDM_free(_block_poly3d->cell_centers_to_compute);

  if (_block_poly3d->_parent_num != NULL) {
    if (_block_poly3d->parent_num_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_parent_num[i]);
      }
    }
    PDM_free(_block_poly3d->_parent_num);
  }

  if (_block_poly3d->_parent_entity_g_num != NULL) {
    if (_block_poly3d->parent_num_owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        PDM_free(_block_poly3d->_parent_entity_g_num[i]);
      }
    }
    PDM_free(_block_poly3d->_parent_entity_g_num);
  }

  PDM_free(_block_poly3d);
  return NULL;
}


/**
 *
 * \brief Update blocks identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static
void
_update_elmt_sections_id
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  int n_section = 0;

  if (pmne->sections_std != NULL) {
    n_section += pmne->n_section_std;
  }

  if (pmne->sections_poly2d != NULL) {
    n_section += pmne->n_section_poly2d;
  }

  if (pmne->sections_poly3d != NULL) {
    n_section += pmne->n_section_poly3d;
  }

  if (pmne->n_section < n_section) {
    PDM_realloc(pmne->sections_id, pmne->sections_id, n_section, int);
  }

  int k = 0;
  if (pmne->sections_std != NULL) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      pmne->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_STD;
    }
  }

  if (pmne->sections_poly2d != NULL) {
    for (int i = 0; i < pmne->n_section_poly2d; i++) {
      pmne->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }

  if (pmne->sections_poly3d != NULL) {
    for (int i = 0; i < pmne->n_section_poly3d; i++) {
      pmne->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }

  pmne->n_section = n_section;
}


/**
 *
 * \brief Build tetrahedron nodal connectivity from faces connectivity
 *
 *   \param[in]  vtx         Vertices coordinates
 *   \param[in]  tria_vtx    Faces connectivity
 *   \param[out] tetra_vtx   Tetrahedron connectivity
 *
 */
static
void
_connec_tetra
(
  const double *vtx_coord,
        int    *tria_vtx,
        int     tetra_vtx[]
)
{

  /* Initialization */
  tetra_vtx[0] = tria_vtx[0];
  tetra_vtx[1] = tria_vtx[1];
  tetra_vtx[2] = tria_vtx[2];

  for (int i = 3; i < 11; i++) {
    if ((tria_vtx[i] != tetra_vtx[0]) &&
        (tria_vtx[i] != tetra_vtx[1]) &&
        (tria_vtx[i] != tetra_vtx[2]))
      tetra_vtx[3] = tria_vtx[i];
  }

  /* Orientation */
  const double *_coords = vtx_coord;
  double v1[3];
  double v2[3];
  double v3[3];
  double n[3];

  for (int i = 0; i < 3; i++) {
    v1[i] = _coords[3*(tetra_vtx[1] - 1) + i] - _coords[3*(tetra_vtx[0] - 1) + i];
    v2[i] = _coords[3*(tetra_vtx[2] - 1) + i] - _coords[3*(tetra_vtx[0] - 1) + i];
    v3[i] = _coords[3*(tetra_vtx[3] - 1) + i] - _coords[3*(tetra_vtx[0] - 1) + i];
  }

  PDM_CROSS_PRODUCT(n, v1, v2);
  double orient = PDM_DOT_PRODUCT(v3, n);

  if (orient < 0) {
    tetra_vtx[0] = tria_vtx[2];
    tetra_vtx[1] = tria_vtx[1];
    tetra_vtx[2] = tria_vtx[0];
  }
}

/**
 *
 * \brief Build prism nodal connectivity from faces connectivity
 *
 *   \param[in]  vtx         Vertices coordinates
 *   \param[in]  tria_vtx    Faces connectivity
 *   \param[in]  quad_vtx    Faces connectivity
 *   \param[out] prism_vtx   Prism connectivity
 *
 */

static
void
_connec_prism
(
  const double *vtx_coord,
        int    *tria_vtx,
        int    *quad_vtx,
        int     prism_vtx[]
)
{

  /* Initialisation */

  for (int i = 0; i < 6; i++)
    prism_vtx[i] = tria_vtx[i];

  /* Orientation des faces */
  const double *_coords = vtx_coord;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++) {
      c[3*i+k] = 0.;
    }
    for (int j = 0; j < 3; j++) {
      int isom = prism_vtx[3*i+j] - 1;
      for (int k = 0; k < 3; k++) {
        c[3*i+k] += _coords[3*isom+k];
      }
    }
    for (int k = 0; k < 3; k++) {
      c[3*i+k] /= 3.;
    }

    for (int k = 0; k < 3; k++) {
      n[3*i+k] = 0.;
    }

    double v1[3];
    double v2[3];
    int isom3 = prism_vtx[3*i+2] - 1 ;
    int isom2 = prism_vtx[3*i+1] - 1;
    int isom1 = prism_vtx[3*i  ] - 1;

    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom2+k] - _coords[3*isom1+k];
      v2[k] = _coords[3*isom3+k] - _coords[3*isom1+k];
    }
    PDM_CROSS_PRODUCT(n + 3*i, v1, v2);
  }

  double cc[3];
  for (int k = 0; k < 3; k++) {
    cc[k] = c[3+k] - c[k];
  }

  double orientation  = PDM_DOT_PRODUCT(cc, n);
  double orientation2 = PDM_DOT_PRODUCT(cc, n+3);

  if (orientation < 0) {
    int tmp = prism_vtx[1];
    prism_vtx[1] = prism_vtx[2];
    prism_vtx[2] = tmp;
  }

  if (orientation2 < 0) {
    int tmp = prism_vtx[4];
    prism_vtx[4] = prism_vtx[5];
    prism_vtx[5] = tmp;
  }

  /* Permutation circulaire */

  int id1 = -1;
  for (int j = 0; j < 12; j++) {
    if (quad_vtx[j] == prism_vtx[0]) {
      id1 = j;
      break;
    }
  }

  int id2 = (id1 / 4) * 4 + (id1 + 1) % 4;
  if ((quad_vtx[id2] == prism_vtx[1]) ||
      (quad_vtx[id2] == prism_vtx[2])) {
    id2 =  (id1 / 4) * 4 + (id1 + 3) % 4;
  }
  int id_deb = -1;
  for (int j = 0; j < 3; j++) {
    if (quad_vtx[id2] == prism_vtx[3+j]) {
      id_deb = j;
      break;
    }
  }

  int tmp[3];
  for (int j = 0; j < 3; j++) {
    tmp[j] = prism_vtx[3+j];
  }

  for (int j = 0; j < 3; j++) {
    int idx = (id_deb + j) % 3;
    prism_vtx[3+j] = tmp[idx];
  }

}


/**
 *
 * \brief Build pyramid nodal connectivity from faces connectivity
 *
 *   \param[in]  vtx         Vertices coordinates
 *   \param[in]  tria_vtx    Faces connectivity
 *   \param[in]  quad_vtx    Faces connectivity
 *   \param[out] pyramid_vtx Pyramid connectivity
 *
 */

static
void
_connec_pyramid
(
  const double  *vtx_coord,
        int     *tria_vtx,
        int     *quad_vtx,
        int      pyramid_vtx[]
)
{

  /* Initialisation */

  pyramid_vtx[0] = quad_vtx[0];
  pyramid_vtx[1] = quad_vtx[1];
  pyramid_vtx[2] = quad_vtx[2];
  pyramid_vtx[3] = quad_vtx[3];

  for (int i = 0; i < 9; i++) {
    if ((tria_vtx[i] != pyramid_vtx[0]) &&
        (tria_vtx[i] != pyramid_vtx[1]) &&
        (tria_vtx[i] != pyramid_vtx[2]) &&
        (tria_vtx[i] != pyramid_vtx[3])) {
      pyramid_vtx[4] = tria_vtx[i];
      break;
    }
  }

  /* Orientation */

  const double *_coords = vtx_coord;

  double c[3];
  double n[3];

  for (int k = 0; k < 3; k++) {
    c[k] = 0.;
  }
  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    for (int k = 0; k < 3; k++) {
      c[k] += _coords[3*isom+k];
    }
  }
  for (int k = 0; k < 3; k++) {
    c[k] *= 0.25;
  }

  for (int k = 0; k < 3; k++) {
    n[k] = 0.;
  }

  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    int suiv = (j+1) % 4;
    int isom_suiv = pyramid_vtx[suiv] - 1;

    double v1[3];
    double v2[3];
    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom     +k] - c[k];
      v2[k] = _coords[3*isom_suiv+k] - c[k];
    }

    PDM_CROSS_PRODUCT(n, v1, v2);

  }

  double cc[3];
  for (int k = 0; k < 3; k++) {
    cc[k] = _coords[3*(pyramid_vtx[3] - 1) + k] - c[k];
  }

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = PDM_DOT_PRODUCT(cc, n);

  if (orientation < 0) {
    int tmp = pyramid_vtx[0];
    pyramid_vtx[0] = pyramid_vtx[3];
    pyramid_vtx[3] = tmp;
    tmp = pyramid_vtx[1];
    pyramid_vtx[1] = pyramid_vtx[2];
    pyramid_vtx[2] = tmp;
  }

}


/**
 *
 * \brief Build hexahedron nodal connectivity from faces connectivity
 *
 *   \param[in]  vtx         Vertices coordinates
 *   \param[in]  quad_vtx    Faces connectivity
 *   \param[out] hexa_vtx    Hexahedron connectivity
 *
 */

static
void
_connec_hexa
(
  const double  *vtx_coord,
        int     *quad_vtx,
        int      hexa_vtx[]
)
{

  /* Initialization */

  hexa_vtx[0] = quad_vtx[0];
  hexa_vtx[1] = quad_vtx[1];
  hexa_vtx[2] = quad_vtx[2];
  hexa_vtx[3] = quad_vtx[3];

  int face_contact[4] = {-1, -1, -1, -1};

  for (int i = 1; i < 6; i++) {
    int cpt = 0;
    for (int j = 0; j < 4; j++) {
      int som_courant = quad_vtx[4*i+j];
      if ((som_courant != hexa_vtx[0]) &&
          (som_courant != hexa_vtx[1]) &&
          (som_courant != hexa_vtx[2]) &&
          (som_courant != hexa_vtx[3]))
        cpt += 1;
    }
    if (cpt == 4) {
      hexa_vtx[4] = quad_vtx[4*i];
      hexa_vtx[5] = quad_vtx[4*i+1];
      hexa_vtx[6] = quad_vtx[4*i+2];
      hexa_vtx[7] = quad_vtx[4*i+3];
    }
    if (cpt == 2) {
      face_contact[0] = quad_vtx[4*i];
      face_contact[1] = quad_vtx[4*i+1];
      face_contact[2] = quad_vtx[4*i+2];
      face_contact[3] = quad_vtx[4*i+3];
    }
  }

  /* Calcul des centres et normales de la base et de la face opposee */

  const double *_coords = vtx_coord;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++) {
      c[3*i+k] = 0.;
    }
    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      for (int k = 0; k < 3; k++) {
        c[3*i+k] += _coords[3*isom+k];
      }
    }
    for (int k = 0; k < 3; k++) {
      c[3*i+k] *= 0.25;
    }

    for (int k = 0; k < 3; k++) {
      n[3*i+k] = 0.;
    }

    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      int suiv = (j+1) % 4;
      int isom_suiv = hexa_vtx[4*i+suiv] - 1;

      double v1[3];
      double v2[3];
      for (int k = 0; k < 3; k++) {
        v1[k] = _coords[3*isom     +k] - c[3*i+k];
        v2[k] = _coords[3*isom_suiv+k] - c[3*i+k];
      }

      PDM_CROSS_PRODUCT(n + 3*i, v1, v2);

    }

  }

  double cc[3];
  for (int k = 0; k < 3; k++) {
    cc[k] = c[3+k] - c[k];
  }

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation  = PDM_DOT_PRODUCT(cc, n);
  double orientation2 = PDM_DOT_PRODUCT(cc, n+3);

  if (orientation < 0) {
    int tmp = hexa_vtx[0];
    hexa_vtx[0] = hexa_vtx[3];
    hexa_vtx[3] = tmp;
    tmp = hexa_vtx[1];
    hexa_vtx[1] = hexa_vtx[2];
    hexa_vtx[2] = tmp;
  }

  if (orientation2 < 0) {
    int tmp = hexa_vtx[4];
    hexa_vtx[4] = hexa_vtx[7];
    hexa_vtx[7] = tmp;
    tmp = hexa_vtx[5];
    hexa_vtx[5] = hexa_vtx[6];
    hexa_vtx[6] = tmp;
  }

  /* Permutation circulaire eventuelle de la face sup */

  int id1 = -1;
  int k1 = -1;
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      if (face_contact[j] == hexa_vtx[k]) {
        id1 = j;
        k1 = k;
        break;
      }
      if (id1 != -1) {
        break;
      }
    }
  }

  if (k1 == -1) {
    PDM_printf("Error connect_hexa : %d %d %d %d %d %d %d %d\n",
               hexa_vtx[0],
               hexa_vtx[1],
               hexa_vtx[2],
               hexa_vtx[3],
               hexa_vtx[4],
               hexa_vtx[5],
               hexa_vtx[6],
               hexa_vtx[7]);

    for (int i10 = 0; i10 < 4; i10++) {
      PDM_printf("   face %d : %d %d %d %d\n", i10+1, quad_vtx[4*i10],
                 quad_vtx[4*i10+1],
                 quad_vtx[4*i10+2],
                 quad_vtx[4*i10+3]);
    }
    abort();

  }

  int id2 = (id1 + 1) % 4;
  int k2 = (k1 + 1) % 4;
  int k3 = (k1 + 3) % 4;

  if ((face_contact[id2] == hexa_vtx[k2]) ||
      (face_contact[id2] == hexa_vtx[k3])) {
    id2 = (id1 + 3) % 4;
  }

  int id_deb = -1;
  for (int j = 0; j < 4; j++) {
    if (face_contact[id2] == hexa_vtx[4+j]) {
      id_deb = (j - k1);
      if (id_deb < 0) {
        id_deb += 4;
      }
      id_deb = id_deb % 4;
      break;
    }
  }

  int tmp[4];
  for (int j = 0; j < 4; j++) {
    tmp[j] = hexa_vtx[4+j];
  }

  for (int j = 0; j < 4; j++) {
    int idx = (id_deb + j) % 4;
    hexa_vtx[4+j] = tmp[idx];
  }
}



static
int
_binary_search
(
  const PDM_l_num_t  elem,
  const PDM_l_num_t  array[],
  const PDM_l_num_t  n,
        PDM_bool_t  *in_array
)
{
  int l = 0;
  int r = n;

  *in_array = PDM_FALSE;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m]) {
      r = m;
    }
    else {
      l = m;
    }
  }

  if (array[l] == elem) {
    *in_array = PDM_TRUE;
    return l;
  }
  else if (array[l] < elem) {
    return l + 1;
  }
  else {
    return l;
  }
}


static
void
_compute_cell_vtx_connectivity
(
  const PDM_l_num_t   n_cell,
  const PDM_l_num_t   n_face,
  const PDM_l_num_t  *face_vtx_idx,
  const PDM_l_num_t  *face_vtx,
  const PDM_l_num_t  *cell_face_idx,
  const PDM_l_num_t  *cell_face,
        PDM_l_num_t **cell_vtx_idx,
        PDM_l_num_t **cell_vtx
)
{
  PDM_UNUSED(n_face);

  const int dbg_enabled = 0;

  PDM_malloc(*cell_vtx_idx, n_cell + 1, int);
  PDM_l_num_t *_cell_vtx_idx = *cell_vtx_idx;

  _cell_vtx_idx[0] = 0;

  size_t s_cell_vtx = 10 * n_cell;
  PDM_malloc(*cell_vtx, s_cell_vtx, PDM_l_num_t);

  PDM_bool_t already_in_cell;
  int pos, i;
  PDM_l_num_t icell, iface, ivtx, id_face, id_vtx;

  /* Loop on cells */
  PDM_l_num_t n_vtx_cell;
  for (icell = 0; icell < n_cell; icell++) {

    PDM_l_num_t *_cell_vtx = *cell_vtx + _cell_vtx_idx[icell];
    n_vtx_cell = 0;

    /* Loop on current cell's faces */
    for (iface = cell_face_idx[icell]; iface < cell_face_idx[icell+1]; iface++) {
      id_face = PDM_ABS (cell_face[iface]) - 1;
      /* Loop on current face's vertices */
      for (ivtx = face_vtx_idx[id_face]; ivtx < face_vtx_idx[id_face+1]; ivtx++) {
        id_vtx = face_vtx[ivtx];

        pos = _binary_search(id_vtx,
                             _cell_vtx,
                             n_vtx_cell,
                            &already_in_cell);

        if (already_in_cell == PDM_TRUE) {
          continue;
        }

        if (n_vtx_cell + _cell_vtx_idx[icell] >= (int) s_cell_vtx) {
          s_cell_vtx = PDM_MAX ((int) (2*s_cell_vtx), n_vtx_cell + _cell_vtx_idx[icell]);
          PDM_realloc(*cell_vtx, *cell_vtx, s_cell_vtx, PDM_l_num_t);
          _cell_vtx = *cell_vtx + _cell_vtx_idx[icell];
        }

        for (i = n_vtx_cell; i > pos; i--) {
          _cell_vtx[i] = _cell_vtx[i-1];
        }
        _cell_vtx[pos] = id_vtx;
        n_vtx_cell++;

      } // End of loop on current face's vertices

    } // End of loop on current cell's faces

    _cell_vtx_idx[icell+1] = _cell_vtx_idx[icell] + n_vtx_cell;

    if (dbg_enabled) {
      printf("cell #%d vtx =", icell);
      for (int j = _cell_vtx_idx[icell]; j < _cell_vtx_idx[icell+1]; j++) {
        printf(" %d", (*cell_vtx)[j]);
      }
      printf("\n");
    }

  } // End of loop on cells

  PDM_realloc(*cell_vtx, *cell_vtx, _cell_vtx_idx[n_cell], PDM_l_num_t);
}

inline
static
PDM_Mesh_nodal_elt_t
_type_cell_3D
(
 const int     n_face_cell,
 const int    *cell_face,
 const int    *face_vtx_idx,
 const int    *face_vtx,
       int     tria_vtx[],
       int     quad_vtx[]
)
{

  int  n_trias = 0;
  int  n_quads = 0;

  if (n_face_cell > 6) {
    return PDM_MESH_NODAL_POLY_3D;
  }

  for (int i = 0; i < n_face_cell; i++) {

    const int face_id = PDM_ABS(cell_face[i]) - 1;
    const int n_som_face = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
    int idx = face_vtx_idx[face_id] ;

    if (n_som_face == 3) {
      int *cell_som_tria_courant = tria_vtx + 3*n_trias;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_tria_courant[j-idx] = face_vtx[j];
      }
      n_trias += 1;
    }
    else if (n_som_face == 4) {
      int *cell_som_quad_courant = quad_vtx + 4*n_quads;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_quad_courant[j-idx] = face_vtx[j];
      }
      n_quads += 1;
    }
    else {
      return PDM_MESH_NODAL_POLY_3D;
    }
  }

  PDM_Mesh_nodal_elt_t cell_type;

  if ((n_quads == 0) && (n_trias == 4)) {
    cell_type = PDM_MESH_NODAL_TETRA4;
  }
  else if (n_quads == 6) {
    cell_type = PDM_MESH_NODAL_HEXA8;
  }
  else if ((n_quads == 1) && (n_trias == 4)) {
    cell_type = PDM_MESH_NODAL_PYRAMID5;
  }
  else if ((n_quads == 3) && (n_trias == 2)) {
    int trias[6];
    n_trias = 0;
    for (int i = 0; i < n_face_cell; i++) {

      const int face_id = PDM_ABS(cell_face[i]) - 1;
      const int ideb = face_vtx_idx[face_id] ;

      const int n_som_face = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];

      if (n_som_face == 3) {
        for (int j = 0; j < 3; j++) {
          trias[3*n_trias+j] = face_vtx[ideb+j];
        }
        n_trias += 1;
      }
      if (n_trias >= 2) {
        break;
      }
    }

    cell_type = PDM_MESH_NODAL_PRISM6;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (trias[i] == trias[3+j]) {
          cell_type = PDM_MESH_NODAL_POLY_3D;
          break;
        }
      }
      if (cell_type == PDM_MESH_NODAL_POLY_3D) {
        break;
      }
    }
  }

  else {
    cell_type = PDM_MESH_NODAL_POLY_3D;
  }

  return cell_type;

}



/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
  const int          mesh_dimension,
  const int          n_part,
  const PDM_MPI_Comm comm
)
{
  PDM_part_mesh_nodal_elmts_t *pmne;
  PDM_malloc(pmne, 1, PDM_part_mesh_nodal_elmts_t);

  pmne->comm             = comm;
  pmne->mesh_dimension   = mesh_dimension;
  pmne->n_part           = n_part;

  pmne->n_elmts          = PDM_array_zeros_int(n_part);

  pmne->n_section        = 0;
  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  pmne->sections_id      = NULL;
  pmne->sections_std     = NULL;
  pmne->sections_poly2d  = NULL;
  pmne->sections_poly3d  = NULL;

  pmne->prepa_blocks             = NULL;
  pmne->num_elmt_parent_to_local = NULL;
  pmne->numabs                   = NULL;

  pmne->ownership_group  = NULL;
  pmne->ownership_numabs = PDM_OWNERSHIP_KEEP;
  pmne->n_group          = 0;
  pmne->n_group_elmt     = NULL;
  pmne->group_elmt       = NULL;
  pmne->group_ln_to_gn   = NULL;

  return pmne;
}


int
PDM_part_mesh_nodal_elmts_add
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const PDM_Mesh_nodal_elt_t         t_elt
)
{
  CHECK_PMNE(pmne)

  int elt_dim = PDM_Mesh_nodal_elt_dim_get(t_elt);
  if (elt_dim != pmne->mesh_dimension) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_part_mesh_nodal_elmts_add = expected = %i and given = %i \n",
               pmne->mesh_dimension, elt_dim);
  }

  int id_section = -1;

  switch (t_elt) {

  case PDM_MESH_NODAL_POINT        :
  case PDM_MESH_NODAL_BAR2         :
  case PDM_MESH_NODAL_TRIA3        :
  case PDM_MESH_NODAL_QUAD4        :
  case PDM_MESH_NODAL_TETRA4       :
  case PDM_MESH_NODAL_PYRAMID5     :
  case PDM_MESH_NODAL_PRISM6       :
  case PDM_MESH_NODAL_HEXA8        :
  case PDM_MESH_NODAL_BARHO        :
  case PDM_MESH_NODAL_TRIAHO       :
  case PDM_MESH_NODAL_BARHO_BEZIER :
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
  case PDM_MESH_NODAL_QUADHO       :
  case PDM_MESH_NODAL_TETRAHO      :
  case PDM_MESH_NODAL_PYRAMIDHO    :
  case PDM_MESH_NODAL_PRISMHO      :
  case PDM_MESH_NODAL_HEXAHO       :
    {
      /* Mise a jour du tableau de stockage */

      pmne->n_section_std++;

      PDM_realloc(pmne->sections_std, pmne->sections_std, pmne->n_section_std, PDM_Mesh_nodal_block_std_t *);

      id_section = pmne->n_section_std-1;

      /* Intialisation du bloc */
      PDM_malloc(pmne->sections_std[id_section], 1, PDM_Mesh_nodal_block_std_t);
      pmne->sections_std[id_section]->t_elt  = t_elt;
      pmne->sections_std[id_section]->n_part = pmne->n_part;

      /* Ownership */
      pmne->sections_std[id_section]->owner               = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_section]->cell_centers_owner  = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_section]->numabs_int_owner    = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_section]->numabs_owner        = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_section]->parent_num_owner    = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_section]->elt_to_entity_owner = PDM_OWNERSHIP_KEEP;

      PDM_malloc(pmne->sections_std[id_section]->n_elt,   pmne->sections_std[id_section]->n_part, int          );
      PDM_malloc(pmne->sections_std[id_section]->_connec, pmne->sections_std[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_std[id_section]->_numabs, pmne->sections_std[id_section]->n_part, PDM_g_num_t *);
      pmne->sections_std[id_section]->numabs_int              = NULL;
      pmne->sections_std[id_section]->_parent_num             = NULL;
      pmne->sections_std[id_section]->_parent_entity_g_num    = NULL;
      pmne->sections_std[id_section]->cell_centers            = NULL;
      pmne->sections_std[id_section]->cell_centers_to_compute = NULL;
      pmne->sections_std[id_section]->order                   = 1;
      pmne->sections_std[id_section]->ho_ordering             = NULL;
      pmne->sections_std[id_section]->_elt_to_entity          = NULL;


      for (int i = 0; i < pmne->sections_std[id_section]->n_part; i++) {
        pmne->sections_std[id_section]->n_elt  [i] = 0;
        pmne->sections_std[id_section]->_connec[i] = NULL;
        pmne->sections_std[id_section]->_numabs[i] = NULL;
      }

      id_section += PDM_BLOCK_ID_BLOCK_STD;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard blocks must be less than %d\n",
                  PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      /* Mise a jour du tableau de stockage */

      pmne->n_section_poly2d++;

      PDM_realloc(pmne->sections_poly2d, pmne->sections_poly2d, pmne->n_section_poly2d, PDM_Mesh_nodal_block_poly2d_t *);

      id_section = pmne->n_section_poly2d-1;

      /* Intialisation du bloc */
      PDM_malloc(pmne->sections_poly2d[id_section], 1, PDM_Mesh_nodal_block_poly2d_t);

      /* Ownership */
      pmne->sections_poly2d[id_section]->owner               = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly2d[id_section]->cell_centers_owner  = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly2d[id_section]->elt_vtx_owner       = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly2d[id_section]->numabs_int_owner    = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly2d[id_section]->numabs_owner        = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly2d[id_section]->parent_num_owner    = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly2d[id_section]->elt_to_entity_owner = PDM_OWNERSHIP_KEEP;

      pmne->sections_poly2d[id_section]->n_part = pmne->n_part;
      PDM_malloc(pmne->sections_poly2d[id_section]->n_elt,       pmne->sections_poly2d[id_section]->n_part, int          );
      PDM_malloc(pmne->sections_poly2d[id_section]->_connec_idx, pmne->sections_poly2d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly2d[id_section]->_connec,     pmne->sections_poly2d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly2d[id_section]->_numabs,     pmne->sections_poly2d[id_section]->n_part, PDM_g_num_t *);
      pmne->sections_poly2d[id_section]->numabs_int              = NULL;
      pmne->sections_poly2d[id_section]->cell_centers            = NULL;
      pmne->sections_poly2d[id_section]->cell_centers_to_compute = NULL;
      pmne->sections_poly2d[id_section]->_parent_num             = NULL;
      pmne->sections_poly2d[id_section]->_parent_entity_g_num    = NULL;
      pmne->sections_poly2d[id_section]->_elt_to_entity          = NULL;

      for (int i = 0; i < pmne->sections_poly2d[id_section]->n_part; i++) {
        pmne->sections_poly2d[id_section]->n_elt      [i] = 0;
        pmne->sections_poly2d[id_section]->_connec_idx[i] = NULL;
        pmne->sections_poly2d[id_section]->_connec    [i] = NULL;
        pmne->sections_poly2d[id_section]->_numabs    [i] = NULL;
      }

      id_section += PDM_BLOCK_ID_BLOCK_POLY2D;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of polygon blocks must be less than %d\n",
                  PDM_BLOCK_ID_BLOCK_POLY3D - PDM_BLOCK_ID_BLOCK_POLY2D);
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      pmne->n_section_poly3d++;

      PDM_realloc(pmne->sections_poly3d, pmne->sections_poly3d, pmne->n_section_poly3d, PDM_Mesh_nodal_block_poly3d_t *);

      id_section = pmne->n_section_poly3d-1;

      /* Intialisation du bloc */

      PDM_malloc(pmne->sections_poly3d[id_section], 1, PDM_Mesh_nodal_block_poly3d_t);
      pmne->sections_poly3d[id_section]->n_part       = pmne->n_part;

      /* Ownership */
      pmne->sections_poly3d[id_section]->owner               = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly3d[id_section]->cell_centers_owner  = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly3d[id_section]->elt_vtx_owner       = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly3d[id_section]->numabs_int_owner    = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly3d[id_section]->numabs_owner        = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly3d[id_section]->parent_num_owner    = PDM_OWNERSHIP_KEEP;
      pmne->sections_poly3d[id_section]->elt_to_entity_owner = PDM_OWNERSHIP_KEEP;

      PDM_malloc(pmne->sections_poly3d[id_section]->n_elt,          pmne->sections_poly3d[id_section]->n_part, int          );
      PDM_malloc(pmne->sections_poly3d[id_section]->n_face,         pmne->sections_poly3d[id_section]->n_part, int          );
      PDM_malloc(pmne->sections_poly3d[id_section]->_facvtx_idx,    pmne->sections_poly3d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_facvtx,        pmne->sections_poly3d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_face_ln_to_gn, pmne->sections_poly3d[id_section]->n_part, PDM_g_num_t *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_cellfac_idx,   pmne->sections_poly3d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_cellfac,       pmne->sections_poly3d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_cellvtx_idx,   pmne->sections_poly3d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_cellvtx,       pmne->sections_poly3d[id_section]->n_part, int         *);
      PDM_malloc(pmne->sections_poly3d[id_section]->_numabs,        pmne->sections_poly3d[id_section]->n_part, PDM_g_num_t *);
      pmne->sections_poly3d[id_section]->numabs_int              = NULL;
      pmne->sections_poly3d[id_section]->cell_centers            = NULL;
      pmne->sections_poly3d[id_section]->cell_centers_to_compute = NULL;
      pmne->sections_poly3d[id_section]->_parent_num             = NULL;
      pmne->sections_poly3d[id_section]->_parent_entity_g_num    = NULL;
      pmne->sections_poly3d[id_section]->_elt_to_entity          = NULL;

      for (int i = 0; i < pmne->sections_poly3d[id_section]->n_part; i++) {
        pmne->sections_poly3d[id_section]->n_elt         [i] = 0;
        pmne->sections_poly3d[id_section]->n_face        [i] = 0;
        pmne->sections_poly3d[id_section]->_facvtx_idx   [i] = NULL;
        pmne->sections_poly3d[id_section]->_facvtx       [i] = NULL;
        pmne->sections_poly3d[id_section]->_face_ln_to_gn[i] = NULL;
        pmne->sections_poly3d[id_section]->_cellfac_idx  [i] = NULL;
        pmne->sections_poly3d[id_section]->_cellfac      [i] = NULL;
        pmne->sections_poly3d[id_section]->_cellvtx_idx  [i] = NULL;
        pmne->sections_poly3d[id_section]->_cellvtx      [i] = NULL;
        pmne->sections_poly3d[id_section]->_numabs       [i] = NULL;
      }

      id_section += PDM_BLOCK_ID_BLOCK_POLY3D;
    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _update_elmt_sections_id(pmne);

  return id_section;
}


void
PDM_part_mesh_nodal_elmts_std_set
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
  const int                          n_elt,
  const int                         *connec,
  const PDM_g_num_t                 *numabs,
  const int                         *parent_num,
  const PDM_g_num_t                 *parent_entity_g_num,
        PDM_ownership_t              owner
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  /* Mapping */
  pmne->n_elmts [id_part] += -block->n_elt[id_part];
  pmne->n_elmts [id_part] += n_elt;
  block->n_elt  [id_part]  = n_elt;
  block->_connec[id_part]  = (int *) connec;
  block->_numabs[id_part]  = (PDM_g_num_t *) numabs;

  if (owner != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->owner            != PDM_OWNERSHIP_USER) block->owner            = owner;
    if (block->numabs_owner     != PDM_OWNERSHIP_USER) block->numabs_owner     = owner;
    if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = owner;
  }

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      PDM_malloc(block->_parent_num, block->n_part, int *);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (int *) parent_num;
  }

  if (parent_entity_g_num != NULL) {
    if (block->_parent_entity_g_num == NULL) {
      PDM_malloc(block->_parent_entity_g_num, block->n_part, PDM_g_num_t *);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_entity_g_num[i] = NULL;
      }
    }
    block->_parent_entity_g_num[id_part] = (PDM_g_num_t *) parent_entity_g_num;
  }

  block->order       = 1;
  block->ho_ordering = NULL;
}


void
PDM_part_mesh_nodal_elmts_std_ho_set
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
  const int                          n_elt,
  const int                         *connec,
  const PDM_g_num_t                 *numabs,
  const int                         *parent_num,
  const PDM_g_num_t                 *parent_entity_g_num,
  const int                          order,
  const char                        *ho_ordering,
        PDM_ownership_t              owner
)
{
  PDM_part_mesh_nodal_elmts_std_set(pmne,
                                    id_section,
                                    id_part,
                                    n_elt,
                                    connec,
                                    numabs,
                                    parent_num,
                                    parent_entity_g_num,
                                    owner);

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

  block->order = order;
  if (block->ho_ordering != NULL) {
    PDM_free(block->ho_ordering);
    block->ho_ordering = NULL;
  }
  if (ho_ordering != NULL) {
    PDM_malloc(block->ho_ordering, strlen(ho_ordering) + 1, char);
    strcpy(block->ho_ordering, ho_ordering);
  }
}


void
PDM_part_mesh_nodal_elmts_section_std_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const int                           id_part,
        int                         **connec,
        PDM_g_num_t                 **numabs,
        int                         **parent_num,
        PDM_g_num_t                 **parent_entity_g_num,
        PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  *connec     = block->_connec[id_part];
  *numabs     = block->_numabs[id_part];
  *parent_num = NULL;
  if(block->_parent_num != NULL) {
    *parent_num = block->_parent_num[id_part];
  }

  *parent_entity_g_num = NULL;
  if(block->_parent_entity_g_num != NULL) {
    *parent_entity_g_num = block->_parent_entity_g_num[id_part];
  }

  // ownership
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->owner            != PDM_OWNERSHIP_USER) block->owner            = ownership;
    if (block->numabs_owner     != PDM_OWNERSHIP_USER) block->numabs_owner     = ownership;
    if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = ownership;
  }
}


void
PDM_part_mesh_nodal_elmts_section_std_ho_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const int                           id_part,
        int                         **connec,
        PDM_g_num_t                 **numabs,
        int                         **parent_num,
        PDM_g_num_t                 **parent_entity_g_num,
        int                          *order,
  const char                        **ho_ordering,
        PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  *connec = block->_connec[id_part];
  *numabs = block->_numabs[id_part];
  if (block->_parent_num != NULL) {
    *parent_num = block->_parent_num[id_part];
  }
  *parent_entity_g_num = NULL;
  if(block->_parent_entity_g_num != NULL) {
    *parent_entity_g_num = block->_parent_entity_g_num[id_part];
  }
  *order       = block->order;
  *ho_ordering = block->ho_ordering;

  // ownership
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->owner            != PDM_OWNERSHIP_USER) block->owner            = ownership;
    if (block->numabs_owner     != PDM_OWNERSHIP_USER) block->numabs_owner     = ownership;
    if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = ownership;
  }
}


void
PDM_part_mesh_nodal_elmts_section_poly2d_set
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
  const int                          n_elt,
  const int                         *connec_idx,
  const int                         *connec,
  const PDM_g_num_t                 *numabs,
  const int                         *parent_num,
        PDM_ownership_t              owner
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  /* Mapping */
  pmne->n_elmts[id_part] += -block->n_elt[id_part];
  pmne->n_elmts[id_part] += n_elt;
  block->n_elt[id_part]       = n_elt;
  block->_connec_idx[id_part] = (int *) connec_idx;
  block->_connec[id_part]     = (int *) connec;
  block->_numabs[id_part]     = (PDM_g_num_t *) numabs;

  // ownership
  if (owner != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->owner            != PDM_OWNERSHIP_USER) block->owner            = owner;
    if (block->numabs_owner     != PDM_OWNERSHIP_USER) block->numabs_owner     = owner;
    if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = owner;
    if (block->elt_vtx_owner    != PDM_OWNERSHIP_USER) block->elt_vtx_owner    = owner;
  }

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      PDM_malloc(block->_parent_num, block->n_part, int *);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (int *) parent_num;
  }
}


void
PDM_part_mesh_nodal_elmts_section_poly3d_set
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
  const int                          n_elt,
  const int                          n_face,
  const int                         *facvtx_idx,
  const int                         *facvtx,
  const PDM_g_num_t                 *face_ln_to_gn,
  const int                         *cellfac_idx,
  const int                         *cellfac,
  const PDM_g_num_t                 *numabs,
  const int                         *parent_num,
  const PDM_g_num_t                 *parent_entity_g_num,
        PDM_ownership_t              owner
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;


  PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  pmne->n_elmts[id_part] += -block->n_elt[id_part];
  pmne->n_elmts[id_part] += n_elt;

  block->n_elt         [id_part] = n_elt;
  block->n_face        [id_part] = n_face;
  block->_facvtx_idx   [id_part] = (int         *) facvtx_idx;
  block->_facvtx       [id_part] = (int         *) facvtx;
  block->_face_ln_to_gn[id_part] = (PDM_g_num_t *) face_ln_to_gn;
  block->_cellfac_idx  [id_part] = (int         *) cellfac_idx;
  block->_cellfac      [id_part] = (int         *) cellfac;
  block->_numabs       [id_part] = (PDM_g_num_t *) numabs;

  // ownership
  if (owner != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->owner            != PDM_OWNERSHIP_USER) block->owner            = owner;
    if (block->numabs_owner     != PDM_OWNERSHIP_USER) block->numabs_owner     = owner;
    if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = owner;
  }

  /* Compute cell-vertex connectivity */
  _compute_cell_vtx_connectivity(n_elt,
                                 n_face,
                                 facvtx_idx,
                                 facvtx,
                                 cellfac_idx,
                                 cellfac,
                                &block->_cellvtx_idx[id_part],
                                &block->_cellvtx    [id_part]);

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      PDM_malloc(block->_parent_num, block->n_part, int *);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (int *) parent_num;
  }

  if (parent_entity_g_num != NULL) {
    if (block->_parent_entity_g_num == NULL) {
      PDM_malloc(block->_parent_entity_g_num, block->n_part, PDM_g_num_t *);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_entity_g_num[i] = NULL;
      }
    }
    block->_parent_entity_g_num[id_part] = (PDM_g_num_t *) parent_entity_g_num;
  }
}


void
PDM_part_mesh_nodal_elmts_section_poly2d_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const int                           id_part,
        int                         **connec_idx,
        int                         **connec,
        PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  *connec_idx = block->_connec_idx[id_part];
  *connec     = block->_connec    [id_part];

  // ownership
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->elt_vtx_owner != PDM_OWNERSHIP_USER) block->elt_vtx_owner = ownership;
  }
}


void
PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const int                           id_part,
        int                         **cell_vtx_idx,
        int                         **cell_vtx,
        PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;


  PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  *cell_vtx_idx = block->_cellvtx_idx[id_part];
  *cell_vtx     = block->_cellvtx    [id_part];

  // ownership
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->elt_vtx_owner != PDM_OWNERSHIP_USER) block->elt_vtx_owner = ownership;
  }
}


void
PDM_part_mesh_nodal_elmts_section_poly3d_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const int                           id_part,
        int                          *n_face,
        PDM_g_num_t                 **face_ln_to_gn,
        int                         **face_vtx_idx,
        int                         **face_vtx,
        PDM_g_num_t                 **numabs,
        int                         **cell_face_idx,
        int                         **cell_face,
        int                         **parent_num,
        PDM_g_num_t                 **parent_entity_g_num,
        PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;


  PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

  CHECK_BLOCK (block)
  CHECK_I_PART(block, id_part)

  *n_face              = block->n_face        [id_part];
  *face_vtx_idx        = block->_facvtx_idx   [id_part];
  *face_vtx            = block->_facvtx       [id_part];
  *cell_face_idx       = block->_cellfac_idx  [id_part];
  *cell_face           = block->_cellfac      [id_part];
  *numabs              = block->_numabs       [id_part];
  *face_ln_to_gn       = block->_face_ln_to_gn[id_part];
  if (block->_parent_num != NULL) {
    *parent_num = block->_parent_num[id_part];
  }
  else {
    *parent_num = NULL;
  }
  if (block->_parent_entity_g_num != NULL) {
    *parent_entity_g_num = block->_parent_entity_g_num[id_part];
  }
  else{
    *parent_entity_g_num = NULL;
  }

  // ownership
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    if (block->owner            != PDM_OWNERSHIP_USER) block->owner            = ownership;
    if (block->numabs_owner     != PDM_OWNERSHIP_USER) block->numabs_owner     = ownership;
    if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = ownership;
  }
}


PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_elmts_section_type_get
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section
)
{
  CHECK_PMNE(pmne)

  PDM_Mesh_nodal_elt_t t_elt;
  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    const PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[id_section];

    CHECK_BLOCK(block)

    t_elt = block->t_elt;
  }
  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {
    t_elt = PDM_MESH_NODAL_POLY_2D;
  }
  else {
    t_elt = PDM_MESH_NODAL_POLY_3D;
  }

  return t_elt;
}


int
PDM_part_mesh_nodal_elmts_section_n_elt_get
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part
)
{
  CHECK_PMNE(pmne)

  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    CHECK_BLOCK(block)
    CHECK_I_PART(block, id_part)

    return block->n_elt[id_part];
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    CHECK_BLOCK(block)
    CHECK_I_PART(block, id_part)

    return block->n_elt[id_part];
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    CHECK_BLOCK(block)
    CHECK_I_PART(block, id_part)

    return block->n_elt[id_part];
  }

}


int
PDM_part_mesh_nodal_elmts_n_section_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  CHECK_PMNE(pmne)
  return pmne->n_section;
}


int *
PDM_part_mesh_nodal_elmts_sections_id_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  CHECK_PMNE(pmne)
  return pmne->sections_id;
}


void
PDM_part_mesh_nodal_elmts_free
(
  PDM_part_mesh_nodal_elmts_t* pmne
)
{
  if (pmne != NULL) {
    PDM_free(pmne->n_elmts);

    /* free standard blocks */
    if (pmne->sections_std != NULL) {
      for (int i = 0; i < pmne->n_section_std; i++) {
        _block_std_free(pmne->sections_std[i]);
      }
      PDM_free(pmne->sections_std);
    }

    /* Free polygon blocks */
    if (pmne->sections_poly2d != NULL) {
      for (int i = 0; i < pmne->n_section_poly2d; i++) {
        _block_poly2d_free(pmne->sections_poly2d[i]);
      }
      PDM_free(pmne->sections_poly2d);
    }

    /* Free polyhedron blocks */
    if (pmne->sections_poly3d != NULL) {
      for (int i = 0; i < pmne->n_section_poly3d; i++) {
        _block_poly3d_free(pmne->sections_poly3d[i]);
      }
      PDM_free(pmne->sections_poly3d);
    }

    PDM_free(pmne->sections_id);

    if (pmne->numabs != NULL) {
      if(pmne->ownership_numabs == PDM_OWNERSHIP_KEEP) {
        for (int i = 0; i < pmne->n_part; i++) {
          PDM_free(pmne->numabs[i]);
        }
        PDM_free(pmne->numabs);
      }
    }

    if (pmne->num_elmt_parent_to_local != NULL) {
      for (int i_part = 0; i_part < pmne->n_part; i_part++) {
        PDM_free(pmne->num_elmt_parent_to_local[i_part]);
      }
      PDM_free(pmne->num_elmt_parent_to_local);
    }

    if (pmne->n_group_elmt != NULL) {

      for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

        for(int i_group = 0; i_group < pmne->n_group; ++i_group) {
          if(pmne->ownership_group[i_part][i_group] == PDM_OWNERSHIP_KEEP) {
            PDM_free(pmne->group_elmt    [i_part][i_group]);
            PDM_free(pmne->group_ln_to_gn[i_part][i_group]);
          }
        }
        PDM_free(pmne->n_group_elmt   [i_part]);
        PDM_free(pmne->group_elmt     [i_part]);
        PDM_free(pmne->group_ln_to_gn [i_part]);
        PDM_free(pmne->ownership_group[i_part]);
      }
      PDM_free(pmne->n_group_elmt   );
      PDM_free(pmne->group_elmt     );
      PDM_free(pmne->group_ln_to_gn );
      PDM_free(pmne->ownership_group);
    }
  }
  PDM_free(pmne);
}


int *
PDM_part_mesh_nodal_elmts_parent_num_get
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
        PDM_ownership_t              ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section;

  int *_parent_num = NULL;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->parent_num_owner != PDM_OWNERSHIP_USER) block->parent_num_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }

  return _parent_num;
}


PDM_g_num_t *
PDM_part_mesh_nodal_elmts_g_num_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
      PDM_ownership_t              ownership
)
{
  CHECK_PMNE(pmne)

  int _id_section;

  PDM_g_num_t *_g_num = NULL;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->numabs_owner != PDM_OWNERSHIP_USER) block->numabs_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_numabs != NULL) {
      _g_num = block->_numabs[id_part];
    }
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->numabs_owner != PDM_OWNERSHIP_USER) block->numabs_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_numabs != NULL) {
      _g_num = block->_numabs[id_part];
    }
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->numabs_owner != PDM_OWNERSHIP_USER) block->numabs_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_numabs != NULL) {
      _g_num = block->_numabs[id_part];
    }
  }

  return _g_num;
}


void
PDM_part_mesh_nodal_elmts_elt_extents_compute
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
  const double                       tolerance,
        double                      *vtx_coord,
        double                      *extents
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  const double eps_extents = 1.e-16;

  PDM_l_num_t *cell_vtx     = NULL;
  PDM_l_num_t *cell_vtx_idx = NULL;
  PDM_l_num_t  n_elt, n_vtx_elt = 0;

  int _id_section;

  double *lagrange_coord = NULL;
  double *bezier_coord   = NULL;
  double *matrix         = NULL;
  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_N_ELEMENT_TYPES;
  int order = -1;

  /* Polyhedra */
  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    n_elt        = block->n_elt       [id_part];
    cell_vtx_idx = block->_cellvtx_idx[id_part];
    cell_vtx     = block->_cellvtx    [id_part];
  }

  /* Polygons */
  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    n_elt        = block->n_elt      [id_part];
    cell_vtx_idx = block->_connec_idx[id_part];
    cell_vtx     = block->_connec    [id_part];
  }

  /* Standard elements */
  else {

    _id_section = id_section;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    CHECK_BLOCK(block)
    CHECK_I_PART(block, id_part)

    n_elt    = block->n_elt  [id_part];
    cell_vtx = block->_connec[id_part];

    order = block->order;

    n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get (block->t_elt, order);

    if (order > 1                                   &&
        block->t_elt != PDM_MESH_NODAL_BARHO_BEZIER &&
        block->t_elt != PDM_MESH_NODAL_TRIAHO_BEZIER) {
      t_elt = block->t_elt;
      int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);
      PDM_malloc(lagrange_coord, n_vtx_elt * 3,               double);
      PDM_malloc(bezier_coord,   n_vtx_elt * 3,               double);
      PDM_malloc(matrix,         n_nodes_quad * n_nodes_quad, double);
    }
  }

  /* Loop on elements */
  int idx = 0;
  for (PDM_l_num_t ielt = 0; ielt < n_elt; ielt++) {

    double *_extents = extents + 6 * ielt;

    for (int idim = 0; idim < 3; idim++) {
      _extents[idim]   =  HUGE_VAL;
      _extents[3+idim] = -HUGE_VAL;
    }


    if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
      idx = cell_vtx_idx[ielt];
      n_vtx_elt = cell_vtx_idx[ielt+1] - cell_vtx_idx[ielt];
    }

    double *coord = NULL;

    if (bezier_coord != NULL) {
      for (int ivtx = 0; ivtx < n_vtx_elt; ivtx++) {
        PDM_l_num_t id_vtx = cell_vtx[idx + ivtx] - 1;
        memcpy(lagrange_coord + 3*ivtx, vtx_coord + 3*id_vtx, sizeof(double) * 3);
      }

      switch (t_elt) {
      case PDM_MESH_NODAL_BARHO:
        PDM_lagrange_to_bezier_bar(order,
                                   lagrange_coord,
                                   bezier_coord,
                                   matrix);
        break;
      case PDM_MESH_NODAL_TRIAHO:
        PDM_lagrange_to_bezier_tria(order,
                                    lagrange_coord,
                                    bezier_coord,
                                    matrix);
        break;
      case PDM_MESH_NODAL_QUADHO:
        PDM_lagrange_to_bezier_quad(order,
                                    lagrange_coord,
                                    bezier_coord,
                                    matrix);
        break;
      case PDM_MESH_NODAL_TETRAHO:
        PDM_lagrange_to_bezier_tetra(order,
                                     lagrange_coord,
                                     bezier_coord,
                                     matrix);
        break;
      case PDM_MESH_NODAL_PYRAMIDHO:
        PDM_lagrange_to_bezier_pyramid(order,
                                       lagrange_coord,
                                       bezier_coord,
                                       matrix);
        break;
      case PDM_MESH_NODAL_PRISMHO:
        PDM_lagrange_to_bezier_prism(order,
                                     lagrange_coord,
                                     bezier_coord,
                                     matrix);
        break;
      case PDM_MESH_NODAL_HEXAHO:
        PDM_lagrange_to_bezier_hexa(order,
                                    lagrange_coord,
                                    bezier_coord,
                                    matrix);
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid elt type %d\n", t_elt);
      }
    }

    for (int ivtx = 0; ivtx < n_vtx_elt; ivtx++) {
      PDM_l_num_t id_vtx = cell_vtx[idx++] - 1;

      if (bezier_coord != NULL) {
        coord = bezier_coord + 3*ivtx;
      }
      else {
        coord = vtx_coord + 3*id_vtx;
      }

      for (int idim = 0; idim < 3; idim++) {
        double x = coord[idim];

        if (x < _extents[idim]) {
          _extents[idim] = x;
        }
        if (x > _extents[3+idim]) {
          _extents[3+idim] = x;
        }
      }
    }

    /* Expand bounding box */
    double l_max = 0.;
    for (int idim = 0; idim < 3; idim++) {
      double x = _extents[3+idim] - _extents[idim];

      if (l_max < x) {
        l_max = x;
      }
    }
    double delta = l_max;
    for (int idim = 0; idim < 3; idim++) {
      double x = _extents[3+idim] - _extents[idim];
      if (x < eps_extents){
        delta = l_max * PDM_MAX(tolerance, eps_extents);
      }
      else {
        delta = l_max * tolerance;
      }
      _extents[idim]   -= delta;
      _extents[3+idim] += delta;
    }

  } // End of loop on elements

  PDM_free(lagrange_coord);
  PDM_free(bezier_coord);
  PDM_free(matrix);
}


void
PDM_part_mesh_nodal_elmts_elt_center_compute
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
  const int                          n_vtx,
        double                      *vtx_coord,
  const PDM_ownership_t              ownership
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  /* Polyhedra */
  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    block->cell_centers_owner = ownership;

    if (block->cell_centers_to_compute == NULL) {
      PDM_malloc(block->cell_centers_to_compute, block->n_part, int);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers_to_compute[i] = 1;
      }
    }

    if (block->cell_centers == NULL) {
      PDM_malloc(block->cell_centers, block->n_part, double *);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers[i] = NULL;
      }
    }
    else if (!block->cell_centers_to_compute[id_part]) {
      return;
    }

    if (block->cell_centers[id_part] == NULL) {
      PDM_malloc(block->cell_centers[id_part], 3*block->n_elt[id_part], double);
    }

    double *volume         = NULL;
    double *charac_length  = NULL;
    int    *is_degenerated = NULL;
    PDM_malloc(volume        , block->n_elt[id_part], double);
    PDM_malloc(charac_length , block->n_elt[id_part], double);
    PDM_malloc(is_degenerated, block->n_elt[id_part], int   );

    PDM_geom_elem_polyhedra_properties(0,
                                       block->n_elt[id_part],
                                       block->n_face[id_part],
                                       block->_facvtx_idx[id_part],
                                       block->_facvtx[id_part],
                                       block->_cellfac_idx[id_part],
                                       block->_cellfac[id_part],
                                       n_vtx,
                                       vtx_coord,
                                       volume,
                                       block->cell_centers[id_part],
                                       charac_length,
                                       is_degenerated);
    PDM_free(volume);
    PDM_free(charac_length);
    PDM_free(is_degenerated);
  }

  /* Polygons */
  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    block->cell_centers_owner = ownership;

    if (block->cell_centers_to_compute == NULL) {
      PDM_malloc(block->cell_centers_to_compute, block->n_part, int);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers_to_compute[i] = 1;
      }
    }


    if (block->cell_centers == NULL) {
      PDM_malloc(block->cell_centers, block->n_part, double *);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers[i] = NULL;
      }
    }
    else if (!block->cell_centers_to_compute[id_part]) {
      return;
    }

    block->cell_centers_to_compute[id_part] = 0;

    if (block->cell_centers[id_part]==NULL) {
      PDM_malloc(block->cell_centers[id_part], 3*block->n_elt[id_part], double);
    }

    double *surface_vector = NULL;
    double *charac_length  = NULL;
    int    *is_degenerated = NULL;
    PDM_malloc(surface_vector, 3 * block->n_elt[id_part], double);
    PDM_malloc(charac_length ,     block->n_elt[id_part], double);
    PDM_malloc(is_degenerated,     block->n_elt[id_part], int   );

    PDM_geom_elem_polygon_properties(block->n_elt[id_part],
                                     block->_connec_idx[id_part],
                                     block->_connec[id_part],
                                     vtx_coord,
                                     surface_vector,
                                     block->cell_centers[id_part],
                                     charac_length,
                                     is_degenerated);
    PDM_free(surface_vector);
    PDM_free(charac_length);
    PDM_free(is_degenerated);
  }

  /* Standard elements */
  else {

    int _id_section = id_section;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    block->cell_centers_owner = ownership;

    if (block->cell_centers_to_compute == NULL) {
      PDM_malloc(block->cell_centers_to_compute, block->n_part, int);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers_to_compute[i] = 1;
      }
    }

    if (block->cell_centers == NULL) {
      PDM_malloc(block->cell_centers, block->n_part, double *);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers[i] = NULL;
      }
    }
    else if (!block->cell_centers_to_compute[id_part]) {
      return;
    }

    if (block->cell_centers[id_part] == NULL) {
      PDM_malloc(block->cell_centers[id_part], 3*block->n_elt[id_part], double);
    }

    double *charac_length  = NULL;
    int    *is_degenerated = NULL;
    PDM_malloc(charac_length , block->n_elt[id_part], double);
    PDM_malloc(is_degenerated, block->n_elt[id_part], int   );

    switch (block->t_elt) {
    case PDM_MESH_NODAL_POINT:
      memcpy(block->cell_centers[id_part],
             vtx_coord,
             sizeof(double)*3*(block->n_elt[id_part]) );
      break;

    case PDM_MESH_NODAL_BAR2:
    {
      double *length;
      PDM_malloc(length, block->n_elt[id_part], double);
      PDM_geom_elem_edges_properties(block->n_elt[id_part],
                                     block->_connec[id_part],
                                     vtx_coord,
                                     length,
                                     block->cell_centers[id_part],
                                     charac_length,
                                     is_degenerated);
      PDM_free(length);
      break;
    }

    case PDM_MESH_NODAL_TRIA3:
    {
      double *surface_vector;
      PDM_malloc(surface_vector, 3*block->n_elt[id_part], double);
      PDM_geom_elem_tria_properties(block->n_elt[id_part],
                                    block->_connec[id_part],
                                    vtx_coord,
                                    surface_vector,
                                    block->cell_centers[id_part],
                                    charac_length,
                                    is_degenerated);
      PDM_free(surface_vector);
      break;
    }
    case PDM_MESH_NODAL_QUAD4:
    {
      double *surface_vector;
      PDM_malloc(surface_vector, 3*block->n_elt[id_part], double);
      PDM_geom_elem_quad_properties(block->n_elt[id_part],
                                    block->_connec[id_part],
                                    vtx_coord,
                                    surface_vector,
                                    block->cell_centers[id_part],
                                    charac_length,
                                    is_degenerated);
      PDM_free(surface_vector);
      break;
    }

    case  PDM_MESH_NODAL_TETRA4:
    {
      double *volume;
      PDM_malloc(volume, block->n_elt[id_part], double);
      PDM_geom_elem_tetra_properties(block->n_elt[id_part],
                                     block->_connec[id_part],
                                     vtx_coord,
                                     volume,
                                     block->cell_centers[id_part],
                                     charac_length,
                                     is_degenerated);
      PDM_free(volume);
      break;
    }

    case PDM_MESH_NODAL_PYRAMID5:
    {
      double *volume;
      PDM_malloc(volume, block->n_elt[id_part], double);
      PDM_geom_elem_pyramid_properties(block->n_elt[id_part],
                                       block->_connec[id_part],
                                       n_vtx,
                                       vtx_coord,
                                       volume,
                                       block->cell_centers[id_part],
                                       charac_length,
                                       is_degenerated);
      PDM_free(volume);
      break;
    }

    case PDM_MESH_NODAL_PRISM6:
    {
      double *volume;
      PDM_malloc(volume, block->n_elt[id_part], double);
      PDM_geom_elem_prism_properties(block->n_elt[id_part],
                                     block->_connec[id_part],
                                     n_vtx,
                                     vtx_coord,
                                     volume,
                                     block->cell_centers[id_part],
                                     charac_length,
                                     is_degenerated);
      PDM_free(volume);
      break;
    }

    case PDM_MESH_NODAL_HEXA8:
    {
      double *volume;
      PDM_malloc(volume, block->n_elt[id_part], double);
      PDM_geom_elem_hexa_properties(block->n_elt[id_part],
                                    block->_connec[id_part],
                                    n_vtx,
                                    vtx_coord,
                                    volume,
                                    block->cell_centers[id_part],
                                    charac_length,
                                    is_degenerated);
      PDM_free(volume);
      break;
    }
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER:
    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
      PDM_error(__FILE__, __LINE__, 0, "Cell center computation not yet implemented for HO elements\n");
    case PDM_MESH_NODAL_POLY_2D:
    case PDM_MESH_NODAL_POLY_3D:
      break;
    default:
      break;
    }//end switch t_elt


    PDM_free(charac_length);
    PDM_free(is_degenerated);
  }
}


const double *
PDM_part_mesh_nodal_elmts_elt_center_get
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
        PDM_ownership_t              ownership
)
{
  CHECK_PMNE(pmne)
  CHECK_I_PART(pmne, id_part)

  double* elt_centers = NULL;
  int _id_section;

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {
    _id_section = id_section;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->cell_centers_owner            != PDM_OWNERSHIP_USER) block->cell_centers_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->cell_centers != NULL) {
      elt_centers = block->cell_centers[id_part];
    }
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->cell_centers_owner            != PDM_OWNERSHIP_USER) block->cell_centers_owner = ownership;
    }

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->cell_centers != NULL) {
      elt_centers = block->cell_centers[id_part];
    }
  }

  else {
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->cell_centers_owner            != PDM_OWNERSHIP_USER) block->cell_centers_owner = ownership;
    }

    if (block->cell_centers != NULL) {
      elt_centers = block->cell_centers[id_part];
    }
  }

  return elt_centers;
}


void
PDM_part_mesh_nodal_elmts_elt_center_reset
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part
)
{
  CHECK_PMNE(pmne)
  CHECK_I_PART(pmne, id_part)

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->cell_centers_to_compute == NULL) {
      PDM_malloc(block->cell_centers_to_compute, block->n_part, int);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers_to_compute[i] = 1;
      }
    }
    else {
      block->cell_centers_to_compute[id_part] = 1;
    }

  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->cell_centers_to_compute == NULL) {
      PDM_malloc(block->cell_centers_to_compute, block->n_part, int);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers_to_compute[i] = 1;
      }
    }
    else {
      block->cell_centers_to_compute[id_part] = 1;
    }

  }

  else {

    int _id_section = id_section;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->cell_centers_to_compute == NULL) {
      PDM_malloc(block->cell_centers_to_compute, block->n_part, int);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers_to_compute[i] = 1;
      }
    }
    else {
      block->cell_centers_to_compute[id_part] = 1;
    }

  }
}


void
PDM_part_mesh_nodal_elmts_reset
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  CHECK_PMNE(pmne)

  if (pmne->sections_std != NULL) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      _block_std_free(pmne->sections_std[i]);
    }
    PDM_free(pmne->sections_std);
  }

  if (pmne->sections_poly2d != NULL) {
    for (int i = 0; i < pmne->n_section_poly2d; i++) {
      _block_poly2d_free(pmne->sections_poly2d[i]);
    }
    PDM_free(pmne->sections_poly2d);
  }

  if (pmne->sections_poly3d != NULL) {
    for (int i = 0; i < pmne->n_section_poly3d; i++) {
      _block_poly3d_free(pmne->sections_poly3d[i]);
    }
    PDM_free(pmne->sections_poly3d);
  }

  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  pmne->sections_std    = NULL;
  pmne->sections_poly2d = NULL;
  pmne->sections_poly3d = NULL;
  PDM_free(pmne->sections_id);
  pmne->n_section    = 0;
  pmne->prepa_blocks = NULL;

  if (pmne->num_elmt_parent_to_local != NULL) {
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      if (pmne->num_elmt_parent_to_local[i_part] != NULL)
        PDM_free(pmne->num_elmt_parent_to_local[i_part]);
    }
    PDM_free(pmne->num_elmt_parent_to_local);
    pmne->num_elmt_parent_to_local = NULL;
  }

  for (int i = 0; i < pmne->n_part; i++) {
    pmne->n_elmts[i] = 0;
  }
}


void
PDM_part_mesh_nodal_elmts_g_num_in_section_compute
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)

  PDM_gen_gnum_t *gnum_gen = PDM_gnum_create (3,
                                              pmne->n_part,
                                              PDM_FALSE,
                                              1e-3,
                                              pmne->comm,
                                              PDM_OWNERSHIP_USER); /* The result is getted and you are owner */

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    CHECK_BLOCK(block)

    block->numabs_int_owner = ownership;

    if (block->numabs_int == NULL) {
      PDM_malloc(block->numabs_int, block->n_part, PDM_g_num_t *);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      PDM_gnum_free(gnum_gen);
      return;
    }

    for (int i = 0; i < block->n_part; i++) {
      PDM_gnum_set_from_parents(gnum_gen, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;


    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    CHECK_BLOCK(block)

    block->numabs_int_owner = ownership;

    if (block->numabs_int == NULL) {
      PDM_malloc(block->numabs_int, block->n_part, PDM_g_num_t *);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      PDM_gnum_free(gnum_gen);
      return;
    }

    for (int i = 0; i < block->n_part; i++) {
      PDM_gnum_set_from_parents(gnum_gen, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  else {

    int _id_section = id_section;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    CHECK_BLOCK(block)

    block->numabs_int_owner = ownership;

    if (block->numabs_int == NULL) {
      PDM_malloc(block->numabs_int, block->n_part, PDM_g_num_t *);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      PDM_gnum_free(gnum_gen);
      return;
    }

    for (int i = 0; i < block->n_part; i++) {
      PDM_gnum_set_from_parents(gnum_gen, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  PDM_gnum_compute (gnum_gen);

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    for (int i = 0; i < block->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get(gnum_gen, i);
    }
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    for (int i = 0; i < block->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get(gnum_gen, i);
    }
  }

  else {

    int _id_section = id_section;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    for (int i = 0; i < block->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get(gnum_gen, i);
    }
  }

  PDM_gnum_free(gnum_gen);
}


int
PDM_part_mesh_nodal_elmts_n_elmts_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_part
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  return pmne->n_elmts[id_part];
}


PDM_g_num_t *
PDM_part_mesh_nodal_elmts_g_num_get_from_part
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_part,
        PDM_ownership_t               ownership
)
{
  // Mandatory
  if (pmne->n_section == 0 || pmne->n_part == 0) {
    return NULL;
  }

  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  if (pmne->numabs == NULL) {
    PDM_malloc(pmne->numabs, pmne->n_part, PDM_g_num_t*);
    // Safe because we early exit if n_section == 0
    int is_not_parent_num = (PDM_part_mesh_nodal_elmts_parent_num_get(pmne, pmne->sections_id[0], 0 /*i_part*/, PDM_OWNERSHIP_KEEP) == NULL);
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      // Check
      for (int i_section = 0; i_section < pmne->n_section; i_section++) {
        int lis_not_parent_num = (PDM_part_mesh_nodal_elmts_parent_num_get(pmne, pmne->sections_id[i_section], i_part, PDM_OWNERSHIP_KEEP) == NULL);
        if(is_not_parent_num != lis_not_parent_num) {
          PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_elmts_g_num_get_from_part have strange mix of parent_num : is_not_parent_num = %i / current_section = %i (lis_not_parent_num=%i) \n",
                    is_not_parent_num, i_section, lis_not_parent_num);
        }
      }
    }

    if (is_not_parent_num) {

      for (int i = 0; i < pmne->n_part; i++) {
        int k = 0;
        PDM_malloc(pmne->numabs[i], pmne->n_elmts[i], PDM_g_num_t);
        for (int i1 = 0; i1 < pmne->n_section_std; i1++) {
          for (int i2 = 0; i2 < pmne->sections_std[i1]->n_elt[i]; i2++) {
            pmne->numabs[i][k++] = pmne->sections_std[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < pmne->n_section_poly2d; i1++) {
          for (int i2 = 0; i2 < pmne->sections_poly2d[i1]->n_elt[i]; i2++) {
            pmne->numabs[i][k++] = pmne->sections_poly2d[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < pmne->n_section_poly3d; i1++) {
          for (int i2 = 0; i2 < pmne->sections_poly3d[i1]->n_elt[i]; i2++) {
            pmne->numabs[i][k++] = pmne->sections_poly3d[i1]->_numabs[i][i2];
          }
        }
      }
    }

    else {
      for (int i = 0; i < pmne->n_part; i++) {
        PDM_malloc(pmne->numabs[i], pmne->n_elmts[i], PDM_g_num_t);
        for (int i1 = 0; i1 < pmne->n_section_std; i1++) {
          for (int i2 = 0; i2 < pmne->sections_std[i1]->n_elt[i]; i2++) {
            pmne->numabs[i][pmne->sections_std[i1]->_parent_num[i][i2]] = pmne->sections_std[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < pmne->n_section_poly2d; i1++) {
          for (int i2 = 0; i2 < pmne->sections_poly2d[i1]->n_elt[i]; i2++) {
            pmne->numabs[i][pmne->sections_poly2d[i1]->_parent_num[i][i2]] = pmne->sections_poly2d[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < pmne->n_section_poly3d; i1++) {
          for (int i2 = 0; i2 < pmne->sections_poly3d[i1]->n_elt[i]; i2++) {
            pmne->numabs[i][pmne->sections_poly3d[i1]->_parent_num[i][i2]] = pmne->sections_poly3d[i1]->_numabs[i][i2];
          }
        }
      }
    }
  }

  // ownership
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    if (pmne->ownership_numabs != PDM_OWNERSHIP_USER) pmne->ownership_numabs = ownership;
  }

  return pmne->numabs[id_part];
}


void
PDM_part_mesh_nodal_elmts_partial_free
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{

  if (pmne->sections_std != NULL) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      _block_std_free_partial(pmne->sections_std[i]);
    }
  }

  if (pmne->sections_poly2d != NULL) {
    for (int i = 0; i < pmne->n_section_poly2d; i++) {
      _block_poly2d_free_partial(pmne->sections_poly2d[i]);
    }
  }

  if (pmne->sections_poly3d != NULL) {
    for (int i = 0; i < pmne->n_section_poly3d; i++) {
      _block_poly3d_free_partial(pmne->sections_poly3d[i]);
    }
  }
}


PDM_g_num_t *
PDM_part_mesh_nodal_elmts_section_g_num_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_section,
  const int                           id_part,
        PDM_ownership_t               ownership
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->numabs_int_owner != PDM_OWNERSHIP_USER) block->numabs_int_owner = ownership;
    }

    return block->numabs_int[id_part];
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->numabs_int_owner != PDM_OWNERSHIP_USER) block->numabs_int_owner = ownership;
    }

    return block->numabs_int[id_part];
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    // ownership
    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      if (block->numabs_int_owner != PDM_OWNERSHIP_USER) block->numabs_int_owner = ownership;
    }

    return block->numabs_int[id_part];
  }
}


int *
PDM_part_mesh_nodal_elmts_num_elmt_parent_to_local_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_part
)
{
  CHECK_PMNE(pmne)
  CHECK_I_PART(pmne, id_part)

  if (pmne->num_elmt_parent_to_local != NULL) {
    return pmne->num_elmt_parent_to_local[id_part];
  }
  else {
    return NULL;
  }
}


void
PDM_part_mesh_elmts_nodal_cell3d_cellface_add
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_part,
  const int                           n_cell,
  const int                           n_face,
  const int                          *face_vtx_idx,
  const int                          *face_vtx,
  const PDM_g_num_t                  *face_ln_to_gn,
  const int                          *cell_face_idx,
  const int                          *cell_face,
  const PDM_g_num_t                  *cell_ln_to_gn,
        PDM_Mesh_nodal_vtx_t        **vtx,
  const PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)
  CHECK_I_PART(pmne, id_part)

  int adjust = 0;
  if (n_cell > 0) {
    if (cell_face_idx[0] == 1) {
      adjust = 1;
    }
  }

  int n_part = 0;

  if (pmne->num_elmt_parent_to_local == NULL) {
    PDM_malloc(pmne->num_elmt_parent_to_local, pmne->n_part, PDM_l_num_t *);
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      pmne->num_elmt_parent_to_local[i_part] = NULL;
    }
  }

  PDM_malloc(pmne->num_elmt_parent_to_local[id_part], n_cell, PDM_l_num_t);
  for (int i = 0; i < n_cell; i++) {
    pmne->num_elmt_parent_to_local[id_part][i] = 0;
  }

  if (pmne->prepa_blocks == NULL) {
    PDM_malloc(pmne->prepa_blocks, 1, PDM_Mesh_nodal_prepa_blocks_t);
    pmne->prepa_blocks->t_add = 1;
    pmne->prepa_blocks->n_tria_proc    = 0;  /* Nb de triangles par proc */
    pmne->prepa_blocks->n_quad_proc    = 0;  /* Nb de quads par proc     */
    pmne->prepa_blocks->n_poly2d_proc  = 0;  /* Nb de poly2d par proc    */
    pmne->prepa_blocks->n_tetra_proc   = 0;  /* Nb de tetra par proc     */
    pmne->prepa_blocks->n_hexa_proc    = 0;  /* Nb d'hexa par proc       */
    pmne->prepa_blocks->n_prism_proc   = 0;  /* Nb de prisme par proc    */
    pmne->prepa_blocks->n_pyramid_proc = 0;  /* Nb de pyramide par proc  */
    pmne->prepa_blocks->n_poly3d_proc  = 0;  /* Nb de poly3d par proc    */

    PDM_malloc(pmne->prepa_blocks->n_cell,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_face,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_tetra,       pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_hexa,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_prism,       pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_pyramid,     pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_poly3d,      pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->face_vtx_idx,  pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->face_vtx,      pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->cell_face_idx, pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->cell_face,     pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->add_etat,      pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->numabs,        pmne->n_part, PDM_g_num_t *);
    PDM_malloc(pmne->prepa_blocks->face_ln_to_gn, pmne->n_part, PDM_g_num_t *);
    for (int i = 0; i < pmne->n_part; i++) {
      pmne->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (pmne->prepa_blocks->t_add != 1) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_elmts_nodal_cell3d_cellface_add : Another type of elements is currently is still in progress \n");
  }

  /* Determination du type de chaque element */

  PDM_l_num_t cell_som_tria[18]; /* 6 triangles max in _type_cell_3D */
  PDM_l_num_t cell_som_quad[24]; /* 6 quadrangles max in _type_cell_3D */
  PDM_l_num_t n_tetra   = 0;
  PDM_l_num_t n_hexa    = 0;
  PDM_l_num_t n_prism   = 0;
  PDM_l_num_t n_pyramid = 0;
  PDM_l_num_t n_poly3d  = 0;

  for (int i = 0; i < n_cell; i++) {

    PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_idx[i+1] - cell_face_idx[i],
                                                   cell_face + cell_face_idx[i] - adjust,
                                                   face_vtx_idx,
                                                   face_vtx,
                                                   cell_som_tria,
                                                   cell_som_quad);
    switch(cell_type) {
    case PDM_MESH_NODAL_TETRA4 :
      n_tetra += 1;
      break;
    case PDM_MESH_NODAL_PYRAMID5 :
      n_pyramid += 1;
      break;
    case PDM_MESH_NODAL_PRISM6 :
      n_prism += 1;
      break;
    case PDM_MESH_NODAL_HEXA8 :
      n_hexa += 1;
      break;
    case PDM_MESH_NODAL_POLY_3D :
      n_poly3d += 1;
      break;
    default :
      break;
    }
  }

  pmne->prepa_blocks->n_tetra_proc          += n_tetra;
  pmne->prepa_blocks->n_hexa_proc           += n_hexa;
  pmne->prepa_blocks->n_prism_proc          += n_prism;
  pmne->prepa_blocks->n_pyramid_proc        += n_pyramid;
  pmne->prepa_blocks->n_poly3d_proc         += n_poly3d;
  pmne->prepa_blocks->n_tetra      [id_part] = n_tetra;
  pmne->prepa_blocks->n_hexa       [id_part] = n_hexa;
  pmne->prepa_blocks->n_prism      [id_part] = n_prism;
  pmne->prepa_blocks->n_pyramid    [id_part] = n_pyramid;
  pmne->prepa_blocks->n_poly3d     [id_part] = n_poly3d;
  pmne->prepa_blocks->face_vtx_idx [id_part] = (PDM_l_num_t *) face_vtx_idx;
  pmne->prepa_blocks->face_vtx     [id_part] = (PDM_l_num_t *) face_vtx;
  pmne->prepa_blocks->cell_face_idx[id_part] = (PDM_l_num_t *) cell_face_idx;
  pmne->prepa_blocks->cell_face    [id_part] = (PDM_l_num_t *) cell_face;
  pmne->prepa_blocks->numabs       [id_part] = (PDM_g_num_t *) cell_ln_to_gn;
  pmne->prepa_blocks->face_ln_to_gn[id_part] = (PDM_g_num_t *) face_ln_to_gn;
  pmne->prepa_blocks->add_etat     [id_part] = 1;
  pmne->prepa_blocks->n_face       [id_part] = n_face;
  pmne->prepa_blocks->n_cell       [id_part] = n_cell;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < pmne->n_part; i++) {
    if (pmne->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (pmne->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[5];
    PDM_l_num_t som_elts[5];

    elts[0] = pmne->prepa_blocks->n_tetra_proc   > 0;
    elts[1] = pmne->prepa_blocks->n_hexa_proc    > 0;
    elts[2] = pmne->prepa_blocks->n_prism_proc   > 0;
    elts[3] = pmne->prepa_blocks->n_pyramid_proc > 0;
    elts[4] = pmne->prepa_blocks->n_poly3d_proc  > 0;

    PDM_MPI_Allreduce(elts, som_elts, 5, PDM_MPI_INT, PDM_MPI_SUM, pmne->comm);

    int id_bloc_tetra4   = -1;
    int id_bloc_hexa8    = -1;
    int id_bloc_prism6   = -1;
    int id_bloc_pyramid5 = -1;
    int id_bloc_poly_3d  = -1;

    if (som_elts[0] > 0) {
      id_bloc_tetra4 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_TETRA4);
    }

    if (som_elts[1] > 0) {
      id_bloc_hexa8 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_HEXA8);
    }

    if (som_elts[2] > 0) {
      id_bloc_prism6 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_PRISM6);
    }

    if (som_elts[3] > 0) {
      id_bloc_pyramid5 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_PYRAMID5);
    }

    if (som_elts[4] > 0) {
      id_bloc_poly_3d = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_POLY_3D);
    }

    /* Determination de la connectivite de chaque element */


    for (int i_part = 0; i_part < pmne->n_part; i_part++) {

      assert(vtx[i_part] != NULL);
      double *vtx_coord = vtx[i_part]->_coords;
      assert(vtx_coord != NULL);

      PDM_l_num_t n_cell_courant                    = pmne->prepa_blocks->n_cell       [i_part];
      PDM_l_num_t *num_cell_parent_to_local_courant = pmne->num_elmt_parent_to_local   [i_part];
      PDM_l_num_t *face_som_idx_courant             = pmne->prepa_blocks->face_vtx_idx [i_part];
      PDM_l_num_t *face_som_courant                 = pmne->prepa_blocks->face_vtx     [i_part];
      PDM_l_num_t *cell_face_idx_courant            = pmne->prepa_blocks->cell_face_idx[i_part];
      PDM_l_num_t *cell_face_courant                = pmne->prepa_blocks->cell_face    [i_part];
      PDM_g_num_t *numabs_courant                   = pmne->prepa_blocks->numabs       [i_part];
      PDM_l_num_t n_face_part                       = pmne->prepa_blocks->n_face       [i_part];

      PDM_l_num_t n_tetra_part   = pmne->prepa_blocks->n_tetra  [i_part];
      PDM_l_num_t n_hexa_part    = pmne->prepa_blocks->n_hexa   [i_part];
      PDM_l_num_t n_prism_part   = pmne->prepa_blocks->n_prism  [i_part];
      PDM_l_num_t n_pyramid_part = pmne->prepa_blocks->n_pyramid[i_part];
      PDM_l_num_t n_poly3d_part  = pmne->prepa_blocks->n_poly3d [i_part];

      PDM_l_num_t *connec_tetra   = NULL;
      PDM_l_num_t *connec_hexa    = NULL;
      PDM_l_num_t *connec_prism   = NULL;
      PDM_l_num_t *connec_pyramid = NULL;

      PDM_g_num_t *numabs_tetra   = NULL;
      PDM_g_num_t *numabs_hexa    = NULL;
      PDM_g_num_t *numabs_prism   = NULL;
      PDM_g_num_t *numabs_pyramid = NULL;
      PDM_g_num_t *numabs_poly3d  = NULL;

      PDM_l_num_t *num_parent_tetra   = NULL;
      PDM_l_num_t *num_parent_hexa    = NULL;
      PDM_l_num_t *num_parent_prism   = NULL;
      PDM_l_num_t *num_parent_pyramid = NULL;
      PDM_l_num_t *num_parent_poly3d  = NULL;

      adjust = 0;
      if (n_cell_courant > 0) {
        if (cell_face_idx_courant[0] == 1) {
          adjust = 1;
        }
      }

      // n_tetra_part > 0
      if (som_elts[0] > 0) {
        PDM_malloc(connec_tetra    , 4 * n_tetra_part, PDM_l_num_t);
        PDM_malloc(numabs_tetra    ,     n_tetra_part, PDM_g_num_t);
        PDM_malloc(num_parent_tetra,     n_tetra_part, PDM_l_num_t);
      }

      // n_hexa_part > 0
      if (som_elts[1] > 0) {
        PDM_malloc(connec_hexa    , 8 * n_hexa_part, PDM_l_num_t);
        PDM_malloc(numabs_hexa    ,     n_hexa_part, PDM_g_num_t);
        PDM_malloc(num_parent_hexa,     n_hexa_part, PDM_l_num_t);
      }

      // n_prism_part > 0
      if (som_elts[2] > 0) {
        PDM_malloc(connec_prism    , 6 * n_prism_part, PDM_l_num_t);
        PDM_malloc(numabs_prism    ,     n_prism_part, PDM_g_num_t);
        PDM_malloc(num_parent_prism,     n_prism_part, PDM_l_num_t);
      }

      // n_pyramid_part > 0
      if (som_elts[3] > 0) {
        PDM_malloc(connec_pyramid    , 5 * n_pyramid_part, PDM_l_num_t);
        PDM_malloc(numabs_pyramid    ,     n_pyramid_part, PDM_g_num_t);
        PDM_malloc(num_parent_pyramid,     n_pyramid_part, PDM_l_num_t);
      }

      // n_poly3d_part > 0
      if (som_elts[4] > 0) {
        PDM_malloc(numabs_poly3d    , n_poly3d_part, PDM_g_num_t);
        PDM_malloc(num_parent_poly3d, n_poly3d_part, PDM_l_num_t);
      }

      PDM_l_num_t *num_parent_tetra_courant   = num_parent_tetra;
      PDM_l_num_t *num_parent_hexa_courant    = num_parent_hexa;
      PDM_l_num_t *num_parent_prism_courant   = num_parent_prism;
      PDM_l_num_t *num_parent_pyramid_courant = num_parent_pyramid;
      PDM_l_num_t *num_parent_poly3d_courant  = num_parent_poly3d;

      PDM_l_num_t *connec_tetra_courant   = connec_tetra;
      PDM_l_num_t *connec_hexa_courant    = connec_hexa;
      PDM_l_num_t *connec_prism_courant   = connec_prism;
      PDM_l_num_t *connec_pyramid_courant = connec_pyramid;

      PDM_g_num_t *numabs_tetra_courant   = numabs_tetra;
      PDM_g_num_t *numabs_hexa_courant    = numabs_hexa;
      PDM_g_num_t *numabs_prism_courant   = numabs_prism;
      PDM_g_num_t *numabs_pyramid_courant = numabs_pyramid;
      PDM_g_num_t *numabs_poly3d_courant  = numabs_poly3d;

      PDM_l_num_t *tag_face_poly3d     = NULL;
      PDM_l_num_t  n_face_poly         = 0;
      PDM_l_num_t *facsom_poly_idx     = NULL;
      PDM_l_num_t *facsom_poly         = NULL;
      PDM_l_num_t *cellfac_poly_idx    = NULL;
      PDM_l_num_t *cellfac_poly        = NULL;
      PDM_l_num_t  l_cellfac_poly      = 0;
      PDM_g_num_t *block_face_ln_to_gn = NULL;

      if (n_poly3d_part > 0) {
        PDM_malloc(tag_face_poly3d, n_face_part, PDM_l_num_t);
        for (int i = 0; i < n_face_part; i++) {
          tag_face_poly3d[i] = -1;
        }
        PDM_malloc(cellfac_poly_idx, n_poly3d_part + 1, PDM_l_num_t);
        cellfac_poly_idx[0] = 0;
      }

      PDM_l_num_t idx_tetra   = 0;
      PDM_l_num_t idx_hexa    = idx_tetra   + n_tetra_part;
      PDM_l_num_t idx_prism   = idx_hexa    + n_hexa_part;
      PDM_l_num_t idx_pyramid = idx_prism   + n_prism_part;
      PDM_l_num_t idx_poly3d  = idx_pyramid + n_pyramid_part;

      n_poly3d_part = 0;
      for (int i = 0; i < n_cell_courant; i++) {
        num_cell_parent_to_local_courant[i] = 0;
        PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_idx_courant[i+1] - cell_face_idx_courant[i],
                                                       cell_face_courant + cell_face_idx_courant[i] - adjust,
                                                       face_som_idx_courant,
                                                       face_som_courant,
                                                       cell_som_tria,
                                                       cell_som_quad);

        switch(cell_type) {
        case PDM_MESH_NODAL_TETRA4 :
          _connec_tetra(vtx_coord,
                        cell_som_tria,
                        connec_tetra_courant);
          *numabs_tetra_courant = numabs_courant[i];
          numabs_tetra_courant += 1;
          connec_tetra_courant += 4;
          *num_parent_tetra_courant = i;
          num_parent_tetra_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tetra++;
          break;
        case PDM_MESH_NODAL_HEXA8 :
          _connec_hexa(vtx_coord,
                       cell_som_quad,
                       connec_hexa_courant);
          *numabs_hexa_courant = numabs_courant[i];
          numabs_hexa_courant += 1;
          connec_hexa_courant += 8;
          *num_parent_hexa_courant = i;
          num_parent_hexa_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_hexa++;
          break;
        case PDM_MESH_NODAL_PRISM6 :
          _connec_prism(vtx_coord,
                        cell_som_tria,
                        cell_som_quad,
                        connec_prism_courant);
          *numabs_prism_courant = numabs_courant[i];
          numabs_prism_courant += 1;
          connec_prism_courant += 6;
          *num_parent_prism_courant = i;
          num_parent_prism_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_prism++;
          break;
        case PDM_MESH_NODAL_PYRAMID5 :
          _connec_pyramid(vtx_coord,
                          cell_som_tria,
                          cell_som_quad,
                          connec_pyramid_courant);
          *numabs_pyramid_courant = numabs_courant[i];
          numabs_pyramid_courant += 1;
          connec_pyramid_courant += 5;
          *num_parent_pyramid_courant = i;
          num_parent_pyramid_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_pyramid++;
          break;
        case PDM_MESH_NODAL_POLY_3D :
          {
            PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - adjust;
            for (int j = 0; j < cell_face_idx_courant[i+1] - cell_face_idx_courant[i]; j++) {
              tag_face_poly3d[PDM_ABS(cell_face_cell[j]) - 1] = 0;
            }
            *numabs_poly3d_courant = numabs_courant[i];
            numabs_poly3d_courant += 1;
            l_cellfac_poly += cell_face_idx_courant[i+1] - cell_face_idx_courant[i];
            cellfac_poly_idx[n_poly3d_part+1] = l_cellfac_poly;
            n_poly3d_part += 1;
            *num_parent_poly3d_courant = i;
            num_parent_poly3d_courant += 1;
            num_cell_parent_to_local_courant[i] = idx_poly3d++;
            break;
          }
        default :
          break;
        }
      }

      if (n_poly3d_part > 0) {
        PDM_malloc(cellfac_poly, l_cellfac_poly, PDM_l_num_t);

        /* Stockage des faces du bloc */

        n_face_poly = 0;
        PDM_l_num_t l_facsom_poly = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] == 0) {
            tag_face_poly3d[i] = n_face_poly++;
            l_facsom_poly += face_som_idx_courant[i+1] - face_som_idx_courant[i];
          }
        }

        PDM_malloc(facsom_poly_idx, n_face_poly + 1, PDM_l_num_t);
        PDM_malloc(facsom_poly    , l_facsom_poly  , PDM_l_num_t);
        if (pmne->prepa_blocks->face_ln_to_gn[i_part] != NULL) {
          PDM_malloc(block_face_ln_to_gn,n_face_poly,PDM_g_num_t);
        }

        facsom_poly_idx[0] = 0;
        PDM_l_num_t idx_facsom_poly = 0;
        PDM_l_num_t idx_facsom = 0;
        n_face_poly = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] >= 0) {
            if (pmne->prepa_blocks->face_ln_to_gn[i_part] != NULL) {
              block_face_ln_to_gn[n_face_poly++] = pmne->prepa_blocks->face_ln_to_gn[i_part][i];
            }
            PDM_l_num_t ideb = face_som_idx_courant[i] - adjust;
            PDM_l_num_t ifin = ideb + face_som_idx_courant[i+1] - face_som_idx_courant[i];
            facsom_poly_idx[idx_facsom+1] = facsom_poly_idx[idx_facsom] + face_som_idx_courant[i+1] - face_som_idx_courant[i];
            idx_facsom += 1;
            for (int j = ideb; j < ifin; j++) {
              facsom_poly[idx_facsom_poly++] = face_som_courant[j];
            }
          }
        }

        /* Remplissage de la structure cellfac_poly */

        l_cellfac_poly = 0;
        for (int i = 0; i < n_cell_courant; i++) {
          PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_idx_courant[i+1] - cell_face_idx_courant[i],
                                                         cell_face_courant + cell_face_idx_courant[i] - adjust,
                                                         face_som_idx_courant,
                                                         face_som_courant,
                                                         cell_som_tria,
                                                         cell_som_quad);

          switch(cell_type) {

          case PDM_MESH_NODAL_POLY_3D :
            {
              PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - adjust;
              for (int j = 0; j < cell_face_idx_courant[i+1]-cell_face_idx_courant[i]; j++) {
                cellfac_poly[l_cellfac_poly++] = tag_face_poly3d[PDM_ABS(cell_face_cell[j]) - 1] + 1;

                if (cell_face_cell[j] < 0) {
                  cellfac_poly[l_cellfac_poly-1] = -cellfac_poly[l_cellfac_poly-1];
                }

              }
              break;
            }
          default:
            break;
          }
        }
        PDM_free(tag_face_poly3d);
      }

      if (som_elts[0] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_tetra4,
                                          i_part,
                                          n_tetra_part,
                                          connec_tetra,
                                          numabs_tetra,
                                          num_parent_tetra,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[1] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_hexa8,
                                          i_part,
                                          n_hexa_part,
                                          connec_hexa,
                                          numabs_hexa,
                                          num_parent_hexa,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[2] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_prism6,
                                          i_part,
                                          n_prism_part,
                                          connec_prism,
                                          numabs_prism,
                                          num_parent_prism,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[3] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_pyramid5,
                                          i_part,
                                          n_pyramid_part,
                                          connec_pyramid,
                                          numabs_pyramid,
                                          num_parent_pyramid,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[4] > 0) {
        PDM_part_mesh_nodal_elmts_section_poly3d_set(pmne,
                                                     id_bloc_poly_3d,
                                                     i_part,
                                                     n_poly3d_part,
                                                     n_face_poly,
                                                     facsom_poly_idx,
                                                     facsom_poly,
                                                     block_face_ln_to_gn,
                                                     cellfac_poly_idx,
                                                     cellfac_poly,
                                                     numabs_poly3d,
                                                     num_parent_poly3d,
                                                     NULL,//parent_entity_g_num,
                                                     ownership);
        // PDM_log_trace_array_int(num_parent_poly3d, n_poly3d_part, "num_parent_poly3d ::");
      }
    }

    if (pmne->prepa_blocks != NULL) {
      PDM_free(pmne->prepa_blocks->n_cell);
      PDM_free(pmne->prepa_blocks->n_face);
      PDM_free(pmne->prepa_blocks->n_tetra);
      PDM_free(pmne->prepa_blocks->n_hexa);
      PDM_free(pmne->prepa_blocks->n_prism);
      PDM_free(pmne->prepa_blocks->n_pyramid);
      PDM_free(pmne->prepa_blocks->n_poly3d);
      PDM_free(pmne->prepa_blocks->face_vtx_idx);
      PDM_free(pmne->prepa_blocks->face_vtx);
      PDM_free(pmne->prepa_blocks->cell_face_idx);
      PDM_free(pmne->prepa_blocks->cell_face);
      PDM_free(pmne->prepa_blocks->add_etat);
      PDM_free(pmne->prepa_blocks->numabs);
      PDM_free(pmne->prepa_blocks->face_ln_to_gn);
      PDM_free(pmne->prepa_blocks);
      pmne->prepa_blocks = NULL;
    }
  }
}


void
PDM_part_mesh_nodal_elmts_face2d_faceedge_add
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_part,
  const int                          n_face,
  const int                          n_edge,
  const int                         *edge_vtx,
  const int                         *face_edge_idx,
  const int                         *face_edge,
  const PDM_g_num_t                 *face_ln_to_gn,
  const int                          n_vtx,
  const PDM_ownership_t              ownership
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  if (pmne->num_elmt_parent_to_local == NULL) {
    PDM_malloc(pmne->num_elmt_parent_to_local, pmne->n_part, PDM_l_num_t *);
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      pmne->num_elmt_parent_to_local[i_part] = NULL;
    }
  }

  PDM_malloc(pmne->num_elmt_parent_to_local[id_part], n_face, PDM_l_num_t);
  for (int i_face = 0; i_face < n_face; i_face++) {
    pmne->num_elmt_parent_to_local[id_part][i_face] = 0;
  }

  if (pmne->prepa_blocks == NULL) {
    PDM_malloc(pmne->prepa_blocks, 1, PDM_Mesh_nodal_prepa_blocks_t);
    pmne->prepa_blocks->t_add         = 2; // From face->edge
    pmne->prepa_blocks->n_tria_proc   = 0;
    pmne->prepa_blocks->n_quad_proc   = 0;
    pmne->prepa_blocks->n_poly2d_proc = 0;
    PDM_malloc(pmne->prepa_blocks->n_cell,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_face,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_vtx,           pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_tria,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_quad,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_poly2d,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->l_connec_poly2d, pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->face_vtx,        pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->cell_face_idx,   pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->cell_face,       pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->add_etat,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->numabs,          pmne->n_part, PDM_g_num_t *);
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      pmne->prepa_blocks->add_etat[i_part] = 0;
    }
  }

  if (pmne->prepa_blocks->t_add != 2) {
    PDM_error(__FILE__, __LINE__, 0,
              "Error in PDM_part_mesh_nodal_elmts_face2d_faceedge_add : prepa_blocks already used for another type of connectivity (%d)\n",
              pmne->prepa_blocks->t_add);
  }

  /* Count number of elements of each type */
  PDM_l_num_t n_tria = 0;
  PDM_l_num_t n_quad = 0;
  PDM_l_num_t n_poly = 0;
  PDM_l_num_t l_connec_poly = 0;

  for (int i_face = 0; i_face < n_face; i_face++) {
    PDM_l_num_t _n_edge = face_edge_idx[i_face+1] - face_edge_idx[i_face];
    if (_n_edge == 3) {
      n_tria++;
    }
    else if (_n_edge == 4) {
      n_quad++;
    }
    else if (_n_edge > 4) {
      n_poly++;
      l_connec_poly += _n_edge;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Invalid 2D element with only %d edge(s)\n", _n_edge);
    }
  }

  /* Setup 'prepa' struct */
  pmne->prepa_blocks->n_tria_proc              += n_tria;
  pmne->prepa_blocks->n_quad_proc              += n_quad;
  pmne->prepa_blocks->n_poly2d_proc            += n_poly;
  pmne->prepa_blocks->add_etat       [id_part] = 1;
  pmne->prepa_blocks->n_cell         [id_part] = n_face;
  pmne->prepa_blocks->n_face         [id_part] = n_edge;
  pmne->prepa_blocks->n_vtx          [id_part] = n_vtx;
  pmne->prepa_blocks->n_tria         [id_part] = n_tria;
  pmne->prepa_blocks->n_quad         [id_part] = n_quad;
  pmne->prepa_blocks->n_poly2d       [id_part] = n_poly;
  pmne->prepa_blocks->l_connec_poly2d[id_part] = l_connec_poly;
  pmne->prepa_blocks->face_vtx       [id_part] = (PDM_l_num_t *) edge_vtx;
  pmne->prepa_blocks->cell_face_idx  [id_part] = (PDM_l_num_t *) face_edge_idx;
  pmne->prepa_blocks->cell_face      [id_part] = (PDM_l_num_t *) face_edge;
  pmne->prepa_blocks->numabs         [id_part] = (PDM_g_num_t *) face_ln_to_gn;

  /* If all parts have not already been set, end here */
  int n_part_set = 0;
  for (int i_part = 0; i_part < pmne->n_part; i_part++) {
    n_part_set += pmne->prepa_blocks->add_etat[i_part];
  }

  if (n_part_set < pmne->n_part) {
    return;
  }

  /* Otherwise, build nodal connectivities */
  int i_have_elts[3];
  i_have_elts[0] = pmne->prepa_blocks->n_tria_proc   > 0;
  i_have_elts[1] = pmne->prepa_blocks->n_quad_proc   > 0;
  i_have_elts[2] = pmne->prepa_blocks->n_poly2d_proc > 0;

  int we_have_elts[3];
  PDM_MPI_Allreduce(i_have_elts, we_have_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, pmne->comm);


  // Add sections
  int id_tria;
  int id_quad;
  int id_poly;
  if (we_have_elts[0] > 0) {
    id_tria = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_TRIA3);
  }
  if (we_have_elts[1] > 0) {
    id_quad = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_QUAD4);
  }
  if (we_have_elts[2] > 0) {
    id_poly = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_POLY_2D);
  }


  /* Determine whether face_edge connectivity is signed */
  int _is_signed = 0;
  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    PDM_l_num_t  _n_face        = pmne->prepa_blocks->n_cell       [i_part];
    PDM_l_num_t *_face_edge_idx = pmne->prepa_blocks->cell_face_idx[i_part];
    PDM_l_num_t *_face_edge     = pmne->prepa_blocks->cell_face    [i_part];

    for (int i = 0; i < _face_edge_idx[_n_face]; i++) {
      if (_face_edge[i] < 0) {
        _is_signed = 1;
        break;
      }
    }
  }

  int is_signed;
  PDM_MPI_Allreduce(&_is_signed, &is_signed, 1, PDM_MPI_INT, PDM_MPI_MAX, pmne->comm);



  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    PDM_l_num_t  _n_face          = pmne->prepa_blocks->n_cell       [i_part];
    PDM_l_num_t *_parent_to_local = pmne->num_elmt_parent_to_local   [i_part];
    PDM_l_num_t *_edge_vtx        = pmne->prepa_blocks->face_vtx     [i_part];
    PDM_l_num_t *_face_edge_idx   = pmne->prepa_blocks->cell_face_idx[i_part];
    PDM_l_num_t *_face_edge       = pmne->prepa_blocks->cell_face    [i_part];
    PDM_g_num_t *_face_ln_to_gn   = pmne->prepa_blocks->numabs       [i_part];

    PDM_l_num_t adjust = 0;
    if (_n_face > 0) {
      adjust = _face_edge_idx[0];
    }

    n_tria = pmne->prepa_blocks->n_tria  [i_part];
    n_quad = pmne->prepa_blocks->n_quad  [i_part];
    n_poly = pmne->prepa_blocks->n_poly2d[i_part];
    l_connec_poly = pmne->prepa_blocks->l_connec_poly2d[i_part];

    PDM_l_num_t *tria_vtx     = NULL;
    PDM_l_num_t *quad_vtx     = NULL;
    PDM_l_num_t *poly_vtx     = NULL;
    PDM_l_num_t *poly_vtx_idx = NULL;

    PDM_g_num_t *tria_ln_to_gn = NULL;
    PDM_g_num_t *quad_ln_to_gn = NULL;
    PDM_g_num_t *poly_ln_to_gn = NULL;

    PDM_l_num_t *tria_to_parent = NULL;
    PDM_l_num_t *quad_to_parent = NULL;
    PDM_l_num_t *poly_to_parent = NULL;

    if (we_have_elts[0] > 0) {
      PDM_malloc(tria_vtx,       n_tria * 3, PDM_l_num_t);
      PDM_malloc(tria_ln_to_gn,  n_tria,     PDM_g_num_t);
      PDM_malloc(tria_to_parent, n_tria,     PDM_l_num_t);
    }

    if (we_have_elts[1] > 0) {
      PDM_malloc(quad_vtx,       n_quad * 4, PDM_l_num_t);
      PDM_malloc(quad_ln_to_gn,  n_quad,     PDM_g_num_t);
      PDM_malloc(quad_to_parent, n_quad,     PDM_l_num_t);
    }

    if (we_have_elts[2] > 0) {
      PDM_malloc(poly_vtx_idx,   n_poly + 1,    PDM_l_num_t);
      PDM_malloc(poly_vtx,       l_connec_poly, PDM_l_num_t);
      PDM_malloc(poly_ln_to_gn,  n_poly,        PDM_g_num_t);
      PDM_malloc(poly_to_parent, n_poly,        PDM_l_num_t);
      poly_vtx_idx[0] = 0;
    }


    // Generate face->vtx connectivity
    int *_face_vtx = NULL;
    for (int i_face = 0; i_face < _n_face; i_face++) {
      _face_edge_idx[i_face+1] -= adjust;
    }
    if (is_signed) {
      PDM_compute_face_vtx_from_face_and_edge(_n_face,
                                              _face_edge_idx,
                                              _face_edge,
                                              _edge_vtx,
                                              &_face_vtx);
    }
    else {
      PDM_compute_face_vtx_from_face_and_edge_unsigned(_n_face,
                                                       _face_edge_idx,
                                                       _face_edge,
                                                       _edge_vtx,
                                                       &_face_vtx);
    }
    for (int i_face = 0; i_face < _n_face; i_face++) {
      _face_edge_idx[i_face+1] += adjust;
    }

    // Sort by element types
    int i_tria = 0;
    int i_quad = 0;
    int i_poly = 0;
    int idx_read = 0;
    for (int i_face = 0; i_face < _n_face; i_face++) {

      int _n_edge = _face_edge_idx[i_face+1] - _face_edge_idx[i_face];

      int *_elt_vtx = NULL;

      if (_n_edge == 3) {
        // Triangle
        tria_to_parent[i_tria] = i_face;
        tria_ln_to_gn [i_tria] = _face_ln_to_gn[i_face];
        _elt_vtx = &tria_vtx[3*i_tria];
        _parent_to_local[i_face] = i_tria++;
      }
      else if (_n_edge == 4) {
        // Quadrangle
        quad_to_parent[i_quad] = i_face;
        quad_ln_to_gn [i_quad] = _face_ln_to_gn[i_face];
        _elt_vtx = &quad_vtx[4*i_quad];
        _parent_to_local[i_face] = i_quad++;
      }
      else {
        // Polygon
        poly_to_parent[i_poly] = i_face;
        poly_ln_to_gn [i_poly] = _face_ln_to_gn[i_face];
        _elt_vtx = &poly_vtx[poly_vtx_idx[i_poly]];
        poly_vtx_idx[i_poly+1] = poly_vtx_idx[i_poly] + _n_edge;
        _parent_to_local[i_face] = i_poly++;
      }

      for (int idx_vtx = 0; idx_vtx < _n_edge; idx_vtx++) {
        _elt_vtx[idx_vtx] = _face_vtx[idx_read++];
      }

    } // End loop on faces
    PDM_free(_face_vtx);

    // Set sections
    if (we_have_elts[0] > 0) {
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_tria,
                                          i_part,
                                          n_tria,
                                          tria_vtx,
                                          tria_ln_to_gn,
                                          tria_to_parent,
                                          NULL,
                                          ownership);
    }
    if (we_have_elts[1] > 0) {
      PDM_part_mesh_nodal_elmts_std_set(pmne,
                                        id_quad,
                                        i_part,
                                        n_quad,
                                        quad_vtx,
                                        quad_ln_to_gn,
                                        quad_to_parent,
                                        NULL,
                                        ownership);
    }
    if (we_have_elts[2] > 0) {
      PDM_part_mesh_nodal_elmts_section_poly2d_set(pmne,
                                                   id_poly,
                                                   i_part,
                                                   n_poly,
                                                   poly_vtx_idx,
                                                   poly_vtx,
                                                   poly_ln_to_gn,
                                                   poly_to_parent,
                                                   ownership);
    }

  } // End loop on parts

  if (pmne->prepa_blocks != NULL) {
    PDM_free(pmne->prepa_blocks->n_cell);
    PDM_free(pmne->prepa_blocks->n_face);
    PDM_free(pmne->prepa_blocks->n_vtx);
    PDM_free(pmne->prepa_blocks->n_tria);
    PDM_free(pmne->prepa_blocks->n_quad);
    PDM_free(pmne->prepa_blocks->n_poly2d);
    PDM_free(pmne->prepa_blocks->l_connec_poly2d);
    PDM_free(pmne->prepa_blocks->face_vtx);
    PDM_free(pmne->prepa_blocks->cell_face_idx);
    PDM_free(pmne->prepa_blocks->cell_face);
    PDM_free(pmne->prepa_blocks->add_etat);
    PDM_free(pmne->prepa_blocks->numabs);
    PDM_free(pmne->prepa_blocks);
  }
}


void
PDM_part_mesh_nodal_elmts_cells_cellvtx_add
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_part,
  const int                           n_cell,
  const int                          *cell_vtx_idx,
  const int                          *cell_vtx,
  const PDM_g_num_t                  *numabs,
  const PDM_ownership_t               ownership
)
{
  CHECK_PMNE(pmne)
  CHECK_I_PART(pmne, id_part)

  int adjust = 0;
  if (n_cell > 0) {
    if (cell_vtx_idx[0] == 1) {
      adjust = 1;
    }
  }


  int n_part = 0;

  if (pmne->num_elmt_parent_to_local == NULL) {
    PDM_malloc(pmne->num_elmt_parent_to_local, pmne->n_part, PDM_l_num_t *);
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      pmne->num_elmt_parent_to_local[i_part] = NULL;
    }
  }

  PDM_malloc(pmne->num_elmt_parent_to_local[id_part], n_cell, PDM_l_num_t);
  for (int i = 0; i < n_cell; i++) {
    pmne->num_elmt_parent_to_local[id_part][i] = 0;
  }

  if (pmne->prepa_blocks == NULL) {
    PDM_malloc(pmne->prepa_blocks, 1, PDM_Mesh_nodal_prepa_blocks_t);
    pmne->prepa_blocks->t_add = 1;
    pmne->prepa_blocks->n_tetra_proc   = 0;  /* Nb de tetra par proc     */
    pmne->prepa_blocks->n_hexa_proc    = 0;  /* Nb d'hexa par proc       */
    pmne->prepa_blocks->n_prism_proc   = 0;  /* Nb de prisme par proc    */
    pmne->prepa_blocks->n_pyramid_proc = 0;  /* Nb de pyramide par proc  */
    pmne->prepa_blocks->n_poly3d_proc  = 0;  /* Nb de poly3d par proc    */
    PDM_malloc(pmne->prepa_blocks->n_cell,       pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_tetra,      pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_hexa,       pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_prism,      pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_pyramid,    pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_poly3d,     pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->cell_vtx_idx, pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->cell_vtx,     pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->add_etat,     pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->numabs,       pmne->n_part, PDM_g_num_t *);
    for (int i = 0; i < pmne->n_part; i++) {
      pmne->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (pmne->prepa_blocks->t_add != 1) {
    PDM_error(__FILE__, __LINE__, 0, "Error in PDM_part_mesh_nodal_elmts_cells_cellvtx_add : Another type of elements is currently is still in progress \n");
    abort();
  }

  /* Determination du type de chaque element */

  PDM_l_num_t n_tetra   = 0;
  PDM_l_num_t n_hexa    = 0;
  PDM_l_num_t n_prism   = 0;
  PDM_l_num_t n_pyramid = 0;
  PDM_l_num_t n_poly3d  = 0;

  for (int i = 0; i < n_cell; i++) {

    PDM_l_num_t n_som_cell = cell_vtx_idx[i+1] - cell_vtx_idx[i];
    if (n_som_cell == 4)
      n_tetra += 1;
    else if (n_som_cell == 5)
      n_pyramid += 1;
    else if (n_som_cell == 6)
      n_prism += 1;
    else if (n_som_cell == 8)
      n_hexa += 1;
    else {
      n_poly3d  += 1;
    }
  }

  pmne->prepa_blocks->n_tetra_proc          += n_tetra;
  pmne->prepa_blocks->n_hexa_proc           += n_hexa;
  pmne->prepa_blocks->n_prism_proc          += n_prism;
  pmne->prepa_blocks->n_pyramid_proc        += n_pyramid;
  pmne->prepa_blocks->n_poly3d_proc         += n_poly3d;
  pmne->prepa_blocks->n_tetra     [id_part] = n_tetra;
  pmne->prepa_blocks->n_hexa      [id_part] = n_hexa;
  pmne->prepa_blocks->n_prism     [id_part] = n_prism;
  pmne->prepa_blocks->n_pyramid   [id_part] = n_pyramid;
  pmne->prepa_blocks->n_poly3d    [id_part] = n_poly3d;
  pmne->prepa_blocks->cell_vtx_idx[id_part] = (PDM_l_num_t *) cell_vtx_idx;
  pmne->prepa_blocks->cell_vtx    [id_part] = (PDM_l_num_t *) cell_vtx;
  pmne->prepa_blocks->numabs      [id_part] = (PDM_g_num_t *) numabs;
  pmne->prepa_blocks->add_etat    [id_part] = 1;
  pmne->prepa_blocks->n_cell      [id_part] = n_cell;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < pmne->n_part; i++) {
    if (pmne->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (pmne->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[5];
    PDM_l_num_t som_elts[5];

    elts[0] = pmne->prepa_blocks->n_tetra_proc   > 0;
    elts[1] = pmne->prepa_blocks->n_hexa_proc    > 0;
    elts[2] = pmne->prepa_blocks->n_prism_proc   > 0;
    elts[3] = pmne->prepa_blocks->n_pyramid_proc > 0;
    elts[4] = pmne->prepa_blocks->n_poly3d_proc  > 0;

    PDM_MPI_Allreduce(elts, som_elts, 5, PDM_MPI_INT, PDM_MPI_SUM, pmne->comm);

    int id_bloc_tetra4   = -1;
    int id_bloc_hexa8    = -1;
    int id_bloc_prism6   = -1;
    int id_bloc_pyramid5 = -1;

    if (som_elts[0] > 0) {
      id_bloc_tetra4 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_TETRA4);
    }

    if (som_elts[1] > 0) {
      id_bloc_hexa8 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_HEXA8);
    }

    if (som_elts[2] > 0) {
      id_bloc_prism6 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_PRISM6);
    }

    if (som_elts[3] > 0) {
      id_bloc_pyramid5 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_PYRAMID5);
    }

    if (som_elts[4] > 0) {
      PDM_error(__FILE__, __LINE__, 0, "Non standard element detected\n");
    }

    /* Determination de la connectivite de chaque element */


    for (int i_part = 0; i_part < pmne->n_part; i_part++) {

      PDM_l_num_t  n_cell_courant                   = pmne->prepa_blocks->n_cell      [i_part];
      PDM_l_num_t *num_cell_parent_to_local_courant = pmne->num_elmt_parent_to_local  [i_part];
      PDM_l_num_t *cell_vtx_idx_courant             = pmne->prepa_blocks->cell_vtx_idx[i_part];
      PDM_l_num_t *cell_vtx_courant                 = pmne->prepa_blocks->cell_vtx    [i_part];
      PDM_g_num_t *numabs_courant                   = pmne->prepa_blocks->numabs      [i_part];

      PDM_l_num_t n_tetra_part   = pmne->prepa_blocks->n_tetra  [i_part];
      PDM_l_num_t n_hexa_part    = pmne->prepa_blocks->n_hexa   [i_part];
      PDM_l_num_t n_prism_part   = pmne->prepa_blocks->n_prism  [i_part];
      PDM_l_num_t n_pyramid_part = pmne->prepa_blocks->n_pyramid[i_part];

      PDM_l_num_t *connec_tetra   = NULL;
      PDM_l_num_t *connec_hexa    = NULL;
      PDM_l_num_t *connec_prism   = NULL;
      PDM_l_num_t *connec_pyramid = NULL;

      PDM_g_num_t *numabs_tetra   = NULL;
      PDM_g_num_t *numabs_hexa    = NULL;
      PDM_g_num_t *numabs_prism   = NULL;
      PDM_g_num_t *numabs_pyramid = NULL;

      PDM_l_num_t *num_parent_tetra   = NULL;
      PDM_l_num_t *num_parent_hexa    = NULL;
      PDM_l_num_t *num_parent_prism   = NULL;
      PDM_l_num_t *num_parent_pyramid = NULL;

      adjust = 0;
      if (n_cell_courant > 0) {
        if (cell_vtx_idx_courant[0] == 1) {
          adjust = 1;
        }
      }

      // n_tetra_part > 0
      if (som_elts[0] > 0) {
        PDM_malloc(connec_tetra    , 4 * n_tetra_part, PDM_l_num_t);
        PDM_malloc(numabs_tetra    ,     n_tetra_part, PDM_g_num_t);
        PDM_malloc(num_parent_tetra,     n_tetra_part, PDM_l_num_t);
      }

      // n_hexa_part > 0
      if (som_elts[1] > 0) {
        PDM_malloc(connec_hexa    , 8 * n_hexa_part, PDM_l_num_t);
        PDM_malloc(numabs_hexa    ,     n_hexa_part, PDM_g_num_t);
        PDM_malloc(num_parent_hexa,     n_hexa_part, PDM_l_num_t);
      }

      // n_prism_part > 0
      if (som_elts[2] > 0) {
        PDM_malloc(connec_prism    , 6 * n_prism_part, PDM_l_num_t);
        PDM_malloc(numabs_prism    ,     n_prism_part, PDM_g_num_t);
        PDM_malloc(num_parent_prism,     n_prism_part, PDM_l_num_t);
      }

      // n_pyramid_part > 0
      if (som_elts[3] > 0) {
        PDM_malloc(connec_pyramid    , 5 * n_pyramid_part, PDM_l_num_t);
        PDM_malloc(numabs_pyramid    ,     n_pyramid_part, PDM_g_num_t);
        PDM_malloc(num_parent_pyramid,     n_pyramid_part, PDM_l_num_t);
      }

      PDM_l_num_t *num_parent_tetra_courant   = num_parent_tetra;
      PDM_l_num_t *num_parent_hexa_courant    = num_parent_hexa;
      PDM_l_num_t *num_parent_prism_courant   = num_parent_prism;
      PDM_l_num_t *num_parent_pyramid_courant = num_parent_pyramid;

      PDM_l_num_t *connec_tetra_courant   = connec_tetra;
      PDM_l_num_t *connec_hexa_courant    = connec_hexa;
      PDM_l_num_t *connec_prism_courant   = connec_prism;
      PDM_l_num_t *connec_pyramid_courant = connec_pyramid;

      PDM_g_num_t *numabs_tetra_courant   = numabs_tetra;
      PDM_g_num_t *numabs_hexa_courant    = numabs_hexa;
      PDM_g_num_t *numabs_prism_courant   = numabs_prism;
      PDM_g_num_t *numabs_pyramid_courant = numabs_pyramid;

      PDM_l_num_t idx_tetra   = 0;
      PDM_l_num_t idx_hexa    = n_tetra_part;
      PDM_l_num_t idx_prism   = idx_hexa + n_hexa_part;
      PDM_l_num_t idx_pyramid = idx_prism + n_prism_part;

      PDM_l_num_t n_som_cell = 0;

      for (int i = 0; i < n_cell_courant; i++) {
        n_som_cell = cell_vtx_idx_courant[i+1] - cell_vtx_idx_courant[i];
        PDM_l_num_t idx_som_cell = cell_vtx_idx_courant[i] - adjust;
        PDM_l_num_t *connec_courant;

        if (n_som_cell == 4) {
          *num_parent_tetra_courant = i;
          num_parent_tetra_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tetra++;
          *numabs_tetra_courant = numabs_courant[i];
          numabs_tetra_courant += 1;
          connec_courant = connec_tetra_courant;
          connec_tetra_courant += n_som_cell;
        }
        else if (n_som_cell == 5) {
          *num_parent_pyramid_courant = i;
          num_parent_pyramid_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_pyramid++;
          *numabs_pyramid_courant = numabs_courant[i];
          numabs_pyramid_courant += 1;
          connec_courant = connec_pyramid_courant;
          connec_pyramid_courant += n_som_cell;
        }
        else if (n_som_cell == 6) {
          *num_parent_prism_courant = i;
          num_parent_prism_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_prism++;
          *numabs_prism_courant = numabs_courant[i];
          numabs_prism_courant += 1;
          connec_courant = connec_prism_courant;
          connec_prism_courant += n_som_cell;
        }
        else {
          *num_parent_hexa_courant = i;
          num_parent_hexa_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_hexa++;
          *numabs_hexa_courant = numabs_courant[i];
          numabs_hexa_courant += 1;
          connec_courant = connec_hexa_courant;
          connec_hexa_courant += n_som_cell;
        }

        /* Remplissage de la connectivite */
        for (int j = 0; j < n_som_cell; j++) {
          connec_courant[j] = cell_vtx_courant[idx_som_cell++];
        }
      }

      if (som_elts[0] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_tetra4,
                                          i_part,
                                          n_tetra_part,
                                          connec_tetra,
                                          numabs_tetra,
                                          num_parent_tetra,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[1] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_hexa8,
                                          i_part,
                                          n_hexa_part,
                                          connec_hexa,
                                          numabs_hexa,
                                          num_parent_hexa,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[2] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_prism6,
                                          i_part,
                                          n_prism_part,
                                          connec_prism,
                                          numabs_prism,
                                          num_parent_prism,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[3] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_pyramid5,
                                          i_part,
                                          n_pyramid_part,
                                          connec_pyramid,
                                          numabs_pyramid,
                                          num_parent_pyramid,
                                          NULL,//parent_entity_g_num,
                                          ownership);
    }

    if (pmne->prepa_blocks != NULL) {
      PDM_free(pmne->prepa_blocks->n_cell);
      PDM_free(pmne->prepa_blocks->n_tetra);
      PDM_free(pmne->prepa_blocks->n_hexa);
      PDM_free(pmne->prepa_blocks->n_prism);
      PDM_free(pmne->prepa_blocks->n_pyramid);
      PDM_free(pmne->prepa_blocks->n_poly3d);
      PDM_free(pmne->prepa_blocks->cell_vtx_idx);
      PDM_free(pmne->prepa_blocks->cell_vtx);
      PDM_free(pmne->prepa_blocks->add_etat);
      PDM_free(pmne->prepa_blocks->numabs);
      PDM_free(pmne->prepa_blocks);
      pmne->prepa_blocks = NULL;
    }
  }
}


void
PDM_part_mesh_nodal_elmts_faces_facevtx_add
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           id_part,
  const int                           n_face,
  const int                          *face_vtx_idx,
  const int                          *face_vtx,
  const PDM_g_num_t                  *numabs,
  const PDM_ownership_t               ownership
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  int adjust = 0;
  if (n_face > 0) {
    if (face_vtx_idx[0] == 1) {
      adjust = 1;
    }
  }

  int n_part = 0;
  if (pmne->num_elmt_parent_to_local == NULL) {
    PDM_malloc(pmne->num_elmt_parent_to_local, pmne->n_part, PDM_l_num_t *);
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      pmne->num_elmt_parent_to_local[i_part] = NULL;
    }
  }

  PDM_malloc(pmne->num_elmt_parent_to_local[id_part], n_face, PDM_l_num_t);
  for (int i = 0; i < n_face; i++) {
    pmne->num_elmt_parent_to_local[id_part][i] = 0;
  }

  if (pmne->prepa_blocks == NULL) {
    PDM_malloc(pmne->prepa_blocks, 1, PDM_Mesh_nodal_prepa_blocks_t);
    pmne->prepa_blocks->t_add = 3;
    pmne->prepa_blocks->n_tria_proc  = 0;   /* Nb de triangles par proc */
    pmne->prepa_blocks->n_quad_proc  = 0;   /* Nb de quads par proc */
    pmne->prepa_blocks->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    PDM_malloc(pmne->prepa_blocks->n_face,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_tria,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_quad,          pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->n_poly2d,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->l_connec_poly2d, pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->face_vtx_idx,    pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->face_vtx,        pmne->n_part, PDM_l_num_t *);
    PDM_malloc(pmne->prepa_blocks->add_etat,        pmne->n_part, PDM_l_num_t  );
    PDM_malloc(pmne->prepa_blocks->numabs,          pmne->n_part, PDM_g_num_t *);
    for (int i = 0; i < pmne->n_part; i++) {
      pmne->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (pmne->prepa_blocks->t_add != 3) {
    PDM_error(__FILE__, __LINE__, 0, "Error in PDM_part_mesh_nodal_elmts_cells_cellvtx_add : Another type of elements is currently is still in progress \n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d  = 0;

  for (int i = 0; i < n_face; i++) {

    PDM_l_num_t n_som_face = face_vtx_idx[i+1] - face_vtx_idx[i];
    if (n_som_face == 3)
      n_tria += 1;
    else if (n_som_face == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += n_som_face;
    }
  }

  pmne->prepa_blocks->n_tria_proc              += n_tria;
  pmne->prepa_blocks->n_quad_proc              += n_quad;
  pmne->prepa_blocks->n_poly2d_proc            += n_poly2d;
  pmne->prepa_blocks->add_etat       [id_part] = 1;
  pmne->prepa_blocks->n_tria         [id_part] = n_tria;
  pmne->prepa_blocks->n_quad         [id_part] = n_quad;
  pmne->prepa_blocks->n_poly2d       [id_part] = n_poly2d;
  pmne->prepa_blocks->l_connec_poly2d[id_part] = l_connec_poly2d;
  pmne->prepa_blocks->face_vtx_idx   [id_part] = (PDM_l_num_t *) face_vtx_idx;
  pmne->prepa_blocks->face_vtx       [id_part] = (PDM_l_num_t *) face_vtx;
  pmne->prepa_blocks->numabs         [id_part] = (PDM_g_num_t *) numabs;
  pmne->prepa_blocks->add_etat       [id_part] = 1;
  pmne->prepa_blocks->n_face         [id_part] = n_face;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < pmne->n_part; i++) {
    if (pmne->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (pmne->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = pmne->prepa_blocks->n_tria_proc > 0;
    elts[1] = pmne->prepa_blocks->n_quad_proc > 0;
    elts[2] = pmne->prepa_blocks->n_poly2d_proc > 0;

    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, pmne->comm);

    int id_bloc_tria3 = -1;
    int id_bloc_quad4 = -1;
    int id_bloc_poly_2d = -1;

    if (som_elts[0] > 0) {
      id_bloc_tria3 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_TRIA3);
    }

    if (som_elts[1] > 0) {
      id_bloc_quad4 = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_QUAD4);
    }

    if (som_elts[2] > 0) {
      id_bloc_poly_2d = PDM_part_mesh_nodal_elmts_add(pmne, PDM_MESH_NODAL_POLY_2D);
    }

    /* Determination de la connectivite de chaque element */

    for (int i_part = 0; i_part < pmne->n_part; i_part++) {

      PDM_l_num_t  n_face_courant                   = pmne->prepa_blocks->n_face      [i_part];
      PDM_l_num_t *num_cell_parent_to_local_courant = pmne->num_elmt_parent_to_local  [i_part];
      PDM_l_num_t *face_som_idx_courant             = pmne->prepa_blocks->face_vtx_idx[i_part];
      PDM_l_num_t *face_som_courant                 = pmne->prepa_blocks->face_vtx    [i_part];
      PDM_g_num_t *numabs_courant                   = pmne->prepa_blocks->numabs      [i_part];

      adjust = 0;
      if (n_face_courant > 0) {
        if (face_som_idx_courant[0] == 1) {
          adjust = 1;
        }
      }

      n_tria          = pmne->prepa_blocks->n_tria         [i_part];
      n_quad          = pmne->prepa_blocks->n_quad         [i_part];
      n_poly2d        = pmne->prepa_blocks->n_poly2d       [i_part];
      l_connec_poly2d = pmne->prepa_blocks->l_connec_poly2d[i_part];

      PDM_l_num_t *connec_tria       = NULL;
      PDM_l_num_t *connec_quad       = NULL;
      PDM_l_num_t *connec_poly2d     = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;

      PDM_g_num_t *numabs_tria   = NULL;
      PDM_g_num_t *numabs_quad   = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;

      PDM_l_num_t *num_parent_tria   = NULL;
      PDM_l_num_t *num_parent_quad   = NULL;
      PDM_l_num_t *num_parent_poly2d = NULL;


      if (som_elts[0] > 0) {
        PDM_malloc(connec_tria    ,  3 * n_tria, PDM_l_num_t);
        PDM_malloc(numabs_tria    ,      n_tria, PDM_g_num_t);
        PDM_malloc(num_parent_tria,      n_tria, PDM_l_num_t);
      }

      if (som_elts[1] > 0) {
        PDM_malloc(connec_quad    , 4 * n_quad, PDM_l_num_t);
        PDM_malloc(numabs_quad    ,     n_quad, PDM_g_num_t);
        PDM_malloc(num_parent_quad,     n_quad, PDM_l_num_t);
      }

      if (som_elts[2] > 0) {
        PDM_malloc(connec_poly2d_idx, n_poly2d + 1, PDM_l_num_t);
        connec_poly2d_idx[0] = 0;
        PDM_malloc(connec_poly2d    , l_connec_poly2d, PDM_l_num_t);
        PDM_malloc(numabs_poly2d    , n_poly2d       , PDM_g_num_t);
        PDM_malloc(num_parent_poly2d, n_poly2d       , PDM_l_num_t);
      }

      PDM_l_num_t *connec_tria_courant       = connec_tria;
      PDM_l_num_t *connec_quad_courant       = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant     = connec_poly2d;

      PDM_l_num_t *num_parent_tria_courant   = num_parent_tria;
      PDM_l_num_t *num_parent_quad_courant   = num_parent_quad;
      PDM_l_num_t *num_parent_poly2d_courant = num_parent_poly2d;

      PDM_g_num_t *numabs_tria_courant   = numabs_tria;
      PDM_g_num_t *numabs_quad_courant   = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      PDM_l_num_t n_som_face = 0;

      for (int i = 0; i < n_face_courant; i++) {
        n_som_face = face_som_idx_courant[i+1] - face_som_idx_courant[i];
        PDM_l_num_t idx_som_face = face_som_idx_courant[i] - adjust;
        PDM_l_num_t *connec_courant;

        if (n_som_face == 3) {
          *num_parent_tria_courant = i;
          num_parent_tria_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_som_face;
        }
        else if (n_som_face == 4) {
          *num_parent_quad_courant = i;
          num_parent_quad_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_som_face;
        }
        else {
          *num_parent_poly2d_courant = i;
          num_parent_poly2d_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) + n_som_face;
          connec_poly2d_idx_courant += 1;
          connec_courant = connec_poly2d_courant;
          connec_poly2d_courant += n_som_face;
        }

        /* Remplissage de la connectivite */

        for (int j = 0; j < n_som_face; j++)
          connec_courant[j] = face_som_courant[idx_som_face++];
      }

      if (som_elts[0] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_tria3,
                                          i_part,
                                          n_tria,
                                          connec_tria,
                                          numabs_tria,
                                          num_parent_tria,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[1] > 0)
        PDM_part_mesh_nodal_elmts_std_set(pmne,
                                          id_bloc_quad4,
                                          i_part,
                                          n_quad,
                                          connec_quad,
                                          numabs_quad,
                                          num_parent_quad,
                                          NULL,//parent_entity_g_num,
                                          ownership);

      if (som_elts[2] > 0)
        PDM_part_mesh_nodal_elmts_section_poly2d_set(pmne,
                                                     id_bloc_poly_2d,
                                                     i_part,
                                                     n_poly2d,
                                                     connec_poly2d_idx,
                                                     connec_poly2d,
                                                     numabs_poly2d,
                                                     num_parent_poly2d,
                                                     ownership);
    }
    if (pmne->prepa_blocks != NULL) {
      PDM_free(pmne->prepa_blocks->n_face);
      PDM_free(pmne->prepa_blocks->n_tria);
      PDM_free(pmne->prepa_blocks->n_quad);
      PDM_free(pmne->prepa_blocks->n_poly2d);
      PDM_free(pmne->prepa_blocks->l_connec_poly2d);
      PDM_free(pmne->prepa_blocks->face_vtx_idx);
      PDM_free(pmne->prepa_blocks->face_vtx);
      PDM_free(pmne->prepa_blocks->add_etat);
      PDM_free(pmne->prepa_blocks->numabs);
      PDM_free(pmne->prepa_blocks);
      pmne->prepa_blocks = NULL;
    }
  }
}


void
PDM_part_mesh_nodal_elmts_extend_to_encompassing_comm
(
  const PDM_MPI_Comm                  comm,
  const int                           n_part,
        PDM_part_mesh_nodal_elmts_t **pmne
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_part_mesh_nodal_elmts_t *_pmne = *pmne;

  int send_buf[3];
  int *recv_buf;
  PDM_malloc(recv_buf,3 * n_rank,int);

  int  n_block        = 0;
  int *blocks_id      = NULL;
  int  mesh_dimension = -1;

  if (_pmne != NULL) {
    assert(n_part == _pmne->n_part);

    n_block   = PDM_part_mesh_nodal_elmts_n_section_get  (_pmne);
    blocks_id = PDM_part_mesh_nodal_elmts_sections_id_get(_pmne);

    mesh_dimension = _pmne->mesh_dimension;
  }
  send_buf[0] = (_pmne == NULL);
  send_buf[1] = n_block;
  send_buf[2] = mesh_dimension;

  PDM_MPI_Allgather(send_buf, 3, PDM_MPI_INT, recv_buf, 3, PDM_MPI_INT, comm);

  // Find lowest rank with non-null mesh
  int master = -1;
  int n_null_rank = 0;
  int *i_null_rank;
  PDM_malloc(i_null_rank, n_rank, int);
  for (int i = 0; i < n_rank; i++) {
    if (recv_buf[3*i] == 1) {
      i_null_rank[n_null_rank++] = i;
    }
    else if (master < 0) {
      master = i;
    }
  }

  assert(master >= 0);


  n_block        = recv_buf[3*master+1];
  mesh_dimension = recv_buf[3*master+2];
  PDM_Mesh_nodal_elt_t *block_type;
  int  *block_order;
  int  *block_len_ho_ordering;
  char *char_buf = NULL;
  PDM_malloc(block_type,            n_block, PDM_Mesh_nodal_elt_t);
  PDM_malloc(block_order,           n_block, int                 );
  PDM_malloc(block_len_ho_ordering, n_block, int                 );
  PDM_free(recv_buf);

  int s_char_buf = 0;
  if (i_rank == master) {
    for (int iblock = 0; iblock < n_block; iblock++) {
      block_type[iblock] = PDM_part_mesh_nodal_elmts_section_type_get(_pmne,
                                                                      blocks_id[iblock]);
      if (PDM_Mesh_nodal_elmt_is_ho(block_type[iblock])) {
        block_order[iblock] = _pmne->sections_std[blocks_id[iblock]]->order;
        if (_pmne->sections_std[blocks_id[iblock]]->ho_ordering != NULL) {
          block_len_ho_ordering[iblock] = strlen(_pmne->sections_std[blocks_id[iblock]]->ho_ordering) + 1;
        }
        else {
          block_len_ho_ordering[iblock] = 0;
        }
      }
      else {
        block_order[iblock] = 1;
        block_len_ho_ordering[iblock] = 0;
      }

      s_char_buf += block_len_ho_ordering[iblock];
    }


    for (int dest = 0; dest < n_null_rank; dest++) {
      PDM_MPI_Send(block_type,            n_block, PDM_MPI_INT, i_null_rank[dest], 1, comm);
      PDM_MPI_Send(block_order,           n_block, PDM_MPI_INT, i_null_rank[dest], 1, comm);
      PDM_MPI_Send(block_len_ho_ordering, n_block, PDM_MPI_INT, i_null_rank[dest], 1, comm);
    }

    if (s_char_buf > 0) {
      PDM_malloc(char_buf, s_char_buf, char);
      int idx = 0;
      for (int iblock = 0; iblock < n_block; iblock++) {
        if (PDM_Mesh_nodal_elmt_is_ho(block_type[iblock])) {
          char *ho_ordering = _pmne->sections_std[blocks_id[iblock]]->ho_ordering;
          for (int i = 0; i < block_len_ho_ordering[iblock]-1; i++) {
            char_buf[idx++] = ho_ordering[i];
          }
          char_buf[idx++] = '\0';
        }
      }

      for (int dest = 0; dest < n_null_rank; dest++) {
        PDM_MPI_Send(char_buf, s_char_buf, PDM_MPI_CHAR, i_null_rank[dest], 1, comm);
      }
    }
  }
  else if (_pmne == NULL) {
    PDM_MPI_Recv(block_type,            n_block, PDM_MPI_INT, master, 1, comm);
    PDM_MPI_Recv(block_order,           n_block, PDM_MPI_INT, master, 1, comm);
    PDM_MPI_Recv(block_len_ho_ordering, n_block, PDM_MPI_INT, master, 1, comm);

    for (int i = 0; i < n_block; i++) {
      s_char_buf += block_len_ho_ordering[i];
    }

    if (s_char_buf > 0) {
      PDM_malloc(char_buf, s_char_buf, char);
      PDM_MPI_Recv(char_buf, s_char_buf, PDM_MPI_CHAR, master, 1, comm);
    }

    /* Create empty part_mesh_nodal_elmts */
    *pmne = PDM_part_mesh_nodal_elmts_create(mesh_dimension, n_part, comm);

    int idx = 0;
    for (int i = 0; i < n_block; i++) {

      PDM_Mesh_nodal_elt_t type = (PDM_Mesh_nodal_elt_t) block_type[i];

      int id_section = PDM_part_mesh_nodal_elmts_add(*pmne, type);

      if (type == PDM_MESH_NODAL_POLY_2D) {
        for (int ipart = 0; ipart < n_part; ipart++) {
          PDM_part_mesh_nodal_elmts_section_poly2d_set(*pmne,
                                                       id_section,
                                                       ipart,
                                                       0,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       PDM_OWNERSHIP_KEEP);
        }
      }
      else if (type == PDM_MESH_NODAL_POLY_3D) {
        for (int ipart = 0; ipart < n_part; ipart++) {
          PDM_part_mesh_nodal_elmts_section_poly3d_set(*pmne,
                                                       id_section,
                                                       ipart,
                                                       0,
                                                       0,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       PDM_OWNERSHIP_KEEP);
        }
      }
      else {
        char *ho_ordering = NULL;
        if (char_buf != NULL) {
          ho_ordering = &char_buf[idx];
        }
        idx += block_len_ho_ordering[i];
        if (n_part == 0) {
          (*pmne)->sections_std[id_section]->order = block_order[i];
          if (ho_ordering != NULL) {
            PDM_malloc((*pmne)->sections_std[id_section]->ho_ordering, block_len_ho_ordering[i], char);
            strcpy((*pmne)->sections_std[id_section]->ho_ordering, ho_ordering);
          }
        }

        for (int ipart = 0; ipart < n_part; ipart++) {
          PDM_part_mesh_nodal_elmts_std_ho_set(*pmne,
                                               id_section,
                                               ipart,
                                               0,
                                               NULL,
                                               NULL,
                                               NULL,
                                               NULL,
                                               block_order[i],
                                               ho_ordering,
                                               PDM_OWNERSHIP_KEEP);
        }
      }

    }

  }
  PDM_free(i_null_rank);
  PDM_free(block_type);
  PDM_free(block_order);
  PDM_free(block_len_ho_ordering);
  PDM_free(char_buf);
}


void
PDM_part_mesh_nodal_elmts_n_group_set
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           n_group
)
{
  CHECK_PMNE(pmne)
  if(pmne->n_group_elmt == NULL) {
    PDM_malloc(pmne->n_group_elmt   , pmne->n_part, int              *);
    PDM_malloc(pmne->group_elmt     , pmne->n_part, int             **);
    PDM_malloc(pmne->group_ln_to_gn , pmne->n_part, PDM_g_num_t     **);
    PDM_malloc(pmne->ownership_group, pmne->n_part, PDM_ownership_t  *);

    for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
      pmne->n_group_elmt  [i_part] = NULL;
      pmne->group_elmt    [i_part] = NULL;
      pmne->group_ln_to_gn[i_part] = NULL;
    }
  }

  pmne->n_group = n_group;
  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
    PDM_malloc(pmne->n_group_elmt   [i_part], n_group, int              );
    PDM_malloc(pmne->group_elmt     [i_part], n_group, int             *);
    PDM_malloc(pmne->group_ln_to_gn [i_part], n_group, PDM_g_num_t     *);
    PDM_malloc(pmne->ownership_group[i_part], n_group, PDM_ownership_t  );

    for(int i_group = 0; i_group < pmne->n_group; ++i_group) {
      pmne->n_group_elmt   [i_part][i_group] = 0;
      pmne->group_elmt     [i_part][i_group] = NULL;
      pmne->group_ln_to_gn [i_part][i_group] = NULL;
      pmne->ownership_group[i_part][i_group] = PDM_OWNERSHIP_KEEP;
    }
  }
}


void
PDM_part_mesh_nodal_elmts_group_set
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           i_part,
  const int                           i_group,
        int                           n_group_elmt,
        int                          *group_elmt,
        PDM_g_num_t                  *group_ln_to_gn,
        PDM_ownership_t               ownership_group
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, i_part)
  CHECK_GROUP (pmne, i_group)

  pmne->n_group_elmt  [i_part][i_group] = n_group_elmt;
  pmne->group_elmt    [i_part][i_group] = group_elmt;
  pmne->group_ln_to_gn[i_part][i_group] = group_ln_to_gn;

  pmne->ownership_group[i_part][i_group] = ownership_group;
}


void
PDM_part_mesh_nodal_elmts_group_get
(
        PDM_part_mesh_nodal_elmts_t   *pmne,
  const int                            i_part,
  const int                            i_group,
        int                           *n_group_elmt,
        int                          **group_elmt,
        PDM_g_num_t                  **group_ln_to_gn,
        PDM_ownership_t                ownership_group
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, i_part)
  CHECK_GROUP (pmne, i_group)

  *n_group_elmt   = pmne->n_group_elmt  [i_part][i_group];
  *group_elmt     = pmne->group_elmt    [i_part][i_group];
  *group_ln_to_gn = pmne->group_ln_to_gn[i_part][i_group];

  if (ownership_group != PDM_OWNERSHIP_BAD_VALUE) {
    pmne->ownership_group[i_part][i_group] = ownership_group;
  }
}


int
PDM_part_mesh_nodal_elmts_n_group_get
(
  PDM_part_mesh_nodal_elmts_t  *pmne
)
{
  CHECK_PMNE(pmne)
  return pmne->n_group;
}


int *
PDM_part_mesh_nodal_elmts_compute_sections_idx
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                     id_part
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  int n_section = pmne->n_section;
  int *section_elmt_idx;
  PDM_malloc(section_elmt_idx, n_section+1, int);
  section_elmt_idx[0] = 0;
  for(int i_section = 0; i_section < n_section; ++i_section) {
    int id_section = pmne->sections_id[i_section];
    int n_elmt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, id_part);
    section_elmt_idx[i_section+1] = section_elmt_idx[i_section] + n_elmt;
  }

  return section_elmt_idx;
}


int
PDM_part_mesh_nodal_elmts_cell_vtx_connect_get
(
        PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                           i_part,
        int                         **cell_vtx_idx,
        int                         **cell_vtx
)
{
  if (pmne == NULL) {
    *cell_vtx_idx = PDM_array_zeros_int(1);
    *cell_vtx     = NULL;
    return 0;
  }

  CHECK_I_PART(pmne, i_part)

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int n_cell = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne,
                                                     i_part);

  *cell_vtx_idx = PDM_array_zeros_int(n_cell + 1);

  int shift = 0;
  for (int isection = 0; isection < n_section; isection++) {

    int id_section = sections_id[isection];

    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            id_section);

    int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                               id_section,
                                                               i_part,
                                                               PDM_OWNERSHIP_BAD_VALUE);

    int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            i_part);
    int *connec_idx;
    int *connec;

    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                   id_section,
                                                   i_part,
                                                   &connec_idx,
                                                   &connec,
                                                   PDM_OWNERSHIP_BAD_VALUE);
    }
    else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
      PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    i_part,
                                                                    &connec_idx,
                                                                    &connec,
                                                                    PDM_OWNERSHIP_BAD_VALUE);
    }

    if (t_elt == PDM_MESH_NODAL_POLY_2D ||
        t_elt == PDM_MESH_NODAL_POLY_3D) {

      if (parent_num != NULL) {
        for (int i = 0; i < n_elt; i++) {
          (*cell_vtx_idx)[parent_num[i]+1] = connec_idx[i+1] - connec_idx[i];
        }
      } else {
        for (int i = 0; i < n_elt; i++) {
          (*cell_vtx_idx)[shift+i+1] = connec_idx[i+1] - connec_idx[i];
        }
      }
    } else {

      PDM_g_num_t *numabs;
      int         *_parent_num;
      PDM_g_num_t *parent_entity_g_num;
      int          order;
      const char  *ho_ordering;
      PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   i_part,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_BAD_VALUE);

      int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

      if (parent_num != NULL) {
        for (int i = 0; i < n_elt; i++) {
          (*cell_vtx_idx)[parent_num[i]+1] = n_vtx_elt;
        }
      } else {
        for (int i = 0; i < n_elt; i++) {
          (*cell_vtx_idx)[shift+i+1] = n_vtx_elt;
        }
      }
    }
    shift += n_elt;
  }

  PDM_array_accumulate_int(*cell_vtx_idx, n_cell+1);


  PDM_malloc(*cell_vtx, (*cell_vtx_idx)[n_cell], int);

  shift = 0;
  for (int isection = 0; isection < n_section; isection++) {

    int id_section = sections_id[isection];

    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            id_section);

    int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                               id_section,
                                                               i_part,
                                                               PDM_OWNERSHIP_BAD_VALUE);

    int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            i_part);
    int *connec_idx;
    int *connec;

    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                   id_section,
                                                   i_part,
                                                   &connec_idx,
                                                   &connec,
                                                   PDM_OWNERSHIP_BAD_VALUE);
    }
    else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
      PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    i_part,
                                                                    &connec_idx,
                                                                    &connec,
                                                                    PDM_OWNERSHIP_BAD_VALUE);
    }

    if (t_elt == PDM_MESH_NODAL_POLY_2D ||
        t_elt == PDM_MESH_NODAL_POLY_3D) {

      if (parent_num != NULL) {
        for (int i = 0; i < n_elt; i++) {
          for (int j = 0; j < connec_idx[i+1] - connec_idx[i]; j++) {
            (*cell_vtx)[(*cell_vtx_idx)[parent_num[i]] + j] = connec[connec_idx[i] + j];
          }
        }
      } else {
        for (int i = 0; i < n_elt; i++) {
          for (int j = 0; j < connec_idx[i+1] - connec_idx[i]; j++) {
            (*cell_vtx)[(*cell_vtx_idx)[shift+i] + j] = connec[connec_idx[i] + j];
          }
        }
      }

    }
    else {

      PDM_g_num_t *numabs;
      int         *_parent_num;
      PDM_g_num_t *parent_entity_g_num;
      int          order;
      const char  *ho_ordering;
      PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   i_part,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_BAD_VALUE);

      int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

      if (parent_num != NULL) {
        for (int i = 0; i < n_elt; i++) {
          for (int j = 0; j < n_vtx_elt; j++) {
            (*cell_vtx)[(*cell_vtx_idx)[parent_num[i]] + j] = connec[n_vtx_elt*i + j];
          }
        }
      } else {
        for (int i = 0; i < n_elt; i++) {
          for (int j = 0; j < n_vtx_elt; j++) {
            (*cell_vtx)[(*cell_vtx_idx)[shift+i] + j] = connec[n_vtx_elt*i + j];
          }
        }
      }
    }

    shift += n_elt;
  }

  return n_cell;
}


void
PDM_part_mesh_nodal_elmts_section_elt_to_entity_set
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
        int                         *elt_to_entity,
        PDM_ownership_t              ownership
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_elt_to_entity == NULL) {
      PDM_malloc(block->_elt_to_entity, pmne->n_part, int *);
    }
    block->_elt_to_entity[id_part] = elt_to_entity;

    block->elt_to_entity_owner = ownership;
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[id_section - PDM_BLOCK_ID_BLOCK_POLY2D];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_elt_to_entity == NULL) {
      PDM_malloc(block->_elt_to_entity, pmne->n_part, int *);
    }
    block->_elt_to_entity[id_part] = elt_to_entity;

    block->elt_to_entity_owner = ownership;
  }

  else {

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[id_section - PDM_BLOCK_ID_BLOCK_POLY3D];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_elt_to_entity == NULL) {
      PDM_malloc(block->_elt_to_entity, pmne->n_part, int *);
    }
    block->_elt_to_entity[id_part] = elt_to_entity;

    block->elt_to_entity_owner = ownership;
  }
}


int *
PDM_part_mesh_nodal_elmts_section_elmt_to_entity_get
(
        PDM_part_mesh_nodal_elmts_t *pmne,
  const int                          id_section,
  const int                          id_part,
        PDM_ownership_t              ownership
)
{
  CHECK_PMNE  (pmne)
  CHECK_I_PART(pmne, id_part)

  int *elt_to_entity = NULL;

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[id_section];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_elt_to_entity != NULL) {
      elt_to_entity = block->_elt_to_entity[id_part];
    }

    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      block->elt_to_entity_owner = ownership;
    }
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[id_section - PDM_BLOCK_ID_BLOCK_POLY2D];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_elt_to_entity != NULL) {
      elt_to_entity = block->_elt_to_entity[id_part];
    }

    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      block->elt_to_entity_owner = ownership;
    }
  }

  else {

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[id_section - PDM_BLOCK_ID_BLOCK_POLY3D];

    CHECK_BLOCK (block)
    CHECK_I_PART(block, id_part)

    if (block->_elt_to_entity != NULL) {
      elt_to_entity = block->_elt_to_entity[id_part];
    }

    if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
      block->elt_to_entity_owner = ownership;
    }
  }

  return elt_to_entity;
}

void
PDM_part_mesh_nodal_elmts_group_to_tag
(
  PDM_part_mesh_nodal_elmts_t   *pmne,
  int                         ***out_tag
)
{
  CHECK_PMNE  (pmne)

  int n_part = pmne->n_part;

  int **tag = NULL;
  PDM_malloc(tag, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_elmt  = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
    int n_group = PDM_part_mesh_nodal_elmts_n_group_get(pmne);

    PDM_malloc(tag[i_part], n_elmt, int);

    for(int i = 0; i < n_elmt; ++i) {
      tag[i_part][i] = -1;
    }

    int n_elmt_tag                  = 0;
    int elt_group_is_multiple       = 0;
    int n_elmt_with_different_group = 0;
    for(int i_group = 0; i_group < n_group; ++i_group) {

      int          n_group_elmt   = 0;
      int         *group_elmt     = 0;
      PDM_g_num_t *group_ln_to_gn = 0;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                          &n_group_elmt,
                                          &group_elmt,
                                          &group_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);

      for(int idx_group = 0; idx_group < n_group_elmt; ++idx_group) {
        int i_elt = group_elmt[idx_group]-1;
        if(tag[i_part][i_elt] == -1) {
          tag[i_part][i_elt] = i_group;
        } else {
          elt_group_is_multiple = 1;
          n_elmt_with_different_group++;
        }
        n_elmt_tag += 1;
      }
    }

    if(n_elmt_tag != n_elmt) {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_part_mesh_nodal_elmts_group_to_tag - All elements from PDM_part_mesh_nodal_elmts are not in groups for dimension %d (n_elmt_tag = %d, n_elmt = %d)\n",
                pmne->mesh_dimension, n_elmt_tag, n_elmt);
    }

    if(elt_group_is_multiple == 1) {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_part_mesh_nodal_elmts_group_to_tag - Several elements are more than one group associated (n_elmt_with_different_group=%i, n_elmt = %d) \n",
                n_elmt_with_different_group, n_elmt);
    }

  }

  *out_tag = tag;
}


void
PDM_part_mesh_nodal_elmts_tag_to_group
(
  PDM_part_mesh_nodal_elmts_t     *pmne,
  int                              n_group,
  int                            **tag
)
{
  CHECK_PMNE  (pmne)

  int n_part = pmne->n_part;

  /*
   * Set n_group in structure (not yet done normaly)
   */
  PDM_part_mesh_nodal_elmts_n_group_set(pmne, n_group);

  /*
   * Loop over part and fill structure
   */
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_elmt  = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
    int* group_elt_n = PDM_array_zeros_int(n_group);
    for(int i_elt = 0; i_elt < n_elmt; ++i_elt) {
      group_elt_n[tag[i_part][i_elt]]++;
    }

    int **group_elmt = NULL;
    PDM_malloc(group_elmt, n_group, int *);
    for(int i_group = 0; i_group < n_group; ++i_group) {
      PDM_malloc(group_elmt[i_group], group_elt_n[i_group], int);
      group_elt_n[i_group] = 0;
    }

    for(int i_elt = 0; i_elt < n_elmt; ++i_elt) {
      int i_group = tag[i_part][i_elt];
      group_elmt[i_group][group_elt_n[i_group]++] = i_elt+1;
    }

    for(int i_group = 0; i_group < n_group; ++i_group) {
      PDM_part_mesh_nodal_elmts_group_set(pmne,
                                          i_part,
                                          i_group,
                                          group_elt_n[i_group],
                                          group_elmt [i_group],
                                          NULL,
                                          PDM_OWNERSHIP_KEEP);
    }

    PDM_free(group_elt_n);
    PDM_free(group_elmt);
  }
}
