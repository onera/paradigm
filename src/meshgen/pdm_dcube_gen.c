#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_priv.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _dcube_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

typedef struct  {
  PDM_MPI_Comm     comm;            /*!< MPI communicator                          */
  PDM_g_num_t      n_vtx_seg;       /*!< Number of vertices in segments            */
  double           length;          /*!< Segment length                            */
  double           zero_x;          /*!< Coordinates of the origin                 */
  double           zero_y;          /*!< Coordinates of the origin                 */
  double           zero_z;          /*!< Coordinates of the origin                 */
  int              n_face_group;    /*!< Number of faces groups                    */
  int              dn_cell;         /*!< Number of cells stored in this process    */
  int              dn_face;         /*!< Number of faces stored in this process    */
  int              dn_vtx;          /*!< Number of vertices stored in this process */
  PDM_g_num_t     *dface_cell;      /*!< Faces from cells connectivity             */
  int             *dface_vtx_idx;   /*!< Faces from vertices connectivity index    */
  PDM_g_num_t     *dface_vtx;       /*!< Faces from vertices connectivity          */
  double          *dvtx_coord;      /*!< Vertices coordinates                      */
  int             *dface_group_idx; /*!< Faces groups index                        */
  PDM_g_num_t     *dface_group;     /*!< Faces groups                              */
  PDM_ownership_t  owner;           /*!< Which have the responsabilities of results*/
} _dcube_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dcubes  = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

static _dcube_t *
_get_from_id
(
 int  id
)
{
  _dcube_t *dcube = (_dcube_t *) PDM_Handles_get (_dcubes, id);

  if (dcube == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "PDM_part_dcube error : Bad dcube identifier\n");
  }

  return dcube;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a distributed cube
 *
 * \param [out]  id             dcube identifier
 * \param [in]   comm           Communicator
 * \param [in]   n_vtx_seg      Number of vertices in segments
 * \param [in]   length         Segment length
 * \param [in]   zero_x         Coordinates of the origin
 * \param [in]   zero_y         Coordinates of the origin
 * \param [in]   zero_z         Coordinates of the origin
 *
 */

void
PDM_dcube_gen_init
(
      int             *id,
      PDM_MPI_Comm     comm,
const PDM_g_num_t      n_vtx_seg,
const double           length,
const double           zero_x,
const double           zero_y,
const double           zero_z,
      PDM_ownership_t  owner
)
{

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  /*
   * Search a dcube free id
   */

  if (_dcubes == NULL) {
    _dcubes = PDM_Handles_create (4);
  }

  _dcube_t *dcube = (_dcube_t *) malloc(sizeof(_dcube_t));

  *id = PDM_Handles_store (_dcubes, dcube);

  /*
   * Build dcube structure
   */

  dcube->comm      = comm;
  dcube->n_vtx_seg = n_vtx_seg;
  dcube->length    = length;
  dcube->zero_x    = zero_x;
  dcube->zero_y    = zero_y;
  dcube->zero_z    = zero_z;
  dcube->owner     = owner;

  PDM_g_num_t n_vtx       = n_vtx_seg * n_vtx_seg * n_vtx_seg;
  PDM_g_num_t n_face_seg  = n_vtx_seg - 1;
  PDM_g_num_t n_face      = 3 * n_face_seg * n_face_seg * n_vtx_seg;
  PDM_g_num_t n_cell      = n_face_seg * n_face_seg * n_face_seg;
  PDM_g_num_t n_face_face = n_face_seg * n_face_seg;
  PDM_g_num_t n_vtx_face  = n_vtx_seg * n_vtx_seg;
  PDM_g_num_t n_face_lim  = 6 * n_face_face;
  double step = length / (double) n_face_seg;
  PDM_g_num_t *distrib_vtx      = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_face     = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_cell     = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_face_lim = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  //
  // Define distribution

  distrib_vtx[0]      = 0;
  distrib_face[0]     = 0;
  distrib_cell[0]     = 0;
  distrib_face_lim[0] = 0;

  PDM_g_num_t step_vtx      = n_vtx / n_rank;
  PDM_g_num_t remainder_vtx = n_vtx % n_rank;

  PDM_g_num_t step_face      = n_face / n_rank;
  PDM_g_num_t remainder_face = n_face % n_rank;

  PDM_g_num_t step_cell      = n_cell / n_rank;
  PDM_g_num_t remainder_cell = n_cell % n_rank;

  PDM_g_num_t step_face_im       = n_face_lim / n_rank;
  PDM_g_num_t remainder_face_lim = n_face_lim % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distrib_vtx[i]     = step_vtx;
    distrib_face[i]    = step_face;
    distrib_cell[i]    = step_cell;
    distrib_face_lim[i] = step_face_im;
    const int i1 = i - 1;
    if (i1 < remainder_vtx)
      distrib_vtx[i]  += 1;
    if (i1 < remainder_face)
      distrib_face[i]  += 1;
    if (i1 < remainder_cell)
      distrib_cell[i]  += 1;
    if (i1 < remainder_face_lim)
      distrib_face_lim[i]  += 1;
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distrib_vtx[i]  += distrib_vtx[i-1];
    distrib_face[i] += distrib_face[i-1];
    distrib_cell[i] += distrib_cell[i-1];
    distrib_face_lim[i] += distrib_face_lim[i-1];
  }

  dcube->n_face_group = 6;
  PDM_g_num_t _dn_cell = distrib_cell[i_rank+1] - distrib_cell[i_rank];
  dcube->dn_cell       = (int) _dn_cell;
  PDM_g_num_t _dn_face = distrib_face[i_rank+1]    - distrib_face[i_rank];
  dcube->dn_face       = (int) _dn_face;
  PDM_g_num_t _dn_vtx  = distrib_vtx[i_rank+1]     - distrib_vtx[i_rank];
  dcube->dn_vtx        = (int) _dn_vtx;
  PDM_g_num_t _dn_face_lim = distrib_face_lim[i_rank+1] - distrib_face_lim[i_rank];
  int dn_face_lim = (int) _dn_face_lim;

  dcube->dface_cell      = (PDM_g_num_t *) malloc(2*(dcube->dn_face) * sizeof(PDM_g_num_t *));
  dcube->dface_vtx_idx   = (int *)          malloc((dcube->dn_face + 1) * sizeof(int *));
  dcube->dface_vtx       = (PDM_g_num_t *) malloc(4*(dcube->dn_face)  * sizeof(PDM_g_num_t *));
  dcube->dvtx_coord      = (double *)       malloc(3*(dcube->dn_vtx)  * sizeof(double *));
  dcube->dface_group_idx = (int *)          malloc((dcube->n_face_group + 1)  * sizeof(int *));
  dcube->dface_group     = (PDM_g_num_t *) malloc(dn_face_lim * sizeof(PDM_g_num_t *));

  PDM_g_num_t  *_dface_cell      = dcube->dface_cell;
  int          *_dface_vtx_idx   = dcube->dface_vtx_idx;
  PDM_g_num_t  *_dface_vtx       = dcube->dface_vtx;
  double       *_dvtx_coord      = dcube->dvtx_coord;
  int          *_dface_group_idx = dcube->dface_group_idx;
  PDM_g_num_t  *_dface_group     = dcube->dface_group;

  _dface_vtx_idx[0] = 0;
  for (int i = 1; i < dcube->dn_face + 1; i++) {
    _dface_vtx_idx[i] = 4 + _dface_vtx_idx[i-1];
  }

  //
  // Coordinates

  const PDM_g_num_t b_vtx_z = distrib_vtx[i_rank] / n_vtx_face;
  const PDM_g_num_t r_vtx_z = distrib_vtx[i_rank] % n_vtx_face;

  const PDM_g_num_t b_vtx_y = r_vtx_z / n_vtx_seg;
  const PDM_g_num_t b_vtx_x = r_vtx_z % n_vtx_seg;

  int i_vtx = 0;
  int cpt   = 0;

  for(PDM_g_num_t k = b_vtx_z; k < n_vtx_seg; k++) {
    PDM_g_num_t _b_vtx_y = 0;
    if (k == b_vtx_z)
      _b_vtx_y = b_vtx_y;
    for(PDM_g_num_t j = _b_vtx_y; j < n_vtx_seg; j++) {
      PDM_g_num_t _b_vtx_x = 0;
      if ((k == b_vtx_z) && (j == b_vtx_y))
        _b_vtx_x = b_vtx_x;
      for(PDM_g_num_t i = _b_vtx_x; i < n_vtx_seg; i++) {
        _dvtx_coord[3 * i_vtx    ] = i * step + zero_x;
        _dvtx_coord[3 * i_vtx + 1] = j * step + zero_y;
        _dvtx_coord[3 * i_vtx + 2] = k * step + zero_z;
        cpt   += 1;
        i_vtx += 1;
        if (cpt == dcube->dn_vtx)
          break;
      }
      if (cpt == dcube->dn_vtx)
        break;
    }
    if (cpt == dcube->dn_vtx)
      break;
  }

  //
  // face_vtx et face_cell

  cpt = 0;

  PDM_g_num_t serie   = n_face / 3;
  PDM_g_num_t i_serie = distrib_face[i_rank] / serie;
  PDM_g_num_t r_serie = distrib_face[i_rank] % serie;

  PDM_g_num_t b1 = 0;
  PDM_g_num_t r1 = 0;

  PDM_g_num_t b2 = 0;
  PDM_g_num_t b3 = 0;

  b1 = r_serie / n_face_face;
  r1 = r_serie % n_face_face;

  b2 = r1 / n_face_seg;
  b3 = r1 % n_face_seg;


  if (i_serie == 0) {
  /* switch (i_serie) { */

  /* case 0 : */

    //
    // Faces zmin -> zmax

    for(PDM_g_num_t k = b1; k < n_vtx_seg; k++) {
      PDM_g_num_t _b2 = 0;
      if (k == b1)
        _b2 = b2;
      for(PDM_g_num_t j = _b2; j < n_face_seg; j++) {
        PDM_g_num_t _b3 = 0;
        if ((k == b1) && (j == b2))
          _b3 = b3;
        for(PDM_g_num_t i = _b3; i < n_face_seg; i++) {
          _dface_vtx[cpt * 4    ] = k * n_vtx_seg * n_vtx_seg + (    j * n_vtx_seg + i + 1);
          _dface_vtx[cpt * 4 + 1] = k * n_vtx_seg * n_vtx_seg + ((j+1) * n_vtx_seg + i + 1);
          _dface_vtx[cpt * 4 + 2] = k * n_vtx_seg * n_vtx_seg + ((j+1) * n_vtx_seg + i + 2);
          _dface_vtx[cpt * 4 + 3] = k * n_vtx_seg * n_vtx_seg + (    j * n_vtx_seg + i + 2);
          if (k == 0) {
            _dface_cell[2*cpt + 0] = j * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] = 0;
          }
          else if (k == n_face_seg) {
            _dface_cell[2*cpt + 0] = (k-1) * n_face_seg * n_face_seg + j * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] = 0;
          }
          else {
            _dface_cell[2*cpt + 0] = (k-1) * n_face_seg * n_face_seg + j * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] =     k * n_face_seg * n_face_seg + j * n_face_seg + i + 1;
          }
          cpt += 1;
          if (cpt == dcube->dn_face)
            break;
        }
        if (cpt == dcube->dn_face)
          break;
      }
      if (cpt == dcube->dn_face)
        break;
    }
    b1 = 0;
    b2 = 0;
    b3 = 0;
  }


  /* if (cpt == dcube->dn_face) */
  /*     break; */

  if ((i_serie == 1) || ((i_serie == 0) &&  (cpt != dcube->dn_face))) {

    //
    // Faces xmin -> xmax

    for(PDM_g_num_t i = b1; i < n_vtx_seg; i++) {
      PDM_g_num_t _b2 = 0;
      if (i == b1)
        _b2 = b2;
      for(PDM_g_num_t k = _b2; k < n_face_seg; k++) {
        PDM_g_num_t _b3 = 0;
        if ((i == b1) && (k == b2))
          _b3 = b3;
        for(PDM_g_num_t j = _b3; j < n_face_seg; j++) {
          _dface_vtx[cpt * 4    ] =     k * n_vtx_seg * n_vtx_seg +     j * n_vtx_seg + i + 1;
          _dface_vtx[cpt * 4 + 1] =     k * n_vtx_seg * n_vtx_seg + (j+1) * n_vtx_seg + i + 1;
          _dface_vtx[cpt * 4 + 2] = (k+1) * n_vtx_seg * n_vtx_seg + (j+1) * n_vtx_seg + i + 1;
          _dface_vtx[cpt * 4 + 3] = (k+1) * n_vtx_seg * n_vtx_seg +     j * n_vtx_seg + i + 1;
          if (i == 0) {
            _dface_cell[2*cpt + 0] = k * n_face_seg * n_face_seg + j * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] = 0;
          }

          else if (i == n_face_seg) {
            _dface_cell[2*cpt + 0] = k * n_face_seg * n_face_seg + j * n_face_seg + i;
            _dface_cell[2*cpt + 1] = 0;
          }

          else {
            _dface_cell[2*cpt + 0] = k * n_face_seg * n_face_seg + j * n_face_seg + i ;
            _dface_cell[2*cpt + 1] = k * n_face_seg * n_face_seg + j * n_face_seg + i + 1;
          }
          cpt += 1;
          if (cpt == dcube->dn_face)
            break;
        }
        if (cpt == dcube->dn_face)
          break;
      }
      if (cpt == dcube->dn_face)
        break;
    }
    b1 = 0;
    b2 = 0;
    b3 = 0;
  }


  /* if (cpt == dcube->dn_face) */
  /*   break; */

  if ((i_serie == 2) || ((i_serie == 1 || i_serie == 0) && (cpt != dcube->dn_face))) {
    /* case 2 : */

    //
    // Faces ymin -> ymax

    for(PDM_g_num_t j = b1; j < n_vtx_seg; j++) {
      PDM_g_num_t _b2 = 0;
      if (j == b1)
        _b2 = b2;
      for(PDM_g_num_t i = _b2; i < n_face_seg; i++) {
        PDM_g_num_t _b3 = 0;
        if ((j == b1) && (i == b2))
          _b3 = b3;
        for(PDM_g_num_t k = _b3; k < n_face_seg; k++) {
          _dface_vtx[cpt * 4    ] =     k * n_vtx_seg * n_vtx_seg + j * n_vtx_seg + i + 1    ;
          _dface_vtx[cpt * 4 + 1] =     k * n_vtx_seg * n_vtx_seg + j * n_vtx_seg + i + 1 + 1;
          _dface_vtx[cpt * 4 + 2] = (k+1) * n_vtx_seg * n_vtx_seg + j * n_vtx_seg + i + 1 + 1;
          _dface_vtx[cpt * 4 + 3] = (k+1) * n_vtx_seg * n_vtx_seg + j * n_vtx_seg + i + 1    ;
          if (j == 0) {
            _dface_cell[2*cpt + 0] = k * n_face_seg * n_face_seg + j * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] = 0;
          }

          else if (j == n_face_seg) {
            _dface_cell[2*cpt + 0] =  k * n_face_seg * n_face_seg + (j-1) * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] = 0;
          }

          else {
            _dface_cell[2*cpt + 0] = k * n_face_seg * n_face_seg + (j-1) * n_face_seg + i + 1;
            _dface_cell[2*cpt + 1] = k * n_face_seg * n_face_seg +     j * n_face_seg + i + 1;
          }
          cpt += 1;
          if (cpt == dcube->dn_face)
            break;
        }
        if (cpt == dcube->dn_face)
          break;
      }
      if (cpt == dcube->dn_face)
        break;
    }
  }

  //
  // Faces limite

  cpt = 0;
  PDM_g_num_t b_face;
  int cpt1 = 0;
  int cpt3 = 0;
  int first_group = 0;

  serie   = n_face_lim / dcube->n_face_group;
  i_serie = distrib_face_lim[i_rank] / serie;
  r_serie = distrib_face_lim[i_rank] % serie;

  for (int i = 0; i < dcube->n_face_group + 1; i++)
    _dface_group_idx[i] = 0;


  //  switch (i_serie) {

  if (i_serie == 0) {

    // case 0 :

    //
    // Faces zmin

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = 0;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_face_seg; j++) {
      for(PDM_g_num_t i = 0; i < n_face_seg; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * n_face_seg + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[1] = cpt - cpt1;

    /* if (cpt == dn_face_lim) */
    /*   break; */
    first_group = 0;
  }


  if ((i_serie == 1) || ((i_serie == 0) && (cpt != dn_face_lim))) {
    //  case 1 :

    //
    // Faces zmax

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_seg * n_face_seg * n_face_seg;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_face_seg; j++) {
      for(PDM_g_num_t i = 0; i < n_face_seg; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * n_face_seg + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[2] = cpt - cpt1;
     first_group = 0;
  }
    /* if (cpt == dn_face_lim) */
    /*   break; */


  if ((i_serie == 2) || (((i_serie == 0) || (i_serie == 1)) && (cpt != dn_face_lim))) {
    //  case 2 :

    //
    // Faces xmin

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_seg * n_face_seg * n_vtx_seg;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_face_seg; j++) {
      for(PDM_g_num_t i = 0; i < n_face_seg; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * n_face_seg + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[3] = cpt - cpt1;
     first_group = 0;
  }
    /* if (cpt == dn_face_lim) */
    /*   break; */


  if ((i_serie == 3) || (((i_serie == 0) || (i_serie == 1)  || (i_serie == 2)) && (cpt != dn_face_lim))) {
    //  case 3 :

    //
    // Faces xmax

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_seg * n_face_seg * (n_vtx_seg + n_face_seg);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_face_seg; j++) {
      for(PDM_g_num_t i = 0; i < n_face_seg; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * n_face_seg + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[4] = cpt - cpt1;
    first_group = 0;
  }
    /* if (cpt == dn_face_lim) */
    /*   break; */


  if ((i_serie == 4) || (((i_serie == 0) || (i_serie == 1)  || (i_serie == 2) || (i_serie == 3)) && (cpt != dn_face_lim))) {
    //  case 4 :

    //
    // Faces ymin

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_seg * n_face_seg * (n_vtx_seg + n_vtx_seg);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_face_seg; j++) {
      for(PDM_g_num_t i = 0; i < n_face_seg; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * n_face_seg + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[5] = cpt - cpt1;
    first_group = 0;
  }

    /* if (cpt == dn_face_lim) */
    /*   break; */


  if ((i_serie == 5) || (((i_serie == 0) || (i_serie == 1)  || (i_serie == 2) || (i_serie == 3) || (i_serie == 4)) && (cpt != dn_face_lim))) {
  /* case 5 : */

    //
    // Faces ymax

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_seg * n_face_seg * (n_vtx_seg + n_vtx_seg + n_face_seg);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_face_seg; j++) {
      for(PDM_g_num_t i = 0; i < n_face_seg; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * n_face_seg + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

  _dface_group_idx[6] = cpt - cpt1;
  first_group = 0;

  }


  for (int i = 1; i < dcube->n_face_group + 1; i++)
    _dface_group_idx[i] += _dface_group_idx[i-1];

  free(distrib_vtx);
  free(distrib_face);
  free(distrib_cell);
  free(distrib_face_lim);

}

void
PROCF (pdm_dcube_gen_init, PDM_DCUBE_GEN_INIT)
(
      int             *id,
const PDM_MPI_Fint    *comm,
const PDM_g_num_t     *n_vtx_seg,
const double          *length,
const double          *zero_x,
const double          *zero_y,
const double          *zero_z,
      PDM_ownership_t *owner
)
{

  PDM_MPI_Fint comm1 = *((PDM_MPI_Fint *) comm);

  PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c(comm1);

  PDM_dcube_gen_init (id,
                      c_comm,
                      *n_vtx_seg,
                      *length,
                      *zero_x,
                      *zero_y,
                      *zero_z,
                      *owner);
}


/**
 *
 * \brief Return distributed cube size
 *
 * \param [in]   id            dcube identifier
 * \param [out]  n_face_group  Number of faces groups
 * \param [out]  dn_cell       Number of cells stored in this process
 * \param [out]  dn_face       Number of faces stored in this process
 * \param [out]  dn_vtx        Number of vertices stored in this process
 * \param [out]  sface_vtx     Length of dface_vtx array
 * \param [out]  sface_group   Length of dface_group array
 *
 */

void
PDM_dcube_gen_dim_get
(
 int                id,
 int                *n_face_group,
 int                *dn_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
)
{
  _dcube_t *dcube = _get_from_id(id);

  *n_face_group = dcube->n_face_group;
  *dn_cell      = dcube->dn_cell;
  *dn_face      = dcube->dn_face;
  *dn_vtx       = dcube->dn_vtx;
  *sface_vtx    = dcube->dface_vtx_idx[dcube->dn_face];
  *sface_group  = dcube->dface_group_idx[dcube->n_face_group];
}


void
PROCF(pdm_dcube_gen_dim_get, PDM_DCUBE_GEN_DIM_GET)
(
 int                *id,
 int                *n_face_group,
 int                *dn_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
)
{
  PDM_dcube_gen_dim_get (*id,
                         n_face_group,
                         dn_cell,
                         dn_face,
                         dn_vtx,
                         sface_vtx,
                         sface_group);

}

/**
 *
 * \brief Return distributed cube data
 *
 * \param [in]  id              dcube identifier
 * \param [out] dface_cell      Faces from cells connectivity (size = 2 * dn_face)
 * \param [out] dface_vtx_idx    Faces from vertices connectivity index (size = dn_face + 1)
 * \param [out] dface_vtx       Faces from vertices connectivity (size = sface_vtx)
 * \param [out] dvtx_coord      Vertices coordinates (size = 3 * dn_vtx)
 * \param [out] dface_group_idx Faces groups index (size = n_face_group + 1)
 * \param [out] dface_group     Faces groups (size = sface_group)
 *
 */

void
PDM_dcube_gen_data_get
(
 int                 id,
 PDM_g_num_t       **dface_cell,
 int               **dface_vtx_idx,
 PDM_g_num_t       **dface_vtx,
 double            **dvtx_coord,
 int               **dface_group_idx,
 PDM_g_num_t       **dface_group
)
{
  _dcube_t *dcube = _get_from_id(id);

  *dface_cell      = dcube->dface_cell;
  *dface_vtx_idx   = dcube->dface_vtx_idx;
  *dface_vtx       = dcube->dface_vtx;
  *dvtx_coord      = dcube->dvtx_coord;
  *dface_group_idx = dcube->dface_group_idx;
  *dface_group     = dcube->dface_group;
}


void
PROCF (pdm_dcube_gen_data_get, PDM_DCUBE_GEN_DATA_GET)
(
 int              *id,
 PDM_g_num_t      *dface_cell,
 int              *dface_vtx_idx,
 PDM_g_num_t      *dface_vtx,
 double           *dvtx_coord,
 int              *dface_group_idx,
 PDM_g_num_t      *dface_group
)
{
  _dcube_t *dcube = _get_from_id(*id);

  for (int i = 0; i < 2 * dcube->dn_face; i++)
    dface_cell[i] = dcube->dface_cell[i];

  for (int i = 0; i < dcube->dn_face + 1 ; i++)
    dface_vtx_idx[i] = dcube->dface_vtx_idx[i];

  for (int i = 0; i < dcube->dface_vtx_idx[dcube->dn_face] ; i++)
    dface_vtx[i] = dcube->dface_vtx[i];

  for (int i = 0; i < 3 * dcube->dn_vtx; i++)
    dvtx_coord[i] = dcube->dvtx_coord[i];

  for (int i = 0; i < dcube->n_face_group + 1; i++)
    dface_group_idx[i] = dcube->dface_group_idx[i];

  for (int i = 0; i < dcube->dface_group_idx[dcube->n_face_group]; i++)
    dface_group[i] = dcube->dface_group[i];
}


/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  id            dcube identifier
 *
 */

void
PDM_dcube_gen_free
(
const int id
)
{
  _dcube_t *dcube = _get_from_id(id);

  if(dcube->owner == PDM_OWNERSHIP_KEEP) {
    if (dcube->dface_cell  != NULL)
      free(dcube->dface_cell);

    if (dcube->dface_vtx_idx  != NULL)
      free(dcube->dface_vtx_idx);

    if (dcube->dface_vtx  != NULL)
      free(dcube->dface_vtx);

    if (dcube->dvtx_coord  != NULL)
      free(dcube->dvtx_coord);

    if (dcube->dface_group_idx  != NULL)
      free(dcube->dface_group_idx);

    if (dcube->dface_group  != NULL)
      free(dcube->dface_group);
  }

  free(dcube);

  PDM_Handles_handle_free (_dcubes, id, PDM_FALSE);

  const int n_dcube = PDM_Handles_n_get (_dcubes);

  if (n_dcube == 0) {
    PDM_Handles_free (_dcubes);

  }
}


void
PROCF (pdm_dcube_gen_free, PDM_DCUBE_GEN_FREE)
(
int *id
)
{
  PDM_dcube_gen_free (*id);
}
