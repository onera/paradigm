#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dcube_nodal_gen_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"


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

PDM_dcube_nodal_t*
PDM_dcube_nodal_gen_init
(
      PDM_MPI_Comm          comm,
const PDM_g_num_t           n_vtx_seg,
const double                length,
const double                zero_x,
const double                zero_y,
const double                zero_z,
      PDM_Mesh_nodal_elt_t  t_elt,
      PDM_ownership_t       owner
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_nodal_t *dcube = (PDM_dcube_nodal_t *) malloc(sizeof(PDM_dcube_nodal_t));

  /*
   * Build dcube structure
   */

  dcube->comm      = comm;
  dcube->n_vtx_seg = n_vtx_seg;
  dcube->length    = length;
  dcube->zero_x    = zero_x;
  dcube->zero_y    = zero_y;
  dcube->zero_z    = zero_z;
  dcube->t_elt     = t_elt;
  dcube->owner     = owner;

  PDM_g_num_t n_vtx           = n_vtx_seg  * n_vtx_seg * n_vtx_seg;
  PDM_g_num_t n_hexa_cell_seg = n_vtx_seg - 1;
  PDM_g_num_t n_hexa_cell     = n_hexa_cell_seg * n_hexa_cell_seg * n_hexa_cell_seg;

  PDM_g_num_t n_quad_seg_face = n_hexa_cell_seg * n_hexa_cell_seg;
  double step = length / (double) n_hexa_cell_seg;

  PDM_g_num_t n_cell          = -1;
  PDM_g_num_t n_elmt_lim      = -1;
  PDM_g_num_t stride_elmt     = -1;
  PDM_g_num_t stride_elmt_lim = -1;

  switch (t_elt) {
    case PDM_MESH_NODAL_TETRA4    :
    {
      n_cell     = n_hexa_cell * 5; // Because 1 Hexa = 5 Tetra
      n_elmt_lim = 6 * 2 * n_quad_seg_face;
      // abort();
      stride_elmt_lim = 3;
      stride_elmt     = 5 * 4;

    }
    break;

    case PDM_MESH_NODAL_PRISM6    :
    {
      n_cell     = n_hexa_cell * 2; // Because 1 Hexa = 2 Prims
      n_elmt_lim = 6 * n_quad_seg_face;
      // abort();
      stride_elmt     = 6;
      stride_elmt_lim = 4;
      stride_elmt = 2 * 6;

    }
    break;
    case PDM_MESH_NODAL_HEXA8    :
    {
      n_cell     = n_hexa_cell;
      n_elmt_lim = 6 * n_quad_seg_face;
      stride_elmt = 8;
      stride_elmt_lim = 4;

    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }


  dcube->dmesh_nodal = PDM_DMesh_nodal_create(dcube->comm,
                                              3,
                                              n_vtx,
                                              n_cell,
                                              -1,   /* n_face */
                                              -1);  /* n_edge */


  PDM_g_num_t *distrib_vtx      = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_cell     = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_face_lim = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  //
  // Define distribution
  distrib_vtx[0]      = 0;
  distrib_cell[0]     = 0;
  distrib_face_lim[0] = 0;

  PDM_g_num_t step_vtx      = n_vtx / n_rank;
  PDM_g_num_t remainder_vtx = n_vtx % n_rank;

  PDM_g_num_t step_cell      = n_cell / n_rank;
  PDM_g_num_t remainder_cell = n_cell % n_rank;

  PDM_g_num_t step_face_im       = n_elmt_lim / n_rank;
  PDM_g_num_t remainder_face_lim = n_elmt_lim % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distrib_vtx[i]     = step_vtx;
    distrib_cell[i]    = step_cell;
    distrib_face_lim[i] = step_face_im;
    const int i1 = i - 1;
    if (i1 < remainder_vtx)
      distrib_vtx[i]  += 1;
    if (i1 < remainder_cell)
      distrib_cell[i]  += 1;
    if (i1 < remainder_face_lim)
      distrib_face_lim[i]  += 1;
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distrib_vtx[i]  += distrib_vtx[i-1];
    distrib_cell[i] += distrib_cell[i-1];
    distrib_face_lim[i] += distrib_face_lim[i-1];
  }

  dcube->n_face_group = 6;
  PDM_g_num_t _dn_cell     = distrib_cell[i_rank+1] - distrib_cell[i_rank];
  dcube->dn_cell           = (int) _dn_cell;
  PDM_g_num_t _dn_vtx      = distrib_vtx[i_rank+1]     - distrib_vtx[i_rank];
  dcube->dn_vtx            = (int) _dn_vtx;
  PDM_g_num_t _dn_elmt_lim = distrib_face_lim[i_rank+1] - distrib_face_lim[i_rank];
  int dn_elmt_lim          = (int) _dn_elmt_lim;

  dcube->dvtx_coord      = (double      *) malloc(              3 * (dcube->dn_vtx          ) * sizeof(double      *));

  // Faux si prisme par exemple
  dcube->delmt_vtx       = (PDM_g_num_t *) malloc(stride_elmt     * (dcube->dn_vtx          ) * sizeof(PDM_g_num_t *));
  dcube->delmt_lim_vtx   = (PDM_g_num_t *) malloc(stride_elmt_lim * (dcube->dn_vtx          ) * sizeof(PDM_g_num_t *));
  dcube->dface_group_idx = (int         *) malloc(                  (dcube->n_face_group + 1) * sizeof(int         *));
  dcube->dface_group     = (PDM_g_num_t *) malloc(                   dn_elmt_lim              * sizeof(PDM_g_num_t *));

  PDM_g_num_t  *_delmt_vtx      = dcube->delmt_vtx;
  PDM_g_num_t  *_delmt_lim_vtx  = dcube->delmt_lim_vtx;

  double       *_dvtx_coord      = dcube->dvtx_coord;
  int          *_dface_group_idx = dcube->dface_group_idx;
  PDM_g_num_t  *_dface_group     = dcube->dface_group;

  /*
   * Generate vertex
   */
  for(int i_vtx = 0; i_vtx < _dn_vtx; ++i_vtx) {

    PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

    PDM_g_num_t idx  = g_vtx;
    PDM_g_num_t indk =  idx                                 / ( n_vtx_seg * n_vtx_seg );
    PDM_g_num_t indj = (idx - indk * n_vtx_seg * n_vtx_seg) / n_vtx_seg;
    PDM_g_num_t indi = (idx - indk * n_vtx_seg * n_vtx_seg - indj * n_vtx_seg);

    printf(" idx = %i -> [%i/%i/%i] \n", idx, indi, indj, indk);

    _dvtx_coord[3 * i_vtx    ] = indi * step + zero_x;
    _dvtx_coord[3 * i_vtx + 1] = indj * step + zero_y;
    _dvtx_coord[3 * i_vtx + 2] = indk * step + zero_z;


  }

  /*
   * Generate elements
   */
  // for(int i_cell = 0; i_cell < _dn_cell/2; ++i_cell) { // HEXA
  // for(int i_cell = 0; i_cell < _dn_cell/2; ++i_cell) { // PRISM
  for(int i_cell = 0; i_cell < _dn_cell/5; ++i_cell) { // TETRA

    /* We need to adapt for each type of elemt to generate */
    PDM_g_num_t g_cell = distrib_cell[i_rank] + i_cell;

    PDM_g_num_t idx  = g_cell;
    PDM_g_num_t indk =  idx                                             / ( n_hexa_cell_seg * n_hexa_cell_seg );
    PDM_g_num_t indj = (idx - indk * n_hexa_cell_seg * n_hexa_cell_seg) / n_hexa_cell_seg;
    PDM_g_num_t indi = (idx - indk * n_hexa_cell_seg * n_hexa_cell_seg - indj * n_hexa_cell_seg);

    if(t_elt == PDM_MESH_NODAL_HEXA8) {
      _delmt_vtx[stride_elmt * i_cell    ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 1] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 2] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 3] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 4] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 5] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 6] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 7] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

    } else if( t_elt == PDM_MESH_NODAL_PRISM6) {

      _delmt_vtx[stride_elmt * i_cell     ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 1 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 2 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 3 ] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 4 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 5 ] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

      _delmt_vtx[stride_elmt * i_cell + 6 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 7 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 8 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 9 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 10] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
      _delmt_vtx[stride_elmt * i_cell + 11] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

    } else if( t_elt == PDM_MESH_NODAL_PYRAMID5) {

      // On a besoin d'un pint au milieu !!!
      _delmt_vtx[stride_elmt * i_cell     ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 1 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 2 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 3 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 4 ] = -1;

      _delmt_vtx[stride_elmt * i_cell + 5 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 6 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 7 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 8 ] = -1;
      _delmt_vtx[stride_elmt * i_cell + 9 ] = -1;

      _delmt_vtx[stride_elmt * i_cell + 10] = -1;
      _delmt_vtx[stride_elmt * i_cell + 11] = -1;
      _delmt_vtx[stride_elmt * i_cell + 12] = -1;
      _delmt_vtx[stride_elmt * i_cell + 13] = -1;
      _delmt_vtx[stride_elmt * i_cell + 14] = -1;

      _delmt_vtx[stride_elmt * i_cell + 15] = -1;
      _delmt_vtx[stride_elmt * i_cell + 16] = -1;
      _delmt_vtx[stride_elmt * i_cell + 17] = -1;
      _delmt_vtx[stride_elmt * i_cell + 18] = -1;
      _delmt_vtx[stride_elmt * i_cell + 19] = -1;

    } else if( t_elt == PDM_MESH_NODAL_TETRA4) {

      PDM_g_num_t n = indi+indj+indk;

      PDM_g_num_t ind1 = (indi     ) + (indj     ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // A (  i,  j,k  )
      PDM_g_num_t ind2 = (indi + 1 ) + (indj     ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // B (i+1,  j,k  )
      PDM_g_num_t ind3 = (indi + 1 ) + (indj + 1 ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // C (i+1,j+1,k  )
      PDM_g_num_t ind4 = (indi     ) + (indj + 1 ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // D (  i,j+1,k  )
      PDM_g_num_t ind5 = (indi     ) + (indj     ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // E (  i,  j,k+1)
      PDM_g_num_t ind6 = (indi + 1 ) + (indj     ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // F (i+1,  j,k+1)
      PDM_g_num_t ind7 = (indi + 1 ) + (indj + 1 ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // G (i+1,j+1,k+1)
      PDM_g_num_t ind8 = (indi     ) + (indj + 1 ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // H (  i,j+1,k+1)
      if( n % 2 == 0) {

        _delmt_vtx[stride_elmt * i_cell     ] = ind1;
        _delmt_vtx[stride_elmt * i_cell + 1 ] = ind2;
        _delmt_vtx[stride_elmt * i_cell + 2 ] = ind4;
        _delmt_vtx[stride_elmt * i_cell + 3 ] = ind5;

        _delmt_vtx[stride_elmt * i_cell + 4 ] = ind2;
        _delmt_vtx[stride_elmt * i_cell + 5 ] = ind3;
        _delmt_vtx[stride_elmt * i_cell + 6 ] = ind4;
        _delmt_vtx[stride_elmt * i_cell + 7 ] = ind7;

        _delmt_vtx[stride_elmt * i_cell + 8 ] = ind4;
        _delmt_vtx[stride_elmt * i_cell + 9 ] = ind5;
        _delmt_vtx[stride_elmt * i_cell + 10] = ind7;
        _delmt_vtx[stride_elmt * i_cell + 11] = ind8;

        _delmt_vtx[stride_elmt * i_cell + 12] = ind2;
        _delmt_vtx[stride_elmt * i_cell + 13] = ind5;
        _delmt_vtx[stride_elmt * i_cell + 14] = ind6;
        _delmt_vtx[stride_elmt * i_cell + 15] = ind7;

        _delmt_vtx[stride_elmt * i_cell + 16] = ind2;
        _delmt_vtx[stride_elmt * i_cell + 17] = ind4;
        _delmt_vtx[stride_elmt * i_cell + 18] = ind5;
        _delmt_vtx[stride_elmt * i_cell + 19] = ind7;

      } else {

        _delmt_vtx[stride_elmt * i_cell     ] = ind1;
        _delmt_vtx[stride_elmt * i_cell + 1 ] = ind3;
        _delmt_vtx[stride_elmt * i_cell + 2 ] = ind4;
        _delmt_vtx[stride_elmt * i_cell + 3 ] = ind8;

        _delmt_vtx[stride_elmt * i_cell + 4 ] = ind1;
        _delmt_vtx[stride_elmt * i_cell + 5 ] = ind2;
        _delmt_vtx[stride_elmt * i_cell + 6 ] = ind3;
        _delmt_vtx[stride_elmt * i_cell + 7 ] = ind6;

        _delmt_vtx[stride_elmt * i_cell + 8 ] = ind3;
        _delmt_vtx[stride_elmt * i_cell + 9 ] = ind6;
        _delmt_vtx[stride_elmt * i_cell + 10] = ind7;
        _delmt_vtx[stride_elmt * i_cell + 11] = ind8;

        _delmt_vtx[stride_elmt * i_cell + 12] = ind6;
        _delmt_vtx[stride_elmt * i_cell + 13] = ind5;
        _delmt_vtx[stride_elmt * i_cell + 14] = ind8;
        _delmt_vtx[stride_elmt * i_cell + 15] = ind1;

        _delmt_vtx[stride_elmt * i_cell + 16] = ind6;
        _delmt_vtx[stride_elmt * i_cell + 17] = ind3;
        _delmt_vtx[stride_elmt * i_cell + 18] = ind1;
        _delmt_vtx[stride_elmt * i_cell + 19] = ind8;

      }
    }
  }



  /*
   * Generate elements lim
   */

  return dcube;

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
PDM_dcube_nodal_gen_dim_get
(
 PDM_dcube_nodal_t  *dcube,
 int                *n_face_group,
 int                *dn_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
)
{
  *n_face_group = dcube->n_face_group;
  *dn_cell      = dcube->dn_cell;
  *dn_face      = -1;
  *dn_vtx       = dcube->dn_vtx;
  *sface_vtx    = 0; // dcube->dface_vtx_idx[dcube->dn_face];
  *sface_group  = 0; // dcube->dface_group_idx[dcube->n_face_group];
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
PDM_dcube_nodal_gen_data_get
(
 PDM_dcube_nodal_t  *pdm_dcube_nodal,
 PDM_g_num_t       **delmt_vtx,
 double            **dvtx_coord,
 int               **dface_group_idx,
 PDM_g_num_t       **dface_group
)
{
  *delmt_vtx       = pdm_dcube_nodal->delmt_vtx;
  *dvtx_coord      = pdm_dcube_nodal->dvtx_coord;
  *dface_group_idx = pdm_dcube_nodal->dface_group_idx;
  *dface_group     = pdm_dcube_nodal->dface_group;
}


PDM_dmesh_nodal_t*
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t  *pdm_dcube_nodal
)
{
  return pdm_dcube_nodal->dmesh_nodal;
}


/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  id            dcube identifier
 *
 */

void
PDM_dcube_nodal_gen_free
(
PDM_dcube_nodal_t        *dcube
)
{
  // if(dcube->owner == PDM_OWNERSHIP_KEEP) {
  //   if (dcube->dface_cell  != NULL)
  //     free(dcube->dface_cell);

  //   if (dcube->dface_vtx_idx  != NULL)
  //     free(dcube->dface_vtx_idx);

  //   if (dcube->dface_vtx  != NULL)
  //     free(dcube->dface_vtx);

  //   if (dcube->dvtx_coord  != NULL)
  //     free(dcube->dvtx_coord);

  //   if (dcube->dface_group_idx  != NULL)
  //     free(dcube->dface_group_idx);

  //   if (dcube->dface_group  != NULL)
  //     free(dcube->dface_group);
  // }

  if(dcube->owner == PDM_OWNERSHIP_KEEP) {
    /* Si l'utilisateur fait le get il doit liberer le dmesh_nodal */
    PDM_DMesh_nodal_free(dcube->dmesh_nodal, 0);
  }

  free(dcube);

}
