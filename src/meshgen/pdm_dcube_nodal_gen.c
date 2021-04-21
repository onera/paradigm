#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dcube_nodal_gen_priv.h"
#include "pdm_mpi.h"
#include "pdm_distrib.h"
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

static
PDM_g_num_t
_get_n_cell_abs
(
PDM_g_num_t          n_hexa,
PDM_Mesh_nodal_elt_t t_elmt
)
{
  switch (t_elmt) {
    case PDM_MESH_NODAL_TETRA4    :
    {
      // Each hexa in split in 5 tetra and boundary
      return n_hexa * 5;
    }
    break;

    case PDM_MESH_NODAL_PRISM6    :
    {
      return n_hexa * 2;
    }
    break;

    case PDM_MESH_NODAL_HEXA8    :
    {
      return n_hexa;
    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }
  return -1;
}

static
void
_generate_tetra_from_hexa
(
 PDM_dcube_nodal_t* dcube_nodal,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  PDM_UNUSED(dcube_nodal);
  PDM_UNUSED(dmesh_nodal);

  // n_cell          = n_g_hexa_cell * 5; // Because 1 Hexa = 5 Tetra
  // n_quad_lim      = 6 * 2 * n_quad_seg_face;
  // stride_quad_lim = 3;
  // stride_elmt     = 5 * 4;



}


static
void
_generate_prism_from_hexa
(
 PDM_dcube_nodal_t* dcube_nodal,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  PDM_UNUSED(dcube_nodal);
  PDM_UNUSED(dmesh_nodal);

  // n_cell          = n_g_hexa_cell * 2; // Because 1 Hexa = 2 Prims
  // n_quad_lim      = n_quad_seg_face * 4 + n_quad_seg_face * 2;
  // stride_elmt     = 6;
  // stride_quad_lim = 4;
  // stride_elmt = 2 * 6;




}

static
void
_generate_hexa_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  PDM_UNUSED(dmesh_nodal);

  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  // n_cell             = n_g_hexa_cell;
  // n_quad_lim         = 6 * n_quad_seg_face;
  int n_vtx_per_elmt     = 8;
  int n_vtx_per_elmt_lim = 4;

  PDM_g_num_t ngh       = dcube->n_g_hexa_cell_seg;
  PDM_g_num_t ngh2x     = ngh * ngh;
  PDM_g_num_t n_vtx_seg = dcube->n_vtx_seg;

  /* Setup volumic part */
  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt * dcube->dn_vtx ) * sizeof(PDM_g_num_t));

  for(int i_cell = 0; i_cell < dcube->dn_hexa_cell; ++i_cell) { // HEXA

    PDM_g_num_t g_cell = dcube->distrib_hexa[i_rank] + i_cell;

    PDM_g_num_t indk = ( g_cell                             ) / ( ngh2x );
    PDM_g_num_t indj = ( g_cell - indk * ngh2x              ) / ( ngh   );
    PDM_g_num_t indi = ( g_cell - indk * ngh2x - indj * ngh );

    delmt_vtx[n_vtx_per_elmt * i_cell    ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 1] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 2] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 3] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 4] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 5] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 6] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 7] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

  }

  int id_hexa = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_HEXA8);
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_hexa,
                                  dcube->dn_hexa_cell,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);

  /*
   * Setup surfacic part : each surface is distributed
   */
  // for(int i_group = 0; i_group < dcube->n_face_group; ++i_group){

  //   PDM_g_num_t* delmt_lim_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));


  //   // Premier numero de face du group
  //   PDM_g_num_t beg_group = i_group * dcube->distrib_quad_lim[i_rank];
  //   for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

  //     /* Pour le remplissage du elmt_group */
  //     PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

  //     // // I min
  //     // int i_min = 1;
  //     // for(int i = 0; i < )
  //   }

  //   int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
  //   PDM_DMesh_nodal_section_std_set(dmesh_nodal,
  //                                   id_quad,
  //                                   dcube->dn_hexa_cell,
  //                                   delmt_lim_vtx,
  //                                   dcube->owner_for_dmesh_nodal);

  //   /*
  //    *  Remplisage des face group (pour l'instant conforme au numero d'element global = aprés les volumiques !!!)
  //    */

  // }

  PDM_g_num_t* delmt_lim_imin_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));
  PDM_g_num_t* delmt_lim_imax_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));
  PDM_g_num_t* delmt_lim_jmin_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));
  PDM_g_num_t* delmt_lim_jmax_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));
  PDM_g_num_t* delmt_lim_kmin_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));
  PDM_g_num_t* delmt_lim_kmax_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_vtx ) * sizeof(PDM_g_num_t));

  /*
   * Generate IMIN plane
   */
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = 0; // Donc maillage en J,K

    delmt_lim_imin_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_imin_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_imin_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_imin_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;

    printf("IMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  /*
   * Generate IMAX plane (Opposite sens of imin to respect face orientation outside )
   */
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ngh; // Donc maillage en J,K

    delmt_lim_imax_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_imax_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_imax_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_imax_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    printf("IMAX :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  /*
   * Generate JMIN plane
   */
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = 0; // Donc maillage en I,K
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );

    delmt_lim_jmin_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_jmin_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_jmin_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_jmin_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  /*
   * Generate GMAX plane (Opposite sens of imin to respect face orientation outside )
   */
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = ngh; // Donc maillage en J,K
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );

    delmt_lim_jmax_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_jmax_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_jmax_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_jmax_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;

    printf("JMAX :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }


  /*
   * Generate KMIN plane
   */
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = 0; // Donc maillage en I,J
    PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

    delmt_lim_kmin_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_kmin_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_kmin_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_kmin_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    printf("kMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  /*
   * Generate KMAX plane (Opposite sens of imin to respect face orientation outside )
   */
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ngh; // Donc maillage en I,J
    PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

    delmt_lim_kmax_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_kmax_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_kmax_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_lim_kmax_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    printf("KMAX :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }


  int id_quad_imin = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
  int id_quad_imax = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
  int id_quad_jmin = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
  int id_quad_jmax = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
  int id_quad_kmin = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
  int id_quad_kmax = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad_imin,
                                  dcube->dn_quad_seq_lim,
                                  delmt_lim_imin_vtx,
                                  dcube->owner_for_dmesh_nodal);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad_imax,
                                  dcube->dn_quad_seq_lim,
                                  delmt_lim_imax_vtx,
                                  dcube->owner_for_dmesh_nodal);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad_jmin,
                                  dcube->dn_quad_seq_lim,
                                  delmt_lim_jmin_vtx,
                                  dcube->owner_for_dmesh_nodal);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad_jmax,
                                  dcube->dn_quad_seq_lim,
                                  delmt_lim_jmax_vtx,
                                  dcube->owner_for_dmesh_nodal);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad_kmin,
                                  dcube->dn_quad_seq_lim,
                                  delmt_lim_kmin_vtx,
                                  dcube->owner_for_dmesh_nodal);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad_kmax,
                                  dcube->dn_quad_seq_lim,
                                  delmt_lim_kmax_vtx,
                                  dcube->owner_for_dmesh_nodal);


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

  // Si l'utilisateur prends les résulats, c'est à travers dmesh_nodal
  dcube->owner_for_dmesh_nodal = owner;

  PDM_g_num_t n_vtx        = n_vtx_seg  * n_vtx_seg * n_vtx_seg;
  dcube->n_g_hexa_cell_seg = n_vtx_seg - 1;
  dcube->n_g_hexa_cell     = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;

  PDM_g_num_t n_quad_seg_face = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;
  double step = length / (double) dcube->n_g_hexa_cell_seg;

  /*
   * Create the dmesh_nodal that hold the resulting mesh
   */
  PDM_g_num_t n_cell_abs = _get_n_cell_abs(dcube->n_g_hexa_cell, t_elt);
  dcube->dmesh_nodal = PDM_DMesh_nodal_create(dcube->comm,
                                              3,
                                              n_vtx,
                                              n_cell_abs,   /* n_cell */
                                              -1,           /* n_face */
                                              -1);          /* n_edge */

  dcube->dn_vtx = PDM_compute_uniform_dn_entity(dcube->comm, n_vtx);

  double* dvtx_coord = (double *) malloc( 3 * (dcube->dn_vtx ) * sizeof(double *));
  PDM_DMesh_nodal_coord_set(dcube->dmesh_nodal,
                            dcube->dn_vtx,
                            dvtx_coord,
                            dcube->owner_for_dmesh_nodal); /* Le responsable de la mémoire est le dmesh_nodal */

  PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dcube->dmesh_nodal);

  PDM_g_num_t _dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  dcube->dn_vtx       = (int) _dn_vtx;

  printf(" _dn_vtx = %i \n", _dn_vtx);
  /*
   * Generate vertex
   */
  for(int i_vtx = 0; i_vtx < _dn_vtx; ++i_vtx) {

    PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

    PDM_g_num_t indk =  g_vtx                                 / ( n_vtx_seg * n_vtx_seg );
    PDM_g_num_t indj = (g_vtx - indk * n_vtx_seg * n_vtx_seg) / n_vtx_seg;
    PDM_g_num_t indi = (g_vtx - indk * n_vtx_seg * n_vtx_seg - indj * n_vtx_seg);

    // printf(" g_vtx = "PDM_FMT_G_NUM" -> ["PDM_FMT_G_NUM"/"PDM_FMT_G_NUM"/"PDM_FMT_G_NUM"] \n", g_vtx, indi, indj, indk);

    dvtx_coord[3 * i_vtx    ] = indi * step + zero_x;
    dvtx_coord[3 * i_vtx + 1] = indj * step + zero_y;
    dvtx_coord[3 * i_vtx + 2] = indk * step + zero_z;

  }

  /*
   * A Prevoir ici le rajout des vtx lié à la decomposition des hexa, ou raffinement (ie HO)
   */


  /*
   * Create the real hexa
   */
  dcube->distrib_hexa = PDM_compute_uniform_entity_distribution(dcube->comm, dcube->n_g_hexa_cell);

  PDM_g_num_t n_g_quad_lim = 6 * n_quad_seg_face;
  dcube->distrib_quad_lim     = PDM_compute_uniform_entity_distribution(dcube->comm, n_g_quad_lim);
  dcube->distrib_quad_seg_lim = PDM_compute_uniform_entity_distribution(dcube->comm, n_quad_seg_face);

  dcube->n_face_group = 6;

  PDM_g_num_t _dn_hexa_cell = dcube->distrib_hexa[i_rank+1] - dcube->distrib_hexa[i_rank];
  dcube->dn_hexa_cell       = (int) _dn_hexa_cell;

  PDM_g_num_t _dn_quad_lim  = dcube->distrib_quad_lim[i_rank+1] - dcube->distrib_quad_lim[i_rank];
  dcube->dn_quad_lim        = (int) _dn_quad_lim;

  PDM_g_num_t _dn_quad_seq_lim  = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  dcube->dn_quad_seq_lim        = (int) _dn_quad_seq_lim;


  switch (t_elt) {
    case PDM_MESH_NODAL_TETRA4    :
    {
      // Each hexa in split in 5 tetra and boundary
      _generate_tetra_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_PRISM6    :
    {
      _generate_prism_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_HEXA8    :
    {
      _generate_hexa_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }

  // PDM_g_num_t n_cell          = -1;
  // PDM_g_num_t n_quad_lim      = -1;
  // PDM_g_num_t stride_elmt     = -1;
  // PDM_g_num_t stride_quad_lim = -1;

  // switch (t_elt) {
  //   case PDM_MESH_NODAL_TETRA4    :
  //   {
  //     // Each hexa in split in 5 tetra and boundary
  //     n_cell          = n_g_hexa_cell * 5; // Because 1 Hexa = 5 Tetra
  //     n_quad_lim      = 6 * 2 * n_quad_seg_face;
  //     stride_quad_lim = 3;
  //     stride_elmt     = 5 * 4;
  //   }
  //   break;

  //   case PDM_MESH_NODAL_PRISM6    :
  //   {
  //   }
  //   break;

  //   case PDM_MESH_NODAL_HEXA8    :
  //   {

  //   }
  //   break;

  //   default :
  //     PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
  //     break;

  // }


  // Faux si prisme par exemple
  // dcube->delmt_vtx       = (PDM_g_num_t *) malloc(stride_elmt     * (dcube->dn_vtx          ) * sizeof(PDM_g_num_t *));
  // dcube->delmt_lim_vtx   = (PDM_g_num_t *) malloc(stride_quad_lim * (dcube->dn_vtx          ) * sizeof(PDM_g_num_t *));
  // dcube->dface_group_idx = (int         *) malloc(                  (dcube->n_face_group + 1) * sizeof(int         *));
  // dcube->dface_group     = (PDM_g_num_t *) malloc(                   dn_elmt_lim              * sizeof(PDM_g_num_t *));

  // PDM_g_num_t  *_delmt_vtx      = dcube->delmt_vtx;
  // PDM_g_num_t  *_dquad_lim_vtx  = dcube->dquad_lim_vtx;

  // int          *_dface_group_idx = dcube->dface_group_idx;
  // PDM_g_num_t  *_dface_group     = dcube->dface_group;

  /*
   * Generate elements
   */
  // for(int i_cell = 0; i_cell < _dn_hexa_cell/2; ++i_cell) { // HEXA
  // for(int i_cell = 0; i_cell < _dn_hexa_cell/2; ++i_cell) { // PRISM
  // for(int i_cell = 0; i_cell < _dn_hexa_cell/5; ++i_cell) { // TETRA

  //   /* We need to adapt for each type of elemt to generate */
  //   PDM_g_num_t g_cell = distrib_hexa[i_rank] + i_cell;

  //   PDM_g_num_t idx  = g_cell;
  //   PDM_g_num_t indk =  idx                                             / ( n_g_hexa_cell_seg * n_g_hexa_cell_seg );
  //   PDM_g_num_t indj = (idx - indk * n_g_hexa_cell_seg * n_g_hexa_cell_seg) / n_g_hexa_cell_seg;
  //   PDM_g_num_t indi = (idx - indk * n_g_hexa_cell_seg * n_g_hexa_cell_seg - indj * n_g_hexa_cell_seg);

  //   if(t_elt == PDM_MESH_NODAL_HEXA8) {
  //     _delmt_vtx[stride_elmt * i_cell    ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 1] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 2] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 3] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 4] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 5] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 6] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 7] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

  //   } else if( t_elt == PDM_MESH_NODAL_PRISM6) {

  //     _delmt_vtx[stride_elmt * i_cell     ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 1 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 2 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 3 ] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 4 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 5 ] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

  //     _delmt_vtx[stride_elmt * i_cell + 6 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 7 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 8 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 9 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 10] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
  //     _delmt_vtx[stride_elmt * i_cell + 11] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

  //   } else if( t_elt == PDM_MESH_NODAL_PYRAMID5) {

  //     // On a besoin d'un pint au milieu !!!
  //     _delmt_vtx[stride_elmt * i_cell     ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 1 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 2 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 3 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 4 ] = -1;

  //     _delmt_vtx[stride_elmt * i_cell + 5 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 6 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 7 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 8 ] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 9 ] = -1;

  //     _delmt_vtx[stride_elmt * i_cell + 10] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 11] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 12] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 13] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 14] = -1;

  //     _delmt_vtx[stride_elmt * i_cell + 15] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 16] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 17] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 18] = -1;
  //     _delmt_vtx[stride_elmt * i_cell + 19] = -1;

  //   } else if( t_elt == PDM_MESH_NODAL_TETRA4) {

  //     PDM_g_num_t n = indi+indj+indk;

  //     PDM_g_num_t ind1 = (indi     ) + (indj     ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // A (  i,  j,k  )
  //     PDM_g_num_t ind2 = (indi + 1 ) + (indj     ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // B (i+1,  j,k  )
  //     PDM_g_num_t ind3 = (indi + 1 ) + (indj + 1 ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // C (i+1,j+1,k  )
  //     PDM_g_num_t ind4 = (indi     ) + (indj + 1 ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // D (  i,j+1,k  )
  //     PDM_g_num_t ind5 = (indi     ) + (indj     ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // E (  i,  j,k+1)
  //     PDM_g_num_t ind6 = (indi + 1 ) + (indj     ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // F (i+1,  j,k+1)
  //     PDM_g_num_t ind7 = (indi + 1 ) + (indj + 1 ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // G (i+1,j+1,k+1)
  //     PDM_g_num_t ind8 = (indi     ) + (indj + 1 ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // H (  i,j+1,k+1)
  //     if( n % 2 == 0) {

  //       _delmt_vtx[stride_elmt * i_cell     ] = ind1;
  //       _delmt_vtx[stride_elmt * i_cell + 1 ] = ind2;
  //       _delmt_vtx[stride_elmt * i_cell + 2 ] = ind4;
  //       _delmt_vtx[stride_elmt * i_cell + 3 ] = ind5;

  //       _delmt_vtx[stride_elmt * i_cell + 4 ] = ind2;
  //       _delmt_vtx[stride_elmt * i_cell + 5 ] = ind3;
  //       _delmt_vtx[stride_elmt * i_cell + 6 ] = ind4;
  //       _delmt_vtx[stride_elmt * i_cell + 7 ] = ind7;

  //       _delmt_vtx[stride_elmt * i_cell + 8 ] = ind4;
  //       _delmt_vtx[stride_elmt * i_cell + 9 ] = ind5;
  //       _delmt_vtx[stride_elmt * i_cell + 10] = ind7;
  //       _delmt_vtx[stride_elmt * i_cell + 11] = ind8;

  //       _delmt_vtx[stride_elmt * i_cell + 12] = ind2;
  //       _delmt_vtx[stride_elmt * i_cell + 13] = ind5;
  //       _delmt_vtx[stride_elmt * i_cell + 14] = ind6;
  //       _delmt_vtx[stride_elmt * i_cell + 15] = ind7;

  //       _delmt_vtx[stride_elmt * i_cell + 16] = ind2;
  //       _delmt_vtx[stride_elmt * i_cell + 17] = ind4;
  //       _delmt_vtx[stride_elmt * i_cell + 18] = ind5;
  //       _delmt_vtx[stride_elmt * i_cell + 19] = ind7;

  //     } else {

  //       _delmt_vtx[stride_elmt * i_cell     ] = ind1;
  //       _delmt_vtx[stride_elmt * i_cell + 1 ] = ind3;
  //       _delmt_vtx[stride_elmt * i_cell + 2 ] = ind4;
  //       _delmt_vtx[stride_elmt * i_cell + 3 ] = ind8;

  //       _delmt_vtx[stride_elmt * i_cell + 4 ] = ind1;
  //       _delmt_vtx[stride_elmt * i_cell + 5 ] = ind2;
  //       _delmt_vtx[stride_elmt * i_cell + 6 ] = ind3;
  //       _delmt_vtx[stride_elmt * i_cell + 7 ] = ind6;

  //       _delmt_vtx[stride_elmt * i_cell + 8 ] = ind3;
  //       _delmt_vtx[stride_elmt * i_cell + 9 ] = ind6;
  //       _delmt_vtx[stride_elmt * i_cell + 10] = ind7;
  //       _delmt_vtx[stride_elmt * i_cell + 11] = ind8;

  //       _delmt_vtx[stride_elmt * i_cell + 12] = ind6;
  //       _delmt_vtx[stride_elmt * i_cell + 13] = ind5;
  //       _delmt_vtx[stride_elmt * i_cell + 14] = ind8;
  //       _delmt_vtx[stride_elmt * i_cell + 15] = ind1;

  //       _delmt_vtx[stride_elmt * i_cell + 16] = ind6;
  //       _delmt_vtx[stride_elmt * i_cell + 17] = ind3;
  //       _delmt_vtx[stride_elmt * i_cell + 18] = ind1;
  //       _delmt_vtx[stride_elmt * i_cell + 19] = ind8;

  //     }
  //   }
  // }



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
 * \param [out]  dn_hexa_cell       Number of cells stored in this process
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
 int                *dn_hexa_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
)
{
  *n_face_group = dcube->n_face_group;
  *dn_hexa_cell = dcube->dn_hexa_cell;
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
 PDM_dcube_nodal_t  *dcube,
 PDM_g_num_t       **delmt_vtx,
 double            **dvtx_coord,
 int               **dface_group_idx,
 PDM_g_num_t       **dface_group
)
{
  *delmt_vtx       = dcube->delmt_vtx;
  *dvtx_coord      = dcube->dvtx_coord;
  *dface_group_idx = dcube->dface_group_idx;
  *dface_group     = dcube->dface_group;
}


PDM_dmesh_nodal_t*
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t  *dcube
)
{
  return dcube->dmesh_nodal;
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
  free(dcube->distrib_hexa);
  free(dcube->distrib_quad_lim);
  free(dcube->distrib_quad_seg_lim);

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
