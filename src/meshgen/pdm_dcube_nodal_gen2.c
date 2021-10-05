#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dcube_nodal_gen2.h"
#include "pdm_dcube_nodal_gen2_priv.h"
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
    case PDM_MESH_NODAL_TRIA3    :
    {
      // Each hexa in split in 5 tetra and boundary
      return n_hexa * 2;
    }
    break;

    case PDM_MESH_NODAL_QUAD4    :
    {
      // Each hexa in split in 5 tetra and boundary
      return n_hexa;
    }
    break;

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


static inline PDM_g_num_t
sub2ind
(
 PDM_dcube_nodal2_t *dcube,
 PDM_g_num_t         icell,
 PDM_g_num_t         jcell,
 PDM_g_num_t         kcell,
 int                 i,
 int                 j,
 int                 k
 )
{
  return 1 + dcube->order*icell + i + (dcube->order*dcube->nx + 1)*(dcube->order*jcell + j + (dcube->order*dcube->ny + 1)*(dcube->order*kcell + k));
}



static
void
_generate_tria_surf
(
 PDM_dcube_nodal2_t* dcube,
 PDM_dmesh_nodal_t*  dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIA3, order);

  int dn_elt = 2 * dcube->dn_quad;

  /* Set up surface part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_quad; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = g / dcube->nx;

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order - j; i++) {
        delt_vtx[idx++] = sub2ind(dcube, indi, indj, 0, i, j, 0);
      }
    }

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order - j; i++) {
        delt_vtx[idx++] = sub2ind(dcube, indi, indj, 0, dcube->order-i, dcube->order-j, 0);
      }
    }

  }


  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_quad_surf
(
 PDM_dcube_nodal2_t* dcube,
 PDM_dmesh_nodal_t*  dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUAD4, order);

  int dn_elt = dcube->dn_quad;

  /* Set up surface part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_quad; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = g / dcube->nx;

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order; i++) {
        delt_vtx[idx++] = sub2ind(dcube, indi, indj, 0, i, j, 0);
      }
    }
  }


  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_QUAD4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}



static
void
_generate_tetra_vol
(
 PDM_dcube_nodal2_t* dcube,
 PDM_dmesh_nodal_t*  dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRA4, order);

  int dn_elt = 5 * dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_hexa; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    if ((indi + indj + indk) % 2 == 0) {
      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi+1, indj, indk+1, -i, j, -k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj+1, indk+1, i, -j, -k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi+1, indj+1, indk, -i, -j, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, order-j-k, i+j, i+k);
          }
        }
      }
    }

    else {
      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i+j+k, j, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, order-j, order-k, i+j+k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, j, order-i, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, i+j+k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i+j, j+k, k+i);
          }
        }
      }
    }

  }

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_TETRA4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}

static
void
_generate_prism_vol
(
 PDM_dcube_nodal2_t* dcube,
 PDM_dmesh_nodal_t*  dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PRISM6, order);

  int dn_elt = 2 * dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_hexa; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, k);
        }
      }
    }

    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, order-i, order-j, k);
        }
      }
    }

  }

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_PRISM6);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_hexa_vol
(
 PDM_dcube_nodal2_t* dcube,
 PDM_dmesh_nodal_t*  dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXA8, order);

  int dn_elt = dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_hexa; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order; i++) {
          delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, k);
        }
      }
    }

  }

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_HEXA8);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_corners
(
 PDM_dcube_nodal2_t *dcube,
 PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int dim = dmesh_nodal->mesh_dimension;

  PDM_g_num_t gn_corner = 1 << dim;

  PDM_g_num_t *distrib_corner = PDM_compute_uniform_entity_distribution(dcube->comm, gn_corner);
  int dn_corner = (int) (distrib_corner[i_rank+1] - distrib_corner[i_rank]);

  PDM_g_num_t *delt_vtx = malloc (sizeof(PDM_g_num_t) * dn_corner);
  for (int icorner = 0; icorner < dn_corner; icorner++) {
    PDM_g_num_t g = distrib_corner[i_rank] + icorner;

    PDM_g_num_t i = g % 2;
    PDM_g_num_t k = g / 4;
    PDM_g_num_t j = (g - i - 4*k) / 2;

    delt_vtx[icorner] = sub2ind(dcube,
                                i*dcube->nx, j*dcube->ny, k*dcube->nz,
                                0, 0, 0);
  }

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->corner, PDM_MESH_NODAL_POINT);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->corner,
                                        id_section,
                                        dn_corner,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);

  dmesh_nodal->corner->n_g_elmts = gn_corner;
}


static
void
_generate_ridges
(
 PDM_dcube_nodal2_t *dcube,
 PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int dim = dmesh_nodal->mesh_dimension;
  int order = dcube->order;

  int n_group = 2*dim;

  int *dgroup_elt_idx = malloc (sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;

  // Group 1: from (0,0,0) to (1,0,0)
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

PDM_dcube_nodal2_t*
PDM_dcube_nodal_gen2_init
(
 PDM_MPI_Comm          comm,
 const PDM_g_num_t     nx,
 const PDM_g_num_t     ny,
 const PDM_g_num_t     nz,
 const double          length,
 const double          zero_x,
 const double          zero_y,
 const double          zero_z,
 PDM_Mesh_nodal_elt_t  t_elt,
 const int             order,
 PDM_ownership_t       owner
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_nodal2_t *dcube = (PDM_dcube_nodal2_t *) malloc(sizeof(PDM_dcube_nodal2_t));

  double t1 = PDM_MPI_Wtime();


  /*
   * Build dcube structure
   */
  dcube->comm   = comm;
  dcube->nx     = nx;
  dcube->ny     = ny;
  dcube->nz     = nz;
  dcube->length = length;
  dcube->zero_x = zero_x;
  dcube->zero_y = zero_y;
  dcube->zero_z = zero_z;
  dcube->t_elt  = t_elt;
  dcube->order  = order;
  dcube->owner  = owner;

  int dim = 3;
  if(t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_QUAD4 ) {
    dim = 2;
  }
  if(t_elt == PDM_MESH_NODAL_POINT || t_elt == PDM_MESH_NODAL_BAR2){
    PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt for PDM_dcube_nodal_gen_init\n");
  }

  PDM_g_num_t n_vtx_x = order*nx + 1;
  PDM_g_num_t n_vtx_y = order*ny + 1;
  PDM_g_num_t n_vtx_z = order*nz + 1;


  PDM_g_num_t gn_vtx  = 0;
  //PDM_g_num_t gn_bar  = 0;
  PDM_g_num_t gn_quad = 0;
  PDM_g_num_t gn_hexa = 0;
  PDM_g_num_t gn_elt  = 0;

  if (dim == 2) {
    gn_vtx  = n_vtx_x * n_vtx_y;
    //gn_bar  = 2*(nx + ny);
    gn_quad = nx * ny;
    gn_elt = gn_quad;
  } else {
    gn_vtx  = n_vtx_x * n_vtx_y * n_vtx_z;
    //gn_bar  = 4*(nx + ny + nz);
    gn_quad = 2*(nx*ny + ny*nz + nz*nx);
    gn_hexa = nx * ny * nz;
    gn_elt = gn_hexa;
  }


  /*
   * Create the dmesh_nodal that hold the resulting mesh
   */
  PDM_g_num_t gn_cell_abs = _get_n_cell_abs(gn_elt, t_elt);
  dcube->dmesh_nodal = PDM_DMesh_nodal_create(dcube->comm,
                                              dim,
                                              gn_vtx,
                                              gn_cell_abs,  /* n_cell */
                                              -1,           /* n_face */
                                              -1);          /* n_edge */


  dcube->dn_vtx = PDM_compute_uniform_dn_entity(dcube->comm, gn_vtx);

  double* dvtx_coord = (double *) malloc( 3 * (dcube->dn_vtx ) * sizeof(double *));
  PDM_DMesh_nodal_coord_set(dcube->dmesh_nodal,
                            dcube->dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP); /* Le responsable de la mémoire est le dmesh_nodal */

  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dcube->dmesh_nodal);
  dcube->dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);

  /*
   * Generate vertices
   */
  double step_x = length / (double) (order*nx);
  double step_y = length / (double) (order*ny);
  double step_z = length / (double) (order*nz);

  if (dim == 2) {
    for (int i_vtx = 0; i_vtx < dcube->dn_vtx; ++i_vtx) {

      PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

      PDM_g_num_t indi = g_vtx % n_vtx_x;
      PDM_g_num_t indj = g_vtx / n_vtx_x;

      dvtx_coord[3 * i_vtx    ] = indi * step_x + zero_x;
      dvtx_coord[3 * i_vtx + 1] = indj * step_y + zero_y;
      dvtx_coord[3 * i_vtx + 2] = zero_z;
    }
  }
  else {
    for (int i_vtx = 0; i_vtx < dcube->dn_vtx; ++i_vtx) {

      PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

      PDM_g_num_t indi = g_vtx % n_vtx_x;
      PDM_g_num_t indj = ((g_vtx - indi) / n_vtx_x) % n_vtx_y;
      PDM_g_num_t indk = g_vtx / (n_vtx_x * n_vtx_y);

      dvtx_coord[3 * i_vtx    ] = indi * step_x + zero_x;
      dvtx_coord[3 * i_vtx + 1] = indj * step_y + zero_y;
      dvtx_coord[3 * i_vtx + 2] = indk * step_z + zero_z;
    }
  }


  //dcube->distrib_bar = PDM_compute_uniform_entity_distribution(dcube->comm, gn_bar);
  //dcube->dn_bar = (int) (dcube->distrib_bar[i_rank+1] - dcube->distrib_bar[i_rank]);

  dcube->distrib_quad = PDM_compute_uniform_entity_distribution(dcube->comm, gn_quad);
  dcube->dn_quad = (int) (dcube->distrib_quad[i_rank+1] - dcube->distrib_quad[i_rank]);
  if (dim == 3) {
    dcube->distrib_hexa = PDM_compute_uniform_entity_distribution(dcube->comm, gn_hexa);
    dcube->dn_hexa = (int) (dcube->distrib_hexa[i_rank+1] - dcube->distrib_hexa[i_rank]);
  }



  PDM_dmesh_nodal_elmts_t* dmn_elmts_vol    = NULL;
  PDM_dmesh_nodal_elmts_t* dmn_elmts_surf   = NULL;
  PDM_dmesh_nodal_elmts_t* dmn_elmts_ridge  = NULL;
  PDM_dmesh_nodal_elmts_t* dmn_elmts_corner = NULL;

  dmn_elmts_corner = PDM_DMesh_nodal_elmts_create(dcube->comm, 0, -1);
  dmn_elmts_ridge  = PDM_DMesh_nodal_elmts_create(dcube->comm, 1, -1);
  dmn_elmts_surf   = PDM_DMesh_nodal_elmts_create(dcube->comm, 2, -1);
  if (dim == 3) {
    dmn_elmts_vol  = PDM_DMesh_nodal_elmts_create(dcube->comm, 3, -1);
  }

  if (dmn_elmts_vol != NULL) {
    PDM_Mesh_nodal_add_dmesh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_vol   );
  }
  if (dmn_elmts_surf != NULL) {
    PDM_Mesh_nodal_add_dmesh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_surf  );
  }
  if (dmn_elmts_ridge != NULL) {
    PDM_Mesh_nodal_add_dmesh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_ridge );
  }
  if (dmn_elmts_corner != NULL) {
    PDM_Mesh_nodal_add_dmesh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_corner);
  }

  switch (t_elt) {
    case PDM_MESH_NODAL_TRIA3:
    {
      _generate_tria_surf (dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_QUAD4:
    {
      _generate_quad_surf (dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_TETRA4:
    {
      _generate_tetra_vol (dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_PRISM6:
    {
      _generate_prism_vol (dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_HEXA8:
    {
      _generate_hexa_vol (dcube, dcube->dmesh_nodal);
    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }

  _generate_corners (dcube, dcube->dmesh_nodal);
  //_generate_ridges (dcube, dcube->dmesh_nodal);
  if (dim == 3) {
    //_generate_surfaces (dcube, dcube->dmesh_nodal);
  }

  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;
  double delta_max;
  double delta_min;

  PDM_MPI_Allreduce (&delta_t,
                     &delta_max,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     dcube->comm);

  PDM_MPI_Allreduce (&delta_t,
                     &delta_min,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MIN,
                     dcube->comm);

  if(i_rank == 0) {
    printf("[%i] PDM_dcube_nodal : duration min/max -> %12.5e %12.5e \n", n_rank, delta_min, delta_max);
  }

  return dcube;
}


PDM_dmesh_nodal_t*
PDM_dcube_nodal_gen2_dmesh_nodal_get
(
 PDM_dcube_nodal2_t  *dcube
)
{
  return dcube->dmesh_nodal;
}
