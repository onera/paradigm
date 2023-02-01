/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_mesh_location.h"
#include "pdm_mesh_location_priv.h"
#include "pdm_point_location.h"
#include "pdm_ho_location.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_doctree.h"
// #include "pdm_mpi_priv.h"

#include "pdm_binary_search.h"
#include "pdm_para_octree.h"
#include "pdm_gnum.h"
#include "pdm_sort.h"
#include "pdm_logging.h"
#include "pdm_writer.h"
#include "pdm_vtk.h"
#include "pdm_extract_part.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_gnum_location.h"
#include "pdm_order.h"
#include "pdm_unique.h"

#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_mesh_nodal_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

//#define NTIMER_MESH_LOCATION 15

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ml_timer_step_t
 *
 */

typedef enum {

  BEGIN                                   = 0,
  BUILD_BOUNDING_BOXES                    = 1,
  STORE_CONNECTIVITY                      = 2,
  EXTRACT_ENTITIES_OF_INTEREST            = 3,
  SEARCH_CANDIDATES                       = 4,
  LOAD_BALANCING                          = 5,
  COMPUTE_ELEMENTARY_LOCATIONS            = 6,
  MERGE_LOCATION_DATA                     = 7,
  TRANSFER_TO_INITIAL_PARTITIONS          = 8,
  FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS = 9,
  END                                     = 10

} _ml_timer_step_t;


struct _pdm_mpi_double_int_t {
  double val;
  int    rank;
};
typedef struct _pdm_mpi_double_int_t PDM_MPI_double_int_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


static
void
end_timer_and_log_from_dt(const char* msg, PDM_MPI_Comm comm, double delta_t){

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_MPI_double_int_t l_info;

  l_info.val  = delta_t;
  l_info.rank = i_rank;

  PDM_MPI_double_int_t g_max_info;
  PDM_MPI_double_int_t g_min_info;


  PDM_MPI_Allreduce (&l_info,
                     &g_max_info,
                     1,
                     PDM_MPI_DOUBLE_INT,
                     PDM_MPI_MAXLOC,
                     comm);

  PDM_MPI_Allreduce (&l_info,
                     &g_min_info,
                     1,
                     PDM_MPI_DOUBLE_INT,
                     PDM_MPI_MINLOC,
                     comm);

  log_trace("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
           n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
  if(i_rank == 0) {
    printf("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
           n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
  }
}


/**
 *
 * \brief Compute point location
 *
 * \param [in]   ml  Pointer to \ref PDM_mesh_location object
 *
 */

static void
_store_cell_vtx
(
  PDM_mesh_location_t  *ml,
  PDM_geometry_kind_t   geom_kind
)
{

  int n_blocks = 0;
  int n_parts  = 0;
  int *blocks_id = NULL;

  if (ml->mesh_nodal != NULL) {
    n_blocks = PDM_part_mesh_nodal_n_section_in_geom_kind_get(ml->mesh_nodal, geom_kind);
    n_parts  = PDM_part_mesh_nodal_n_part_get(ml->mesh_nodal);
    blocks_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(ml->mesh_nodal, geom_kind);
  }

  int **n_vtx_per_elt = malloc (sizeof(int *) * n_parts);

  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
                                                geom_kind,
                                                ipart);

    n_vtx_per_elt[ipart] = malloc (sizeof(int) * n_elt);
    // ml->cell_vtx_idx[ipart] = malloc (sizeof(int) * (n_elt+1));
    // ml->cell_vtx_idx[ipart][0] = 0;
  }

  PDM_g_num_t n_g_cell = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
                                                geom_kind,
                                                ipart);
    int ielt = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      // int id_block = blocks_id[iblock];
      int id_section_in_geom_kind = blocks_id[iblock];
      int i_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(ml->mesh_nodal,
                                                                        geom_kind,
                                                                        id_section_in_geom_kind);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(ml->mesh_nodal,
                                                                            i_section);

      int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(ml->mesh_nodal,
                                                                   i_section,
                                                                   ipart);

      PDM_g_num_t *g_num = PDM_part_mesh_nodal_g_num_get(ml->mesh_nodal,
                                                         i_section,
                                                         ipart);

      int n_elt_block = PDM_part_mesh_nodal_section_n_elt_get(ml->mesh_nodal,
                                                              i_section,
                                                              ipart);

      for (int i = 0; i < n_elt_block; i++) {
        n_g_cell = PDM_MAX (n_g_cell, g_num[i]);
      }
      int n_vtx_elt = 0;
      switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
        n_vtx_elt = 1;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_BAR2:
        n_vtx_elt = 2;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_TRIA3:
        n_vtx_elt = 3;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
        n_vtx_elt = 4;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_PYRAMID5:
        n_vtx_elt = 5;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_PRISM6:
        n_vtx_elt = 6;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_HEXA8:
        n_vtx_elt = 8;
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][parent_num[i]] = n_vtx_elt;
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_per_elt[ipart][ielt++] = n_vtx_elt;
          }
        }
        break;
      case PDM_MESH_NODAL_POLY_2D:
        {
          int *connec_idx;
          int *connec;
          PDM_part_mesh_nodal_section_poly2d_get(ml->mesh_nodal,
                                                 i_section,
                                                 ipart,
                                                 &connec_idx,
                                                 &connec);
          if (parent_num != NULL) {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][parent_num[i]] = connec_idx[i+1] - connec_idx[i];
            }
          }
          else {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][ielt++] = connec_idx[i+1] - connec_idx[i];
            }
          }
          break;
        }
      case PDM_MESH_NODAL_POLY_3D:
        {
          int *connec_idx;
          int *connec;
          PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get(ml->mesh_nodal,
                                                                  i_section,
                                                                  ipart,
                                                                  &connec_idx,
                                                                  &connec);
          if (parent_num != NULL) {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][parent_num[i]] = connec_idx[i+1] - connec_idx[i];
            }
          }
          else {
            for (int i = 0; i < n_elt_block; i++) {
              n_vtx_per_elt[ipart][ielt++] = connec_idx[i+1] - connec_idx[i];
            }
          }
          break;
        }
      default :
        PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad Element type\n");
      }
    }

    ml->cell_vtx_idx[ipart] = malloc (sizeof(int) * (n_elt+1));
    ml->cell_vtx_idx[ipart][0] = 0;
    for (int i = 0; i < n_elt; i++) {
      ml->cell_vtx_idx[ipart][i+1] = ml->cell_vtx_idx[ipart][i] + n_vtx_per_elt[ipart][i];
    }
    ml->cell_vtx[ipart] = malloc (sizeof(int) * ml->cell_vtx_idx[ipart][n_elt]);
  }

  PDM_g_num_t _n_g_cell = 0;
  PDM_MPI_Allreduce (&n_g_cell, &_n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ml->comm);
  n_g_cell = _n_g_cell;

  for (int ipart = 0; ipart < n_parts; ipart++) {
    int ielt = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      // int id_block = blocks_id[iblock];
      int id_section_in_geom_kind = blocks_id[iblock];
      int i_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(ml->mesh_nodal,
                                                                        geom_kind,
                                                                        id_section_in_geom_kind);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(ml->mesh_nodal,
                                                                            i_section);

      int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(ml->mesh_nodal,
                                                                   i_section,
                                                                   ipart);

      int n_elt_block = PDM_part_mesh_nodal_section_n_elt_get(ml->mesh_nodal,
                                                              i_section,
                                                              ipart);


      int n_vtx_elt = 0;
      switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
      case PDM_MESH_NODAL_BAR2:
      case PDM_MESH_NODAL_TRIA3:
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
      case PDM_MESH_NODAL_PYRAMID5:
      case PDM_MESH_NODAL_PRISM6:
      case PDM_MESH_NODAL_HEXA8: {
        int *connec;
        PDM_g_num_t *numabs;
        int         *_parent_num;
        PDM_g_num_t *parent_entity_g_num;
        PDM_part_mesh_nodal_section_std_get(ml->mesh_nodal,
                                            i_section,
                                            ipart,
                                            &connec,
                                            &numabs,
                                            &_parent_num,
                                            &parent_entity_g_num);
        if (parent_num != NULL) {
          int idx2 = 0;
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx2++];
            }
          }
        }
        else {
          int idx2 = 0;
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx2++];
            }
          }
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_2D: {
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_section_poly2d_get(ml->mesh_nodal,
                                               i_section,
                                               ipart,
                                               &connec_idx,
                                               &connec);
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_3D:{
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get(ml->mesh_nodal,
                                                                i_section,
                                                                ipart,
                                                                &connec_idx,
                                                                &connec);
        if (parent_num != NULL) {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][parent_num[i]];
            int idx = ml->cell_vtx_idx[ipart][parent_num[i]];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        else {
          for (int i = 0; i < n_elt_block; i++) {
            n_vtx_elt = n_vtx_per_elt[ipart][ielt];
            int idx = ml->cell_vtx_idx[ipart][ielt++];
            int idx1 = connec_idx[i];
            for (int j = 0; j < n_vtx_elt; j++) {
              ml->cell_vtx[ipart][idx+j] = connec[idx1+j];
            }
          }
        }
        break;
      }
      default :
        PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad Element type\n");
      }
    }

    // int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
    //                                             geom_kind,
    //                                             ipart);
    // PDM_log_trace_array_int(n_vtx_per_elt[ipart],
    //                         n_elt,
    //                         "n_vtx_per_elt : ");
    // PDM_log_trace_connectivity_int(ml->cell_vtx_idx[ipart],
    //                                ml->cell_vtx[ipart],
    //                                n_elt,
    //                                "ml->cell_vtx : ");
  }

  for (int ipart = 0; ipart < n_parts; ipart++) {
    free (n_vtx_per_elt[ipart]);
  }
  free (n_vtx_per_elt);
}



static void _pmesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                          hack
 )
{
  if (hack) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      free(pmne->sections_std[i]);
    }
    free(pmne->sections_id);
  }
  free(pmne->sections_std   );
  free(pmne->sections_poly2d);
  free(pmne->sections_poly3d);
  free(pmne->n_elmts);
  free(pmne);
}


/* Fonction provisoire, en attente de fusion mesh_nodal/part_mesh_nodal */
static PDM_part_mesh_nodal_elmts_t *
_mesh_nodal_to_pmesh_nodal_elmts
(
 PDM_Mesh_nodal_t *mesh_nodal,
 PDM_MPI_Comm     location_comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank(location_comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(location_comm, &n_rank);

  PDM_part_mesh_nodal_elmts_t *pmne = NULL;

  int  n_block   = 0;
  int  n_part    = 0;
  int *blocks_id = NULL;

  int is_null = 1;

  if (mesh_nodal != NULL) {
    is_null   = 0;
    n_block   = PDM_Mesh_nodal_n_blocks_get (mesh_nodal);
    n_part    = PDM_Mesh_nodal_n_part_get   (mesh_nodal);
    blocks_id = PDM_Mesh_nodal_blocks_id_get(mesh_nodal);
  }


  /* Infer mesh dimension from mesh_nodal */
  int has_dim[4] = {0, 0, 0, 0};
  for (int iblock = 0; iblock < n_block; iblock++) {
    int id_block = blocks_id[iblock];

    PDM_Mesh_nodal_elt_t t_elt = PDM_Mesh_nodal_block_type_get(mesh_nodal,
                                                               id_block);

    int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);
    has_dim[dim] = 1;
  }
  // PDM_log_trace_array_int(has_dim, 4, "has_dim : ");

  int mesh_dimension = 0;
  for (int i = 3; i >= 0; i--) {
    if (has_dim[i] != 0) {
      for (int j = 0; j < i; j++) {
        assert(has_dim[j] == 0);
      }

      mesh_dimension = i;
    }
  }


  int send_buf[3] = {is_null, n_block, mesh_dimension};
  int *recv_buf = malloc(sizeof(int) * 3 * n_rank);

  PDM_MPI_Allgather(send_buf, 3, PDM_MPI_INT, recv_buf, 3, PDM_MPI_INT, location_comm);


  // Find smallest rank with non-null mesh_nodal
  int master = -1;
  int n_null_rank = 0;
  int *i_null_rank = malloc(sizeof(int) * n_rank);
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
  int *block_type = malloc(sizeof(int) * n_block);
  mesh_dimension = recv_buf[3*master+2];
  free(recv_buf);

  if (i_rank == master) {

    for (int iblock = 0; iblock < n_block; iblock++) {
      int id_block = blocks_id[iblock];

      block_type[iblock] = (int) PDM_Mesh_nodal_block_type_get(mesh_nodal,
                                                               id_block);
    }


    for (int dest = 0; dest < n_null_rank; dest++) {
      PDM_MPI_Send(block_type, n_block, PDM_MPI_INT, i_null_rank[dest], 1, location_comm);
    }
  }
  else if (is_null) {
    PDM_MPI_Recv(block_type, n_block, PDM_MPI_INT, master, 1, location_comm);
  }
  free(i_null_rank);


  pmne = PDM_part_mesh_nodal_elmts_create(mesh_dimension,
                                          n_part,
                                          location_comm);

  pmne->n_section        = n_block;
  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  if (is_null) {
    pmne->sections_id = malloc(sizeof(int) * n_block);
    for (int i = 0; i < n_block; i++) {
      if (block_type[i] == PDM_MESH_NODAL_POLY_2D) {
        pmne->sections_id[i] = i + PDM_BLOCK_ID_BLOCK_POLY2D;
      }
      else if (block_type[i] == PDM_MESH_NODAL_POLY_3D) {
        pmne->sections_id[i] = i + PDM_BLOCK_ID_BLOCK_POLY3D;
      }
      else {
        pmne->sections_id[i] = i;
      }
    }
  }
  else {
    pmne->sections_id = blocks_id;
  }



  for (int iblock = 0; iblock < n_block; iblock++) {

    int id_block = pmne->sections_id[iblock];

    if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {
      pmne->n_section_std++;
    }
    else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {
      id_block -= PDM_BLOCK_ID_BLOCK_POLY2D;
      pmne->n_section_poly2d++;
    }
    else {
      id_block -= PDM_BLOCK_ID_BLOCK_POLY3D;
      pmne->n_section_poly3d++;
    }

  }


  pmne->sections_std    = malloc(sizeof(PDM_Mesh_nodal_block_std_t    *) * pmne->n_section_std   );
  pmne->sections_poly2d = malloc(sizeof(PDM_Mesh_nodal_block_poly2d_t *) * pmne->n_section_poly2d);
  pmne->sections_poly3d = malloc(sizeof(PDM_Mesh_nodal_block_poly3d_t *) * pmne->n_section_poly3d);

  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  if (is_null) {

    // Add empty nodal blocks with appropriate types

    for (int isection = 0; isection < n_block; isection++) {
      PDM_Mesh_nodal_elt_t type = (PDM_Mesh_nodal_elt_t) block_type[isection];

      if (type == PDM_MESH_NODAL_POLY_2D) {
        pmne->sections_poly2d[pmne->n_section_poly2d++] = NULL;
      }
      else if (type == PDM_MESH_NODAL_POLY_3D) {
        pmne->sections_poly3d[pmne->n_section_poly3d++] = NULL;
      }
      else {
        pmne->sections_std[pmne->n_section_std++] = malloc(sizeof(PDM_Mesh_nodal_block_std_t));
        PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[pmne->n_section_std-1];
        block->t_elt = type;
        block->n_part = n_part;
        block->order = 1;
        block->ho_ordering = NULL;
      }

    }


    for (int i = 0; i < n_part; i++) {
      pmne->n_elmts[i] = 0;
    }
  }

  else {

    for (int iblock = 0; iblock < n_block; iblock++) {

      int id_block = blocks_id[iblock];

      if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {
        pmne->sections_std[pmne->n_section_std++] = mesh_nodal->blocks_std[id_block];
      }
      else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {
        id_block -= PDM_BLOCK_ID_BLOCK_POLY2D;
        pmne->sections_poly2d[pmne->n_section_poly2d++] = mesh_nodal->blocks_poly2d[id_block];
      }
      else {
        id_block -= PDM_BLOCK_ID_BLOCK_POLY3D;
        pmne->sections_poly3d[pmne->n_section_poly3d++] = mesh_nodal->blocks_poly3d[id_block];
      }

    }

    for (int i = 0; i < n_part; i++) {
      pmne->n_elmts[i] = PDM_Mesh_nodal_n_cell_get(mesh_nodal,
                                                   i);
    }
  }

  free(block_type);


  return pmne;
}



/**
 * NOTES:
 * - _points_in_element_t : separate point clouds?
 * - origin frame for elements (one per cloud?)
 *   -> _store_cell_vtx_selected (one per cloud for PDM_mesh_location_cell_vertex_get?)
 */

 static void
_dump_point_cloud
(
 const char     *name,
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *pn_pts,
 double        **ppts_coord,
 PDM_g_num_t   **ppts_ln_to_gn,
 const int       n_field,
 const char    **field_name,
 double       ***field_value
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t **g_num = ppts_ln_to_gn;

  if (1) {
    PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,
                                               n_part,
                                               PDM_FALSE,
                                               1.,
                                               comm,
                                               PDM_OWNERSHIP_USER);

    for (int ipart = 0; ipart < n_part; ipart++) {
      PDM_gnum_set_from_parents(gen_gnum,
                                ipart,
                                pn_pts[ipart],
                                ppts_ln_to_gn[ipart]);
    }

    PDM_gnum_compute(gen_gnum);

    g_num = malloc(sizeof(PDM_g_num_t *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      g_num[ipart] = PDM_gnum_get(gen_gnum,
                                  ipart);
    }

    PDM_gnum_free(gen_gnum);
  }

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var[n_field];
  for (int i = 0; i < n_field; i++) {
    id_var[i] = PDM_writer_var_create(wrt,
                                      PDM_WRITER_OFF,
                                      PDM_WRITER_VAR_SCALAR,
                                      PDM_WRITER_VAR_ELEMENTS,
                                      field_name[i]);
  }

  PDM_writer_step_beg(wrt, 0.);

  int id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom,
                                          PDM_WRITER_POINT,
                                          PDM_OWNERSHIP_USER);

  PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t  *) * n_part);
  int **connec = malloc(sizeof(int *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    // PDM_log_trace_array_long(ppts_ln_to_gn[ipart],
    //                          pn_pts[ipart],
    //                          "ppts_ln_to_gn[ipart] : ");

    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              ipart,
                              pn_pts[ipart],
                              ppts_coord[ipart],
                              g_num[ipart],
                              PDM_OWNERSHIP_USER);

    connec[ipart] = malloc(sizeof(int) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      connec[ipart][i] = i+1;
    }
    PDM_writer_geom_bloc_std_set(wrt,
                                 id_geom,
                                 id_block,
                                 ipart,
                                 pn_pts[ipart],
                                 connec[ipart],
                                 g_num[ipart]);

    val_num_part[ipart] = malloc(sizeof(PDM_real_t) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      val_num_part[ipart][i] = n_part*i_rank + ipart;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       ipart,
                       val_num_part[ipart]);

    for (int i = 0; i < n_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var[i],
                         id_geom,
                         ipart,
                         field_value[i][ipart]);
    }
  }

  PDM_writer_geom_write(wrt,
                        id_geom);


  PDM_writer_var_write(wrt,
                       id_var_num_part);
  PDM_writer_var_free(wrt,
                      id_var_num_part);

  for (int i = 0; i < n_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var[i]);
    PDM_writer_var_free(wrt,
                        id_var[i]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(val_num_part[ipart]);
    free(connec[ipart]);
    if (g_num != ppts_ln_to_gn) {
      free(g_num[ipart]);
    }
  }
  free(val_num_part);
  free(connec);
  if (g_num != ppts_ln_to_gn) {
    free(g_num);
  }

  PDM_writer_step_end(wrt);

  PDM_writer_geom_free(wrt,
                        id_geom);

  PDM_writer_free(wrt);
}



static void
_dump_pmne
(
      PDM_MPI_Comm                  comm,
const char                         *prefix,
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           n_part,
      PDM_g_num_t                 **pelt_ln_to_gn,
      int                          *pn_vtx,
      double                      **pvtx_coord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  char filename[999];

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int ipart = 0; ipart < n_part; ipart++) {

    for (int isection = 0; isection < n_section; isection++) {

      sprintf(filename, "%s_part_%d_section_%d_%3.3d.vtk",
              prefix, ipart, isection, i_rank);

      int id_section = sections_id[isection];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                             id_section);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                             id_section,
                                                             ipart);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart);

      // if (parent_num != NULL) {
      //   log_trace("section %d (type %d), ", isection, t_elt);
      //   PDM_log_trace_array_int(parent_num, n_elt, "parent_num : ");
      // }


      PDM_g_num_t *gnum = malloc(sizeof(PDM_g_num_t) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        int icell = i;
        if (parent_num != NULL) {
          icell = parent_num[i];
        }
        gnum[i] = pelt_ln_to_gn[ipart][icell];
      }


      if (t_elt == PDM_MESH_NODAL_POLY_2D) {

        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec);

        PDM_vtk_write_polydata(filename,
                               pn_vtx[ipart],
                               pvtx_coord[ipart],
                               NULL,
                               n_elt,
                               connec_idx,
                               connec,
                               gnum,
                               NULL);

      }

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {

        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        int *cell_face_idx;
        int *cell_face;
        PDM_g_num_t *face_ln_to_gn    = NULL;
        int         *_parent_num      = NULL;
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &elt_ln_to_gn,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_elt_g_num);

        // PDM_log_trace_connectivity_int(cell_face_idx,
        //                                cell_face,
        //                                n_elt,
        //                                "cell_face : ");

        PDM_vtk_write_polydata(filename,
                               pn_vtx[ipart],
                               pvtx_coord[ipart],
                               NULL,
                               n_face,
                               face_vtx_idx,
                               face_vtx,
                               face_ln_to_gn,
                               NULL);

      }

      else {
        int         *elmt_vtx                 = NULL;
        int         *_parent_num              = NULL;
        PDM_g_num_t *numabs                   = NULL;
        PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                id_section,
                                                ipart,
                                                &elmt_vtx,
                                                &numabs,
                                                &_parent_num,
                                                &parent_entitity_ln_to_gn);

        const char *field_name[] = {"parent_num"};
        PDM_vtk_write_std_elements(filename,
                                   pn_vtx[ipart],
                                   pvtx_coord[ipart],
                                   NULL,
                                   t_elt,
                                   n_elt,
                                   elmt_vtx,
                                   gnum,
                                   0,
                                   field_name,
                    (const int **) &parent_num);
      }

      free(gnum);
    }
  }
}



/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute the location of point clouds inta a mesh
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 * \param [in]   owner          OwnerShip
 *
 * \return     Pointer to \ref PDM_mesh_location object
 *
 */

PDM_mesh_location_t*
PDM_mesh_location_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
)
{

  PDM_mesh_location_t *ml = (PDM_mesh_location_t *) malloc(sizeof(PDM_mesh_location_t));

  ml->n_point_cloud = n_point_cloud;
  ml->comm = comm;
  ml->mesh_nature = mesh_nature;

  ml->shared_nodal = 0;
  ml->mesh_nodal   = NULL;
  ml->_mesh_nodal  = NULL;

  ml->mesh_dimension = -1;

  ml->use_user_extract = 0;
  ml->is_elmt_select_by_user = NULL;

  ml->point_clouds =
    (_point_cloud_t*) malloc (sizeof(_point_cloud_t) * n_point_cloud);

  ml->ptp = malloc(sizeof(PDM_part_to_part_t *) * n_point_cloud);

  ml->ptp_ownership = malloc(sizeof(PDM_ownership_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    ml->point_clouds[i].n_part = -1;
    ml->point_clouds[i].n_points = NULL;
    ml->point_clouds[i].coords = NULL;
    ml->point_clouds[i].gnum = NULL;
    ml->point_clouds[i].location = NULL;
    ml->point_clouds[i].uvw = NULL;
    ml->point_clouds[i].weights = NULL;
    ml->point_clouds[i].weights_idx = NULL;
    ml->point_clouds[i].projected_coords = NULL;
    ml->point_clouds[i].n_located = NULL;
    ml->point_clouds[i].n_un_located = NULL;
    ml->point_clouds[i].located = NULL;
    ml->point_clouds[i].un_located = NULL;

    ml->ptp[i] = NULL;
    ml->ptp_ownership[i] = PDM_OWNERSHIP_KEEP;
  }

  ml->points_in_elements = NULL;

  ml->reverse_result = 1;
  ml->tolerance = 0.;

  ml->method = PDM_MESH_LOCATION_OCTREE;
  // ml->method = PDM_MESH_LOCATION_DBBTREE;

  ml->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER_MESH_LOCATION; i++) {
    ml->times_elapsed[i] = 0.;
    ml->times_cpu[i]     = 0.;
    ml->times_cpu_u[i]   = 0.;
    ml->times_cpu_s[i]   = 0.;
  }

  ml->owner = owner;

  ml->tag_unlocated_get = 0;
  ml->tag_located_get = 0;
  ml->tag_point_location_get = 0;
  ml->tag_points_in_elt_get = 0;
  ml->tag_cell_vtx_get = 0;

  return ml;

}



/**
 *
 * \brief Get the number of located points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The number of located points
 *
 */

int
PDM_mesh_location_n_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  return pcloud->n_located[i_part];
}


/**
 *
 * \brief Get the number of unlocated points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The number of unlocated points
 *
 */
int
PDM_mesh_location_n_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  return pcloud->n_un_located[i_part];
}


/**
 *
 * \brief Get the list of unlocated points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The list of unlocated points
 *
 */
int *
PDM_mesh_location_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  ml->tag_unlocated_get = 1;
  return pcloud->un_located[i_part];
}


/**
 *
 * \brief Get the list of located points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The list of located points
 *
 */
int *
PDM_mesh_location_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);
  ml->tag_located_get = 1;
  return pcloud->located[i_part];
}

/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_location_n_part_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  n_part
)
{
  ml->point_clouds[i_point_cloud].n_part = n_part;
  ml->point_clouds[i_point_cloud].n_points =
    realloc(ml->point_clouds[i_point_cloud].n_points, n_part * sizeof(int));
  ml->point_clouds[i_point_cloud].coords =
    realloc(ml->point_clouds[i_point_cloud].coords,
            n_part * sizeof(double *));
  ml->point_clouds[i_point_cloud].gnum =
    realloc(ml->point_clouds[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));

  for (int i = 0; i < n_part; i++) {
    ml->point_clouds[i_point_cloud].n_points[i] = -1;
    ml->point_clouds[i_point_cloud].coords[i] = NULL;
    ml->point_clouds[i_point_cloud].gnum[i] = NULL;
  }

}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */
void
PDM_mesh_location_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
)
{
  ml->point_clouds[i_point_cloud].n_points[i_part] = n_points;
  ml->point_clouds[i_point_cloud].coords[i_part] = coords;
  ml->point_clouds[i_point_cloud].gnum[i_part] = gnum;

}


/**
 *
 * \brief Get a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  n_points        Number of points
 * \param [out]  coords          Point coordinates
 * \param [out]  gnum            Point global number
 *
 */
void
PDM_mesh_location_cloud_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       int                  *n_points,
       double              **coords,
       PDM_g_num_t         **gnum
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *n_points        = pcloud->n_points[i_part];
  *coords          = pcloud->coords[i_part];
  *gnum            = pcloud->gnum[i_part];
}



/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Pointer to \ref PDM_mesh_location object
 * \param [in]   mesh_nodal  Mesh nodal Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_shared_nodal_mesh_set
(
 PDM_mesh_location_t   *ml,
 PDM_part_mesh_nodal_t *mesh_nodal
)
{
  ml->mesh_nodal = mesh_nodal;
  ml->shared_nodal = 1;

  ml->mesh_dimension = mesh_nodal->mesh_dimension;

  int n_part = 0;
  if (mesh_nodal != NULL) {
    n_part = PDM_part_mesh_nodal_n_part_get(mesh_nodal);
  }

  ml->cell_vtx_idx = malloc(sizeof(PDM_l_num_t *) * n_part);
  ml->cell_vtx     = malloc(sizeof(PDM_l_num_t *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    ml->cell_vtx_idx[i_part] = NULL;
    ml->cell_vtx    [i_part] = NULL;
  }

}


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   ml             Pointer to \ref PDM_mesh_location object
 * \param [in]   n_part         Number of partitions
 *
 */

void
PDM_mesh_location_mesh_global_data_set
(
       PDM_mesh_location_t *ml,
 const int                  n_part
)
{

  assert (ml->shared_nodal == 0);

  if (ml->shared_nodal == 0) {
    if(ml->mesh_nodal != NULL) {
      PDM_part_mesh_nodal_free(ml->mesh_nodal);
    }
    int mesh_dimension = 3;//?
    ml->mesh_nodal = PDM_part_mesh_nodal_create(mesh_dimension, n_part, ml->comm);
  }

  if(ml->shared_nodal == 0) {
    ml->face_vtx_n   = malloc(sizeof(PDM_l_num_t *) * n_part);
    ml->cell_face_n  = malloc(sizeof(PDM_l_num_t *) * n_part);
  }

  ml->cell_vtx_idx = malloc(sizeof(PDM_l_num_t *) * n_part);
  ml->cell_vtx     = malloc(sizeof(PDM_l_num_t *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    if(ml->shared_nodal == 0) {
      ml->face_vtx_n  [i_part] = NULL;
      ml->cell_face_n [i_part] = NULL;
    }
    ml->cell_vtx_idx[i_part] = NULL;
    ml->cell_vtx    [i_part] = NULL;
  }

  ml->is_elmt_select_by_user = malloc(sizeof(int *) * n_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    ml->is_elmt_select_by_user[i_part] = NULL;
  }

}


/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */
void
PDM_mesh_location_part_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_face_idx,
 const int                 *cell_face,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_face,
 const int                 *face_vtx_idx,
 const int                 *face_vtx,
 const PDM_g_num_t         *face_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
)
{
  ml->mesh_dimension = 3;

  PDM_UNUSED(face_ln_to_gn);

  /*
   * Creation de mesh nodal
   */

  PDM_part_mesh_nodal_coord_set(ml->mesh_nodal,
                                i_part,
                                n_vtx,
                                coords,
                                vtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);



  // ml->face_vtx_n[i_part]  = malloc (sizeof(PDM_l_num_t) * n_face);
  // ml->cell_face_n[i_part] = malloc (sizeof(PDM_l_num_t) * n_cell);

  // for (int i = 0; i < n_face; i++) {
  //   ml->face_vtx_n[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
  // }

  // for (int i = 0; i < n_cell; i++) {
  //   ml->cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
  // }

  PDM_part_mesh_nodal_cell3d_cellface_add(ml->mesh_nodal,
                                          i_part,
                                          n_cell,
                                          n_face,
                                          face_vtx_idx,
                                          // ml->face_vtx_n[i_part],
                                          face_vtx,
                                          face_ln_to_gn,
                                          cell_face_idx,
                                          // ml->cell_face_n[i_part],
                                          cell_face,
                                          cell_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);
}

/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id                     Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                 Partition to define
 * \param [in]   n_cell                 Number of cells
 * \param [in]   is_elmt_select_by_user Flag to determine if user want or no to extract current cell
 *
 */
void
PDM_mesh_location_user_extract_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                 *is_elmt_select_by_user
)
{
  ml->is_elmt_select_by_user[i_part] = (int *) is_elmt_select_by_user;
}

/**
 *
 * \brief Set a part of a mesh (2d version)
 *
 * \param [in]   id            Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_edge_idx Index in the cell -> edge connectivity
 * \param [in]   cell_edge     cell -> edge connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_edge        Number of edges
 * \param [in]   edge_vtx_idx  Index in the edge -> vertex connectivity
 * \param [in]   edge_vtx      edge -> vertex connectivity
 * \param [in]   edge_ln_to_gn Local edge numbering to global edge numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_mesh_location_part_set_2d
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_edge_idx,
 const int                 *cell_edge,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_edge,
 const int                 *edge_vtx_idx,
 const int                 *edge_vtx,
 const PDM_g_num_t         *edge_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
)
{
  ml->mesh_dimension = 2;

  PDM_UNUSED (edge_ln_to_gn);

  /*
   * Creation de mesh nodal
   */

  PDM_part_mesh_nodal_coord_set(ml->mesh_nodal,
                                i_part,
                                n_vtx,
                                coords,
                                vtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);

  PDM_UNUSED(edge_vtx_idx);

  // ml->face_vtx_n[i_part]  = malloc (sizeof(PDM_l_num_t) * n_edge);
  // ml->cell_face_n[i_part] = malloc (sizeof(PDM_l_num_t) * n_cell);

  // PDM_l_num_t *edge_vtx_nb  = ml->face_vtx_n[i_part];
  // PDM_l_num_t *cell_edge_nb = ml->cell_face_n[i_part];

  // for (int i = 0; i < n_edge; i++) {
  //   edge_vtx_nb[i] = edge_vtx_idx[i+1] - edge_vtx_idx[i];
  // }

  // for (int i = 0; i < n_cell; i++) {
  //   cell_edge_nb[i] = cell_edge_idx[i+1] - cell_edge_idx[i];
  // }


  PDM_part_mesh_nodal_face2d_faceedge_add(ml->mesh_nodal,
                                          i_part,
                                          n_cell,
                                          n_edge,
                                          // edge_vtx_idx,
                                          // edge_vtx_nb,
                                          edge_vtx,
                                          cell_edge_idx,
                                          // cell_edge_nb,
                                          cell_edge,
                                          cell_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);
}





/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   tol             Tolerance
 *
 */
void
PDM_mesh_location_tolerance_set
(
       PDM_mesh_location_t *ml,
 const double               tol
)
{

  ml->tolerance = tol;
}


/**
 *
 * \brief Set the method for computing location
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   method          Method
 *
 */
void
PDM_mesh_location_method_set
(
       PDM_mesh_location_t        *ml,
 const PDM_mesh_location_method_t  method
)
{

  ml->method = method;
}


/**
 *
 * \brief Get point location
 *
 * \param [in]   id                    Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  n_points              Number of points in point cloud
 * \param [out]  coord                 Coordinates of points in point cloud
 * \param [out]  location              The global number of the closest element if the point is located,
 *                                     -1 otherwise
 *
 */

void
PDM_mesh_location_point_location_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       PDM_g_num_t         **location,
       double              **dist2,
       double              **projected_coord
)
{

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  _point_cloud_t *pcloud = ml->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *location        = pcloud->location[i_part];
  // TODO :Leak in python
  // *weights_idx     = pcloud->weights_idx[i_part];
  // *weights         = pcloud->weights[i_part];
  *projected_coord = pcloud->projected_coords[i_part];
  *dist2           = pcloud->dist2[i_part];

  ml->tag_point_location_get = 1;
}


/**
 *
 * \brief get cell vertex connectivity
 *
 * \param [in]   id                    Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  cell_vtx_idx          Index in (size = n_elt + 1)
 * \param [out]  cell_vtx              Cell vertex connectivity
 *
 */
void
PDM_mesh_location_cell_vertex_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
       int                 **cell_vtx_idx,
       int                 **cell_vtx
)
{

  if (!ml->reverse_result) {
    PDM_error(__FILE__, __LINE__, 0,
      "PDM_mesh_location_cell_vertex_get : Reverse results computation is disable, "
      "call PDM_mesh_location_reverse_results_enable to enable its\n");
  }

  *cell_vtx_idx = ml->cell_vtx_idx[i_part];
  *cell_vtx     = ml->cell_vtx[i_part];

  ml->tag_cell_vtx_get = 1;

}


/**
 *
 * \brief Get point list located in elements
 *
 * \param [in]   id                      Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                  Index of partition of the mesh
 * \param [in]   i_point_cloud           Index of cloud
 * \param [out]  elt_pts_inside_idx      Points index (size = n_elt + 1)
 * \param [out]  points_gnum             Points global number
 * \param [out]  points_coords           Points coordinates
 * \param [out]  points_uvw              Points parametric coordinates in elements
 * \param [out]  points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
 * \param [out]  points_weights          Interpolation weights
 * \param [out]  points_dist2            Distance element-points (dist < 0 if the point is inside)
 * \param [out]  points_projected_coords Point projection on element if the point is outside
 *
 */

void
PDM_mesh_location_points_in_elt_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
 const int                   i_point_cloud,
       int                 **elt_pts_inside_idx,
       PDM_g_num_t         **points_gnum,
       double              **points_coords,
       double              **points_uvw,
       int                 **points_weights_idx,
       double              **points_weights,
       double              **points_dist2,
       double              **points_projected_coords
)
{

  if (!ml->reverse_result) {
    PDM_error(__FILE__, __LINE__, 0,
      "PDM_mesh_location_points_in_elt_get : Reverse results computation is disable, "
      "call PDM_mesh_location_reverse_results_enable to enable its\n");
  }

  assert (ml->point_clouds != NULL);
  assert (i_point_cloud < ml->n_point_cloud);

  assert (ml->points_in_elements != NULL);
  assert (i_part < ml->points_in_elements->n_part);

  _points_in_element_t *_points_in_elements = ml->points_in_elements + i_point_cloud;

  *elt_pts_inside_idx      = _points_in_elements->pts_inside_idx[i_part];
  *points_gnum             = _points_in_elements->gnum[i_part];
  *points_coords           = _points_in_elements->coords[i_part];
  *points_uvw              = _points_in_elements->uvw[i_part];
  *points_weights_idx      = _points_in_elements->weights_idx[i_part];
  *points_weights          = _points_in_elements->weights[i_part];
  *points_dist2            = _points_in_elements->dist2[i_part];
  *points_projected_coords = _points_in_elements->projected_coords[i_part];

  ml->tag_points_in_elt_get = 1;
}


/**
 *
 * \brief Free a mesh location structure
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */
void
PDM_mesh_location_free
(
 PDM_mesh_location_t  *ml
)
{

  /* Free point clouds */

  if (ml->point_clouds != NULL) {
    for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
      _point_cloud_t *pcloud = ml->point_clouds + icloud;

      if (ml->reverse_result && ml->points_in_elements != NULL) {
        _points_in_element_t *_points_in_elements = ml->points_in_elements + icloud;
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_points_in_elt_get)) {
          for (int i_part = 0; i_part < _points_in_elements->n_part; ++i_part) {
            free (_points_in_elements->pts_inside_idx[i_part]);
            free (_points_in_elements->gnum[i_part]);
            free (_points_in_elements->uvw[i_part]);
            free (_points_in_elements->coords[i_part]);
            free (_points_in_elements->projected_coords[i_part]);
            free (_points_in_elements->weights_idx[i_part]);
            free (_points_in_elements->weights[i_part]);
            free (_points_in_elements->dist2[i_part]);
          }
        }
        free (_points_in_elements->pts_inside_idx);
        free (_points_in_elements->n_elts);
        free (_points_in_elements->gnum);
        free (_points_in_elements->uvw);
        free (_points_in_elements->coords);
        free (_points_in_elements->projected_coords);
        free (_points_in_elements->weights_idx);
        free (_points_in_elements->weights);
        free (_points_in_elements->dist2);
        // free (_points_in_elements);
      }

      if (pcloud->n_points != NULL) {
        free (pcloud->n_points);
      }

      if (pcloud->coords != NULL) {
        free (pcloud->coords);
      }

      if (pcloud->gnum != NULL) {
        free (pcloud->gnum);
      }

      if (pcloud->location != NULL) {
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_point_location_get)) {
          for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
            if (pcloud->location[ipart] != NULL) {
              free (pcloud->location[ipart]);
            }
            if (pcloud->dist2[ipart] != NULL) {
              free (pcloud->dist2[ipart]);
            }
            if (pcloud->projected_coords[ipart] != NULL) {
              free (pcloud->projected_coords[ipart]);
            }
          }
        }
        free (pcloud->location);
        free (pcloud->dist2);
        free (pcloud->projected_coords);
      }

      if (pcloud->uvw != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->uvw[ipart] != NULL) {
            free (pcloud->uvw[ipart]);
          }
        }
        free (pcloud->uvw);
      }

      if (pcloud->weights_idx != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->weights_idx[ipart] != NULL) {
            free (pcloud->weights_idx[ipart]);
          }
        }
        free (pcloud->weights_idx);
      }

      if (pcloud->weights != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->weights[ipart] != NULL) {
            free (pcloud->weights[ipart]);
          }
        }
        free (pcloud->weights);
      }


      if (pcloud->n_located != NULL) {
        free (pcloud->n_located);
      }

      if (pcloud->n_un_located != NULL) {
        free (pcloud->n_un_located);
      }

      if (pcloud->located != NULL) {

        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_located_get)) {
          for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
            if (pcloud->located[ipart] != NULL) {
              free (pcloud->located[ipart]);
            }
          }
        }
        free (pcloud->located);
        pcloud->located = NULL;
      }

      if (pcloud->un_located != NULL) {
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_unlocated_get)) {
          for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
            if (pcloud->un_located[ipart] != NULL) {
              free (pcloud->un_located[ipart]);
            }
          }
        }
        free (pcloud->un_located);
        pcloud->un_located = NULL;
      }

    }
    if (ml->points_in_elements != NULL) {
      free (ml->points_in_elements);
    }
    ml->points_in_elements = NULL;

    free (ml->point_clouds);
    ml->point_clouds = NULL;
  }

  /* Free mesh nodal */
  //PDM_part_mesh_nodal_partial_free(ml->mesh_nodal);?

  if (ml->mesh_nodal != NULL) {
    int _n_part = PDM_part_mesh_nodal_n_part_get(ml->mesh_nodal);

    if (ml->cell_vtx_idx != NULL) {
      for (int i = 0; i< _n_part; i++) {
        if(( ml->owner == PDM_OWNERSHIP_KEEP ) ||
           ( ml->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !ml->tag_cell_vtx_get)) {
          if(ml->cell_vtx_idx[i] != NULL) {
            free(ml->cell_vtx[i]);
            free(ml->cell_vtx_idx[i]);
          }
        }
      }
      free(ml->cell_vtx);
      free(ml->cell_vtx_idx);
    }

    ml->cell_vtx_idx = NULL;
    ml->cell_vtx = NULL;

    if (!ml->shared_nodal) {

      PDM_part_mesh_nodal_free(ml->mesh_nodal);

      if(ml->cell_face_n != NULL){
        for (int i = 0; i< _n_part; i++) {
          if(ml->cell_face_n[i] != NULL) {
            free(ml->cell_face_n[i]);
          }
        }
        free (ml->cell_face_n);
      }

      if(ml->face_vtx_n != NULL){
        for (int i = 0; i< _n_part; i++) {
          if(ml->face_vtx_n[i] != NULL) {
            free(ml->face_vtx_n[i]);
          }
        }
        free (ml->face_vtx_n);
      }
    }
  }

  if(ml->is_elmt_select_by_user != NULL) {
    free(ml->is_elmt_select_by_user);
  }


  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
    if (ml->ptp_ownership[icloud] == PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(ml->ptp[icloud]);
      ml->ptp[icloud] = NULL;
    }
  }
  free(ml->ptp);
  free(ml->ptp_ownership);

  PDM_timer_free(ml->timer);

  free(ml);
}

/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 *
 */
void
PDM_mesh_location_dump_times
(
PDM_mesh_location_t *ml
)
{

  double t1 = ml->times_elapsed[END] - ml->times_elapsed[BEGIN];
  double t2 = ml->times_cpu[END] - ml->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

  double t_elaps_max[NTIMER_MESH_LOCATION];
  PDM_MPI_Allreduce (ml->times_elapsed,
                     t_elaps_max,
                     NTIMER_MESH_LOCATION,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     ml->comm);

  double t_cpu_max[NTIMER_MESH_LOCATION];
  PDM_MPI_Allreduce (ml->times_cpu,
                     t_cpu_max, NTIMER_MESH_LOCATION,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     ml->comm);

  int rank;
  PDM_MPI_Comm_rank (ml->comm, &rank);

  if (rank == 0) {

    PDM_printf( "mesh_location timer : all (elapsed and cpu) :                                   "
                " %12.5es %12.5es\n",
                t1max, t2max);

    PDM_printf( "mesh_location timer : build bounding boxes (elapsed and cpu)                    "
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BOUNDING_BOXES],
                t_cpu_max[BUILD_BOUNDING_BOXES]);

    PDM_printf( "mesh_location timer : store connectivity (elapsed and cpu) :                    "
                " %12.5es %12.5es\n",
                t_elaps_max[STORE_CONNECTIVITY],
                t_cpu_max[STORE_CONNECTIVITY]);

    PDM_printf( "mesh_location timer : extract entities of interest (elapsed and cpu) :          "
                " %12.5es %12.5es\n",
                t_elaps_max[EXTRACT_ENTITIES_OF_INTEREST],
                t_cpu_max[EXTRACT_ENTITIES_OF_INTEREST]);

    PDM_printf( "mesh_location timer : build trees + search candidates (elapsed and cpu) :       "
                " %12.5es %12.5es\n",
                t_elaps_max[SEARCH_CANDIDATES],
                t_cpu_max[SEARCH_CANDIDATES]);

    PDM_printf( "mesh_location timer : load balancing (elapsed and cpu) :                        "
                " %12.5es %12.5es\n",
                t_elaps_max[LOAD_BALANCING],
                t_cpu_max[LOAD_BALANCING]);

    PDM_printf( "mesh_location timer : compute elementary locations (elapsed and cpu) :          "
                " %12.5es %12.5es\n",
                t_elaps_max[COMPUTE_ELEMENTARY_LOCATIONS],
                t_cpu_max[COMPUTE_ELEMENTARY_LOCATIONS]);

    PDM_printf( "mesh_location timer : merge location data (elapsed and cpu) :                   "
                " %12.5es %12.5es\n",
                t_elaps_max[MERGE_LOCATION_DATA],
                t_cpu_max[MERGE_LOCATION_DATA]);

    PDM_printf( "mesh_location timer : transfer to initial partitions (elapsed and cpu) :        "
                " %12.5es %12.5es\n",
                t_elaps_max[TRANSFER_TO_INITIAL_PARTITIONS],
                t_cpu_max[TRANSFER_TO_INITIAL_PARTITIONS]);

    PDM_printf( "mesh_location timer : finalize transfer to initial partition (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS],
                t_cpu_max[FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS]);
  }
}


/**
 *
 * \brief Compute point location optim
 *
 * \param [in]   id  Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_compute
(
 PDM_mesh_location_t        *ml
)
{
  /*
   *  Parameters
   */
  const double newton_tolerance = 1e-6;
  float extraction_threshold = 0.5; // max size ratio between extracted and original meshes

  const int full_async  = 0;
  const int dbg_enabled = 0;
  const int log_timer   = 0;
  const int dim = 3;

  /* Octree parameters */
  const int octree_depth_max             = 31;
  const int octree_points_in_leaf_max    = 1;
  const int octree_build_leaf_neighbours = 0;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (ml->comm, &i_rank);
  PDM_MPI_Comm_size (ml->comm, &n_rank);


  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  PDM_MPI_Barrier (ml->comm);
  ml->times_elapsed[BEGIN] = PDM_timer_elapsed (ml->timer);
  ml->times_cpu    [BEGIN] = PDM_timer_cpu     (ml->timer);
  ml->times_cpu_u  [BEGIN] = PDM_timer_cpu_user(ml->timer);
  ml->times_cpu_s  [BEGIN] = PDM_timer_cpu_sys (ml->timer);

  b_t_elapsed = ml->times_elapsed[BEGIN];
  b_t_cpu     = ml->times_cpu    [BEGIN];
  b_t_cpu_u   = ml->times_cpu_u  [BEGIN];
  b_t_cpu_s   = ml->times_cpu_s  [BEGIN];
  PDM_timer_resume(ml->timer);


  ml->points_in_elements = malloc(sizeof(_points_in_element_t) * ml->n_point_cloud);

  /*
   *  Compute global extents of source mesh
   *  -------------------------------------
   */


  /* Create dummy part_mesh_nodal */

  // PDM_part_mesh_nodal_elmts_t *pmne = _mesh_nodal_to_pmesh_nodal_elmts(ml->mesh_nodal, ml->comm);

  int n_block = 0;
  int n_part  = 0;
  int *blocks_id = NULL;

  /* Infer geometry kind from part mesh_nodal */
  int mesh_dimension = 0;
  PDM_MPI_Allreduce(&ml->mesh_dimension, &mesh_dimension, 1, PDM_MPI_INT, PDM_MPI_MIN, ml->comm);

  PDM_geometry_kind_t geom_kind;
  switch (mesh_dimension) {
  case 3:
    geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  case 2:
    geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 1:
    geom_kind = PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 0:
    geom_kind = PDM_GEOMETRY_KIND_CORNER;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid mesh_dimension %d\n", mesh_dimension);
  }

  PDM_part_mesh_nodal_elmts_t *pmne = NULL;

  if (dbg_enabled) {
    log_trace("geom_kind = %d, mesh_dimension = %d\n", (int) geom_kind, mesh_dimension);
  }

  if (ml->mesh_nodal == NULL) {
    /* Create empy part_mesh_nodal_elmts */
    // TODO
  }
  else {
    switch (geom_kind) {
    case PDM_GEOMETRY_KIND_VOLUMIC:
      pmne = ml->mesh_nodal->volumic;
      break;
    case PDM_GEOMETRY_KIND_SURFACIC:
      pmne = ml->mesh_nodal->surfacic;
      break;
    case PDM_GEOMETRY_KIND_RIDGE:
      pmne = ml->mesh_nodal->ridge;
      break;
    case PDM_GEOMETRY_KIND_CORNER:
      pmne = ml->mesh_nodal->corner;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid geometry kind\n", (int) geom_kind);
    }

    n_part = PDM_part_mesh_nodal_n_part_get(ml->mesh_nodal);

    n_block   = PDM_part_mesh_nodal_n_section_in_geom_kind_get  (ml->mesh_nodal, geom_kind);
    blocks_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(ml->mesh_nodal, geom_kind);
  }


  int *req_pts_proj_coord = PDM_array_const_int(ml->n_point_cloud, -1);
  int *req_pts_dist2      = PDM_array_const_int(ml->n_point_cloud, -1);


  /* Build the bounding boxes of all mesh elements
    (concatenate nodal blocks for each part) */
  int                   *pn_elt      = malloc(sizeof(int                   ) * n_part);
  PDM_g_num_t          **elt_g_num   = malloc(sizeof(PDM_g_num_t          *) * n_part);
  double               **elt_extents = malloc(sizeof(double               *) * n_part);
  PDM_Mesh_nodal_elt_t **elt_type    = malloc(sizeof(PDM_Mesh_nodal_elt_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
                                                geom_kind,
                                                ipart);

    pn_elt     [ipart] = n_elt;
    elt_g_num  [ipart] = malloc(sizeof(PDM_g_num_t         ) * n_elt);
    elt_extents[ipart] = malloc(sizeof(double              ) * n_elt * 6);
    elt_type   [ipart] = malloc(sizeof(PDM_Mesh_nodal_elt_t) * n_elt);
    int idx = -1;
    for (int iblock = 0; iblock < n_block; iblock++) {
      // int id_block = blocks_id[iblock];
      int id_section_in_geom_kind = blocks_id[iblock];
      int i_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(ml->mesh_nodal,
                                                                        geom_kind,
                                                                        id_section_in_geom_kind);
      int n_elt_in_block = PDM_part_mesh_nodal_section_n_elt_get(ml->mesh_nodal,
                                                                 i_section,
                                                                 ipart);
      double *_extents = malloc(sizeof(double) * n_elt_in_block * 6);
      PDM_part_mesh_nodal_section_elt_extents_compute(ml->mesh_nodal,
                                                      i_section,
                                                      ipart,
                                                      ml->tolerance,
                                                      _extents);


      PDM_g_num_t *_elt_g_num = PDM_part_mesh_nodal_g_num_get(ml->mesh_nodal,
                                                              i_section,
                                                              ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(ml->mesh_nodal,
                                                                            i_section);

      int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(ml->mesh_nodal,
                                                                   i_section,
                                                                   ipart);
      // log_trace("iblock %d/%d : %d elt / %d\n", iblock, n_block, n_elt_in_block, n_elt);
      // if (parent_num != NULL) {
      //   PDM_log_trace_array_int(parent_num, n_elt_in_block, "parent_num : ");
      // }

      for (int ielt = 0; ielt < n_elt_in_block; ielt++) {
        //idx = ielt;
        if (parent_num != NULL) {
          idx = parent_num[ielt];
        }
        else {
          idx++;
        }
        elt_g_num[ipart][idx] = _elt_g_num[ielt];
        elt_type [ipart][idx] = t_elt;
        memcpy(elt_extents[ipart] + 6*idx, _extents + 6*ielt, sizeof(double)*6);
      }
      free(_extents);
    }
  }

  if (dbg_enabled) {
    // int     *pn_vtx     = malloc(sizeof(int     ) * n_part);
    // double **pvtx_coord = malloc(sizeof(double *) * n_part);
    // for (int ipart = 0; ipart < n_part; ipart++) {
    //   pvtx_coord[ipart] = (double *) PDM_part_mesh_nodal_vtx_coord_get(ml->mesh_nodal,
    //                                                                    ipart);

    //   pn_vtx[ipart] = PDM_part_mesh_nodal_n_vtx_get(ml->mesh_nodal,
    //                                                 ipart);
    // }
    // _dump_pmne(ml->comm,
    //            "init_pmne",
    //            pmne,
    //            n_part,
    //            elt_g_num,
    //            pn_vtx,
    //            pvtx_coord);
    // free(pn_vtx    );
    // free(pvtx_coord);
    PDM_part_mesh_nodal_dump_vtk(ml->mesh_nodal,
                                 geom_kind,
                                 "init_pmne");
  }

  if (1 && dbg_enabled) {
    PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_BIN,
                                          PDM_WRITER_TOPO_CST,
                                          PDM_WRITER_OFF,
                                          "mesh_location__source",
                                          "source_mesh",
                                          ml->comm,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);

    int id_geom = PDM_writer_geom_create_from_mesh_nodal(wrt,
                                                         "source_mesh_geom",
                                                         ml->mesh_nodal);

    int id_var_num_part = PDM_writer_var_create(wrt,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "num_part");

    int id_var_gnum = PDM_writer_var_create(wrt,
                                            PDM_WRITER_OFF,
                                            PDM_WRITER_VAR_SCALAR,
                                            PDM_WRITER_VAR_ELEMENTS,
                                            "elt_gnum");

    PDM_writer_step_beg(wrt, 0.);

    PDM_writer_geom_write(wrt,
                          id_geom);

    PDM_real_t **val_num_part = malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_gnum     = malloc (sizeof(PDM_real_t *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
                                                  geom_kind,
                                                  ipart);
      val_num_part[ipart] = malloc(sizeof(PDM_real_t) * n_elt);
      val_gnum    [ipart] = malloc(sizeof(PDM_real_t) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        val_num_part[ipart][i] = (PDM_real_t) (i_rank*n_part + ipart);
        val_gnum    [ipart][i] = (PDM_real_t) elt_g_num[ipart][i];
      }

      PDM_writer_var_set(wrt,
                         id_var_num_part,
                         id_geom,
                         ipart,
                         val_num_part[ipart]);

      PDM_writer_var_set(wrt,
                         id_var_gnum,
                         id_geom,
                         ipart,
                         val_gnum[ipart]);
    }
    PDM_writer_var_write(wrt,
                         id_var_num_part);
    PDM_writer_var_free(wrt,
                        id_var_num_part);

    PDM_writer_var_write(wrt,
                         id_var_gnum);
    PDM_writer_var_free(wrt,
                        id_var_gnum);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(val_num_part[ipart]);
      free(val_gnum    [ipart]);
    }
    free(val_num_part);
    free(val_gnum);

    PDM_writer_step_end(wrt);
    PDM_writer_free(wrt);
  }

  double l_mesh_extents[6] = {  HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                               -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal, geom_kind, ipart);
    for (int i = 0; i < n_elt; i++) {
      for (int j = 0; j < 3; j++) {
        l_mesh_extents[j  ] = PDM_MIN(l_mesh_extents[j  ], elt_extents[ipart][6*i+j  ]);
        l_mesh_extents[j+3] = PDM_MAX(l_mesh_extents[j+3], elt_extents[ipart][6*i+j+3]);
      }
    }
  }

  double g_mesh_extents[6];
  PDM_MPI_Allreduce(l_mesh_extents,   g_mesh_extents,   3,
                    PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
  PDM_MPI_Allreduce(l_mesh_extents+3, g_mesh_extents+3, 3,
                    PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

  if (dbg_enabled && i_rank == 0) {
    printf("g_mesh_extents = %f %f %f / %f %f %f\n",
           g_mesh_extents[0], g_mesh_extents[1], g_mesh_extents[2],
           g_mesh_extents[3], g_mesh_extents[4], g_mesh_extents[5]);
  }

  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed (ml->timer);
  e_t_cpu     = PDM_timer_cpu     (ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

  ml->times_elapsed[BUILD_BOUNDING_BOXES] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu    [BUILD_BOUNDING_BOXES] += e_t_cpu     - b_t_cpu;
  ml->times_cpu_u  [BUILD_BOUNDING_BOXES] += e_t_cpu_u   - b_t_cpu_u;
  ml->times_cpu_s  [BUILD_BOUNDING_BOXES] += e_t_cpu_s   - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  /*
  *  Store cell vertex connectivity
  *  ------------------------------
  */

  _store_cell_vtx(ml, geom_kind);

  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed (ml->timer);
  e_t_cpu     = PDM_timer_cpu     (ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

  ml->times_elapsed[STORE_CONNECTIVITY] += e_t_elapsed - b_t_elapsed;
  ml->times_cpu    [STORE_CONNECTIVITY] += e_t_cpu     - b_t_cpu;
  ml->times_cpu_u  [STORE_CONNECTIVITY] += e_t_cpu_u   - b_t_cpu_u;
  ml->times_cpu_s  [STORE_CONNECTIVITY] += e_t_cpu_s   - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  double dt_search_all = 0.;
  double dt_extract_all = 0.;

  double dt_step1_all = 0.;
  double dt_step2_all = 0.;
  double dt_step3_all = 0.;
  double dt_step4_all = 0.;
  double dt_step5_all = 0.;
  double dt_step6_all = 0.;

  // PDM_mpi_comm_kind_t comm_kind = PDM_MPI_COMM_KIND_P2P;
  PDM_mpi_comm_kind_t comm_kind = PDM_MPI_COMM_KIND_COLLECTIVE;

  /* Big loop on point clouds */
  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    b_t_elapsed = PDM_timer_elapsed(ml->timer);
    b_t_cpu     = PDM_timer_cpu(ml->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);
    PDM_timer_resume(ml->timer);

    if (dbg_enabled) {
      log_trace("Point cloud %d\n", icloud);
    }

    _point_cloud_t *pcloud = ml->point_clouds + icloud;

    if (dbg_enabled) {
      char name[999];
      sprintf(name, "mesh_location_point_cloud_%d", icloud);
      _dump_point_cloud(name,
                        ml->comm,
                        pcloud->n_part,
                        pcloud->n_points,
                        pcloud->coords,
                        pcloud->gnum,
                        0,
                        NULL,
                        NULL);
    }

    /*
     *  Extract points that intersect the source mesh global extents
     *  Brute force (could be accelerated using octree)
     */
    int  *n_select_pts     = malloc(sizeof(int  ) * pcloud->n_part);
    int **select_pts_l_num = malloc(sizeof(int *) * pcloud->n_part);

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {

      n_select_pts    [ipart] = 0;
      select_pts_l_num[ipart] = malloc (sizeof(int) * pcloud->n_points[ipart]);

      for (int i = 0; i < pcloud->n_points[ipart]; i++) {
        int inside = 1;
        for (int j = 0; j < 3; j++) {
          if (pcloud->coords[ipart][3*i+j] < g_mesh_extents[j  ] ||
              pcloud->coords[ipart][3*i+j] > g_mesh_extents[j+3]) {
            inside = 0;
            break;
          }
        }

        if (inside) {
          select_pts_l_num[ipart][n_select_pts[ipart]++] = i;
        }
      } // End of loop on current parition's points

      select_pts_l_num[ipart] = realloc(select_pts_l_num[ipart],
                                        sizeof(int) * n_select_pts[ipart]);

    } // End of loop on current point cloud's partitions

    /* Compute global number of selected points */
    PDM_g_num_t l_n_pts[2] = {0, 0};
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      l_n_pts[0] += pcloud->n_points[ipart];
      l_n_pts[1] += n_select_pts    [ipart];
    }

    PDM_g_num_t g_n_pts[2];
    PDM_MPI_Allreduce(l_n_pts, g_n_pts, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);


    if (g_n_pts[1] == 0) {
      /* We extracted zero points, take a shortcut */
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free(select_pts_l_num[ipart]);
      }
      free(n_select_pts);
      free(select_pts_l_num);


      /* Result from mesh PoV */
      if (ml->reverse_result) {
        _points_in_element_t *pts_in_elt = ml->points_in_elements + icloud;

        pts_in_elt->n_part = n_part;
        pts_in_elt->n_elts = malloc(sizeof(int) * n_part);
        memcpy(pts_in_elt->n_elts, pn_elt, sizeof(int) * n_part);

        pts_in_elt->pts_inside_idx   = malloc(sizeof(int         *) * n_part);
        pts_in_elt->gnum             = malloc(sizeof(PDM_g_num_t *) * n_part);
        pts_in_elt->coords           = malloc(sizeof(double      *) * n_part);
        pts_in_elt->uvw              = malloc(sizeof(double      *) * n_part);
        pts_in_elt->projected_coords = malloc(sizeof(double      *) * n_part);
        pts_in_elt->weights_idx      = malloc(sizeof(int         *) * n_part);
        pts_in_elt->weights          = malloc(sizeof(double      *) * n_part);
        pts_in_elt->dist2            = malloc(sizeof(double      *) * n_part);

        int **pts_in_elt_triplet_idx = malloc(sizeof(int *) * n_part);
        int **pts_in_elt_triplet     = malloc(sizeof(int *) * n_part);
        for (int ipart = 0; ipart < n_part; ipart++) {
          pts_in_elt->pts_inside_idx  [ipart] = PDM_array_zeros_int(pn_elt[ipart] + 1);
          pts_in_elt->gnum            [ipart] = malloc(sizeof(PDM_g_num_t) * 0);
          pts_in_elt->coords          [ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->uvw             [ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->projected_coords[ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->weights_idx     [ipart] = PDM_array_zeros_int(1);
          pts_in_elt->weights         [ipart] = malloc(sizeof(double     ) * 0);
          pts_in_elt->dist2           [ipart] = malloc(sizeof(double     ) * 0);

          pts_in_elt_triplet_idx[ipart] = PDM_array_zeros_int(1);
          pts_in_elt_triplet    [ipart] = malloc(sizeof(int) * 0);
        }


        ml->ptp[icloud] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) elt_g_num,
                                                                    (const int          *) pn_elt,
                                                                    n_part,
                                                                    (const int          *) pcloud->n_points,
                                                                    pcloud->n_part,
                                                                    (const int         **) pts_in_elt->pts_inside_idx, // size = n_elt
                                                                    (const int         **) pts_in_elt_triplet_idx,     // size = pts_inside_idx[n_elt]
                                                                    (const int         **) pts_in_elt_triplet,
                                                                    ml->comm);
        for (int ipart = 0; ipart < n_part; ipart++) {
          free(pts_in_elt_triplet_idx[ipart]);
          free(pts_in_elt_triplet    [ipart]);
        }
        free(pts_in_elt_triplet_idx);
        free(pts_in_elt_triplet    );
      }

      /* Result from point cloud PoV */
      pcloud->n_located        = malloc(sizeof(int          ) * pcloud->n_part);
      pcloud->n_un_located     = malloc(sizeof(int          ) * pcloud->n_part);
      pcloud->located          = malloc(sizeof(int         *) * pcloud->n_part);
      pcloud->un_located       = malloc(sizeof(int         *) * pcloud->n_part);
      pcloud->location         = malloc(sizeof(PDM_g_num_t *) * pcloud->n_part);
      pcloud->projected_coords = malloc(sizeof(double      *) * pcloud->n_part);
      pcloud->dist2            = malloc(sizeof(double      *) * pcloud->n_part);
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        pcloud->n_located   [ipart] = 0;
        pcloud->n_un_located[ipart] = pcloud->n_points[ipart];

        pcloud->located   [ipart] = malloc(sizeof(int) * pcloud->n_located   [ipart]);
        pcloud->un_located[ipart] = malloc(sizeof(int) * pcloud->n_un_located[ipart]);

        for (int i = 0; i < pcloud->n_points[ipart]; i++) {
          pcloud->un_located[ipart][i] = i+1;
        }

        pcloud->location        [ipart] = malloc(sizeof(PDM_g_num_t) * 0);
        pcloud->projected_coords[ipart] = malloc(sizeof(double     ) * 0);
        pcloud->dist2           [ipart] = malloc(sizeof(double     ) * 0);
      }


      continue; // move on the next point cloud
    }



    int use_extracted_pts = (g_n_pts[1] < extraction_threshold * g_n_pts[0]);

    PDM_g_num_t **select_pts_g_num_user    = NULL;
    double      **select_pts_coord         = NULL;
    int         **select_pts_init_location = NULL;
    if (use_extracted_pts) {
      if (dbg_enabled) {
        log_trace("point cloud extraction %d / %d\n", l_n_pts[1], l_n_pts[0]);
      }

      select_pts_g_num_user    = malloc(pcloud->n_part * sizeof(PDM_g_num_t * ));
      select_pts_coord         = malloc(pcloud->n_part * sizeof(double      * ));
      select_pts_init_location = malloc(pcloud->n_part * sizeof(int         * ));

      // Just extract gnum
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        select_pts_g_num_user   [ipart] = malloc(    n_select_pts[ipart] * sizeof(PDM_g_num_t));
        select_pts_coord        [ipart] = malloc(3 * n_select_pts[ipart] * sizeof(double     ));
        select_pts_init_location[ipart] = malloc(3 * n_select_pts[ipart] * sizeof(int        ));

        for (int i = 0; i < n_select_pts[ipart]; i++) {
          int j = select_pts_l_num[ipart][i];

          select_pts_init_location[ipart][3*i  ] = i_rank;
          select_pts_init_location[ipart][3*i+1] = ipart;
          select_pts_init_location[ipart][3*i+2] = j;

          select_pts_g_num_user[ipart][i] = pcloud->gnum[ipart][j];
          for (int k = 0; k < 3; k++) {
            select_pts_coord[ipart][3*i + k] = pcloud->coords[ipart][3*j + k];
          }
        }
      }
    }
    else {
      /* We keep the whole point cloud */
      if (dbg_enabled) {
        log_trace("no point cloud extraction\n");
      }
      free(n_select_pts);
      n_select_pts            = pcloud->n_points;
      select_pts_g_num_user   = pcloud->gnum;
      select_pts_coord        = pcloud->coords;

      select_pts_init_location = malloc(pcloud->n_part * sizeof(int         * ));
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        select_pts_init_location[ipart] = malloc(3 * n_select_pts[ipart] * sizeof(int        ));
        for (int i = 0; i < n_select_pts[ipart]; i++) {
          select_pts_init_location[ipart][3*i  ] = i_rank;
          select_pts_init_location[ipart][3*i+1] = ipart;
          select_pts_init_location[ipart][3*i+2] = i;
        }
      }
    }

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(select_pts_l_num[ipart]);
    }
    free(select_pts_l_num);


    /*
     *  Compute global extents of extracted point cloud
     */
    double l_pts_extents[6] = {
      HUGE_VAL, HUGE_VAL, HUGE_VAL,
      -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
    };

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      for (int i = 0; i < n_select_pts[ipart]; i++) {
        for (int j = 0; j < 3; j++) {
          l_pts_extents[j  ] = PDM_MIN(l_pts_extents[j  ], select_pts_coord[ipart][3*i + j]);
          l_pts_extents[j+3] = PDM_MAX(l_pts_extents[j+3], select_pts_coord[ipart][3*i + j]);
        }
      }
    }

    double g_pts_extents[6];
    PDM_MPI_Allreduce(l_pts_extents,   g_pts_extents,   3,
                      PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
    PDM_MPI_Allreduce(l_pts_extents+3, g_pts_extents+3, 3,
                      PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

    if (dbg_enabled && i_rank == 0) {
      printf("  g_pts_extents = %f %f %f / %f %f %f\n",
             g_pts_extents[0], g_pts_extents[1], g_pts_extents[2],
             g_pts_extents[3], g_pts_extents[4], g_pts_extents[5]);
    }



    /*
     *  Select elements whose bounding box intersect
     *  the global extents of the extracted point clouds
     *  Brute force (could be accelerated using bbtree)
     */
    int  *n_select_elt     = malloc(sizeof(int  ) * n_part);
    int **select_elt_l_num = malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
                                                  geom_kind,
                                                  ipart);

      n_select_elt[ipart] = 0;
      select_elt_l_num[ipart] = malloc(sizeof(int) * n_elt);

      for (int ielt = 0; ielt < n_elt; ielt++) {

        double *box_min = elt_extents[ipart] + 6*ielt;
        double *box_max = box_min + 3;

        int intersect = 1;
        for (int j = 0; j < 3; j++) {
          if (box_min[j] > g_pts_extents[j+3] ||
              box_max[j] < g_pts_extents[j]) {
            intersect = 0;
            break;
          }
        }

        // if (intersect &&
        //     ml->is_elmt_select_by_user[ipart]       !=  NULL &&
        //     ml->is_elmt_select_by_user[ipart][ielt] == 1) {
        //   select_elt_l_num[ipart][n_select_elt[ipart]++] = ielt;
        // }
        if (intersect) {
          select_elt_l_num[ipart][n_select_elt[ipart]++] = ielt;
        }

      } // End of loop on current part's boxes

      select_elt_l_num[ipart] = realloc(select_elt_l_num[ipart],
                                        sizeof(int) * n_select_elt[ipart]);

    } // End of loop on mesh parts

    /* Compute global number of selected elements */
    PDM_g_num_t l_n_elt[2] = {0, 0};
    for (int ipart = 0; ipart < n_part; ipart++) {
      int n_elt = PDM_part_mesh_nodal_n_elmts_get(ml->mesh_nodal,
                                                  geom_kind,
                                                  ipart);
      l_n_elt[0] += n_elt;
      l_n_elt[1] += n_select_elt[ipart];
    }

    PDM_g_num_t g_n_elt[2];
    PDM_MPI_Allreduce(l_n_elt, g_n_elt, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, ml->comm);

    int use_extracted_mesh = (g_n_elt[1] < extraction_threshold * g_n_elt[0]);

    int         **select_elt_init_location_user = NULL;
    double      **select_elt_extents            = NULL;
    PDM_g_num_t **select_elt_g_num_user         = NULL;

    if (use_extracted_mesh) {
      if (dbg_enabled) {
        // log_trace("mesh extraction %d / %d\n", l_n_elt[1], l_n_elt[0]);
        log_trace("mesh extraction %d / %d\n", g_n_elt[1], g_n_elt[0]);
      }

      select_elt_init_location_user = malloc(sizeof(int         *) * n_part);
      select_elt_extents            = malloc(sizeof(double      *) * n_part);
      select_elt_g_num_user         = malloc(sizeof(PDM_g_num_t *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        select_elt_init_location_user[ipart] = malloc(sizeof(int        ) * n_select_elt[ipart] * 3);
        select_elt_extents           [ipart] = malloc(sizeof(double     ) * n_select_elt[ipart] * 6);
        select_elt_g_num_user        [ipart] = malloc(sizeof(PDM_g_num_t) * n_select_elt[ipart]);
        for (int i = 0; i < n_select_elt[ipart]; i++) {
          int elt_id = select_elt_l_num[ipart][i];

          memcpy(select_elt_extents[ipart] + 6*i,
                 elt_extents       [ipart] + 6*elt_id,
                 sizeof(double) * 6);

          select_elt_g_num_user[ipart][i] = elt_g_num[ipart][elt_id];

          select_elt_init_location_user[ipart][3*i  ] = i_rank;
          select_elt_init_location_user[ipart][3*i+1] = ipart;
          select_elt_init_location_user[ipart][3*i+2] = elt_id;

        }
      }
    }
    else {
      if (dbg_enabled) {
        log_trace("no mesh extraction\n");
      }
      free(n_select_elt);
      n_select_elt          = pn_elt;
      select_elt_extents    = elt_extents;
      select_elt_g_num_user = elt_g_num;

      select_elt_init_location_user = malloc(sizeof(int*) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        select_elt_init_location_user[ipart] = malloc(sizeof(int) * n_select_elt[ipart] * 3);
        for (int i = 0; i < n_select_elt[ipart]; i++) {
          select_elt_init_location_user[ipart][3*i  ] = i_rank;
          select_elt_init_location_user[ipart][3*i+1] = ipart;
          select_elt_init_location_user[ipart][3*i+2] = i;
        }
      }
    }

    /*
     *  Redistribute evenly the selected points (ptb_geom)
     */
    int **weight      = malloc(sizeof(int *) * pcloud->n_part);
    int **pstride_one = malloc(sizeof(int *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      weight     [ipart] = PDM_array_const_int(n_select_pts[ipart], 1);
      pstride_one[ipart] = PDM_array_const_int(n_select_pts[ipart], 1);
    }

    /* Use same extents for Hilbert encoding of pts and boxes?? */
    PDM_part_to_block_t *ptb_pts = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                 PDM_PART_TO_BLOCK_POST_MERGE, // TODO: Merge
                                                                 1.,
                                                                 PDM_PART_GEOM_MORTON,
                                                                 select_pts_coord,
                                                                 select_pts_g_num_user,
                                                                 weight,
                                                                 n_select_pts,
                                                                 pcloud->n_part,
                                                                 ml->comm);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(weight[ipart]);
    }
    free(weight);

    int dn_pts = PDM_part_to_block_n_elt_block_get(ptb_pts);
    PDM_g_num_t *dpts_g_num_user = PDM_part_to_block_block_gnum_get(ptb_pts);


    PDM_g_num_t *distrib_pts = PDM_part_to_block_distrib_index_get(ptb_pts);
    PDM_g_num_t *dpts_g_num_geom = malloc(dn_pts * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_pts; ++i) {
      dpts_g_num_geom[i] = distrib_pts[i_rank] + i + 1;
    }

    if (dbg_enabled) {
      PDM_log_trace_array_long(dpts_g_num_user, dn_pts, "dpts_g_num_user        : ");
    }

    /* Exchange coordinates (do this with abstract distrib?) */
    int    *blk_coord_n       = NULL;
    double *tmp_blk_pts_coord = NULL;
    int request_pts_coord = -1;
    PDM_part_to_block_iexch(ptb_pts,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            3 * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            pstride_one,
                  (void **) select_pts_coord,
                            &blk_coord_n,
                  (void **) &tmp_blk_pts_coord,
                            &request_pts_coord);

    /*
     * Exchange init_location of pts --> Is almost the selected with proc / part info
     */
    int *dpts_init_location_pts_n = NULL;
    int *dpts_init_location_pts   = NULL;
    PDM_part_to_block_exch(ptb_pts,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           pstride_one,
                 (void **) select_pts_init_location,
                           &dpts_init_location_pts_n,
                 (void **) &dpts_init_location_pts);

    PDM_part_to_block_iexch_wait(ptb_pts, request_pts_coord);

    double *dpts_coord = malloc(3 * dn_pts * sizeof(double));
    int idx_read  = 0;
    for(int i = 0; i < dn_pts; ++i) {
      dpts_coord[3*i  ] = tmp_blk_pts_coord[3*idx_read  ];
      dpts_coord[3*i+1] = tmp_blk_pts_coord[3*idx_read+1];
      dpts_coord[3*i+2] = tmp_blk_pts_coord[3*idx_read+2];

      idx_read += blk_coord_n[i];
    }
    free(blk_coord_n);
    free(tmp_blk_pts_coord);

    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(pstride_one[ipart]);
    }
    free(pstride_one);

    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "mesh_location_dpts_%d_%3.3d.vtk", icloud, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                dn_pts,
                                dpts_coord,
                                dpts_g_num_user,
                                NULL);
    }

    /*
     *  Redistribute evenly the selected boxes (ptb_geom)
     */
    double **select_box_center = malloc(sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      select_box_center[ipart] = malloc(sizeof(double) * n_select_elt[ipart] * 3);
      for (int i = 0; i < n_select_elt[ipart]; i++) {
        for (int j = 0; j < 3; j++) {
          select_box_center[ipart][3*i+j] = 0.5*(select_elt_extents[ipart][6*i+j  ] +
                                                 select_elt_extents[ipart][6*i+j+3]);
        }
      }
    }

    weight = malloc(sizeof(int *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      weight[ipart] = PDM_array_const_int(n_select_elt[ipart], 1);
    }

    if(0){
      char filename[999];
      sprintf(filename, "mesh_location_extract_boxes_%d_%3.3d.vtk", icloud, i_rank);
      PDM_vtk_write_boxes(filename,
                          n_select_elt[0],
                          select_elt_extents[0],
                          select_elt_g_num_user[0]);
    }

    PDM_part_to_block_t *ptb_elt = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                 PDM_PART_TO_BLOCK_POST_CLEANUP, // A faire en merge
                                                                 1.,
                                                                 PDM_PART_GEOM_MORTON,
                                                                 select_box_center,
                                                                 select_elt_g_num_user,
                                                                 weight,
                                                                 n_select_elt,
                                                                 n_part,
                                                                 ml->comm);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free(select_box_center[ipart]);
      free(weight[ipart]);
    }
    free(select_box_center);
    free(weight);


    int dn_elt1 = PDM_part_to_block_n_elt_block_get(ptb_elt);
    PDM_g_num_t *delt_g_num_user = PDM_part_to_block_block_gnum_get(ptb_elt);
    if (dbg_enabled) {
      PDM_log_trace_array_long(delt_g_num_user, dn_elt1, "delt_g_num_user : ");
    }

    PDM_g_num_t *distrib_elt1 = PDM_part_to_block_distrib_index_get(ptb_elt);
    PDM_g_num_t *delmt_g_num_geom = malloc(dn_elt1 * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_elt1; ++i) {
      delmt_g_num_geom[i] = distrib_elt1[i_rank] + i + 1;
    }


    /* Exchange extents (do this with abstract distrib?) */
    double *delt_extents1 = NULL;
    int request_elt_extents = -1;
    PDM_part_to_block_iexch(ptb_elt,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            6 * sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) select_elt_extents,
                            NULL,
                  (void **) &delt_extents1,
                            &request_elt_extents);

    /*
     * TODO - Exchange box_init_location en merge  !!! Attention on doit merger les blocks mais garder le non mergé pour le renvoie
     */
    double *delt_init_location_user = NULL;
    int request_elt_init_location = -1;
    PDM_part_to_block_iexch(ptb_elt,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            3 * sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) select_elt_init_location_user,
                            NULL,
                  (void **) &delt_init_location_user,
                            &request_elt_init_location);
    PDM_part_to_block_iexch_wait(ptb_elt, request_elt_extents);
    PDM_part_to_block_iexch_wait(ptb_elt, request_elt_init_location);


    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "mesh_location_dboxes_%d_%3.3d.vtk", icloud, i_rank);
      PDM_vtk_write_boxes(filename,
                          dn_elt1,
                          delt_extents1,
                          delt_g_num_user);
    }

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed (ml->timer);
    e_t_cpu     = PDM_timer_cpu     (ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

    ml->times_elapsed[EXTRACT_ENTITIES_OF_INTEREST] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [EXTRACT_ENTITIES_OF_INTEREST] += e_t_cpu     - b_t_cpu;
    ml->times_cpu_u  [EXTRACT_ENTITIES_OF_INTEREST] += e_t_cpu_u   - b_t_cpu_u;
    ml->times_cpu_s  [EXTRACT_ENTITIES_OF_INTEREST] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /*
     *  Location : search candidates
     *  ----------------------------
     */

    /*
     * Get points inside bounding boxes of elements
     */

    int         *box_pts_idx   = NULL;
    PDM_g_num_t *box_pts_g_num = NULL;
    double      *box_pts_coord = NULL;

    int dn_elt2 = 0;
    PDM_g_num_t *delt_g_num_geom2 = NULL;
    int         *delt_pts_n2      = NULL;
    PDM_g_num_t *delt_pts_g_num_geom  = NULL;
    double      *delt_pts_coord2  = NULL;

    char *env_var_oct = getenv("OCTREE_SHARED");
    int use_shared_tree = 0;
    if (env_var_oct != NULL) {
      use_shared_tree = atoi(env_var_oct);
    }

    switch (ml->method) {
      case PDM_MESH_LOCATION_OCTREE: {
        /* Create octree structure */
        PDM_para_octree_t *octree = PDM_para_octree_create(1,
                                                           octree_depth_max,
                                                           octree_points_in_leaf_max,
                                                           octree_build_leaf_neighbours,
                                                           ml->comm);

        /* Set octree point cloud */
        PDM_para_octree_point_cloud_set(octree,
                                        0,
                                        dn_pts,
                                        dpts_coord,
                                        dpts_g_num_geom);

        /* Build parallel octree */
        // PDM_MPI_Barrier(ml->comm);

        if (use_shared_tree == 0) {
          PDM_para_octree_build(octree, NULL);
        }
        else {
          PDM_para_octree_build_shared(octree, NULL);
        }

        if (use_shared_tree == 0) {
          PDM_part_to_block_t *ptb_pib = NULL;
          PDM_para_octree_points_inside_boxes_block_frame(octree,
                                                          dn_elt1,
                                                          delt_extents1,
                                                          delmt_g_num_geom,
                                                          &ptb_pib,
                                                          &delt_pts_n2,
                                                          &delt_pts_g_num_geom,
                                                          &delt_pts_coord2);

          dn_elt2 = PDM_part_to_block_n_elt_block_get(ptb_pib);
          PDM_g_num_t *_g_num = PDM_part_to_block_block_gnum_get(ptb_pib);

          delt_g_num_geom2 = malloc(sizeof(PDM_g_num_t) * dn_elt2);
          memcpy(delt_g_num_geom2, _g_num, sizeof(PDM_g_num_t) * dn_elt2);

          PDM_part_to_block_free(ptb_pib);

        }
        else {
          PDM_part_to_block_t *ptb_pib = NULL;
          PDM_para_octree_points_inside_boxes_shared_block_frame(octree,
                                                                 dn_elt1,
                                                                 delt_extents1,
                                                                 delmt_g_num_geom,
                                                                 &ptb_pib,
                                                                 &delt_pts_n2,
                                                                 &delt_pts_g_num_geom,
                                                                 &delt_pts_coord2);
          dn_elt2 = PDM_part_to_block_n_elt_block_get(ptb_pib);
          PDM_g_num_t *_g_num = PDM_part_to_block_block_gnum_get(ptb_pib);

          delt_g_num_geom2 = malloc(sizeof(PDM_g_num_t) * dn_elt2);
          memcpy(delt_g_num_geom2, _g_num, sizeof(PDM_g_num_t) * dn_elt2);

          PDM_part_to_block_free(ptb_pib);
        }

        PDM_para_octree_free(octree);
        break;
      }
      case PDM_MESH_LOCATION_DOCTREE: {

        // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_KDTREE;
        PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_OCTREE;
        PDM_doctree_t *doct = PDM_doctree_create(ml->comm,
                                                 3,
                                                 1,
                                                 NULL, // global_extents
                                                 local_tree_kind);

        /* pass abstract distrib of pts to doctree? */
        int *init_location_pts = NULL;
        if (0) {
          init_location_pts = malloc(3 * dn_pts * sizeof(int));
          for(int i = 0; i < dn_pts; ++i) {
            init_location_pts[3*i  ] = i_rank;
            init_location_pts[3*i+1] = 0;
            init_location_pts[3*i+2] = i;
          }
        }
        PDM_doctree_point_set(doct,
                              0,
                              dn_pts,
                              init_location_pts,
                              dpts_g_num_geom,
                              dpts_coord);

        /* pass abstract distrib of boxes to doctree? */
        /* or get init location from ad_elt? */
        int *init_location_box = malloc(3 * dn_elt1 * sizeof(int));
        for(int i = 0; i < dn_elt1; ++i) {
          init_location_box[3*i  ] = i_rank;
          init_location_box[3*i+1] = 0;
          init_location_box[3*i+2] = i;
        }

        PDM_doctree_solicitation_set(doct,
                                     PDM_TREE_SOLICITATION_BOXES_POINTS,
                                     1,
                                     &dn_elt1,
                                     &init_location_box,
                                     &delmt_g_num_geom,
                                     &delt_extents1);

        PDM_doctree_build(doct);

        PDM_doctree_results_in_block_frame_get(doct,
                                               &dn_elt2,
                                               &delt_g_num_geom2,
                                               &delt_pts_n2,
                                               &delt_pts_g_num_geom,
                                               &delt_pts_coord2,
                                               PDM_OWNERSHIP_USER);


        PDM_doctree_dump_times(doct);
        PDM_doctree_free(doct);
        if (init_location_pts != NULL) {
          free(init_location_pts);
        }
        free(init_location_box);
        break;
      }
      case PDM_MESH_LOCATION_DBBTREE: {
        /* Compute local extents */
        double l_extents[6] = {
          HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
          -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
        };

        for (int ipart = 0; ipart < n_part; ipart++) {
          for (int i = 0; i < n_select_elt[ipart]; i++) {
            for (int j = 0; j < 3; j++) {
              l_extents[j  ] = PDM_MIN(l_extents[j  ], select_elt_extents[ipart][6*i+j]);
              l_extents[j+3] = PDM_MAX(l_extents[j+3], select_elt_extents[ipart][6*i+3+j]);
            }
          }
        }

        /* Compute global extents */
        double g_extents[6];
        PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, ml->comm);
        PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, ml->comm);

        /* Break symmetry */
        double max_range = 1e-12;
        for (int i = 0; i < 3; i++) {
          max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
        }
        for (int i = 0; i < 3; i++) {
          g_extents[i]   -= max_range * 1.1e-3;
          g_extents[i+3] += max_range * 1.0e-3;
        }

        PDM_dbbtree_t *dbbt = PDM_dbbtree_create(ml->comm, dim, g_extents);

        PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                       1,
                                                       &dn_elt1,
                                (const double      **) &delt_extents1,
                                (const PDM_g_num_t **) &delmt_g_num_geom);


        if (use_shared_tree == 0) {

          PDM_part_to_block_t *ptb_pib = NULL;
          PDM_dbbtree_points_inside_boxes_block_frame(dbbt,
                                                      dn_pts,
                                                      dpts_g_num_geom,
                                                      dpts_coord,
                                                      &ptb_pib,
                                                      &delt_pts_n2,
                                                      &delt_pts_g_num_geom,
                                                      &delt_pts_coord2,
                                                      0);

          dn_elt2 = PDM_part_to_block_n_elt_block_get(ptb_pib);
          PDM_g_num_t *_g_num = PDM_part_to_block_block_gnum_get(ptb_pib);

          delt_g_num_geom2 = malloc(sizeof(PDM_g_num_t) * dn_elt2);
          memcpy(delt_g_num_geom2, _g_num, sizeof(PDM_g_num_t) * dn_elt2);

          PDM_part_to_block_free(ptb_pib);

        }
        else {
          abort();
          // PDM_MPI_Barrier (ml->comm);
          PDM_dbbtree_points_inside_boxes_shared(dbbt,
                                                 dn_pts,
                                                 dpts_g_num_geom,
                                                 dpts_coord,
                                                 dn_elt1,
                                                 delmt_g_num_geom, // Attention faire une distribution part_to_bloc_geom dans le cas octree
                                                 &box_pts_idx,
                                                 &box_pts_g_num,
                                                 &box_pts_coord,
                                                 0);
        }

        PDM_dbbtree_free(dbbt);
        PDM_box_set_destroy(&box_set);
        break;
      }
      default: {
        PDM_error(__FILE__, __LINE__, 0,
                  "PDM_mesh_location : unknown location method %d\n", (int) ml->method);
      }
    }
    free(delt_extents1);
    free(dpts_coord);
    free(delmt_g_num_geom);
    free(dpts_g_num_geom);


    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed (ml->timer);
    e_t_cpu     = PDM_timer_cpu     (ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

    ml->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [SEARCH_CANDIDATES] += e_t_cpu     - b_t_cpu;
    ml->times_cpu_u  [SEARCH_CANDIDATES] += e_t_cpu_u   - b_t_cpu_u;
    ml->times_cpu_s  [SEARCH_CANDIDATES] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);


    if (dbg_enabled) {
      log_trace("before compression\n");
      PDM_log_trace_array_long(delt_g_num_geom2, dn_elt2, "delt_g_num_geom2 : ");
    }

    /* TODO: maybe take into account element type to yield better weights ?? */
    if (dbg_enabled) {
      int *_idx = PDM_array_new_idx_from_sizes_int(delt_pts_n2, dn_elt2);
      PDM_log_trace_connectivity_long(_idx,
                                      delt_pts_g_num_geom,
                                      dn_elt2,
                                      "delt_pts_g_num_geom : ");
      free(_idx);
    }

    /* Compress block (remove elements with zero candidate point) -> frame (2) */
    int tmp_dn_elt2 = 0;
    for (int i = 0; i < dn_elt2; i++) {
      if (delt_pts_n2[i] > 0) {
        delt_pts_n2[tmp_dn_elt2] = delt_pts_n2[i];
        delt_g_num_geom2[tmp_dn_elt2] = delt_g_num_geom2[i];
        tmp_dn_elt2++;
      }
    }
    dn_elt2 = tmp_dn_elt2;

    delt_g_num_geom2 = realloc(delt_g_num_geom2, sizeof(PDM_g_num_t) * dn_elt2);
    delt_pts_n2      = realloc(delt_pts_n2     , sizeof(int        ) * dn_elt2);

    if (dbg_enabled) {
      log_trace("after compression\n");
      PDM_log_trace_array_long(delt_g_num_geom2, dn_elt2, "delt_g_num_geom2 : ");
      int total_weight = 0;
      for (int i = 0; i < dn_elt2; i++) {
        total_weight += delt_pts_n2[i];
      }
      log_trace("total weight = %d\n", total_weight);
    }

    int *delt_pts_idx2 = PDM_array_new_idx_from_sizes_int(delt_pts_n2, dn_elt2);
    free(delt_pts_n2);

    /*
     * Comment faire le lien avec le extract part ?
     *    delt_g_num_geom2 = Numero geometrique mais requilibré fonction de la soliciation
     *    Je pense que c'est presque un block_to_block ...
     */

    PDM_block_to_part_t* btp_elmt_geom_to_elmt_user = PDM_block_to_part_create(distrib_elt1,
                                                        (const PDM_g_num_t **) &delt_g_num_geom2,
                                                                               &dn_elt2,
                                                                               1,
                                                                               ml->comm);

    /*
     * Exchange gnum
     */
    PDM_g_num_t **tmp_delt_parent_g_num2 = NULL;
    int stride_one = 1;
    PDM_block_to_part_exch(btp_elmt_geom_to_elmt_user,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           delt_g_num_user,
                           NULL,
            (void ***)    &tmp_delt_parent_g_num2);
    PDM_g_num_t *delt_parent_g_num2 = tmp_delt_parent_g_num2[0];
    free(tmp_delt_parent_g_num2);

    // TODO : Adpat when we merge elmt and take : delt_init_location_user_unified
    int **tmp_delt_init_location2 = NULL;
    PDM_block_to_part_exch(btp_elmt_geom_to_elmt_user,
                           3 * sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           delt_init_location_user,
                           NULL,
         (void ***)        &tmp_delt_init_location2);
    int *delt_init_location2 = tmp_delt_init_location2[0];
    free(tmp_delt_init_location2);
    free(delt_init_location_user);

    PDM_block_to_part_free(btp_elmt_geom_to_elmt_user);

    if(dbg_enabled) {
      PDM_log_trace_array_long(delt_parent_g_num2, dn_elt2, "delt_parent_g_num2 ::");
      PDM_log_trace_array_int (delt_init_location2, 3 * dn_elt2, "delt_init_location2 ::");
    }

    /* Extract partition associated to current frame (2) */
    if (dbg_enabled) {
      log_trace(">> PDM_extract_part_create\n");
    }
    // TODO: proper get
    // int mesh_dimension = pmne->mesh_dimension;


    PDM_extract_part_t *extrp = PDM_extract_part_create(mesh_dimension,
                                                        n_part,
                                                        1,                                 // n_part_out
                                                        PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                        PDM_SPLIT_DUAL_WITH_PTSCOTCH,      // unused in this case
                                                        PDM_FALSE,                         // compute_child_gnum
                                                        PDM_OWNERSHIP_KEEP,
                                                        ml->comm);

    PDM_extract_part_part_nodal_set(extrp, pmne);

    /* Set vtx_coord */
    for (int ipart = 0; ipart < n_part; ipart++) {
      const double *pvtx_coord = PDM_part_mesh_nodal_vtx_coord_get(ml->mesh_nodal,
                                                                   ipart);

      const int pn_vtx = PDM_part_mesh_nodal_n_vtx_get(ml->mesh_nodal,
                                                       ipart);
      PDM_g_num_t *pvtx_ln_to_gn = (PDM_g_num_t *) PDM_part_mesh_nodal_vtx_g_num_get(ml->mesh_nodal,
                                                                                     ipart);


      int n_cell = 0;
      int n_face = 0;
      int n_edge = 0;
      PDM_g_num_t *cell_g_num = NULL;
      PDM_g_num_t *face_g_num = NULL;
      PDM_g_num_t *edge_g_num = NULL;
      switch (mesh_dimension) {
        case 3: {
          n_cell = pn_elt[ipart];
          cell_g_num = elt_g_num[ipart];
          break;
        }
        case 2: {
          n_face = pn_elt[ipart];
          face_g_num = elt_g_num[ipart];
          break;
        }
        case 1: {
          n_edge = pn_elt[ipart];
          edge_g_num = elt_g_num[ipart];
          break;
        }
        default:
        PDM_error(__FILE__, __LINE__, 0, "incorrect mesh_dimension %d\n", mesh_dimension);
      }

      PDM_extract_part_part_set(extrp,
                                ipart,
                                n_cell,
                                n_face,
                                n_edge,
                                pn_vtx,
                                NULL, // pcell_face_idx[ipart],
                                NULL, // pcell_face[ipart],
                                NULL, // pface_edge_idx[ipart],
                                NULL, // pface_edge[ipart],
                                NULL, // pedge_vtx[ipart],
                                NULL, // pface_vtx_idx[ipart],
                                NULL, // pface_vtx[ipart],
                                cell_g_num,
                                face_g_num,
                                edge_g_num,
                                pvtx_ln_to_gn,
                     (double *) pvtx_coord);
    }


    PDM_mesh_entities_t entity_type = PDM_MESH_ENTITY_CELL;
    switch (mesh_dimension) {
      case 3: {
        entity_type = PDM_MESH_ENTITY_CELL;
        break;
      }
      case 2: {
        entity_type = PDM_MESH_ENTITY_FACE;
        break;
      }
      case 1: {
        entity_type = PDM_MESH_ENTITY_EDGE;
        break;
      }
      default:
      PDM_error(__FILE__, __LINE__, 0, "incorrect mesh_dimension %d\n", mesh_dimension);
    }
    PDM_extract_part_target_set(extrp,
                                0,
                                dn_elt2,
                                delt_parent_g_num2,
                                delt_init_location2);

    if (dbg_enabled) {
      log_trace(">> PDM_extract_part_compute\n");
    }
    PDM_extract_part_compute(extrp);
    if (dbg_enabled) {
      log_trace("Yeah :D\n");
    }

    PDM_part_mesh_nodal_elmts_t *extract_pmne = NULL;
    PDM_extract_part_part_mesh_nodal_get(extrp,
                                         &extract_pmne,
                                         PDM_OWNERSHIP_USER);

    PDM_part_to_part_t *ptp_elt = NULL;
    PDM_extract_part_part_to_part_get(extrp,
                                      entity_type,
                                      &ptp_elt,
                                      PDM_OWNERSHIP_USER);

    int          pextract_n_elt        = 0;
    int          pextract_n_vtx        = 0;
    double      *pextract_vtx_coord    = NULL;

    pextract_n_elt = PDM_extract_part_n_entity_get(extrp,
                                                   0,
                                                   entity_type);
    assert(pextract_n_elt == dn_elt2);

    pextract_n_vtx = PDM_extract_part_n_entity_get(extrp,
                                                   0,
                                                   PDM_MESH_ENTITY_VERTEX);

    PDM_extract_part_vtx_coord_get(extrp,
                                   0,
                                   &pextract_vtx_coord,
                                   PDM_OWNERSHIP_KEEP);

    if (dbg_enabled) {
      _dump_pmne(ml->comm,
                 "extract_pmne",
                 extract_pmne,
                 1,
                 &delt_g_num_geom2,
                 &pextract_n_vtx,
                 &pextract_vtx_coord);
    }


    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed (ml->timer);
    e_t_cpu     = PDM_timer_cpu     (ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys (ml->timer);

    ml->times_elapsed[LOAD_BALANCING] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [LOAD_BALANCING] += e_t_cpu     - b_t_cpu;
    ml->times_cpu_u  [LOAD_BALANCING] += e_t_cpu_u   - b_t_cpu_u;
    ml->times_cpu_s  [LOAD_BALANCING] += e_t_cpu_s   - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);
    // PDM_MPI_Barrier (ml->comm);

    /* Perform elementary point locations */
    double **pelt_pts_distance2   = NULL;
    double **pelt_pts_proj_coord2 = NULL;
    int    **pelt_pts_weight_idx2 = NULL;
    double **pelt_pts_weight2     = NULL;
    double **pelt_pts_uvw2        = NULL;

    PDM_point_location_nodal(extract_pmne,
                              1,
            (const double **) &pextract_vtx_coord,
            (const int    **) &delt_pts_idx2,
            (const double **) &delt_pts_coord2,
                              newton_tolerance,
                              &pelt_pts_distance2,
                              &pelt_pts_proj_coord2,
                              &pelt_pts_weight_idx2,
                              &pelt_pts_weight2,
                              &pelt_pts_uvw2);

    double *delt_pts_distance2   = pelt_pts_distance2  [0];
    double *delt_pts_proj_coord2 = pelt_pts_proj_coord2[0];
    int    *delt_pts_weight_idx2 = pelt_pts_weight_idx2[0];
    double *delt_pts_weight2     = pelt_pts_weight2    [0];
    double *delt_pts_uvw2        = pelt_pts_uvw2       [0];
    free(pelt_pts_distance2  );
    free(pelt_pts_proj_coord2);
    free(pelt_pts_weight_idx2);
    free(pelt_pts_weight2    );
    free(pelt_pts_uvw2       );

    PDM_part_mesh_nodal_elmts_free(extract_pmne);
    PDM_extract_part_free(extrp);
    free(delt_init_location2);

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[COMPUTE_ELEMENTARY_LOCATIONS] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [COMPUTE_ELEMENTARY_LOCATIONS] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [COMPUTE_ELEMENTARY_LOCATIONS] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [COMPUTE_ELEMENTARY_LOCATIONS] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);
    // PDM_MPI_Barrier (ml->comm);

    /*
     *  Pass in block frame for points to select closest element
     */
    int n_pts2 = delt_pts_idx2[dn_elt2];
    int         *part_stride = PDM_array_const_int(n_pts2, 1);
    PDM_g_num_t *part_elt_id = malloc(sizeof(PDM_g_num_t) * n_pts2);
    // double      *part_weight = malloc(sizeof(double     ) * n_pts2);
    for (int ielt = 0; ielt < dn_elt2; ielt++) {
      for (int i = delt_pts_idx2[ielt]; i < delt_pts_idx2[ielt+1]; i++) {
        part_elt_id[i] = delt_g_num_geom2[ielt];
        // part_weight[i] = 1.;
      }
    }

    PDM_g_num_t *pts_ln_to_gn = malloc(sizeof(PDM_g_num_t) * n_pts2);
    int *pts_unique_order = malloc(sizeof(int) * n_pts2);
    memcpy(pts_ln_to_gn, delt_pts_g_num_geom, sizeof(PDM_g_num_t) * n_pts2);
    int n_pts_unique = PDM_inplace_unique_long2(pts_ln_to_gn,
                                                pts_unique_order,
                                                0,
                                                n_pts2-1);
    pts_ln_to_gn = realloc(pts_ln_to_gn, sizeof(PDM_g_num_t) * n_pts_unique);
    if (dbg_enabled) {
      log_trace("%d unique pts / %d\n", n_pts_unique, n_pts2);
    }

    PDM_g_num_t *local_pts_elt_g_num = malloc(sizeof(PDM_g_num_t) * n_pts_unique);
    double      *local_pts_elt_dist2 = malloc(sizeof(double     ) * n_pts_unique);
    double      *part_weight         = malloc(sizeof(double     ) * n_pts_unique);
    for (int i = 0; i < n_pts_unique; i++) {
      local_pts_elt_dist2[i] = HUGE_VAL;
      part_weight[i] = 1.;
    }

    for (int ielt = 0; ielt < dn_elt2; ielt++) {
      for (int i = delt_pts_idx2[ielt]; i < delt_pts_idx2[ielt+1]; i++) {
        int pt_id = pts_unique_order[i];
        if (delt_pts_distance2[i] < local_pts_elt_dist2[pt_id]) {
          local_pts_elt_g_num[pt_id] = delt_g_num_geom2[ielt];
          local_pts_elt_dist2[pt_id] = delt_pts_distance2[i];
        }
      }
    }


    if (dbg_enabled) {
      for (int i = 0; i < n_pts2; i++) {
        log_trace("pt "PDM_FMT_G_NUM" (%f %f %f) : elt "PDM_FMT_G_NUM", dist = %e\n",
                  delt_pts_g_num_geom[i],
                  delt_pts_coord2[3*i], delt_pts_coord2[3*i+1], delt_pts_coord2[3*i+2],
                  part_elt_id[i], delt_pts_distance2[i]);
      }
    }

    if (dbg_enabled) {
      for (int i = 0; i < n_pts_unique; i++) {
        log_trace("  pt "PDM_FMT_G_NUM" : elt "PDM_FMT_G_NUM", dist = %e\n",
                  pts_ln_to_gn[i], local_pts_elt_g_num[i], local_pts_elt_dist2[i]);
      }
    }

    // PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
    //                                                     PDM_PART_TO_BLOCK_POST_MERGE,
    //                                                     1.,
    //                                                     &delt_pts_g_num_geom,
    //                                                     &part_weight,
    //                                                     &n_pts2,
    //                                                     1,
    //                                                     ml->comm);
    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &pts_ln_to_gn,
                                                        &part_weight,
                                                        &n_pts_unique,
                                                        1,
                                                        ml->comm);
    free(part_weight);

    int    *block_pts_elt_n     = NULL;
    double *block_pts_elt_dist2 = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 // (void **) &delt_pts_distance2,
                 (void **) &local_pts_elt_dist2,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_dist2);
    free(block_pts_elt_n);
    free(local_pts_elt_dist2);

    PDM_g_num_t *block_pts_elt_id = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &part_stride,
                 // (void **) &part_elt_id,
                 (void **) &local_pts_elt_g_num,
                           &block_pts_elt_n,
                 (void **) &block_pts_elt_id);
    free(part_elt_id);
    free(part_stride);
    free(local_pts_elt_g_num);


    /* Pick closest elt for each point in current block */
    int block_n_pts = PDM_part_to_block_n_elt_block_get(ptb);
    PDM_g_num_t *block_pts_g_num = PDM_part_to_block_block_gnum_get(ptb);

    if (dbg_enabled) {
      int k = 0;
      log_trace("block_pts_elt_id (pre) :\n");
      for (int i = 0; i < block_n_pts; i++) {
        log_trace(" "PDM_FMT_G_NUM" [%d] -> ", block_pts_g_num[i], i);
        for (int j = 0; j < block_pts_elt_n[i]; j++) {
          log_trace("(%d, %f) ", block_pts_elt_id[k], block_pts_elt_dist2[k]);
          k++;
        }
        log_trace("\n");
      }
    }

    int idx = 0;
    for (int i = 0; i < block_n_pts; i++) {
      double min_dist2 = HUGE_VAL;
      for (int j = 0; j < block_pts_elt_n[i]; j++) {
        if (block_pts_elt_dist2[idx] < min_dist2) {
          min_dist2           = block_pts_elt_dist2[idx];
          block_pts_elt_id[i] = block_pts_elt_id   [idx];
        }
        idx++;
      }
    }

    if (dbg_enabled) {
      int *block_pts_elt_idx = PDM_array_new_idx_from_sizes_int(block_pts_elt_n,
                                                                block_n_pts);
      PDM_log_trace_connectivity_long(block_pts_elt_idx,
                                      block_pts_elt_id,
                                      block_n_pts,
                                      "block_pts_elt_id (post) : ");
      free(block_pts_elt_idx);
    }
    free(block_pts_elt_dist2);


    /* Send back to elements */
    PDM_g_num_t **tmp_part_elt_id = NULL;
    PDM_part_to_block_reverse_exch(ptb,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_CST_INTERLACED,
                                   1,
                                   block_pts_elt_n,
                        (void   *) block_pts_elt_id,
                                   NULL,
                        (void ***) &tmp_part_elt_id);
    part_elt_id = tmp_part_elt_id[0];
    free(tmp_part_elt_id);
    free(block_pts_elt_n);
    free(block_pts_elt_id);
    PDM_part_to_block_free(ptb);

    // if (dbg_enabled) {
    //   PDM_log_trace_connectivity_long(delt_pts_idx2,
    //                                   part_elt_id,
    //                                   dn_elt2,
    //                                   "part_elt_id : ");
    // }

    /* Compress (get rid of false positives) */
    int         *final_elt_pts_n          = PDM_array_zeros_int(dn_elt2);
    PDM_g_num_t *final_elt_pts_g_num_geom = malloc(sizeof(PDM_g_num_t) * n_pts2                      );
    double      *final_elt_pts_coord      = malloc(sizeof(double     ) * n_pts2 * 3                  );
    double      *final_elt_pts_distance   = malloc(sizeof(double     ) * n_pts2                      );
    double      *final_elt_pts_proj_coord = malloc(sizeof(double     ) * n_pts2 * 3                  );
    int         *final_elt_pts_weight_idx = malloc(sizeof(int        ) * (n_pts2+1)                  );
    double      *final_elt_pts_weight     = malloc(sizeof(double     ) * delt_pts_weight_idx2[n_pts2]);
    double      *final_elt_pts_uvw        = malloc(sizeof(double     ) * n_pts2 * 3                  );

    final_elt_pts_weight_idx[0] = 0;
    idx = 0;
    for (int ielt = 0; ielt < dn_elt2; ielt++) {

      int idx_pt = delt_pts_idx2[ielt];
      int n_pts  = delt_pts_idx2[ielt+1] - idx_pt;
      // PDM_g_num_t *_elt_id       = part_elt_id            + idx_pt;
      PDM_g_num_t *_parent_g_num = delt_pts_g_num_geom    + idx_pt;
      double      *_dist2        = delt_pts_distance2     + idx_pt;
      double      *_coord        = delt_pts_coord2        + idx_pt*3;
      double      *_proj         = delt_pts_proj_coord2   + idx_pt*3;
      int         *_weight_idx   = delt_pts_weight_idx2   + idx_pt;
      double      *_uvw          = delt_pts_uvw2          + idx_pt*3;

      for (int i = 0; i < n_pts; i++) {

        // if (_elt_id[i] == delt_g_num_geom2[ielt]) {
        int pt_id = pts_unique_order[idx_pt+i];
        if (part_elt_id[pt_id] == delt_g_num_geom2[ielt]) {

          final_elt_pts_n[ielt]++;
          final_elt_pts_g_num_geom[idx] = _parent_g_num[i];
          final_elt_pts_distance  [idx] = _dist2[i];
          memcpy(final_elt_pts_coord      + 3*idx, _coord + 3*i, sizeof(double)*3);
          memcpy(final_elt_pts_proj_coord + 3*idx, _proj  + 3*i, sizeof(double)*3);

          int _weight_n = _weight_idx[i+1] - _weight_idx[i];
          double *_weight = delt_pts_weight2 + _weight_idx[i];
          final_elt_pts_weight_idx[idx+1] = final_elt_pts_weight_idx[idx] + _weight_n;
          memcpy(final_elt_pts_weight + final_elt_pts_weight_idx[idx],
                 _weight,
                 sizeof(double) * _weight_n);

          memcpy(final_elt_pts_uvw + 3*idx, _uvw + 3*i, sizeof(double)*3);

          idx++;

        }
      } // End of loop on current elt's pts
    } // End of loop on elts in frame 2
    free(delt_pts_idx2         );
    free(delt_pts_coord2       );
    free(delt_pts_distance2    );
    free(delt_pts_proj_coord2  );
    free(delt_pts_weight_idx2  );
    free(delt_pts_weight2      );
    free(delt_pts_uvw2         );
    free(part_elt_id           );
    free(delt_g_num_geom2);
    free(delt_pts_g_num_geom);

    free(pts_unique_order);
    free(pts_ln_to_gn);

    int final_n_pts = idx;
    final_elt_pts_g_num_geom = realloc(final_elt_pts_g_num_geom, sizeof(PDM_g_num_t) * final_n_pts);
    final_elt_pts_coord      = realloc(final_elt_pts_coord     , sizeof(double     ) * final_n_pts*3);
    final_elt_pts_distance   = realloc(final_elt_pts_distance  , sizeof(double     ) * final_n_pts);
    final_elt_pts_proj_coord = realloc(final_elt_pts_proj_coord, sizeof(double     ) * final_n_pts*3);
    final_elt_pts_weight_idx = realloc(final_elt_pts_weight_idx, sizeof(int        ) * (final_n_pts+1));
    final_elt_pts_weight     = realloc(final_elt_pts_weight    , sizeof(double     ) * final_elt_pts_weight_idx[final_n_pts]);
    final_elt_pts_uvw        = realloc(final_elt_pts_uvw       , sizeof(double     ) * final_n_pts*3);

    /*
     * Update results in user frame
     */
    PDM_block_to_part_t* btp_pts_gnum_geom_to_user = PDM_block_to_part_create(distrib_pts,
                                                      (const PDM_g_num_t **)  &final_elt_pts_g_num_geom,
                                                                              &final_n_pts,
                                                                              1,
                                                                              ml->comm);
    /*
     * Exchange gnum
     */
    PDM_g_num_t **tmp_final_elt_pts_g_num = NULL;
    PDM_block_to_part_exch(btp_pts_gnum_geom_to_user,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           dpts_g_num_user,
                           NULL,
            (void ***)    &tmp_final_elt_pts_g_num);
    PDM_g_num_t *final_elt_pts_g_num = tmp_final_elt_pts_g_num[0];
    free(tmp_final_elt_pts_g_num);
    free(final_elt_pts_g_num_geom); // No longer used

    /*
     * Exchange init_location
     */
    int **tmp_final_elt_pts_triplet_n = NULL;
    int **tmp_final_elt_pts_triplet   = NULL;
    PDM_block_to_part_exch(btp_pts_gnum_geom_to_user,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           dpts_init_location_pts_n,
                           dpts_init_location_pts,
                          &tmp_final_elt_pts_triplet_n,
            (void ***)    &tmp_final_elt_pts_triplet);
    int *final_elt_pts_triplet_n = tmp_final_elt_pts_triplet_n[0];
    int *final_elt_pts_triplet   = tmp_final_elt_pts_triplet  [0];
    free(tmp_final_elt_pts_triplet_n);
    free(tmp_final_elt_pts_triplet);

    free(dpts_init_location_pts_n);
    free(dpts_init_location_pts);

    idx_read = 0;
    int *final_elt_pts_triplet_idx = malloc((dn_elt2+1) * sizeof(int));
    final_elt_pts_triplet_idx[0] = 0;
    for(int i = 0; i < dn_elt2; ++i) {
      final_elt_pts_triplet_idx[i+1] = final_elt_pts_triplet_idx[i];
      int n_pts_in_elmt = final_elt_pts_n[i];
      for(int j = 0; j < n_pts_in_elmt; ++j) {
        final_elt_pts_triplet_idx[i+1] += final_elt_pts_triplet_n[idx_read++];
      }
    }
    assert(idx_read == final_n_pts);

    PDM_block_to_part_free(btp_pts_gnum_geom_to_user);

    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[MERGE_LOCATION_DATA] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [MERGE_LOCATION_DATA] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [MERGE_LOCATION_DATA] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [MERGE_LOCATION_DATA] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    // PDM_MPI_Barrier (ml->comm);

    /*
     *  Transfer location data from elt (current frame) to elt (user frame)
     */
    int   request_pts_in_elt_n   = -1;
    int   request_pts_g_num      = -1;
    int   request_pts_triplet    = -1;
    int   request_pts_triplet_n  = -1;
    int   request_pts_dist2      = -1;
          request_pts_coord      = -1;
    int   request_pts_proj_coord = -1;
    int   request_pts_uvw        = -1;
    int   request_pts_weight     = -1;
    int **pts_in_elt_n           = NULL;
    int **stride_pts             = NULL;
    int **stride_pts_triplet     = NULL;
    int **stride_pts_weight      = NULL;

    int *final_elt_pts_triplet_stride = NULL;
    int **pts_in_elt_triplet_n   = NULL;
    int **pts_in_elt_triplet     = NULL;

    _points_in_element_t *pts_in_elt = NULL;

    if (ml->reverse_result) {
      pts_in_elt = ml->points_in_elements + icloud;

      pts_in_elt->n_part = n_part;  // foireux?
      pts_in_elt->n_elts = malloc(sizeof(int) * n_part);
      memcpy(pts_in_elt->n_elts, pn_elt, sizeof(int) * n_part);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_CST_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(int),
                             NULL,
            (const void  **) &final_elt_pts_n,
                             NULL,
            (      void ***) &pts_in_elt_n,
                             &request_pts_in_elt_n);
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(PDM_g_num_t),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_g_num,
                             &stride_pts,
            (      void ***) &pts_in_elt->gnum,
                             &request_pts_g_num);

      /* Echanger triplets des pts pour construire le ptp elt_user <-> pts_user */
      final_elt_pts_triplet_stride = malloc(sizeof(int) * dn_elt2);
      for (int i = 0; i < dn_elt2; i++) {
        final_elt_pts_triplet_stride[i] = final_elt_pts_triplet_idx[i+1] - final_elt_pts_triplet_idx[i];
      }

      // Exchange all triplets
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(int),
            (const int   **) &final_elt_pts_triplet_stride,
            (const void  **) &final_elt_pts_triplet,
                             &stride_pts_triplet,
            (      void ***) &pts_in_elt_triplet,
                             &request_pts_triplet);

      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_g_num);


      // Exchange number of triplet per point
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(int),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_triplet_n,
                             &stride_pts,
            (      void ***) &pts_in_elt_triplet_n,
                             &request_pts_triplet_n);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_distance,
                             &stride_pts,
            (      void ***) &pts_in_elt->dist2,
                             &request_pts_dist2);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_coord,
                             &stride_pts,
            (      void ***) &pts_in_elt->coords,
                             &request_pts_coord);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_proj_coord,
                             &stride_pts,
            (      void ***) &pts_in_elt->projected_coords,
                             &request_pts_proj_coord);

      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             3*sizeof(double),
            (const int   **) &final_elt_pts_n,
            (const void  **) &final_elt_pts_uvw,
                             &stride_pts,
            (      void ***) &pts_in_elt->uvw,
                             &request_pts_uvw);
    }
    // free(final_elt_pts_triplet_n);

    // end_timer_and_log("Step 2 : ", ml->comm, t1_step2);
    // PDM_MPI_Barrier (ml->comm);

    int *final_elt_pts_idx = PDM_array_new_idx_from_sizes_int(final_elt_pts_n, dn_elt2);
    // free(final_elt_pts_n);

    int *elt_pts_weight_stride = NULL;
    if (ml->reverse_result) {
      elt_pts_weight_stride = malloc(sizeof(int) * dn_elt2);
      for (int ielt = 0; ielt < dn_elt2; ielt++) {
        int head = final_elt_pts_idx[ielt  ];
        int tail = final_elt_pts_idx[ielt+1];

        elt_pts_weight_stride[ielt] = final_elt_pts_weight_idx[tail] - final_elt_pts_weight_idx[head];
      }
      if (dbg_enabled) {
        PDM_log_trace_array_int(elt_pts_weight_stride, dn_elt2, "elt_pts_weight_stride : ");
      }
    }

    if (ml->reverse_result) {
      PDM_part_to_part_iexch(ptp_elt,
                             comm_kind,
                             PDM_STRIDE_VAR_INTERLACED,
                             PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                             1,
                             sizeof(double),
            (const int   **) &elt_pts_weight_stride,
            (const void  **) &final_elt_pts_weight,
                             &stride_pts_weight,
            (      void ***) &pts_in_elt->weights,
                             &request_pts_weight);


      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_in_elt_n);


      int  *n_ref_elt = NULL;
      int **ref_elt   = NULL;
      PDM_part_to_part_ref_lnum2_get(ptp_elt,
                                     &n_ref_elt,
                                     &ref_elt);

      int  *n_unref_elt = NULL;
      int **unref_elt   = NULL;
      PDM_part_to_part_unref_lnum2_get(ptp_elt,
                                       &n_unref_elt,
                                       &unref_elt);

      pts_in_elt->pts_inside_idx = malloc(sizeof(int *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        pts_in_elt->pts_inside_idx[ipart] = malloc(sizeof(int) * (pn_elt[ipart]+1));
        pts_in_elt->pts_inside_idx[ipart][0] = 0;

        for (int i = 0; i < n_ref_elt[ipart]; i++) {
          pts_in_elt->pts_inside_idx[ipart][ref_elt[ipart][i]] = pts_in_elt_n[ipart][i];
        }
        free(pts_in_elt_n[ipart]);

        for (int i = 0; i < n_unref_elt[ipart]; i++) {
          pts_in_elt->pts_inside_idx[ipart][unref_elt[ipart][i]] = 0;
        }

        for (int i = 0; i < pn_elt[ipart]; i++) {
          pts_in_elt->pts_inside_idx[ipart][i+1] += pts_in_elt->pts_inside_idx[ipart][i];
        }
      }
      free(pts_in_elt_n);

      pts_in_elt->weights_idx = malloc(sizeof(int *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        pts_in_elt->weights_idx[ipart] = malloc(sizeof(int) * (pts_in_elt->pts_inside_idx[ipart][pn_elt[ipart]] + 1));
        pts_in_elt->weights_idx[ipart][0] = 0;

        for (int ielt = 0; ielt < pn_elt[ipart]; ielt++) {
          int n_vtx = ml->cell_vtx_idx[ipart][ielt+1] - ml->cell_vtx_idx[ipart][ielt];
          for (int idx_pt = pts_in_elt->pts_inside_idx[ipart][ielt]; idx_pt < pts_in_elt->pts_inside_idx[ipart][ielt+1]; idx_pt++) {
            pts_in_elt->weights_idx[ipart][idx_pt+1] = pts_in_elt->weights_idx[ipart][idx_pt] + n_vtx;
          }
        }
      }

      if (!full_async) {
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_dist2);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_coord);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_proj_coord);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_weight);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_uvw);
        // for (int ipart = 0; ipart < n_part; ipart++) {
        //   free(stride_pts_weight[ipart]);
        // }
        // free(stride_pts_weight);
      }

      free(elt_pts_weight_stride);
    }

    /*
     *  Create ptp to exchange data between elt and points (both in user frame)
     */

    if (ml->reverse_result) {

      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_triplet);
      free(final_elt_pts_triplet_stride);
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_triplet_n);
      free(final_elt_pts_triplet_n);

      int  *n_ref_elt = NULL;
      int **ref_elt   = NULL;
      PDM_part_to_part_ref_lnum2_get(ptp_elt,
                                     &n_ref_elt,
                                     &ref_elt);

      int  *n_unref_elt = NULL;
      int **unref_elt   = NULL;
      PDM_part_to_part_unref_lnum2_get(ptp_elt,
                                       &n_unref_elt,
                                       &unref_elt);

      int **pts_in_elt_triplet_idx = malloc(sizeof(int *) * n_part);
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(stride_pts_triplet[ipart]);

        int _n_pts = pts_in_elt->pts_inside_idx[ipart][pn_elt[ipart]];

        // PDM_log_trace_array_int(pts_in_elt_triplet_n[ipart],
        //                         _n_pts,
        //                         "pts_in_elt_triplet_n : ");

        pts_in_elt_triplet_idx[ipart] = malloc(sizeof(int) * (_n_pts + 1));
        pts_in_elt_triplet_idx[ipart][0] = 0;
        int max_triplet_n = 0;
        for (int i = 0; i < _n_pts; i++) {
          max_triplet_n = PDM_MAX(max_triplet_n, pts_in_elt_triplet_n[ipart][i]);

          pts_in_elt_triplet_idx[ipart][i+1] = pts_in_elt_triplet_idx[ipart][i] +
          pts_in_elt_triplet_n[ipart][i]*3;
        }

        if (0) {
          int *_idx = malloc(sizeof(int) * (_n_pts + 1));
          for (int i = 0; i <= _n_pts; i++) {
            _idx[i] = pts_in_elt_triplet_idx[ipart][i]/3;
          }

          PDM_log_trace_array_long(pts_in_elt->gnum[ipart],
                                   _n_pts,
                                   "pts_in_elt_gnum : ");

          PDM_log_trace_graph_nuplet_int(_idx,
                                         pts_in_elt_triplet[ipart],
                                         3,
                                         _n_pts,
                                         "pts_in_elt_triplet      : ");
          free(_idx);
        }



        /* Lexicographic sort each point's triplets */
        int *order = malloc(sizeof(int) * max_triplet_n);
        for (int i = 0; i < _n_pts; i++) {
          PDM_order_lnum_s(pts_in_elt_triplet[ipart] + pts_in_elt_triplet_idx[ipart][i],
                           3,
                           order,
                           pts_in_elt_triplet_n[ipart][i]);
        }
        free(order);

        free(pts_in_elt_triplet_n[ipart]);

      }
      free(stride_pts_triplet);
      free(pts_in_elt_triplet_n);
      free(final_elt_pts_triplet);
      free(final_elt_pts_triplet_idx);
      free(final_elt_pts_idx);
      free(final_elt_pts_g_num);

      // PDM_part_to_part_free(ptp_elt);

      // log_trace("PDM_part_to_part_create_from_num2_triplet start \n");
      ml->ptp[icloud] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) elt_g_num,
                                                                  (const int          *) pn_elt,
                                                                  n_part,
                                                                  (const int          *) pcloud->n_points,
                                                                  pcloud->n_part,
                                                                  (const int         **) pts_in_elt->pts_inside_idx, // size = n_elt
                                                                  (const int         **) pts_in_elt_triplet_idx,     // size = pts_inside_idx[n_elt]
                                                                  (const int         **) pts_in_elt_triplet,
                                                                  ml->comm);

      for (int ipart = 0; ipart < n_part; ipart++) {
        free(pts_in_elt_triplet    [ipart]);
        free(pts_in_elt_triplet_idx[ipart]);
      }
      free(pts_in_elt_triplet);
      free(pts_in_elt_triplet_idx);

    } // end ml->reverse_result
    // TO DO: move further to encompass all ml->ptp[icloud]

    // for (int ipart = 0; ipart < n_part; ipart++) {
    //   free(stride_pts[ipart]);
    // }
    // free(stride_pts);

    // end_timer_and_log("Step 4 : ", ml->comm, t1_step4);
    // PDM_MPI_Barrier (ml->comm);

    // TODO: ownership on located/unlocated??
    int  *n_located = NULL;
    int **located   = NULL;
    PDM_part_to_part_ref_lnum2_get(ml->ptp[icloud],
                                   &n_located,
                                   &located);

    int  *n_unlocated = NULL;
    int **unlocated   = NULL;
    PDM_part_to_part_unref_lnum2_get(ml->ptp[icloud],
                                     &n_unlocated,
                                     &unlocated);

    pcloud->n_located    = malloc(sizeof(int  ) * pcloud->n_part);
    pcloud->n_un_located = malloc(sizeof(int  ) * pcloud->n_part);
    pcloud->located      = malloc(sizeof(int *) * pcloud->n_part);
    pcloud->un_located   = malloc(sizeof(int *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->n_located   [ipart] = n_located  [ipart];
      pcloud->n_un_located[ipart] = n_unlocated[ipart];

      pcloud->located   [ipart] = malloc(sizeof(int) * n_located  [ipart]);
      pcloud->un_located[ipart] = malloc(sizeof(int) * n_unlocated[ipart]);

      memcpy(pcloud->located   [ipart], located  [ipart], sizeof(int) * n_located  [ipart]);
      memcpy(pcloud->un_located[ipart], unlocated[ipart], sizeof(int) * n_unlocated[ipart]);
    }

    // TODO: ownership on gnum1_come_from(_idx)??
    int         **gnum1_come_from_idx = NULL;
    PDM_g_num_t **gnum1_come_from     = NULL;
    PDM_part_to_part_gnum1_come_from_get(ml->ptp[icloud],
                                         &gnum1_come_from_idx,
                                         &gnum1_come_from);

    pcloud->location = malloc(sizeof(PDM_g_num_t *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart] = malloc(sizeof(PDM_g_num_t) * n_located[ipart]);
      for (int i = 0; i < n_located[ipart]; i++) {

        if (gnum1_come_from_idx[ipart][i+1] != gnum1_come_from_idx[ipart][i] + 1) {
          int ipt = pcloud->located[ipart][i] - 1;
          log_trace("point "PDM_FMT_G_NUM" has locations ",
                    pcloud->gnum[ipart][ipt]);
          PDM_log_trace_array_long(gnum1_come_from[ipart] + gnum1_come_from_idx[ipart][i],
                                   gnum1_come_from_idx[ipart][i+1] - gnum1_come_from_idx[ipart][i],
                                   "");
        }
        assert(gnum1_come_from_idx[ipart][i+1] == gnum1_come_from_idx[ipart][i] + 1);

        pcloud->location[ipart][i] = gnum1_come_from[ipart][i];
      }

      if (dbg_enabled) {
        PDM_log_trace_array_long(pcloud->location[ipart], n_located[ipart], "pcloud->location[ipart] : ");
      }
    }


    if (dbg_enabled) {
      double **is_located = malloc(sizeof(double) * pcloud->n_part);
      double **location   = malloc(sizeof(double) * pcloud->n_part);
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        is_located[ipart] = malloc(sizeof(double) * pcloud->n_points[ipart]);
        location  [ipart] = malloc(sizeof(double) * pcloud->n_points[ipart]);
        for (int i = 0; i < pcloud->n_points[ipart]; i++) {
          is_located[ipart][i] = -1;
        }
        for (int i = 0; i < n_located[ipart]; i++) {
          is_located[ipart][located[ipart][i]-1] = 1;
          location  [ipart][located[ipart][i]-1] = (double) pcloud->location[ipart][i];
        }
        for (int i = 0; i < n_unlocated[ipart]; i++) {
          is_located[ipart][unlocated[ipart][i]-1] = 0;
          location  [ipart][unlocated[ipart][i]-1] = -1;
        }
      }

      const char  *field_name[]   = {"is_located", "location"};
      double     **field_value[2] = {is_located, location};

      char name[999];
      sprintf(name, "mesh_location_point_cloud_%d_loc", icloud);
      _dump_point_cloud(name,
                        ml->comm,
                        pcloud->n_part,
                        pcloud->n_points,
                        pcloud->coords,
                        pcloud->gnum,
                        2,
                        field_name,
                        field_value);

      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free(is_located[ipart]);
        free(location  [ipart]);
      }
      free(is_located);
      free(location);
    }

    if (full_async) {
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_proj_coord);
    }
    /* Exchange other fields (dist2, uvw, weights(_idx), proj_coord) */
    /* Really necessary from the points' PoV ?? */
    int **tmp_projected_coords_n = NULL;
    PDM_part_to_part_iexch(ml->ptp[icloud],
                           comm_kind,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           3*sizeof(double),
                           NULL,
          (const void  **) pts_in_elt->projected_coords,
                           &tmp_projected_coords_n,
          (      void ***) &pcloud->projected_coords,
                           &req_pts_proj_coord[icloud]);

    // int request_pts_dist2;
    if (full_async) {
      PDM_part_to_part_iexch_wait(ptp_elt, request_pts_dist2);
    }

    int **tmp_pts_dist2_n = NULL;
    PDM_part_to_part_iexch(ml->ptp[icloud],
                           comm_kind,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           sizeof(double),
                           NULL,
          (const void  **) pts_in_elt->dist2,
                           &tmp_pts_dist2_n,
          (      void ***) &pcloud->dist2,
                           &req_pts_dist2[icloud]);

    // PDM_part_to_part_iexch_wait(ml->ptp[icloud], request_pts_proj_coord);
    // PDM_part_to_part_iexch_wait(ml->ptp[icloud], request_pts_dist2);

    if (ptp_elt != NULL) {
      if (full_async) {
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_uvw);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_weight);
        PDM_part_to_part_iexch_wait(ptp_elt, request_pts_coord);
      }
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(stride_pts_weight[ipart]);
        free(stride_pts       [ipart]);
      }
      free(stride_pts_weight);
      free(stride_pts       );
      PDM_part_to_part_free(ptp_elt);
    }


    PDM_MPI_Barrier (ml->comm);
    PDM_timer_hang_on(ml->timer);
    e_t_elapsed = PDM_timer_elapsed(ml->timer);
    e_t_cpu     = PDM_timer_cpu(ml->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

    ml->times_elapsed[TRANSFER_TO_INITIAL_PARTITIONS] += e_t_elapsed - b_t_elapsed;
    ml->times_cpu    [TRANSFER_TO_INITIAL_PARTITIONS] += e_t_cpu - b_t_cpu;
    ml->times_cpu_u  [TRANSFER_TO_INITIAL_PARTITIONS] += e_t_cpu_u - b_t_cpu_u;
    ml->times_cpu_s  [TRANSFER_TO_INITIAL_PARTITIONS] += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(ml->timer);

    /* Free memory */

    /* See if stuff can be freed earlier.... */
    free(delt_parent_g_num2      );
    free(final_elt_pts_coord     );
    free(final_elt_pts_distance  );
    free(final_elt_pts_proj_coord);
    free(final_elt_pts_weight_idx);
    free(final_elt_pts_weight    );
    free(final_elt_pts_uvw       );
    free(final_elt_pts_n         );


    if (use_extracted_pts) {
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
        free(select_pts_g_num_user[ipart]);
        free(select_pts_coord     [ipart]);
      }
      free(select_pts_g_num_user);
      free(select_pts_coord);
      free(n_select_pts);
    }
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      free(select_pts_init_location[ipart]);
    }
    free(select_pts_init_location);


    if (use_extracted_mesh) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(select_elt_extents   [ipart]);
        free(select_elt_g_num_user[ipart]);
      }
      free(select_elt_extents);
      free(select_elt_g_num_user);
      free(n_select_elt);
    }
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(select_elt_l_num[ipart]);
      free(select_elt_init_location_user[ipart]);
    }
    free(select_elt_l_num);
    free(select_elt_init_location_user);

    PDM_part_to_block_free(ptb_pts);
    PDM_part_to_block_free(ptb_elt);

  } /* End icloud */

  // _pmesh_nodal_elmts_free(pmne,
  //                         (ml->mesh_nodal == NULL));
  if (ml->mesh_nodal == NULL) {
    PDM_part_mesh_nodal_elmts_free(pmne);
  }

  for (int icloud = 0; icloud < ml->n_point_cloud; icloud++) {
    if (ml->ptp[icloud] != NULL) {
      if (req_pts_proj_coord[icloud] != -1) {
        // log_trace("out wait req_pts_proj_coord = %i \n", req_pts_proj_coord[icloud]);
        PDM_part_to_part_iexch_wait(ml->ptp[icloud], req_pts_proj_coord[icloud]);
        // log_trace("out wait req_pts_proj_coord END = %i \n", req_pts_proj_coord[icloud]);
      }

      if (req_pts_dist2[icloud] != -1) {
        // log_trace("out wait req_pts_dist2 = %i \n", req_pts_dist2[icloud]);
        PDM_part_to_part_iexch_wait(ml->ptp[icloud], req_pts_dist2[icloud]);
        // log_trace("out wait req_pts_dist2 END = %i \n", req_pts_dist2[icloud]);
      }
    }
  }

  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  e_t_elapsed = PDM_timer_elapsed(ml->timer);
  e_t_cpu     = PDM_timer_cpu(ml->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(ml->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(ml->timer);

  ml->times_elapsed[FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS] = e_t_elapsed - b_t_elapsed;
  ml->times_cpu    [FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS] = e_t_cpu - b_t_cpu;
  ml->times_cpu_u  [FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS] = e_t_cpu_u - b_t_cpu_u;
  ml->times_cpu_s  [FINALIZE_TRANSFER_TO_INITIAL_PARTITIONS] = e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(ml->timer);

  free(req_pts_proj_coord);
  free(req_pts_dist2);

  if (log_timer) {
    log_trace("dt_extract_all = %12.5e / dt_search_all = %12.5e \n", dt_extract_all, dt_search_all);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(elt_extents[ipart]);
    free(elt_g_num  [ipart]);
    free(elt_type   [ipart]);
  }
  free(elt_extents);
  free(elt_g_num);
  free(elt_type);
  free(pn_elt);

  PDM_MPI_Barrier (ml->comm);
  PDM_timer_hang_on(ml->timer);
  ml->times_elapsed[END] = PDM_timer_elapsed (ml->timer);
  ml->times_cpu    [END] = PDM_timer_cpu     (ml->timer);
  ml->times_cpu_u  [END] = PDM_timer_cpu_user(ml->timer);
  ml->times_cpu_s  [END] = PDM_timer_cpu_sys (ml->timer);
  PDM_timer_resume(ml->timer);

  if (log_timer) {
    end_timer_and_log_from_dt("Step 1 : ", ml->comm, dt_step1_all);
    end_timer_and_log_from_dt("Step 2 : ", ml->comm, dt_step2_all);
    end_timer_and_log_from_dt("Step 3 : ", ml->comm, dt_step3_all);
    end_timer_and_log_from_dt("Step 4 : ", ml->comm, dt_step4_all);
    end_timer_and_log_from_dt("Step 5 : ", ml->comm, dt_step5_all);
    end_timer_and_log_from_dt("Step 6 : ", ml->comm, dt_step6_all);
  }

  // exit(1);
}

/**
 *
 * \brief Get the number of cells
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 * \param [in]  i_part   Index of partition of the mesh
 *
 * \return Number of cells
 */

int
PDM_mesh_location_n_cell_get
(
       PDM_mesh_location_t *ml,
 const int                  i_part
)
{
  return ml->points_in_elements[0].n_elts[i_part];
}


PDM_part_mesh_nodal_t*
PDM_mesh_location_mesh_nodal_get
(
 PDM_mesh_location_t *ml
)
{
  return ml->mesh_nodal;
}


/**
 * Enable reverse results computation (To call PDM_mesh_location_points_in_elt_get)
 */

void
PDM_mesh_location_reverse_results_enable
(
PDM_mesh_location_t *ml
)
{
  printf("WARNING : PDM_mesh_location_reverse_results_enable is deprecated (always enabled)\n");
  ml->reverse_result = 1;
}



/**
 * \brief Get part_to_part object to exchange data between
 * the source mesh and a target point cloud (both in user frame)
 *
 * \param [in ] ml         Pointer to \ref PDM_mesh_location_t object
 * \param [in ] icloud     Point cloud ID
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_mesh_location_part_to_part_get
(
       PDM_mesh_location_t  *ml,
 const int                   icloud,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership
 )
{
  assert(icloud < ml->n_point_cloud);

  *ptp = ml->ptp[icloud];
  ml->ptp_ownership[icloud] = ownership;
}


#ifdef	__cplusplus
}
#endif
