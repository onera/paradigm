#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_geom_elem.h"
#include "pdm_part_coarse_mesh.h"
#include "pdm_part_coarse_mesh_priv.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"

#include "pdm_part.h"
#include "pdm_part_renum.h"
#include "pdm_order.h"
#include "pdm_mpi.h"

#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_ext_wrapper.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

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


/**
 * \def _PDM_part_MIN(a,b)
 * Computes the minimum of \a x and \a y.
 *
 */

#define _PDM_part_MIN(a,b) ((a) > (b) ? (b) : (a))

/**
 * \def _PDM_part_MAX(a,b)
 * Computes the maximum of \a x and \a y.
 *
 */

#define _PDM_part_MAX(a,b) ((a) < (b) ? (b) : (a))

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

static PDM_Handles_t *_cm  = NULL;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * Storage of face renumbering methods
 */

static PDM_Handles_t *_coarse_mesh_methods = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/



/**
 *
 * \brief Return an initialized coarse part object
 *
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   n_part             Number of partitions
 * \param [in]   n_total_part            Total number of partitions
 * \param [in]   n_face_group        Number of boundaries
 * \param [in]   have_cell_tag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_face_tag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtx_tag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cell_weight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_face_weight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_face_group    Presence des tableaux de groupes de faces
 */

static _coarse_mesh_t *
_coarse_mesh_create
(
 const PDM_MPI_Comm  comm,
 const char         *method,
 const char         *renum_cell_method,
 const char         *renum_face_method,
 const int           n_property_cell,
 const int          *renum_properties_cell,
 const int           n_property_face,
 const int          *renum_properties_face,
 const int           n_part,
 const int           n_total_part,
 const int           n_face_group,
 const int           have_cell_tag,
 const int           have_face_tag,
 const int           have_vtx_tag,
 const int           have_cell_weight,
 const int           have_face_weight,
 const int           have_face_group

 )
{
   _coarse_mesh_t *cm = (_coarse_mesh_t *) malloc(sizeof(_coarse_mesh_t));

   cm->n_part = n_part;
   cm->comm  = comm;

   int _method = PDM_coarse_mesh_method_idx_get(method);

   if (_method == -1) {
     PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown coarse mesh method\n", method);
   }

   cm->method = _method;

   /* Reordering */
   _method = PDM_part_renum_method_cell_idx_get(renum_cell_method);

   if (_method == -1) {
     PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering cell method\n", renum_cell_method);
   }

   cm->renum_cell_method = _method;

   _method = PDM_part_renum_method_face_idx_get(renum_face_method);

   if (_method == -1) {
     PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering face method\n", renum_face_method);
   }
   cm->renum_face_method = _method;

   cm->n_property_cell        = n_property_cell;
   cm->renum_properties_cell  = renum_properties_cell;
   cm->n_property_face        = n_property_face;
   cm->renum_properties_face  = renum_properties_face;


   cm->n_total_part = n_total_part;

   cm->n_face_group = n_face_group;

   cm->have_cell_tag    = have_cell_tag;
   cm->have_face_tag    = have_face_tag;
   cm->have_vtx_tag     = have_vtx_tag;
   cm->have_cell_weight = have_cell_weight;
   cm->have_face_weight = have_face_weight;
   cm->have_face_group  = have_face_group;

   cm->part_ini = malloc(sizeof(_part_t *) * n_part); //On déclare un tableau de partitions

   cm->part_res = malloc(sizeof(_coarse_part_t *) * n_part);

   cm->specific_data = NULL;
   cm->specific_func = NULL;

   for (int i = 0; i < n_part; i++) {
     cm->part_ini[i] = _part_create();

     cm->part_res[i] = _coarse_part_create();

     cm->part_res[i]->part->n_face_group = cm->n_face_group;

   }

   return cm;
}

/**
 *
 * \brief Perform the coarse mesh from the SCOTCH graph method
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_coarse_from_scotch
(
_coarse_mesh_t* cm,
const int       i_part,
int            *n_coarse_cell_computed,
int            *cell_cell_idx,
int            *cell_cell,
int            *cell_part
)
{
#ifdef PDM_HAVE_PTSCOTCH

      int check = 0;
      int n_part  = cm->part_res[i_part]->nCoarseCellWanted;

      PDM_SCOTCH_part (cm->part_ini[i_part]->n_cell,
                       cell_cell_idx,
                       cell_cell,
                       (int *) cm->part_ini[i_part]->cellWeight,
                       (int *) cm->part_ini[i_part]->faceWeight,
                       check,
                       n_part,
                       cell_part);

      (*n_coarse_cell_computed) = n_part;

      if (0 == 1) {
        PDM_printf("\nContent of cell_part\n");
        for(int i = 0; i < cm->part_ini[i_part]->n_cell ; i++) {
          PDM_printf(" %d ", cell_part[i]);
        }
        PDM_printf("\n");
      }

#else
      PDM_printf("PDM_part error : Scotch unavailable\n");
      exit(1);
#endif
}


/**
 *
 * \brief Perform the coarse mesh from the SCOTCH graph method
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_coarse_from_metis
(
_coarse_mesh_t* cm,
const int       i_part,
int            *n_coarse_cell_computed,
int            *cell_cell_idx,
int            *cell_cell,
int            *cell_part
)
{
#ifdef PDM_HAVE_PARMETIS
      _part_t * part_ini       = cm->part_ini[i_part];
      _coarse_part_t *part_res = cm->part_res[i_part];

      int *cellWeight = (int *) part_ini->cellWeight;
      int *faceWeight = (int *) part_ini->faceWeight;

      int n_part  = part_res->nCoarseCellWanted;
      //            Define Metis properties

      //          int flag_weights = 0; //0 = False -> weights are unused

      int flag_weights = 1; //0 = False -> weights are unused

      int ncon = 1; //The number of balancing constraints

      int *vwgt = cellWeight; //Weights of the vertices of the graph (NULL if unused)

      int *adjwgt = faceWeight; //Weights of the edges of the graph (NULL if unused)

      if (flag_weights != 0) {
        double *tpwgts = (double *) malloc(ncon * n_part * sizeof(double));
        for (int i = 0; i < ncon * n_part; i++){
          tpwgts[i] = (double) (1./n_part);
        }
      }

      double *tpwgts = NULL;

      if (flag_weights != 0) {
        double *ubvec = (double *) malloc(ncon * sizeof(double));
        for (int i = 0; i < ncon; i++) {
          ubvec[i] = 1.05;
        }
      }

      double *ubvec = NULL;

      //TO ADD: USE OF ADJWGT IN AN IF STATEMENT

      //This value is solely a memory space to be filled by METIS

      int edgecut;

      if (n_part < 8) {

        PDM_printf("\n \t\t\t\t PDM_METIS_PartGraphRecursive\n");
        PDM_METIS_PartGraphRecursive (&(part_ini->n_cell),
                                      &ncon,
                                      cell_cell_idx,
                                      cell_cell,
                                      vwgt,
                                      adjwgt,
                                      &n_part,
                                      tpwgts,
                                      ubvec,
                                      &edgecut,
                                      cell_part);
      }

      else {

        PDM_printf("\n \t\t\t\tPDM_METIS_PartGraphKway \n");
        PDM_METIS_PartGraphKway (&(part_ini->n_cell),
                                 &ncon,
                                 cell_cell_idx,
                                 cell_cell,
                                 vwgt,
                                 adjwgt,
                                 &n_part,
                                 tpwgts,
                                 ubvec,
                                 &edgecut,
                                 cell_part);
      }

      if (1 == 0) {
        PDM_printf("\n Contenu de cell_part : \n");
        for (int i = 0; i < part_ini->n_cell; i++) {
          PDM_printf(" %d ", cell_part[i]);
        }
        PDM_printf("\n");
      }

      if (flag_weights != 0) {
          free(ubvec);
          free(tpwgts);
          free(adjwgt);
      }

      (*n_coarse_cell_computed) = n_part;

#else
      PDM_printf("PDM_part error : METIS unavailable\n");
      exit(1);

#endif
}

/**
 *
 * \brief Return coarse mesh object from its identifier
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

static _coarse_mesh_t *
_get_from_id
(
 int  cmId
)
{
  _coarse_mesh_t *cm = (_coarse_mesh_t *) PDM_Handles_get (_cm, cmId);

  if (cm == NULL) {
    PDM_printf("PPART error : Bad cm identifier\n");
    exit(1);
  }

  return cm;
}


/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */

static void
_quickSort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}


/**
 *
 * \brief Builds dual graph face cell connectivity
 *
 * \param [inout] part_ini                 Part object - fine mesh partition
 *
 * \param [inout] cell_cell_idxCompressed    Array of indexes of the dual graph
 * \param [inout] cell_cellCompressed       Dual graph
 *
 */

static void
_dual_graph_from_face_cell
(
  _part_t        *part_ini,
  int           **cell_cell_idxCompressed,
  int           **cell_cellCompressed
)
{
  //cell_cellN: array of counters of the numbers of connectivities
  //cell_cell: dual graph to be built
  //cell_cell_idx: array of indexes of the dual graph (same as cell_face_idx)

  int *cell_cellN = (int *) malloc(part_ini->n_cell * sizeof(int));
  for (int i = 0; i < part_ini->n_cell; i++) {
    cell_cellN[i] = 0;
  }

  int *cell_cell = (int *) malloc(part_ini->cell_face_idx[part_ini->n_cell] * sizeof(int));
  for (int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
    cell_cell[i] = -1;
  }

  int *cell_cell_idx = (int *) malloc((part_ini->n_cell + 1) * sizeof(int));
  for(int i = 0; i < part_ini->n_cell + 1; i++) {
    cell_cell_idx[i] = part_ini->cell_face_idx[i];
  }

  for (int i = 0; i < part_ini->n_face; i++) {
    int iCell1 = PDM_ABS (part_ini->face_cell[2*i    ]);
    int iCell2 = PDM_ABS (part_ini->face_cell[2*i + 1]);
    //Only the non-boundary faces are stored
    if (iCell2 > 0) {
      int idx1 = cell_cell_idx[iCell1-1] + cell_cellN[iCell1-1];
      cell_cell[idx1] = iCell2;
      cell_cellN[iCell1-1] += 1;

      int idx2 = cell_cell_idx[iCell2-1] + cell_cellN[iCell2-1];
      cell_cell[idx2] = iCell1;
      cell_cellN[iCell2-1] += 1;
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cell_cellN after looping over cell_face: ");
    for(int i = 0; i < part_ini->n_cell; i++) {
      PDM_printf(" %d ", cell_cellN[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cell_cell after looping over cell_face: ");
    for(int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
      PDM_printf(" %d ", cell_cell[i]);
    }
    PDM_printf("\n");
  }

  //cell_cell_idx is rebuilt
  *cell_cell_idxCompressed = malloc((part_ini->n_cell + 1) * sizeof(int));

  (*cell_cell_idxCompressed)[0] = 0;
  for(int i = 0; i < part_ini->n_cell; i++) {
    (*cell_cell_idxCompressed)[i + 1] = (*cell_cell_idxCompressed)[i] + cell_cellN[i];
  }

  //We compress the dual graph since cell_cell_idx was built from cell_face_idx
  //We have then n_face elements in cell_cell whereas it needs to be composed of n_cell elements

  //    PDM_printf("(*cell_cell_idxCompressed)[part_ini->n_cell] : %d \n", (*cell_cell_idxCompressed)[part_ini->n_cell]);
  *cell_cellCompressed = malloc((*cell_cell_idxCompressed)[part_ini->n_cell] * sizeof(int));

  int cpt_cell_cellCompressed = 0;
  for(int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
    //        PDM_printf("I am testing a value for the %d time! \n", i);

    //We have an information to store when a neighboring cell exists
    if(cell_cell[i] > -1){
      //            PDM_printf("I am storing a value for the %d time! \n", i);
      //We add a -1 to have the graph vertices numbered from 0 to n (C numbering)
      (*cell_cellCompressed)[cpt_cell_cellCompressed++] = cell_cell[i] - 1;
      //            PDM_printf("Valeur stockee : %d \n ", (*cell_cellCompressed)[cpt_cell_cellCompressed - 1]);
    }
  }

  if( 0 == 1) {
    PDM_printf("Content of cell_cellCompressed after compression and renumbering: ");
    for(int i = 0; i < (*cell_cell_idxCompressed)[part_ini->n_cell]; i++) {
      PDM_printf(" %d ", (*cell_cellCompressed)[i]);
    }
    PDM_printf("\n");
  }

  /* Free temporary arrays*/

  free(cell_cellN);
  free(cell_cell);
  free(cell_cell_idx);

  //Remove duplicate cells of the dual graph
  //We use the following scheme:
  //We loop over the indexes for the whole array to subdivide it into subarrays
  //We sort locally each subarray (determined thanks to cell_cell_idxCompressed)
  //We loop over each subarray
  //We store the first value of each subarray anyway
  //We store each non-duplicated value and increment the writing index
  //We update the index array at each iteration

  int idx_write = 0;
  int tabIdxTemp = 0;

  for (int i = 0; i < part_ini->n_cell; i++) {
    _quickSort_int((*cell_cellCompressed), tabIdxTemp, (*cell_cell_idxCompressed)[i + 1] - 1);

    int last_value = -1;

    for (int j = tabIdxTemp; j < (*cell_cell_idxCompressed)[i + 1]; j++) {
      //We need to have a local index (between 0 and n_face)
      //If the value is different from the previous one (higher than is the same as different since the array is sorted)

      if(last_value != (*cell_cellCompressed)[j]) {
        (*cell_cellCompressed)[idx_write++] = (*cell_cellCompressed)[j];
        last_value = (*cell_cellCompressed)[j];
      }
    }

    if (0 == 1) {
      PDM_printf("\n Contenu de cell_cellCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < (*cell_cell_idxCompressed)[part_ini->n_cell]; i1++) {
        PDM_printf(" %d ", (*cell_cellCompressed)[i1]);
      }
      PDM_printf("\n");
    }

    tabIdxTemp = (*cell_cell_idxCompressed)[i + 1];
    (*cell_cell_idxCompressed)[i + 1] = idx_write;

    if (0 == 1) {
      PDM_printf("\n Contenu de cell_cell_idxCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < part_ini->n_cell + 1; i1++) {
        PDM_printf(" %d ", (*cell_cell_idxCompressed)[i1]);
      }
      PDM_printf("\n");
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cell_cell_idxCompressed after compression: ");
    for(int i1 = 0; i1 < part_ini->n_cell + 1; i1++) {
      PDM_printf(" %d ", (*cell_cell_idxCompressed)[i1]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cell_cellCompressed after compression: ");
    for(int i1 = 0; i1 < (*cell_cell_idxCompressed)[part_ini->n_cell]; i1++) {
      PDM_printf(" %d ", (*cell_cellCompressed)[i1]);
    }
    PDM_printf("\n");
  }

  //We reallocate the memory in case of duplicated values removed
  //The new array size is idx_write (stored in (*cell_cell_idxCompressed)[part_ini->n_cell])
  *cell_cellCompressed = realloc(*cell_cellCompressed,
                                (*cell_cell_idxCompressed)[part_ini->n_cell] * sizeof(int));

}

/**
 *
 * \brief Splits the graph
 *
 * \param [in]  method       Method to be used: choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]  n_part        Number of partitions
 * \param [in]  part_ini     Part object fine mesh
 *
 * \param [in]  cell_cell                  Dual graph (size : cell_cell_idx[n_cell])
 * \param [in]  cell_cell_idx               Array of indexes of the dual graph (size : n_cell + 1)
 * \param [in]  cellWeight         Cell weight (size = n_cell)
 * \param [in]  faceWeight         Face weight (size = n_face)
 *
 * \param [inout] cell_part  Cell partitioning (size : n_cell)
 *
 */

static void
_split
(
_coarse_mesh_t *cm,
const int       i_part,
int            *n_coarse_cell_computed,
int            *cell_cell_idx,
int            *cell_cell,
int            **cell_part)
{

  _part_t * part_ini       = cm->part_ini[i_part];

  int method = cm->method;

  *cell_part = (int *) malloc(part_ini->n_cell * sizeof(int));

  for (int i = 0; i < part_ini->n_cell; i++){
    (*cell_part)[i] = 0;
  }

  const _coarse_mesh_method_t *method_ptr = (const _coarse_mesh_method_t *)
                                      PDM_Handles_get (_coarse_mesh_methods, method);

  PDM_coarse_mesh_fct_t fct = method_ptr->fct;

  if (fct != NULL) {
    (fct) (cm,
           i_part,
           n_coarse_cell_computed,
           cell_cell_idx,
           cell_cell,
           *cell_part);
  }

}

/**
 *
 * \brief Obtains the cells per processor from the partition of each cell
 *
 * \param [in]  nCoarseCell  Number of partitions ( = number of coarse cells)
 * \param [in]  n_cell        Number of cells before refining
 * \param [in]  cell_part     Cell partitioning (size : n_cell)
 *
 * \param [inout] partCellIdx  Array of indexes of the partitioning array (size : nCoarseCell + 1)
 * \param [inout] partCell     Partitioning array (size : partCellIdx[nCoarseCell] = n_cell)
 *
 */

static void
_partCell_from_cell_part
(
 int            nCoarseCell,
 int            n_cell,
 int           *cell_part,
 int          **partCellIdx,
 int          **partCell
)
{
  //Allocation of an array to count the number of cells per partition
  int * cptCellsPerPartitions = (int *) malloc(nCoarseCell * sizeof(int));
  for (int i = 0; i < nCoarseCell; i++){
    cptCellsPerPartitions[i] = 0;
  }

  for (int i = 0; i < n_cell; i++){
    int color = cell_part[i]; //A color is a number of partition (output of Metis or Scotch)
    cptCellsPerPartitions[color]++;
  }

  if(0 == 1) {
    PDM_printf("\n Contenu de cptCellsPerPartitions : \n");
    for(int i = 0; i < nCoarseCell; i++) {
      PDM_printf(" %d ", cptCellsPerPartitions[i]);
    }
    PDM_printf("\n");
  }

  //Allocation of an array for counter indexes
  *partCellIdx = (int *) malloc((nCoarseCell + 1) * sizeof(int));
  (*partCellIdx)[0] = 0;
  for (int i = 0; i < nCoarseCell; i++){
    (*partCellIdx)[i + 1] = (*partCellIdx)[i] + cptCellsPerPartitions[i];
  }

  if (0 == 1) {
    PDM_printf("\n Contenu de partCellIdx : \n");
    for(int i = 0; i < nCoarseCell + 1; i++) {
      PDM_printf(" %d ", (*partCellIdx)[i]);
    }
    PDM_printf("\n");
  }

  *partCell = (int *) malloc((*partCellIdx)[nCoarseCell] * sizeof(int));

  //cptCellsPerPartitions is reused for building partCell
  for (int i = 0; i < nCoarseCell; i++){
    cptCellsPerPartitions[i] = 0;
  }

  //We store each cell in partCell by means of (*partCellIdx)
  for (int i = 0; i < n_cell; i++){
    int color = cell_part[i]; //A color is a number of partition (output of Metis or Scotch)
    int idx = (*partCellIdx)[color] + cptCellsPerPartitions[color];
    (*partCell)[idx] = i;
    cptCellsPerPartitions[color]++;
  }

  if (0 == 1) {
    PDM_printf("\nContenu de partCell \n");
    for (int i = 0; i < nCoarseCell; i++){
      PDM_printf("Valeur de i + 1 : %d \n", i + 1);
      for (int j = (*partCellIdx)[i]; j < (*partCellIdx)[i + 1]; j++){
        PDM_printf("%d " ,(*partCell)[j]);
      }
      PDM_printf("\n");
    }
    PDM_printf("\n");
  }

  //Free
  free(cptCellsPerPartitions);

}

/**
 *
 * \brief Checks the neighboring cell of each studied cell by filling an array of tags
 *        A tag is set to -1 by default and set to numberGlobalPartition if the studied cell has valid neighbors
 *
 * \param [in]  cellNumber                Number of the cell studied
 * \param [in]  cell_part                  Cell partitioning (size : n_cell)
 * \param [in]  cell_cell                  Dual graph (size : cell_cell_idx[n_cell])
 * \param [in]  cell_cell_idx               Array of indexes of the dual graph (size : n_cell + 1)
 * \param [in]  numberGlobalPartition     Number of the coarse cell to be applied
 *
 * \param [inout] cellCoarseCell          Cell partitioning (size : n_cell) (partitions are equal to coarse cells for a good partitioning)
 * \param [inout] cptCellConnectedLocal   Number of cells that have been tagged
 */

static void
_fill_Neighboring
(
 int            cellNumber,
 int           *cell_part,
 int           *cell_cell,
 int           *cell_cell_idx,
 int          **cellCoarseCell,
 int            numberGlobalPartition,
 int           *cptCellConnectedLocal
)
{
  /*
   * If the studied cell is not tagged, it is tagged
   */

  if ((*cellCoarseCell)[cellNumber] == -1) {
    (*cellCoarseCell)[cellNumber] = numberGlobalPartition;
    (*cptCellConnectedLocal)++;
  }

  /*
   * If the cell is tagged, I break
   */

  else {
    return;
  }

  /*
   * Loop over the partition cells (k is a neighbor cell number)
   */

  for (int k = cell_cell_idx[cellNumber]; k < cell_cell_idx[cellNumber + 1]; k++) {

    /*
     * if the neighboring cell is part of the same partition
     */

    if (cell_part[cellNumber] == cell_part[cell_cell[k]]) {
      _fill_Neighboring (cell_cell[k],
                         cell_part,
                         cell_cell,
                         cell_cell_idx,
                         cellCoarseCell,
                         numberGlobalPartition,
                         cptCellConnectedLocal);
    }
  }
}


/**
 *
 * \brief Checks that the partitioning is fully connected
 * Otherwise, the cells that are not connected are removed for their initial partition to be part of a new one
 * partCell and partCellIdx are updated
 *
 * \param [in]  nCoarseCellChecked  Number of partitions checked ( >= number of coarse cells wanted by the user)
 * \param [in]  n_cell               Number of cells before refining
 * \param [in]  cell_part            Cell partitioning (size : n_cell) *
 * \param [in]  cell_cell            Dual graph (size : cell_cell_idx[n_cell])
 * \param [in]  cell_cell_idx         Array of indexes of the dual graph (size : n_cell + 1)
 * \param [in]  partCell            Partitioning array (size : partCellIdx[nCoarseCell] = n_cell)
 * \param [in]  partCellIdx         Array of indexes of the partitions (size : nCoarseCellWanted + 1)
 *
 * \param [inout] coarse_cell_cell_idx  Array of indexes of the connected partitions (size : nCoarseCellChecked + 1)
 * \param [inout] coarse_cell_cell     Partitioning array (size : coarse_cell_cell_idx[nCoarseCellChecked])
 * \param [inout] cellCoarseCell     Cell partitioning with coarse cells (size : nCoarseCellChecked)
 */

static void
_adapt_Connectedness
(
 int           *nCoarseCellChecked,
 int            n_cell,
 int           *cell_part,
 int          **cellCoarseCell,
 int           *cell_cell,
 int           *cell_cell_idx,
 int           *partCell,
 int           *partCellIdx,
 int          **coarse_cell_cell,
 int          **coarse_cell_cell_idx
)
{
  int numberGlobalPartition = 1;

  *cellCoarseCell = (int *) malloc(n_cell * sizeof(int));

  for (int i = 0; i < n_cell; i++) {
    (*cellCoarseCell)[i] = -1;
  }

  /*
   *  We store the initial number of coarse cells wanted by the user
   */

  int nCoarseCellWanted = (*nCoarseCellChecked);

  if (0 == 1) {
    PDM_printf("Valeur finale de (*nCoarseCellChecked) : %d %d\n", (*nCoarseCellChecked), n_cell);

    PDM_printf("partCell : \n");
    for (int i = 0; i < nCoarseCellWanted; i++) {
      for (int j = partCellIdx[i]; j < partCellIdx[i+1]; j++) {
        PDM_printf(" %d ", partCell[j]);
      }
      PDM_printf("\n");
    }

    PDM_printf("cell_cell : \n");
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_cell_idx[i]; j < cell_cell_idx[i+1]; j++) {
        PDM_printf(" %d ", cell_cell[j]);
      }
      PDM_printf("\n");
    }

  }


  /*
   * cellNumber will be replaced by the first cell of the first partition at the beginning of the loop
   */

  int cellNumber = -1;

  /*
   * Loop over the partitions (i is a partition number)
   */

  for (int i = 0; i < nCoarseCellWanted; i++) {
    /*
     * We study the first cell of the partition
     */

    cellNumber = partCell[partCellIdx[i]];
    int cptCellConnectedLocal = 0;

    /*
     * We tag all the neighboring cells of the cell cellNumber of the partition
     */

    _fill_Neighboring(cellNumber, cell_part, cell_cell, cell_cell_idx, &(*cellCoarseCell), numberGlobalPartition, &cptCellConnectedLocal);

    numberGlobalPartition++;

    /*
     * If the size of array indexes is too low
     */

    int n_cellLocal = partCellIdx[i + 1] - partCellIdx[i];
    int n_cellLocalRemaining = n_cellLocal - cptCellConnectedLocal;

    /*
     * If the partition has not been fully looped over, we will have to create an extra coarse cell
     */

    if (cptCellConnectedLocal < n_cellLocal) {
      /*
       * We reinitialize cptCellConnectedLocal since we have a new coarse cell
       */
      cptCellConnectedLocal = 0;
    }

    /*
     *  As long as the partition has not been fully looped over, we call the recursive function
     */

    while (cptCellConnectedLocal < n_cellLocalRemaining) {
      for (int j = 1; j < n_cellLocal; j++) {
        cellNumber = partCell[partCellIdx[i] + j];
        if ((*cellCoarseCell)[cellNumber] == -1) {
          break;
        }
      }

      _fill_Neighboring(cellNumber, cell_part, cell_cell, cell_cell_idx, &(*cellCoarseCell), numberGlobalPartition, &cptCellConnectedLocal);

      numberGlobalPartition++;

      /*
       * If the size of array indexes is too low
       */

    }

  }

  (*nCoarseCellChecked) = numberGlobalPartition - 1;

  /*
   * Size of *coarse_cell_cell_idx may be dynamic
   */

  *coarse_cell_cell_idx = malloc(((*nCoarseCellChecked) + 1) * sizeof(int));
  for (int i = 0; i < (*nCoarseCellChecked) + 1; i++) {
    (*coarse_cell_cell_idx)[i] = 0;
  }

  for (int i = 0; i < n_cell; i++) {
    (*coarse_cell_cell_idx)[(*cellCoarseCell)[i]]++;
  }

  for (int i = 0; i < (*nCoarseCellChecked); i++) {
    (*coarse_cell_cell_idx)[i+1] += (*coarse_cell_cell_idx)[i];
  }

  if (0 == 1) {
    PDM_printf("Valeur finale de (*nCoarseCellChecked) : %d %d\n", (*nCoarseCellChecked), n_cell);

    PDM_printf("Affichage de *coarse_cell_cell_idx");
    for (int i = 0; i < (*nCoarseCellChecked) + 1; i++) {
      PDM_printf(" %d ", (*coarse_cell_cell_idx)[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cellCoarseCell: ");
    for (int i = 0; i < n_cell; i++) {
      PDM_printf(" %d ", (*cellCoarseCell)[i]);
    }
    PDM_printf("\n");
  }

  /*
   * Creation of coarse_cell_cell from cellCoarseCell and cellCoarseCellIdx
   */

  *coarse_cell_cell = (int *) malloc(n_cell * sizeof(int));

  int * cptCellsPerPartitions = (int *) malloc((*nCoarseCellChecked) * sizeof(int));
  for (int i = 0; i < (*nCoarseCellChecked); i++){
    cptCellsPerPartitions[i] = 0;
  }

  /*
   * We store each cell in partCell by means of (*partCellIdx)
   */

  for (int i = 0; i < n_cell; i++){
    int color = (*cellCoarseCell)[i] - 1; //A color is a number of partition (output of Metis or Scotch)
    int idx = (*coarse_cell_cell_idx)[color] + cptCellsPerPartitions[color];
    (*coarse_cell_cell)[idx] = i + 1;
    cptCellsPerPartitions[color]++;
  }

  free(cptCellsPerPartitions);
}

/**
 *
 * \brief Builds the array faceCoarseCell with all the inner faces removed
 *
 * \param [in]  n_faceChecked                Number of faces after refining ( <= n_face)
 * \param [in]  face_cell                    Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [in]  cellCoarseCell              Cell partitioning with coarse cells (size : nCoarseCellChecked)
 *
 * \param [inout] faceCoarseCell            Face to coarse cell connectivity  (size = 2 * n_faceChecked, numbering : 1 to n)
 * \param [inout] fineFaceToCoarseFace      Fine face - coarse face connectivity (size = n_face)
 * \param [inout] coarse_face_to_fine_face      Coarse face - fine face connectivity (size = n_faceChecked)
 *
 */

static void
_build_faceCoarseCell
(
 int           *n_faceChecked,
 int           *face_cell,
 int           *cellCoarseCell,
 int          **faceCoarseCell,
 int          **fineFaceToCoarseFace,
 int          **coarse_face_to_fine_face
)
{
  /*
   * n_faceChecked = n_face at the beginning of the method
   */

  int n_face = (*n_faceChecked);

  /*
   * Fine face - coarse face connectivity (size = n_face)
   */

  *fineFaceToCoarseFace = malloc(n_face * sizeof(int));

  int *faceCellTemp = (int *) malloc(2 * n_face * sizeof(int));

  for (int i = 0; i < n_face; i++) {
    (*fineFaceToCoarseFace)[i] = -1;
  }

  for (int i = 0; i < 2 * n_face; i++) {
    faceCellTemp[i] = 0;
  }

  /*
   * Loop over face_cell. i = number of face, face_cell[i] = cell number
   * We fill faceCellTemp with the coarse cells associated
   */

  for (int i = 0; i < 2 * n_face; i++) {

    /*
     * If we have a boarding face, we have a -1 -> nothing to do
     * If we have a "real" neighboring cell, we store its coarse cell
     */

    if (face_cell[i] != 0) {
      faceCellTemp[i] = cellCoarseCell[PDM_ABS (face_cell[i]) - 1];
    }
  }

  if (0 == 1) {
    PDM_printf("Content of faceCellTemp: |");
    for (int i = 0; i < 2 * n_face; i++) {
      PDM_printf(" %d ", faceCellTemp[i]);
      if (i % 2 == 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }

  *faceCoarseCell = (int *) malloc(2 * n_face * sizeof(int));

  /*
   * Loop over faceCellTemp which is to be compressed. i = face number
   */

  int idx = 0;
  for (int i = 0; i < n_face; i++) {
    int iCell1 = PDM_ABS (faceCellTemp[2 * i    ]);
    int iCell2 = PDM_ABS (faceCellTemp[2 * i + 1]);

    /*
     * If a face is surrounded by the same coarse cell, it is not stored
     */
    if (iCell1 != iCell2) {
      (*faceCoarseCell)[2 * idx]     = iCell1;
      (*faceCoarseCell)[2 * idx + 1] = iCell2;

      (*fineFaceToCoarseFace)[i] = idx + 1;
      idx++;
    }
  }

  (*n_faceChecked) = idx;
  /*
   * realloc of the correct size
   */

  *faceCoarseCell = realloc((*faceCoarseCell), 2 * (*n_faceChecked) * sizeof(int));

  /*
   * Fine face - coarse face connectivity (size = n_face)
   */

  *coarse_face_to_fine_face = malloc(n_face * sizeof(int));

  int idx_coarse_face_to_fine_face = 0;

  /*
   *  Loop over fineFaceToCoarseFace
   */

  for (int i = 0; i < n_face; i++) {

    /*
     * If the fine face has not been removed, I store it
     */

    if((*fineFaceToCoarseFace)[i] != -1) {
      (*coarse_face_to_fine_face)[idx_coarse_face_to_fine_face++] = i + 1;
    }
  }

  /*
   * At the end of the loop, idx_coarse_face_to_fine_face must be equal to n_faceChecked
   */

  assert(idx_coarse_face_to_fine_face == (*n_faceChecked));

  *coarse_face_to_fine_face = realloc((*coarse_face_to_fine_face), idx_coarse_face_to_fine_face * sizeof(int));

  if(0 == 1) {
    PDM_printf("Valeur finale de (*n_faceChecked) : %d \n", (*n_faceChecked));

    PDM_printf("Final content of faceCoarseCell: |");
    for (int i = 0; i < 2 * (*n_faceChecked); i++) {
      PDM_printf(" %d ", (*faceCoarseCell)[i]);
      if (i % 2 == 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");

    PDM_printf("Final content of fineFaceToCoarseFace: \n");
    for (int i = 0; i < n_face; i++) {
      PDM_printf(" %d ", (*fineFaceToCoarseFace)[i]);
    }
    PDM_printf("\n");


    PDM_printf("Affichage final de (*coarse_face_to_fine_face) \n");
    for (int i = 0; i < (*n_faceChecked); i++) {
      PDM_printf(" %d ", (*coarse_face_to_fine_face)[i]);
    }
    PDM_printf("\n");
  }

  free(faceCellTemp);
}

/**
 *
 * \brief Obtains the faces per coarse cell from the coarse cell of each face
 *
 * \param [in]  nCoarseCellChecked    Number of coarse cells after the connectedness check
 * \param [in]  n_faceChecked          Number of faces obtained after the creation of faceCoarseCell
 * \param [in]  faceCoarseCell        Face to coarse cell connectivity  (size = 2 * n_faceChecked, numbering : 1 to n)

 * \param [inout] coarsecell_face_idx   Array of indexes of the coarse cell to face connectivity (size = nCoarseCellChecked + 1, numbering : 1 to n)
 * \param [inout] coarsecell_face      Coarse cell to face connectivity  (size = coarsecell_face_idx[nCoarseCellChecked], numbering : 1 to n)
 *
 */

static void
_coarsecell_face_from_faceCoarseCell
(
 int            nCoarseCellChecked,
 int            n_faceChecked,
 int           *faceCoarseCell,
 int          **coarsecell_face_idx,
 int          **coarsecell_face
)
{
  /*
   *  Allocation of an array to count the number of faces per coarse cell
   */

  int * cptFacesPerCoarseCell = (int *) malloc(nCoarseCellChecked * sizeof(int));
  for (int i = 0; i < nCoarseCellChecked; i++) {
    cptFacesPerCoarseCell[i] = 0;
  }

  /*
   * Loop over faceCoarseCell. i = number of face
   */

  for (int i = 0; i < n_faceChecked; i++) {
    int coarseCell1 = faceCoarseCell[2 * i    ];
    int coarseCell2 = faceCoarseCell[2 * i + 1];
    cptFacesPerCoarseCell[coarseCell1 - 1]++;

    /*
     * If coarseCell2 != -1, it is not a boarder cell
     * A non-boarder cell touches two coarse cells
     */

    if(coarseCell2 != 0) {
      cptFacesPerCoarseCell[coarseCell2 - 1]++;
    }
  }

  if(0 == 1) {
    PDM_printf("\n Contenu de cptFacesPerCoarseCell : \n");
    for(int i = 0; i < nCoarseCellChecked; i++) {
      PDM_printf(" %d ", cptFacesPerCoarseCell[i]);
    }
    PDM_printf("\n");
  }

  /*
   * Allocation of an array for counter indexes
   */

  *coarsecell_face_idx = (int *)malloc((nCoarseCellChecked + 1) * sizeof(int));
  (*coarsecell_face_idx)[0] = 0;
  for (int i = 0; i < nCoarseCellChecked; i++) {
    (*coarsecell_face_idx)[i + 1] = (*coarsecell_face_idx)[i] + cptFacesPerCoarseCell[i];
  }

  *coarsecell_face = (int *) malloc((*coarsecell_face_idx)[nCoarseCellChecked] * sizeof(int));

  /*
   *  cptFacesPerCoarseCell is reused for building coarsecell_face
   */

  for (int i = 0; i < nCoarseCellChecked; i++){
    cptFacesPerCoarseCell[i] = 0;
  }

  /*
   * We store each face in coarsecell_face by means of (*coarsecell_face_idx)
   * Loop over faceCoarseCell. i = number of face
   */

  for (int i = 0; i < n_faceChecked; i++) {
    int coarseCell1 = faceCoarseCell[2 * i];
    int coarseCell2 = faceCoarseCell[2 * i + 1];

    int idx1 = (*coarsecell_face_idx)[coarseCell1 - 1] + cptFacesPerCoarseCell[coarseCell1 - 1];
    int idx2 = -1;

    /*
     * If the face is not on the boarder, we store it
     */

    if (coarseCell2 != 0) {
      idx2 = (*coarsecell_face_idx)[coarseCell2 - 1] + cptFacesPerCoarseCell[coarseCell2 - 1];
    }

    (*coarsecell_face)[idx1] = i + 1;
    cptFacesPerCoarseCell[coarseCell1 - 1]++;

    /*
     * If idx2 is higher than -1, it means that the face is not on the boarder
     */

    if (idx2 > -1) {
      (*coarsecell_face)[idx2] = i + 1;
      cptFacesPerCoarseCell[coarseCell2 - 1]++;
    }
  }

  if(0 == 1) {
    PDM_printf("Contenu de (*coarsecell_face) \n");
    for (int i = 0; i < (*coarsecell_face_idx)[nCoarseCellChecked]; i++) {
      PDM_printf(" %d ", (*coarsecell_face)[i]);
      if (i % (*coarsecell_face_idx)[1] == (*coarsecell_face_idx)[1] - 1) {
        PDM_printf("|");
      }
    }

    PDM_printf("\n Contenu de (*coarsecell_face_idx) : \n");
    for (int i = 0; i < nCoarseCellChecked + 1; i++) {
      PDM_printf(" %d ", (*coarsecell_face_idx)[i]);
    }
    PDM_printf("\n");
  }

  free(cptFacesPerCoarseCell);
}

/**
 *
 * \brief Builds the array face_vtx with all the inner vertices removed
 *
 * \param [in] n_face                 Number of faces before refining
 * \param [in] n_faceChecked          Number of faces after refining ( <= n_face)
 * \param [in] n_vtx                  Number of vertices before refining
 * \param [in] fineFaceToCoarseFace  Fine face - coarse face connectivity (size = n_face)
 *
 * \param [inout] face_vtx_idx         Face vertex connectivity index (final size = n_faceChecked + 1)
 * \param [inout] face_vtx            Face vertex connectivity (final size = face_vtx_idx[n_faceChecked])
 * \param [inout] n_vtxChecked        Number of vertices before refining becoming the number of vertices after refining
 * \param [inout] fineVtxToCoarseVtx Fine vertex - coarse vertex connectivity (size = n_vtx)
 * \param [inout] coarse_vtx_to_fine_vtx Coarse vertex - fine vertex connectivity (size = n_vtxChecked)
 *
 */

static void
_build_face_vtx
(
 int            n_face,
 int            n_faceChecked,
 int            n_vtx,
 int           *fineFaceToCoarseFace,
 int          **face_vtx_idx,
 int          **face_vtx,
 int           *n_vtxChecked,
 int          **fineVtxToCoarseVtx,
 int          **coarse_vtx_to_fine_vtx
)
{
  int idx_write_face_vtx = 0;

  (*face_vtx_idx)[0] = 0;
  int idx_write_face_vtx_idx = 1;

  /*
   * Loop over the old face_vtx_idx, i = face number
   */

  for (int i = 0; i < n_face; i++) {
    //Loop over the old face_vtx, j = vertex number
    if (fineFaceToCoarseFace[i] != - 1) {

      for (int j = (*face_vtx_idx)[i]; j < (*face_vtx_idx)[i + 1]; j++) {
        //If the face studied has been removed, I skip it
        int vtx = (*face_vtx)[j];
        (*face_vtx)[idx_write_face_vtx++] = vtx;
      }

      (*face_vtx_idx)[idx_write_face_vtx_idx] = (*face_vtx_idx)[i + 1] - (*face_vtx_idx)[i] + (*face_vtx_idx)[idx_write_face_vtx_idx - 1];
      idx_write_face_vtx_idx++;
    }
  }

  *face_vtx_idx = realloc((*face_vtx_idx), (n_faceChecked + 1) * sizeof(int));
  *face_vtx = realloc((*face_vtx), (*face_vtx_idx)[n_faceChecked] * sizeof(int));

  if (0 == 1) {
    PDM_printf("Valeur de (*face_vtx_idx)[n_faceChecked] : %d \n", (*face_vtx_idx)[n_faceChecked]);

    for (int i = 0; i < n_faceChecked; i++) {
      for (int j = (*face_vtx_idx)[i]; j < (*face_vtx_idx)[i + 1]; j++) {
        //If the face studied has been removed, I skip it
        int vtx = (*face_vtx)[j];
        // A supprimer
        for (int j1 = (*face_vtx_idx)[i]; j1 < (*face_vtx_idx)[i + 1]; j1++) {
          //If the face studied has been removed, I skip it
          int vtx1 = (*face_vtx)[j1];
          if (j != j1 && vtx == vtx1) {
            PDM_printf("Error multiple vertex in a face\n");
            abort();
          }
        }
      }
    }
    PDM_printf("\n");
  }

  /*
   * Creation of a correspondence table coarse vertex to fine vertex
   */

  *coarse_vtx_to_fine_vtx = malloc((*face_vtx_idx)[n_faceChecked] * sizeof(int));

  int idx_write_coarse_vtx_to_fine_vtx = 0;

  /*
   * It is a copy of face_vtx at first
   * Then, it is sorted
   * All the double vertices from the sorted array are removed
   * We have our correspondence table
   */

  for (int i = 0; i < (*face_vtx_idx)[n_faceChecked]; i++) {
    (*coarse_vtx_to_fine_vtx)[i] = (*face_vtx)[i];
  }

  _quickSort_int((*coarse_vtx_to_fine_vtx), 0, (*face_vtx_idx)[n_faceChecked] - 1);

  int last_value = -1;

  /*
   * Loop over (*coarse_vtx_to_fine_vtx)
   * Each vertex is stored only once
   */

  for (int i = 0; i < (*face_vtx_idx)[n_faceChecked]; i++) {
    if (last_value != (*coarse_vtx_to_fine_vtx)[i]) {
      (*coarse_vtx_to_fine_vtx)[idx_write_coarse_vtx_to_fine_vtx++] = (*coarse_vtx_to_fine_vtx)[i];
      last_value = (*coarse_vtx_to_fine_vtx)[i];
    }
  }

  (*n_vtxChecked) = idx_write_coarse_vtx_to_fine_vtx;

  (*coarse_vtx_to_fine_vtx) = realloc((*coarse_vtx_to_fine_vtx), (*n_vtxChecked) * sizeof(int));

  if (0 == 1) {
    PDM_printf("\nFinal content of coarse_vtx_to_fine_vtx: ");
    for (int i = 0; i < (*n_vtxChecked); i++) {
      PDM_printf(" %d ", (*coarse_vtx_to_fine_vtx)[i]);
    }
    PDM_printf("\n");
  }

  /*
   * Creation of a correspondence table fine vertex to coarse vertex
   */

  *fineVtxToCoarseVtx = malloc(n_vtx * sizeof(int));

  for (int i = 0; i < n_vtx; i++) {
    (*fineVtxToCoarseVtx)[i] = -1;
  }

  /*
   * Loop over (*coarse_vtx_to_fine_vtx)
   */

  for (int i = 0; i < (*n_vtxChecked); i++) {
    int fineVtx = (*coarse_vtx_to_fine_vtx)[i];
    //        PDM_printf("Valeur de fineVtx : %d \n", fineVtx);
    (*fineVtxToCoarseVtx)[fineVtx - 1] = i + 1;
  }

  if(0 == 1) {
    PDM_printf("Content of fineVtxToCoarseVtx: ");
    for (int i = 0; i < n_vtx; i++) {
      PDM_printf(" %d ", (*fineVtxToCoarseVtx)[i]);
    }
    PDM_printf("\n");
  }

  //Loop over face_vtx to re-number face_vtx thanks to (*coarse_vtx_to_fine_vtx)
  for (int i = 0; i < (*face_vtx_idx)[n_faceChecked]; i++) {
    (*face_vtx)[i] = (*fineVtxToCoarseVtx)[(*face_vtx)[i] - 1];
  }

  if (0 == 1) {
    PDM_printf("Valeur de (*n_vtxChecked) : %d \n", (*n_vtxChecked));
    PDM_printf("Valeur de idx_write_face_vtx : %d \n", idx_write_face_vtx);

    PDM_printf("Final content of face_vtx_idx: ");
    for(int i = 0; i < n_faceChecked + 1; i++) {
      PDM_printf(" %d ", (*face_vtx_idx)[i]);
    }
    PDM_printf("\n");

    PDM_printf("Final content of face_vtx: |");
    for (int i = 0; i < (*face_vtx_idx)[n_faceChecked]; i++) {
      PDM_printf(" %d ", (*face_vtx)[i]);
      if (i % (*face_vtx_idx)[1] == (*face_vtx_idx)[1] - 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }

}

/**
 *
 * \brief Builds the array vtx with all the coordinates of the inner vertices removed
 *
 * \param [in]  n_vtx                 Number of vertices before refining
 * \param [in]  n_vtxChecked          Number of vertices after refining
 * \param [in]  fineVtxToCoarseVtx   Fine vertex - coarse vertex connectivity (size = n_vtx)
 *
 * \param [inout] vtx                Vertex coordinates (size = n_vtxChecked)
 *
 */

static void
_build_vtx
(
 int            n_vtx,
 int            n_vtxChecked,
 int           *fineVtxToCoarseVtx,
 double       **vtx
)
{
  //If no vertex has been removed, nothing to do!
  if (n_vtx == n_vtxChecked) {
      return;
  }

  int idx_write = 0;

  //Loop over fineVtxToCoarseVtx, i = index of a vertex number (vertex number - 1)
  for (int i = 0; i < n_vtx; i++) {
    //We store each vertex that has not been removed
    if (fineVtxToCoarseVtx[i] != -1) {
      double coord1 = (*vtx)[3 * i    ];
      double coord2 = (*vtx)[3 * i + 1];
      double coord3 = (*vtx)[3 * i + 2];

      (*vtx)[idx_write++] = coord1;
      (*vtx)[idx_write++] = coord2;
      (*vtx)[idx_write++] = coord3;
    }
  }

  //Reallocation of vtx at the suitable size
  *vtx = realloc((*vtx), 3 * n_vtxChecked * sizeof(double));

  assert(3 * n_vtxChecked == idx_write);

  if(0 == 1) {
    PDM_printf("Contenu final de vtx\n");
    for (int i = 0; i < 3 * n_vtxChecked; i++) {
      PDM_printf(" %.1f ", (*vtx)[i]);
      if (i % 3 == 2) {
          PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Updates the array cell_tag in an array called coarsecell_tag
 *
 * \param [in] nCoarseCellChecked Number of partitions checked ( >= number of coarse cells wanted by the user)
 * \param [in] coarse_cell_cell_idx  Array of indexes of the connected partitions (size : nCoarseCellChecked + 1)
 * \param [in] coarse_cell_cell     Partitioning array (size : coarse_cell_cell_idx[nCoarseCellChecked])
 * \param [in] cell_tag            Cell tag (size = n_cell)
 *
 * \param [inout] coarsecell_tag   Tag coarse cell connectivity index (size = nCoarseCellChecked)
 *
 */

static void
_build_coarsecell_tag
(
 int            nCoarseCellChecked,
 int           *coarse_cell_cell_idx,
 int           *coarse_cell_cell,
 int           *cell_tag,
 int          **coarsecell_tag
)
{
  if (cell_tag == NULL) {
      return;
  }

  *coarsecell_tag = (int *) malloc(nCoarseCellChecked * sizeof(int));

  //Loop over coarse_cell_cell_idx, i = index of coarse cell
  for (int i = 0; i < nCoarseCellChecked; i++) {
    //This should be the tag of all the fine cells of the coarse cell j
    int tag = cell_tag[coarse_cell_cell[coarse_cell_cell_idx[i]] - 1];

    //Loop over coarse_cell_cell, j = coarse cell number
    for (int j = coarse_cell_cell_idx[i]; j < coarse_cell_cell_idx[i + 1]; j++) {
      //If any fine cell does not have the same tag as the previous one, the cell_tag array is incorrect
      if (cell_tag[coarse_cell_cell[j] - 1] != tag) {
        PDM_printf("Incorrect cell_tag array provided!\n");
        PDM_printf("Please check the fine cell %d\n", j + 1);
        PDM_printf("A default tag of 0 will be written in the coarse cell %d\n", i + 1);
        tag = 0;
        break;
      }
      tag = cell_tag[coarse_cell_cell[j] - 1];
    }
    (*coarsecell_tag)[i] = tag;
  }

  if(0 == 1) {
    PDM_printf("Affichage de (*coarsecell_tag)\n");
    for (int i = 0; i < nCoarseCellChecked; i++) {
      PDM_printf(" %d ", (*coarsecell_tag)[i]);
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Updates the array face_tag
 *
 * \param [in]    n_faceChecked          Number of faces after refining ( <= n_face)
 * \param [in]    coarse_face_to_fine_face  Coarse face - fine face connectivity (size = n_faceChecked)
 *
 * \param [inout] face_tag               Tag face connectivity index (size = n_face at first and n_faceChecked at the end)
 *
 */

static void
_build_face_tag
(
 int            n_faceChecked,
 int           *coarse_face_to_fine_face,
 int          **face_tag
)
{
  if(*face_tag == NULL) {
    return;
  }

  //Loop over coarse_face_to_fine_face, i = number of a face after refinement
  for (int i = 0; i < n_faceChecked; i++) {
    (*face_tag)[i] = (*face_tag)[coarse_face_to_fine_face[i] - 1];
  }

  (*face_tag) = realloc((*face_tag), n_faceChecked * sizeof(int));

  if(0 == 1) {
    PDM_printf("Contenu de (*face_tag)\n");
    for (int i = 0; i < n_faceChecked; i++) {
        PDM_printf(" %d ", (*face_tag)[i]);
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Updates the array vtx_tag
 *
 * \param [in]    n_vtxChecked          Number of vertices after refining
 * \param [in]    coarse_vtx_to_fine_vtx   Coarse vertex - fine vertex connectivity (size = n_vtxChecked)
 *
 * \param [inout] vtx_tag               Tag vertex connectivity index (size = n_vtx at first and n_vtxChecked at the end)
 *
 */

static void
_build_vtx_tag
(
 int            n_vtxChecked,
 int           *coarse_vtx_to_fine_vtx,
 int          **vtx_tag
)
{
  if(*vtx_tag == NULL) {
    return;
  }

  //Loop over coarse_face_to_fine_face, i = number of a face after refinement
  for (int i = 0; i < n_vtxChecked; i++) {
    (*vtx_tag)[i] = (*vtx_tag)[coarse_vtx_to_fine_vtx[i] - 1];
  }

  (*vtx_tag) = realloc((*vtx_tag), n_vtxChecked * sizeof(int));

  if(0 == 1) {
    PDM_printf("Contenu de (*vtx_tag)\n");
    for (int i = 0; i < n_vtxChecked; i++) {
      PDM_printf(" %d ", (*vtx_tag)[i]);
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Updates the array face_group by renumbering the faces and removing the removed faces
 *
 * \param [in] n_face_group            Number of groups of faces
 * \param [in] fineFaceToCoarseFace  Fine face - coarse face connectivity (size = n_face)
 *
 * \param [inout] face_group          Face group index (size = face_group_idx[n_face_group])
 * \param [inout] face_group_idx       Face group index (size = n_face_group + 1)
 *
 */

static void
_build_faceGroup
(
 int            n_face_group,
 int          **face_group,
 int          **face_group_idx,
 int          **coarse_face_group_to_fine_face_group
)
{

  if(*face_group == NULL || *face_group_idx == NULL || n_face_group == 0) {
    return;
  }
  *coarse_face_group_to_fine_face_group = malloc((*face_group_idx)[n_face_group] * sizeof(int));

  //Renumbering of partGroup from the fine numbering to the coarse one
  //Loop over face_group, i = face number
//  for (int i = 0; i < (*face_group_idx)[n_face_group]; i++) {
//      (*face_group)[i] = fineFaceToCoarseFace[(*face_group)[i] - 1];
//  }

  if (0 == 1) {
    PDM_printf("Content of face_group after renumbering: |");
    for (int i = 0; i < (*face_group_idx)[n_face_group]; i++) {
      PDM_printf(" %d ", (*face_group)[i]);
    }
    PDM_printf("\n");
  }

  int idx = 0;

  //Counter of faces per group
  int *cptFacesPerGroup = malloc(n_face_group * sizeof(int));

  for (int i = 0; i < n_face_group; i++) {
    cptFacesPerGroup[i] = 0;
  }

  //face_group_idx is rebuilt
  //Loop over face_group_idx, i = group number
  for (int i = 0; i < n_face_group; i++) {
    int startNumberingFace = 1;
    //Loop over face_group, j = face number
    for (int j = (*face_group_idx)[i]; j < (*face_group_idx)[i + 1]; j++) {
      //If we do not have a -1, the face has not been removed and is saved
      if ((*face_group)[j] != -1) {
        cptFacesPerGroup[i]++;
        (*face_group)[idx] = (*face_group)[j];
        (*coarse_face_group_to_fine_face_group)[idx++] = startNumberingFace;
      }
      startNumberingFace++;
    }
  }

  if (0 == 1) {
    PDM_printf("Contenu de cptFacesPerGroup apres remplissage\n");
    for (int i = 0; i < n_face_group; i++) {
      PDM_printf(" %d ", cptFacesPerGroup[i]);
    }
    PDM_printf("\n");
  }

  //Update of face_group_idx
  (*face_group_idx)[0] = 0;

  //Loop over cptFacesPerGroup, i = group number
  for (int i = 0; i < n_face_group; i++) {
    (*face_group_idx)[i + 1] = (*face_group_idx)[i] + cptFacesPerGroup[i];
  }

  (*face_group) = realloc((*face_group), (*face_group_idx)[n_face_group] * sizeof(int));
  (*coarse_face_group_to_fine_face_group) = realloc((*coarse_face_group_to_fine_face_group), (*face_group_idx)[n_face_group] * sizeof(int));

  if (0 == 1) {
    PDM_printf("Final content of face_group_idx: ");
    for (int i = 0; i < n_face_group + 1; i++) {
      PDM_printf(" %d ", (*face_group_idx)[i]);
    }
    PDM_printf("\n");

    PDM_printf("Final content of face_group: |");
    for (int i = 0; i < (*face_group_idx)[n_face_group]; i++){
      PDM_printf(" %d ", (*face_group)[i]);
    }
    PDM_printf("\n");

    PDM_printf("Final content of coarse_face_group_to_fine_face_group: ");
    for (int i = 0; i < (*face_group_idx)[n_face_group]; i++) {
      PDM_printf(" %d ", (*coarse_face_group_to_fine_face_group)[i]);
    }
    PDM_printf("\n\n");
  }
  free(cptFacesPerGroup);
}

/**
 *
 * \brief SetUp data array fo r coarse mesh
 *
 * \param [out] cgId              Coarse grid identifier
 *
 * \param [in]  i_part              Partition identifier
 * \param [in]  nCoarseCellWanted  Number of cells in the coarse grid wanted by the user
 * \param [in]  n_cell              Number of cells
 * \param [in]  n_face              Number of faces
 * \param [in]  n_vtx               Number of vertices
 * \param [in]  n_face_group         Number of face groups
 * \param [in]  n_face_part_bound     Number of partitioning boundary faces
 * \param [in]  cell_face_idx        Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [in]  cell_face           Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                             numbering : 1 to n)
 * \param [in]  cell_tag            Cell tag (size = n_cell)
 * \param [in]  cellWeight         Cell weight (size = n_cell)
 * \param [in]  faceWeight         Face weight (size = n_face)
 * \param [in]  cell_ln_to_gn         Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [in]  face_cell           Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [in]  face_vtx_idx         Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [in]  face_vtx            Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [in]  face_tag            Face tag (size = n_face)
 * \param [in]  face_ln_to_gn         Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [in]  vtxCoord           Vertex coordinates (size = 3 * nVertex)
 * \param [in]  vtx_tag             Vertex tag (size = nVertex)
 * \param [in]  vtx_ln_to_gn          Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [in]  face_group_idx       Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [in]  face_group          faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [in]  face_group_ln_to_gn    Faces global numbering for each group
 *                                  (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [in]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 * \param [in]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 * \param [in]  face_part_bound      Partitioning boundary faces (size = 4 * n_face_part_bound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number
 *                                          in the connected partition (numbering :1 to n)
 */

static void
_coarse_grid_mesh_input
(
 _coarse_mesh_t     *cm,
 const int           i_part,
 const int           nCoarseCellWanted,
 const int           n_cell,
 const int           n_face,
 const int           n_vtx,
 const int           n_face_group,
 const int           n_face_part_bound,
 const int          *cell_face_idx,
 const int          *cell_face,
 const int          *cell_tag,
 const int          *cellWeight,
 const int          *faceWeight,
 const PDM_g_num_t  *cell_ln_to_gn,
 const int          *face_cell,
 const int          *face_vtx_idx,
 const int          *face_vtx,
 const int          *face_tag,
 const PDM_g_num_t  *face_ln_to_gn,
 const double       *vtxCoord,
 const int          *vtx_tag,
 const PDM_g_num_t  *vtx_ln_to_gn,
 const int          *face_group_idx,
 const int          *face_group,
 const PDM_g_num_t  *face_group_ln_to_gn,
 const int          *face_part_bound_proc_idx,
 const int          *face_part_bound_part_idx,
 const int          *face_part_bound
)
{
  _part_t * part_ini = cm->part_ini[i_part];
  _coarse_part_t *part_res = cm->part_res[i_part];

  // const int *_cellWeight = cellWeight;
  // const int *_faceWeight = faceWeight;

  part_ini->cellWeight = cellWeight;
  part_ini->faceWeight = faceWeight;

  part_ini->n_vtx = n_vtx;
  part_ini->n_cell = n_cell;
  part_ini->n_face = n_face;
  part_ini->n_face_group = n_face_group;
  part_ini->n_face_part_bound = n_face_part_bound;
  part_ini->cell_face_idx = (int *) cell_face_idx;
  part_ini->cell_face = (int *) cell_face;
  part_ini->cell_tag = (int *) cell_tag;
  part_ini->face_cell = (int *) face_cell;
  part_ini->face_vtx_idx = (int *) face_vtx_idx;
  part_ini->face_vtx = (int *) face_vtx;
  part_ini->face_tag = (int *) face_tag;
  part_ini->face_part_bound_proc_idx = (int *) face_part_bound_proc_idx;
  part_ini->face_part_bound_part_idx = (int *) face_part_bound_part_idx;
  part_ini->face_part_bound = (int *) face_part_bound;
  part_ini->face_group_idx = (int *) face_group_idx;
  part_ini->face_group = (int *) face_group;
  part_ini->vtx = (double *) vtxCoord;
  part_ini->vtx_tag = (int *) vtx_tag;

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
  part_ini->cell_ln_to_gn = (PDM_g_num_t *) cell_ln_to_gn;
  part_ini->face_ln_to_gn = (PDM_g_num_t *) face_ln_to_gn;
  part_ini->face_group_ln_to_gn = (PDM_g_num_t *) face_group_ln_to_gn;
  part_ini->vtx_ln_to_gn = (PDM_g_num_t *) vtx_ln_to_gn;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

  part_res->nCoarseCellWanted = nCoarseCellWanted;

}



/**
 *
 * \brief Build a coarse grid prealably setUp with _coarse_grid_mesh_input
 *
 * \param [out] cgId              Coarse grid identifier
 *
 * \param [in]  i_part              Partition identifier
 * \param [in]  nCoarseCellWanted  Number of cells in the coarse grid wanted by the user
 */

static void
_coarse_grid_compute
(
 _coarse_mesh_t     *cm,
 const int           i_part
)
{

  _part_t * part_ini = cm->part_ini[i_part];
  _coarse_part_t *part_res = cm->part_res[i_part];



  cm->timer = PDM_timer_create();
  for (int i = 0; i < 18; i++) {
    cm->times_elapsed[i] = 0.;
    cm->times_cpu[i] = 0.;
    cm->times_cpu_u[i] = 0.;
    cm->times_cpu_s[i] = 0.;
  }

  PDM_timer_resume(cm->timer);

  int *dualGraphIdx = NULL;
  int *dualGraph    = NULL;

  _dual_graph_from_face_cell(part_ini,
                             (int **) &dualGraphIdx,
                             (int **) &dualGraph);

  int itime = 1;
  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Call Metis or Scotch to get the cell_part array
  //cell_part must be allocated before proceeding (the initialization is part of the split method)

  PDM_timer_resume(cm->timer);

  int *cell_part = NULL;

  int n_coarse_cell_computed;

  _split( cm,
          i_part,
         &n_coarse_cell_computed,
         dualGraphIdx,
         dualGraph,
         (int **) &cell_part);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);

  /* Assign size of multigrid */
  part_res->part->n_cell = n_coarse_cell_computed;

  //  From the cell_part array, get the partCell
  int *partCellIdx = NULL;
  int *partCell = NULL;

  // _partCell_from_cell_part(nCoarseCellWanted,
  _partCell_from_cell_part(part_res->part->n_cell,
                          part_ini->n_cell,
                          cell_part,
                          (int **) &partCellIdx,
                          (int **) &partCell);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Check that all partitions are correctly connected

  PDM_timer_resume(cm->timer);

  int *cellCoarseCell = NULL;

  // part_res->part->n_cell = nCoarseCellWanted;

  _adapt_Connectedness(&(part_res->part->n_cell),
                       part_ini->n_cell,
                       cell_part,
                       (int **) &cellCoarseCell,
                       dualGraph,
                       dualGraphIdx,
                       partCell,
                       partCellIdx,
                       (int **) &(part_res->coarse_cell_cell),
                       (int **) &(part_res->coarse_cell_cell_idx));

  free(partCellIdx);
  free(partCell);

  free(cell_part);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Compress the face_cell array to create the faceCoarseCell array

  PDM_timer_resume(cm->timer);
  int *fineFaceToCoarseFace = NULL;

  //Temporary storage of the data of part_ini
  part_res->part->face_cell = part_ini->face_cell;
  part_res->part->n_face = part_ini->n_face;

  _build_faceCoarseCell(&(part_res->part->n_face),
                        part_ini->face_cell,
                        cellCoarseCell,
                        &(part_res->part->face_cell),
                        (int **) &fineFaceToCoarseFace,
                        &(part_res->coarse_face_to_fine_face));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Updates the face_group_idx and face_group arrays
  PDM_timer_resume(cm->timer);

  part_res->part->face_group_idx = NULL;
  part_res->part->face_group = NULL;
  part_res->part->face_group_ln_to_gn = NULL;

  if (part_ini->n_face_group > 0) {
    part_res->part->face_group_idx = malloc((part_ini->n_face_group + 1) * sizeof(int));
    for (int i = 0; i < (part_ini->n_face_group + 1); i++) {
      part_res->part->face_group_idx[i] = part_ini->face_group_idx[i];
    }
    part_res->part->face_group = malloc(part_res->part->face_group_idx[part_ini->n_face_group] * sizeof(int));
    for (int i = 0; i < part_ini->face_group_idx[part_ini->n_face_group]; i++) {
      part_res->part->face_group[i] = fineFaceToCoarseFace[part_ini->face_group[i] - 1];
    }
  }

  _build_faceGroup(part_ini->n_face_group,
                   &(part_res->part->face_group),
                   &(part_res->part->face_group_idx),
                   &(part_res->coarse_face_group_to_fine_face_group));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Compress the cell_face_idx and cell_face arrays

  PDM_timer_resume(cm->timer);

  _coarsecell_face_from_faceCoarseCell(part_res->part->n_cell,
                                      part_res->part->n_face,
                                      part_res->part->face_cell,
                                      &(part_res->part->cell_face_idx),
                                      &(part_res->part->cell_face));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Compress the face_vtx_idx and face_vtx arrays

  PDM_timer_resume(cm->timer);

  part_res->part->n_vtx = part_ini->n_vtx;

  part_res->part->face_vtx_idx = malloc((part_ini->n_face + 1) * sizeof(int));
  for (int i = 0; i < (part_ini->n_face + 1); i++) {
    part_res->part->face_vtx_idx[i] = part_ini->face_vtx_idx[i];
  }

  part_res->part->face_vtx = malloc(part_res->part->face_vtx_idx[part_ini->n_face] * sizeof(int));
  for (int i = 0; i < part_res->part->face_vtx_idx[part_ini->n_face]; i++) {
    part_res->part->face_vtx[i] = part_ini->face_vtx[i];
  }

  int *fineVtxToCoarseVtx = NULL;

  _build_face_vtx(part_ini->n_face,
                 part_res->part->n_face,
                 part_ini->n_vtx,
                 fineFaceToCoarseFace,
                 &(part_res->part->face_vtx_idx),
                 &(part_res->part->face_vtx),
                 &(part_res->part->n_vtx),
                 (int **) &fineVtxToCoarseVtx,
                 (int **) &(part_res->coarse_vtx_to_fine_vtx));


  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //  Compress the vtxCoord array

  PDM_timer_resume(cm->timer);

  part_res->part->vtx = malloc(3 * part_ini->n_vtx * sizeof(double));
  for (int i = 0; i < 3 * part_ini->n_vtx; i++) {
    part_res->part->vtx[i] = part_ini->vtx[i];
  }

  _build_vtx(part_ini->n_vtx,
             part_res->part->n_vtx,
             fineVtxToCoarseVtx,
             &(part_res->part->vtx));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  //Update the tag arrays

  PDM_timer_resume(cm->timer);

  _build_coarsecell_tag(part_res->part->n_cell,
                       part_res->coarse_cell_cell_idx,
                       part_res->coarse_cell_cell,
                       part_ini->cell_tag,
                       &(part_res->part->cell_tag));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);

  if (part_ini->face_tag != NULL) {
    part_res->part->face_tag = malloc(part_ini->n_face * sizeof(int));
    for (int i = 0; i < part_ini->n_face; i++) {
      part_res->part->face_tag[i] = part_ini->face_tag[i];
    }
  }
  else {
    part_res->part->face_tag = NULL;
  }

  _build_face_tag(part_res->part->n_face,
                 part_res->coarse_face_to_fine_face,
                 &(part_res->part->face_tag));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);
  if (part_ini->vtx_tag != NULL) {
    part_res->part->vtx_tag = malloc(part_ini->n_vtx * sizeof(int));
    for (int i = 0; i < part_ini->n_vtx; i++) {
      part_res->part->vtx_tag[i] = part_ini->vtx_tag[i];
    }
  }
  else {
    part_res->part->vtx_tag = NULL;
  }

  _build_vtx_tag(part_res->part->n_vtx,
                part_res->coarse_vtx_to_fine_vtx,
                &(part_res->part->vtx_tag));

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  /* reordering cells and faces */
  // Conflict with New renumbering for OpenMP/Better vecto - Add part_to_block option
  // HEAD ---------------------
  // part_res->part->new_to_old_order_cell = (int *) malloc (sizeof(int) * part_res->part->n_cell);
  // for (int i = 0; i < part_res->part->n_cell; i++){
  //   part_res->part->new_to_old_order_cell[i] = i;
  // }
  // part_res->part->new_to_old_order_face = (int *) malloc (sizeof(int) * part_res->part->n_face);
  // for (int i = 0; i < part_res->part->n_face; i++){
  //   part_res->part->new_to_old_order_face[i] = i;
  // }
  // Conflict with New renumbering for OpenMP/Better vecto - Add part_to_block option VOID

  free(cellCoarseCell);

  free(dualGraphIdx);
  free(dualGraph);

  free(fineFaceToCoarseFace);

  free(fineVtxToCoarseVtx);
}

/**
 *
 * \brief Updates the cell_ln_to_gn array for the coarse mesh are save it in the coarse mesh partition
 *
 * \param [in]   cm                 Coarse mesh
 */

static void
_build_coarsecell_ln_to_gn
(
 _coarse_mesh_t * cm
 )
{
  //Calculation of the number of cells on the processor
  PDM_g_num_t n_cellProc = 0;

  //Loop over the partition numbers, i = partition number
  for (int i = 0; i < cm->n_part; i++) {
    n_cellProc += cm->part_res[i]->part->n_cell;
  }

  //    PDM_printf("\nValeur de n_cellProc : %d \n", n_cellProc);

  //Global numbering of the cells
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&n_cellProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

  //Index to position the local cells
  beg_NumAbs -= n_cellProc;

  int idx_write = 0;

  //Loop over the partition numbers, i = partition number
  for (int i = 0; i < cm->n_part; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    int n_cell = cm->part_res[i]->part->n_cell;
    cmp->cell_ln_to_gn = (PDM_g_num_t *) malloc(n_cell * sizeof(PDM_g_num_t));
    //Loop over the partition cells, j = cell number
    for (int j = 0; j < n_cell; j++) {
      cmp->cell_ln_to_gn[j] = beg_NumAbs + idx_write + 1;
      idx_write++;
    }

  }

  if(0 == 1) {
    for (int i_part = 0; i_part < cm->n_part; i_part++) {
      PDM_printf("\nContenu de cm->part_res[%d]->part->cell_ln_to_gn\n", i_part);
      for (int j = 0; j < cm->part_res[i_part]->part->n_cell; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cm->part_res[i_part]->part->cell_ln_to_gn[j]);
      }
      PDM_printf("\n\n");
    }
  }
}

/**
 *
 * \brief Updates the face_ln_to_gn array for the coarse mesh are save it in the coarse mesh partition
 *
 * \param [in]   cm                 Coarse mesh
 *
 */

static void
_build_face_ln_to_gn
(
_coarse_mesh_t * cm
)
{
  PDM_g_num_t **face_ln_to_gn_part = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));
  int *n_facePart = (int *) malloc(cm->n_part * sizeof(int));

  for (int i = 0; i < cm->n_part; i++) {
    face_ln_to_gn_part[i] = cm->part_ini[i]->face_ln_to_gn;
    n_facePart[i] = cm->part_ini[i]->n_face;
  }

  if(0 == 1) {
    PDM_printf("Contenu de face_ln_to_gn_part\n");
    for (int i = 0; i < cm->n_part; i++) {
      for (int j = 0; j < n_facePart[i]; j++) {
         PDM_printf(" "PDM_FMT_G_NUM" ", face_ln_to_gn_part[i][j]);
      }
    PDM_printf("\n");
    }

  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                     1.,
                                                     (PDM_g_num_t **) face_ln_to_gn_part,
                                                       NULL,
                                                     n_facePart,
                                                     cm->n_part,
                                                     cm->comm);

  PDM_g_num_t **face_ln_to_gnTag = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));

  int idx_write = 0;

  for (int i = 0; i < cm->n_part; i++) {
    idx_write = 0;
    face_ln_to_gnTag[i] = (PDM_g_num_t *) malloc(cm->part_ini[i]->n_face * sizeof(PDM_g_num_t));
      //Loop over coarse_face_to_fine_face, i = index of coarse_face_to_fine_face (from 0 to cm->part_res[i_part]->part->n_face)

    for (int j = 0; j < cm->part_ini[i]->n_face; j++) {
      face_ln_to_gnTag[i][j] = -1;
    }

    for (int j = 0; j < cm->part_res[i]->part->n_face; j++) {
          //If the vertex studied is the same as in coarse_face_to_fine_face, it is to be stored
     int k =  cm->part_res[i]->coarse_face_to_fine_face[j] - 1;
     face_ln_to_gnTag[i][k] = 0;
    }
  }

  if(0 == 1) {
    PDM_printf("Contenu de face_ln_to_gnTag\n");
    for (int i = 0; i < cm->n_part; i++) {
      for (int j = 0; j < cm->part_res[i]->part->n_face; j++) {
          PDM_printf(" "PDM_FMT_G_NUM" ", face_ln_to_gnTag[i][j]);
      }
    }
    PDM_printf("\n");
  }

  PDM_g_num_t *b_tIntersects = NULL;
  int *b_stride_one = NULL;
  int *part_stride = NULL;

  PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            1,
                            &part_stride,
                            (void **) face_ln_to_gnTag,
                            &b_stride_one,
                            (void **) &b_tIntersects);

  //Calculation of the number of faces on the processor
  PDM_g_num_t n_faceProc = 0;

  int size_block = PDM_part_to_block_n_elt_block_get(ptb);
  for (int i = 0; i < size_block; i++) {
    //If the face has not been removed
     if(b_tIntersects[i] == 0) {
         n_faceProc++;
     }
  }

  //Global numbering of the faces
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&n_faceProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

  //Index to position the local vertices
  beg_NumAbs -= n_faceProc;

  idx_write = 0;

  //Loop over the partition numbers, i = partition number

  for (int i = 0; i < size_block; i++) {
    //If the vertex has not been removed
    if(b_tIntersects[i] == 0) {
        b_tIntersects[i] = beg_NumAbs + (idx_write++) + 1;

    }
    else {
        b_tIntersects[i] = -1;
    }
  }

  PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb);

  PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                       (const PDM_g_num_t **) face_ln_to_gn_part,
                                                       n_facePart,
                                                       cm->n_part,
                                                       cm->comm);

  PDM_g_num_t  **face_ln_to_gnFine = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));

  for (int i = 0; i < cm->n_part; i++) {
    face_ln_to_gnFine[i] = (PDM_g_num_t *) malloc(cm->part_ini[i]->n_face * sizeof(PDM_g_num_t));
  }

  int stride_one = 1;

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &stride_one,
                          (void *) b_tIntersects,
                          &part_stride,
                          (void **) face_ln_to_gnFine);

  if(0 == 1) {
    PDM_printf("\nContenu de face_ln_to_gnFine\n");
    for (int i = 0; i < cm->n_part; i++) {
      //Loop over the partition faces, j = face number
      for (int j = 0; j < cm->part_ini[i]->n_face; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", face_ln_to_gnFine[i][j]);
      }
    }
    PDM_printf("\n");
  }

  for (int i = 0; i < cm->n_part; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    int n_face = cm->part_res[i]->part->n_face;
    cmp->face_ln_to_gn = (PDM_g_num_t *) malloc(n_face * sizeof(PDM_g_num_t));
    for (int j = 0; j < n_face; j++) {
      cmp->face_ln_to_gn[j] = (PDM_g_num_t) face_ln_to_gnFine[i][cm->part_res[i]->coarse_face_to_fine_face[j] - 1];
    }
  }

  if(0 == 1) {
    PDM_printf("\nContenu de face_ln_to_gn de la structure\n");
    for (int i = 0; i < cm->n_part; i++) {
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cm->part_res[i]->part->n_face; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cm->part_res[i]->part->face_ln_to_gn[j]);
      }
    }
    PDM_printf("\n");
  }

  free(face_ln_to_gn_part);
  free(n_facePart);

  for (int i = 0; i < cm->n_part; i++) {
    free(face_ln_to_gnTag[i]);
  }
  free(face_ln_to_gnTag);

  for (int i = 0; i < cm->n_part; i++) {
    free(face_ln_to_gnFine[i]);
  }
  free(face_ln_to_gnFine);

  free (b_stride_one);
  free (part_stride);
  free (b_tIntersects);

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
}

/**
 *
 * \brief Updates the vtx_ln_to_gn array for the coarse mesh are save it in the coarse mesh partition
 *
 * \param [in]   cm                 Coarse mesh
 *
 */

static void
_build_vtx_ln_to_gn
(
_coarse_mesh_t * cm
)
{
  PDM_g_num_t **vtx_ln_to_gn_part = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));
  int *n_vtxPart = (int *) malloc(cm->n_part * sizeof(int));

  for (int i = 0; i < cm->n_part; i++) {
    vtx_ln_to_gn_part[i] = cm->part_ini[i]->vtx_ln_to_gn;
    n_vtxPart[i] = cm->part_ini[i]->n_vtx;
  }

  if(0 == 1) {
    PDM_printf("Contenu de vtx_ln_to_gn_part\n");
    for (int i = 0; i < cm->n_part; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", *(vtx_ln_to_gn_part[i]));
    }
    PDM_printf("\n");

    PDM_printf("Contenu de n_vtxPart\n");
    for (int i = 0; i < cm->n_part; i++) {
      PDM_printf(" %d ", n_vtxPart[i]);
    }
    PDM_printf("\n");
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                     1.,
                                                     (PDM_g_num_t **) vtx_ln_to_gn_part,
                                                       NULL,
                                                     n_vtxPart,
                                                     cm->n_part,
                                                     cm->comm);

  PDM_g_num_t **vtx_ln_to_gnTag = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));

  int idx_write = 0;

  for (int i = 0; i < cm->n_part; i++) {
    int nFineVtx = cm->part_ini[i]->n_vtx;
    int nCoarseVtx = cm->part_res[i]->part->n_vtx;
    idx_write = 0;
    vtx_ln_to_gnTag[i] = (PDM_g_num_t *) malloc(nFineVtx * sizeof(PDM_g_num_t));
    //Loop over coarse_face_to_fine_face, i = index of coarse_vtx_to_fine_vtx (from 0 to cm->part_res[i_part]->part->n_vtx)
    for (int j = 0; j < nFineVtx; j++) {
      vtx_ln_to_gnTag[i][j] = -1;
    }

    for (int j = 0; j < nCoarseVtx; j++) {
        //If the vertex studied is the same as in coarse_vtx_to_fine_vtx, it is to be stored
      int k = cm->part_res[i]->coarse_vtx_to_fine_vtx[j] - 1;
      vtx_ln_to_gnTag[i][k] = 0;
    }
  }

  if (0 == 1) {
    PDM_printf("Contenu de vtx_ln_to_gnTag\n");
    for (int i = 0; i < cm->n_part; i++) {
      for (int j = 0; j < cm->part_res[i]->part->n_vtx; j++) {
        PDM_printf(" "PDM_FMT_G_NUM, vtx_ln_to_gnTag[i][j]);
      }
    }
    PDM_printf("\n");
  }

  PDM_g_num_t *b_tIntersects = NULL;
  int *b_stride_one = NULL;
  int *part_stride = NULL;

  PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            1,
                            &part_stride,
                            (void **) vtx_ln_to_gnTag,
                            &b_stride_one,
                            (void **) &b_tIntersects);

  //Calculation of the number of vertices on the processor
  PDM_g_num_t n_vtxProc = 0;

  int size_block = PDM_part_to_block_n_elt_block_get(ptb);
  for (int i = 0; i < size_block; i++) {
    //If the vertex has not been removed
    if(b_tIntersects[i] == 0) {
      n_vtxProc++;
    }
  }

  //Global numbering of the vertices
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&n_vtxProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

  //Index to position the local vertices
  beg_NumAbs -= n_vtxProc;

  idx_write = 0;

  //Loop over the partition numbers, i = partition number
  for (int i = 0; i < size_block; i++) {
    //If the vertex has not been removed
    if (b_tIntersects[i] == 0) {
      b_tIntersects[i] = beg_NumAbs + (idx_write++) + 1;
    }
    else {
      b_tIntersects[i] = -1;
    }
  }

  PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb);

  PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                       (const PDM_g_num_t **) vtx_ln_to_gn_part,
                                                       n_vtxPart,
                                                       cm->n_part,
                                                       cm->comm);

  PDM_g_num_t  **vtx_ln_to_gnFine = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));

  for (int i = 0; i < cm->n_part; i++) {
    vtx_ln_to_gnFine[i] = (PDM_g_num_t *) malloc(cm->part_ini[i]->n_vtx * sizeof(PDM_g_num_t));
  }

  int stride_one = 1;

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &stride_one,
                          (void *) b_tIntersects,
                          &part_stride,
                          (void **) vtx_ln_to_gnFine);

  if(0 == 1) {
    PDM_printf("\nContenu de vtx_ln_to_gnFine\n");
    for (int i = 0; i < cm->n_part; i++) {
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cm->part_ini[i]->n_vtx; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", vtx_ln_to_gnFine[i][j]);
      }
    }
    PDM_printf("\n");
  }

  for (int i = 0; i < cm->n_part; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    int n_vtx = cm->part_res[i]->part->n_vtx;
    cmp->vtx_ln_to_gn = (PDM_g_num_t *) malloc(n_vtx * sizeof(PDM_g_num_t));
    for (int j = 0; j < n_vtx; j++) {
      cmp->vtx_ln_to_gn[j] = (PDM_g_num_t) vtx_ln_to_gnFine[i][cm->part_res[i]->coarse_vtx_to_fine_vtx[j] - 1];
    }
  }

  if(0 == 1) {
    PDM_printf("\nContenu de vtx_ln_to_gn de la structure\n");
    for (int i = 0; i < cm->n_part; i++) {
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cm->part_res[i]->part->n_vtx; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cm->part_res[i]->part->vtx_ln_to_gn[j]);
      }
    }
    PDM_printf("\n");
  }

  free(vtx_ln_to_gn_part);
  free(n_vtxPart);

  for (int i = 0; i < cm->n_part; i++) {
    free(vtx_ln_to_gnTag[i]);
  }
  free(vtx_ln_to_gnTag);

  for (int i = 0; i < cm->n_part; i++) {
    free(vtx_ln_to_gnFine[i]);
  }
  free(vtx_ln_to_gnFine);

  free (b_stride_one);
  free (part_stride);
  free (b_tIntersects);

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
}

/**
 *
 * \brief Updates the face_group_ln_to_gn array for the coarse mesh are save it in the coarse mesh partition
 *
 * \param [in]   cm                 Coarse mesh
 *
 */

static void
_build_faceGroupLNToGN
(
_coarse_mesh_t * cm
)
{
  for (int i = 0; i < cm->n_part; i++) {
    //Si un des face_group_idx est NULL, on n'utilise pas les groupes
    //On quitte donc la boucle et la fonction !
    if(cm->part_ini[i]->face_group_idx == NULL)  {
      return;
    }
  }

  int **face_ln_to_gnTag = (int **) malloc(cm->n_part * sizeof(int *));

  int idx_write = 0;

  for (int i = 0; i < cm->n_part; i++) {
    idx_write = 0;
    face_ln_to_gnTag[i] = (int *) malloc(cm->part_ini[i]->n_face * sizeof(int));
    int n_face = cm->part_res[i]->part->n_face;
    for (int j = 0; j < cm->part_ini[i]->n_face; j++) {
      face_ln_to_gnTag[i][j] = -1;
    }
    //Loop over coarse_face_to_fine_face, i = index of coarse_face_to_fine_face (from 0 to cm->part_res[i_part]->part->n_face)
    for (int j = 0; j < n_face; j++) {
      //If the face studied is the same as in coarse_face_to_fine_face, it is to be stored
      int k =  cm->part_res[i]->coarse_face_to_fine_face[j] - 1;
      face_ln_to_gnTag[i][k] = 0;
    }
  }

  if (0 == 1) {
    PDM_printf("Contenu de face_ln_to_gnTag\n");
    for (int i = 0; i < cm->n_part; i++) {
      for (int j = 0; j < cm->part_ini[i]->n_face; j++) {
        PDM_printf(" %d ", face_ln_to_gnTag[i][j]);
      }
    }
    PDM_printf("\n");
  }

  int **faceGroupLNToGNTag = (int **) malloc(cm->n_part * sizeof(int *));

  idx_write = 0;

  int _rank;

  PDM_MPI_Comm_rank (cm->comm, &_rank);

  for (int i = 0; i < cm->n_part; i++) {
    _part_t *cmp_coarse = cm->part_res[i]->part;
    _part_t *cmp_fine = cm->part_ini[i];
    fflush(stdout);
    faceGroupLNToGNTag[i] = (int *) malloc(cmp_coarse->face_group_idx[cm->n_face_group] * sizeof(int));

    //Loop over face_group_idx, i = index of group of faces (from 0 to cm->part_res[i_part]->part->n_face_group)
    for (int j = 0; j < cmp_fine->face_group_idx[cm->n_face_group]; j++) {
      faceGroupLNToGNTag[i][j] = face_ln_to_gnTag[i][cmp_fine->face_group[j] - 1];
    }
  }

  if(0 == 1) {
    PDM_printf("Contenu de faceGroupLNToGNTag\n");
    for (int i = 0; i < cm->n_part; i++) {
      for (int j = 0; j < cm->part_res[i]->part->face_group_idx[cm->n_face_group]; j++) {
        PDM_printf(" %d ", faceGroupLNToGNTag[i][j]);
      }
    }
    PDM_printf("\n");
  }

  for (int i = 0; i < cm->n_part; i++) {
    free(face_ln_to_gnTag[i]);
  }
  free(face_ln_to_gnTag);

  for (int i = 0; i < cm->n_part; i++) {
    _part_t *cmp = cm->part_res[i]->part;
    cmp->face_group_ln_to_gn = (PDM_g_num_t *) malloc(cmp->face_group_idx[cm->n_face_group] * sizeof(PDM_g_num_t));
  }

  for(int iGroup = 0; iGroup < cm->n_face_group; iGroup++) {

    PDM_g_num_t **faceGroupLNToGn_part = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));
    int *n_face_groupPart = (int *) malloc(cm->n_part * sizeof(int));

    for (int i = 0; i < cm->n_part; i++) {
      _part_t *cmp = cm->part_ini[i];
      faceGroupLNToGn_part[i] = &(cmp->face_group_ln_to_gn[cmp->face_group_idx[iGroup]]);
      n_face_groupPart[i] = cmp->face_group_idx[iGroup + 1] - cmp->face_group_idx[iGroup];
    }

    if(0 == 1) {
      PDM_printf("Contenu de faceGroupLNToGn_part\n");
      for (int i = 0; i < cm->n_part; i++) {
        int n_faceCurrentGroup = cm->part_ini[i]->face_group_idx[iGroup+1] - cm->part_ini[i]->face_group_idx[iGroup];
        for (int j = 0; j < n_faceCurrentGroup; j++)
          PDM_printf(" "PDM_FMT_G_NUM" ", faceGroupLNToGn_part[i][j]);
        PDM_printf("\n");
      }

      PDM_printf("Contenu de n_face_groupPart\n");
      for (int i = 0; i < cm->n_part; i++) {
        PDM_printf(" %d ", n_face_groupPart[i]);
      }
      PDM_printf("\n");
    }

    int rank;
    PDM_MPI_Comm_rank(cm->comm, &rank);

    PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                         PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                         1.,
                                                         (PDM_g_num_t **) faceGroupLNToGn_part,
                                                         NULL,
                                                         n_face_groupPart,
                                                         cm->n_part,
                                                         cm->comm);

    PDM_g_num_t *b_tIntersects = NULL;
    int *b_stride_one = NULL;
    int *part_stride = NULL;

    PDM_g_num_t **faceGroupLNToGNTagGroup = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));

    for (int i = 0; i < cm->n_part; i++) {
      _part_t *cmp = cm->part_res[i]->part;
      int n_facePerGroup = cmp->face_group_idx[iGroup + 1] - cmp->face_group_idx[iGroup];
      faceGroupLNToGNTagGroup[i] = (PDM_g_num_t *) malloc(n_facePerGroup * sizeof(PDM_g_num_t));

      idx_write = 0;
      //Copy of the sub-array face_group_ln_to_gn for each group
      for (int j = cmp->face_group_idx[iGroup]; j < cmp->face_group_idx[iGroup + 1]; j++) {
        faceGroupLNToGNTagGroup[i][idx_write++] = faceGroupLNToGNTag[i][j];
      }
    }

    if (0 == 1) {
      PDM_printf("Contenu de faceGroupLNToGNTagGroup\n");
      for (int i = 0; i < cm->n_part; i++) {
        int n_facePerGroup = cm->part_res[i]->part->face_group_idx[iGroup + 1] - cm->part_res[i]->part->face_group_idx[iGroup];
        for (int j = 0; j < n_facePerGroup; j++) {
          PDM_printf(" %d ", faceGroupLNToGNTag[i][j]);
        }
      }
      PDM_printf("\n");
    }

    PDM_part_to_block_exch (ptb,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            1,
                            &part_stride,
                            (void **) faceGroupLNToGNTagGroup,
                            &b_stride_one,
                            (void **) &b_tIntersects);

    //Calculation of the number of faces for all the groups on the processor
    PDM_g_num_t n_face_groupProc = 0;

    int size_block = PDM_part_to_block_n_elt_block_get(ptb);

    for (int i = 0; i < size_block; i++) {
      //If the face of the group has not been removed
      if (b_tIntersects[i] == 0) {
        n_face_groupProc++;
      }

    }

    //Global numbering of the vertices
    PDM_g_num_t beg_NumAbs;

    PDM_MPI_Scan(&n_face_groupProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cm->comm);

    //Index to position the local vertices
    beg_NumAbs -= n_face_groupProc;

    idx_write = 0;

    //Loop over the partition numbers, i = partition number

    for (int i = 0; i < size_block; i++) {
      //If the vertex has not been removed
      if(b_tIntersects[i] == 0) {
        b_tIntersects[i] = beg_NumAbs + (idx_write++) + 1;

      }
      else {
        b_tIntersects[i] = -1;
      }

    }

    PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb);

    //        PDM_printf("assert : [%d] %d %d\n", rank, size_block, blockDistribIdx[rank+1] - blockDistribIdx[rank] );
    //        assert(blockDistribIdx[rank+1] - blockDistribIdx[rank] ==  size_block);

    PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                         (const PDM_g_num_t **) faceGroupLNToGn_part,
                                                         n_face_groupPart,
                                                         cm->n_part,
                                                         cm->comm);

    PDM_g_num_t  **faceGroupLNToGNFine = (PDM_g_num_t **) malloc(cm->n_part * sizeof(PDM_g_num_t *));

    for (int i = 0; i < cm->n_part; i++) {
      _part_t *cmp = cm->part_ini[i];
      int n_facePerGroup = cmp->face_group_idx[iGroup + 1] - cmp->face_group_idx[iGroup];
      faceGroupLNToGNFine[i] = (PDM_g_num_t *) malloc(n_facePerGroup * sizeof(PDM_g_num_t));
    }

    int stride_one = 1;

    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &stride_one,
                            (void *) b_tIntersects,
                            &part_stride,
                            (void **) faceGroupLNToGNFine);

    if(0 == 1) {
      PDM_printf("\nContenu de faceGroupLNToGNFine\n");
      for (int i = 0; i < cm->n_part; i++) {
        _part_t *cmp = cm->part_ini[i];
        //Loop over the partition vertices, j = vertex number
        int n_facePerGroup = cmp->face_group_idx[iGroup + 1] - cmp->face_group_idx[iGroup];
        for (int j = 0; j < n_facePerGroup; j++) {
          PDM_printf(" "PDM_FMT_G_NUM" ", faceGroupLNToGNFine[i][j]);
        }
      }
      PDM_printf("\n");
    }

    for (int i = 0; i < cm->n_part; i++) {
      _part_t *cmp = cm->part_res[i]->part;int n_facePerGroupCoarse = cmp->face_group_idx[iGroup + 1] - cmp->face_group_idx[iGroup];
      for (int j = 0; j < n_facePerGroupCoarse; j++) {
        int idxConcatenation = cmp->face_group_idx[iGroup];
        cmp->face_group_ln_to_gn[idxConcatenation + j] =
                faceGroupLNToGNFine[i][cm->part_res[i]->coarse_face_group_to_fine_face_group[idxConcatenation + j] - 1];
      }
    }

    free(faceGroupLNToGn_part);
    free(n_face_groupPart);

    for (int i = 0; i < cm->n_part; i++) {
      free(faceGroupLNToGNTagGroup[i]);
    }
    free(faceGroupLNToGNTagGroup);

    for (int i = 0; i < cm->n_part; i++) {
      free(faceGroupLNToGNFine[i]);
    }
    free(faceGroupLNToGNFine);

    free (b_stride_one);
    free (part_stride);
    free (b_tIntersects);

    PDM_part_to_block_free(ptb);
    PDM_block_to_part_free(btp);

  }

  for (int i = 0; i < cm->n_part; i++) {
    free(faceGroupLNToGNTag[i]);
  }
  free(faceGroupLNToGNTag);

  if (0 == 1) {
    for (int i = 0; i < cm->n_part; i++) {
      PDM_printf("\nContenu de face_group_ln_to_gn de la structure %d\n", i);
      _part_t *cmp = cm->part_res[i]->part;
      //Loop over the partition vertices, j = vertex number
      for (int j = 0; j < cmp->face_group_idx[cm->n_face_group]; j++) {
        PDM_printf(" "PDM_FMT_G_NUM" ", cmp->face_group_ln_to_gn[j]);
      }
    }
    PDM_printf("\n");
  }

}

/**
 *
 * \brief Updates the face_part_bound array for the coarse mesh are save it in the coarse mesh partition
 *
 * \param [in]   cm                 Coarse mesh
 *
 */

static void
_build_facePartBound
(
_coarse_mesh_t * cm
)
{
  PDM_printf("_build_facePartBound \n");
  //Number of processors
  int n_proc;
  PDM_MPI_Comm_size(cm->comm, &n_proc);

  //Loop over the partitions of part_ini and part_res
  for (int i_part = 0; i_part < cm->n_part; i_part++) {
    _part_t *cmp_fine = cm->part_ini[i_part];
    _part_t *cmp_coarse = cm->part_res[i_part]->part;

    //Copy of face_part_bound_part_idx for the coarse mesh
    int *coarseFacePartBoundPartIdx = (int *) malloc((cm->n_total_part + 1) * sizeof(int));
    for (int i = 0; i < cm->n_total_part + 1; i++) {
      coarseFacePartBoundPartIdx[i] = cmp_fine->face_part_bound_part_idx[i];
    }
    cmp_coarse->face_part_bound_part_idx = coarseFacePartBoundPartIdx;

    //Copy of face_part_bound_proc_idx for the coarse mesh
    int *coarseFacePartBoundProcIdx = (int *) malloc((n_proc + 1) * sizeof(int));
    for (int i = 0; i < n_proc + 1; i++) {
      coarseFacePartBoundProcIdx[i] = cmp_fine->face_part_bound_proc_idx[i];
    }
    cmp_coarse->face_part_bound_proc_idx = coarseFacePartBoundProcIdx;
  }

  int **fineFaceToCoarseFace = (int **) malloc(cm->n_part * sizeof(int *));

  for (int i_part = 0; i_part < cm->n_part; i_part++) {
    _part_t *cmp_fine = cm->part_ini[i_part];
    _part_t *cmp_coarse = cm->part_res[i_part]->part;

    //Creation of fineFaceToCoarseFace
    fineFaceToCoarseFace[i_part] = (int *) malloc(cmp_fine->n_face * sizeof(int));

    //Initialization to -1
    for (int i = 0; i < cmp_fine->n_face; i++) {
      fineFaceToCoarseFace[i_part][i] = -1;
    }

    //Loop over coarse_face_to_fine_face
    for (int i = 0; i < cmp_coarse->n_face; i++) {
      int fineFace = cm->part_res[i_part]->coarse_face_to_fine_face[i] - 1;
      fineFaceToCoarseFace[i_part][fineFace] = i + 1;
    }

    if(0 == 1) {
      PDM_printf("Final content of fineFaceToCoarseFace[%d]: \n",i_part);
      PDM_printf("Valeur de cm->part_ini[i_part]->n_face : %d \n", cm->part_ini[i_part]->n_face);
      for (int i = 0; i < cm->part_ini[i_part]->n_face; i++) {
        PDM_printf(" %d ", fineFaceToCoarseFace[i_part][i]);
      }
      PDM_printf("\n");
      PDM_printf("------------------------------------------\n\n");
    }

  }

  int *sendIdx = malloc(sizeof(int) * n_proc);
  int *sendN = malloc(sizeof(int) * n_proc);
  PDM_g_num_t *sendBuff = NULL;

  int *recvIdx = malloc(sizeof(int) * n_proc);
  int *recvN = malloc(sizeof(int) * n_proc);
  PDM_g_num_t *recvBuff = NULL;

  for (int i = 0; i < n_proc; i++) {
    sendN[i] = 0;
  }

  int n_t_send = 0;
  for (int i = 0; i < cm->n_part; i++)  {
    _part_t *cmp = cm->part_ini[i];
    n_t_send += cmp->n_face_part_bound;

    for (int j = 0; j < cmp->n_face_part_bound; j++) {
      int iProc    = cmp->face_part_bound[4*j + 1];
      /* int i_part    = cmp->face_part_bound[4*j + 2]; */
      /* int iFacDist = cmp->face_part_bound[4*j + 3]; */
      sendN[iProc] += 1;
    }
  }

  int n_t_recv = 0;
  for (int i = 0; i < cm->n_part; i++)  {
    _part_t *cmp_coarse = cm->part_res[i]->part;
    int n_face_part_bound = cm->part_ini[i]->n_face_part_bound;
    cmp_coarse->n_face_part_bound = n_face_part_bound;

    n_t_recv += cmp_coarse->n_face_part_bound;
  }

  sendIdx[0] = 0;
  for (int i = 1; i < n_proc; i++) {
    sendIdx[i] = sendIdx[i - 1] + sendN[i - 1];
    sendN[i - 1] = 0;
  }

  sendN[n_proc - 1] = 0;

  sendBuff = malloc(sizeof(PDM_g_num_t) * n_t_send * 3);
  recvBuff = malloc(sizeof(PDM_g_num_t) * n_t_recv * 3);

  int **iFaceLocToIPartBound = (int **) malloc(cm->n_part * sizeof(int *));

  for (int i = 0; i < cm->n_part; i++) {
    _part_t *cmp = cm->part_ini[i];
    //Creation of iFaceLocToIPartBound
    iFaceLocToIPartBound[i] = (int *) malloc(cmp->n_face * sizeof(int));

    //Initialization to -1
    for (int cpt = 0; cpt < cmp->n_face; cpt++) {
      iFaceLocToIPartBound[i][cpt] = -1;
    }

    if(0 == 1) {
      PDM_printf("Valeur de n_face_part_bound : %d \n", cmp->n_face_part_bound);

      PDM_printf("Contenu de face_part_bound initial de la partition %d\n", i);
      for (int cpt = 0; cpt < 4* cmp->n_face_part_bound; cpt++) {
        if (cpt % 4 == 0)
          PDM_printf("|");
        PDM_printf(" %d ", cmp->face_part_bound[cpt]);
      }
      PDM_printf("\n");
    }

    for (int j = 0; j < cmp->n_face_part_bound; j++) {
      int iFacLoc  = cmp->face_part_bound[4*j    ];
      int iProc    = cmp->face_part_bound[4*j + 1];
      int i_part    = cmp->face_part_bound[4*j + 2];
      int iFacDist = cmp->face_part_bound[4*j + 3];

      iFaceLocToIPartBound[i][iFacLoc - 1] = j;

      int id = sendIdx[iProc] + sendN[iProc];

      ++sendN[iProc];
//
      sendBuff[3*id    ] = i_part;
      sendBuff[3*id + 1] = iFacDist;
      sendBuff[3*id + 2] = fineFaceToCoarseFace[i][iFacLoc - 1];

    }

    if (0 == 1) {
      PDM_printf("Contenu de iPartBoundToIFacLoc \n");
      for (int j = 0; j < cmp->n_face_part_bound; j++) { //Modif : olp->nLinkedFace = cmp->n_face_part_bound
        PDM_printf(" %d ", cmp->face_part_bound[4*j    ]);
      }
      PDM_printf("\n");

      PDM_printf("Final content of iFaceLocToIPartBound[i]: \n");
      for (int cpt = 0; cpt < cm->part_ini[i]->n_face; cpt++) {
        PDM_printf(" %d ", iFaceLocToIPartBound[i][cpt]);
      }
      PDM_printf("\n");
    }

  }

  for (int i = 0; i < n_proc; i++) {
    sendN[i] *= 3;
    sendIdx[i] *= 3;
  }

  PDM_MPI_Alltoall (sendN, 1, PDM_MPI_INT,
                recvN, 1, PDM_MPI_INT,
                cm->comm);
  recvIdx[0] = 0;

  for (int i = 1; i < n_proc; i++)  {
    recvIdx[i] = recvIdx[i - 1] + recvN[i - 1];
  }

  PDM_MPI_Alltoallv(sendBuff, sendN, sendIdx, PDM__PDM_MPI_G_NUM,
                recvBuff, recvN, recvIdx, PDM__PDM_MPI_G_NUM, cm->comm);

  //Loop over the partitions
  for (int i_part = 0; i_part < cm->n_part; i_part++) {
    _part_t *cmp = cm->part_ini[i_part];
    _part_t *cmp_coarse = cm->part_res[i_part]->part;
    //Memory allocation of coarseFacePartBound (linked to part_res after)
    int *coarseFacePartBound = (int  *) malloc(4 * cmp->n_face_part_bound * sizeof(int));

    //Copy of the face_part_bound of part_ini
    for (int i = 0; i < 4 * cmp->n_face_part_bound; i++) {
      coarseFacePartBound[i] = cmp->face_part_bound[i];
    }
    cmp_coarse->face_part_bound = coarseFacePartBound;
  }

  //Loop over recvBuff
  for (int i = 0; i < (recvIdx[n_proc - 1] + recvN[n_proc - 1]) / 3; i++) {
    int iPartLoc       = (int) recvBuff[3 * i    ];
    int iFacLocFine    = (int) recvBuff[3 * i + 1];
    int iFacDistCoarse = (int) recvBuff[3 * i + 2];

    int posFacePartBoundCoarse = iFaceLocToIPartBound[iPartLoc - 1][iFacLocFine - 1];

    //Update of the data about faces in the face_part_bound of part_res
    _part_t *cmp = cm->part_res[iPartLoc - 1]->part;
    cmp->face_part_bound[4 * posFacePartBoundCoarse    ] = fineFaceToCoarseFace[iPartLoc - 1][iFacLocFine - 1];
    cmp->face_part_bound[4 * posFacePartBoundCoarse + 3] = iFacDistCoarse;
  }

  if(0 == 1) {
    for (int i_part = 0; i_part < cm->n_part; i_part++) {
      PDM_printf("\nContent of face_part_bound of part_res[%d]\n", i_part);
      PDM_printf("Valeur de cm->part_res[%d]->part->n_face : %d\n",i_part, cm->part_res[i_part]->part->n_face);
      for (int i = 0; i < 4 * cm->part_res[i_part]->part->n_face_part_bound; i++) {
        if (i % 4 == 0)
          PDM_printf("|");
        PDM_printf(" %d ", cm->part_res[i_part]->part->face_part_bound[i]);
      }
      PDM_printf("\n");
    }
    PDM_printf("\n");
  }

  for (int i_part = 0; i_part < cm->n_part; i_part++) { //Modif : n_partB => cm->n_part

    free(fineFaceToCoarseFace[i_part]);
    free(iFaceLocToIPartBound[i_part]);
  }

  free(fineFaceToCoarseFace);
  free(iFaceLocToIPartBound);

  free(sendN);
  free(sendIdx);
  free(sendBuff);

  free(recvN);
  free(recvIdx);
  free(recvBuff);
}

/**
 *
 * \brief Displays all the arrays of a partition of type _part_t
 *
 * \param [in]  n_part        Number of partitions to define on this process
 * \param [in]  n_total_part       Total number of partitions
 * \param [in]  n_face_group   Number of boundaries
 *
 */

static void
_part_display
(
  _part_t *part,
  int n_part,
  int n_total_part,
  int n_face_group
)
{
  if (part == NULL) {
    PDM_printf("Incorrect part to display\n");
    return;
  }

  PDM_printf("Value of n_vtx : %d \n", part->n_vtx);
  PDM_printf("Value of n_cell : %d \n", part->n_cell);
  PDM_printf("Value of n_face : %d \n", part->n_face);
  PDM_printf("Value of n_face_part_bound : %d \n", part->n_face_part_bound);

  if (part->cell_face_idx != NULL) {
    PDM_printf("\nContent of cell_face_idx\n");
    for(int i = 0; i < part->n_cell + 1; i++) {
      PDM_printf(" %d ", part->cell_face_idx[i]);
    }
    PDM_printf("\n");
  }

  if (part->gcell_face != NULL) {
    PDM_printf("\nContent of gcell_face\n");
    for(int i = 0; i < part->cell_face_idx[part->n_cell]; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->gcell_face[i]);
    }
    PDM_printf("\n");
  }

  if (part->cell_face != NULL) {
    PDM_printf("\nContent of cell_face\n");
    for(int i = 0; i < part->cell_face_idx[part->n_cell]; i++) {
      PDM_printf(" %d ", part->cell_face[i]);
      if (i % (part->cell_face_idx)[1] == (part->cell_face_idx)[1] - 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }

  if (part->cell_ln_to_gn != NULL) {
    PDM_printf("\nContent of cell_ln_to_gn\n");
    for(int i = 0; i < part->n_cell; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->cell_ln_to_gn[i]);
    }
    PDM_printf("\n");
  }

  if (part->cell_tag != NULL) {
    PDM_printf("\nContent of cell_tag\n");
    for(int i = 0; i < part->n_cell; i++) {
      PDM_printf(" %d ", part->cell_tag[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_cell != NULL) {
    PDM_printf("\nContent of face_cell\n");
    for(int i = 0; i < 2 * part->n_face; i++) {
      PDM_printf(" %d ", part->face_cell[i]);
      if (i % 2 == 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }

  if (part->face_vtx_idx != NULL) {
    PDM_printf("\nContent of face_vtx_idx\n");
    for(int i = 0; i < part->n_face + 1; i++) {
      PDM_printf(" %d ", part->face_vtx_idx[i]);
    }
    PDM_printf("\n");
  }

  if (part->gface_vtx != NULL) {
    PDM_printf("\nContent of gface_vtx\n");
    for(int i = 0; i < part->face_vtx_idx[part->n_face]; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->gface_vtx[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_vtx != NULL) {
    PDM_printf("\nContent of face_vtx\n");
    for(int i = 0; i < part->face_vtx_idx[part->n_face]; i++) {
      PDM_printf(" %d ", part->face_vtx[i]);
      if (i % (part->face_vtx_idx)[1] == (part->face_vtx_idx)[1] - 1) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }

  if (part->face_ln_to_gn != NULL) {
    PDM_printf("\nContent of face_ln_to_gn\n");
    for(int i = 0; i < part->n_face; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->face_ln_to_gn[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_tag != NULL) {
    PDM_printf("\nContent of face_tag\n");
    for(int i = 0; i < part->n_face; i++) {
      PDM_printf(" %d ", (part->face_tag)[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_part_bound_part_idx != NULL) {
    PDM_printf("\nContent of face_part_bound_part_idx\n");
    for(int i = 0; i < n_part + 1; i++) {
      PDM_printf(" %d ", part->face_part_bound_part_idx[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_part_bound_proc_idx != NULL) {
    PDM_printf("\nContent of face_part_bound_proc_idx\n");
    for(int i = 0; i < n_total_part + 1; i++) {
      PDM_printf(" %d ", part->face_part_bound_proc_idx[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_part_bound != NULL) {
    PDM_printf("\nContent of face_part_bound\n");
    for(int i = 0; i < 4 * part->n_face_part_bound; i++) {
      PDM_printf(" %d ", part->face_part_bound[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_group_idx != NULL) {
    PDM_printf("\nContent of face_group_idx\n");
    for(int i = 0; i < n_face_group + 1; i++) {
      PDM_printf(" %d ", part->face_group_idx[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_group != NULL) {
    PDM_printf("\nContent of face_group\n");
    for(int i = 0; i < part->face_group_idx[n_face_group]; i++) {
      PDM_printf(" %d ", part->face_group[i]);
    }
    PDM_printf("\n");
  }

  if (part->face_group_ln_to_gn != NULL) {
    PDM_printf("\nContent of face_group_ln_to_gn\n");
    for(int i = 0; i < part->face_group_idx[n_face_group]; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->face_group_ln_to_gn[i]);
    }
    PDM_printf("\n");
  }

  if (part->vtx != NULL) {
    PDM_printf("\nContent of vtx\n");
    for(int i = 0; i < 3 * part->n_vtx; i++) {
      PDM_printf(" %.1f ", part->vtx[i]);
      if (i % 3 == 2) {
        PDM_printf("|");
      }
    }
    PDM_printf("\n");
  }

  if (part->vtx_ln_to_gn != NULL) {
    PDM_printf("\nContent of vtx_ln_to_gn\n");
    for(int i = 0; i < part->n_vtx; i++) {
      PDM_printf(" "PDM_FMT_G_NUM" ", part->vtx_ln_to_gn[i]);
    }
    PDM_printf("\n");
  }


  if (part->vtx_tag != NULL) {
    PDM_printf("\nContent of vtx_tag\n");
    for(int i = 0; i < part->n_vtx; i++) {
      PDM_printf(" %d ", part->vtx_tag[i]);
    }
    PDM_printf("\n");
  }

}

/**
 *
 * \brief Displays all the arrays of a coarse partition of type _coarse_part_t
 *
 * \param [in]  n_part        Number of partitions to define on this process
 * \param [in]  n_total_part       Total number of partitions
 * \param [in]  n_face_group   Number of boundaries
 *
 */

static void
_coarse_part_display
(
  _coarse_part_t *coarse_part,
  int n_part,
  int n_total_part,
  int n_face_group
)
{
  if (coarse_part == NULL) {
    PDM_printf("Incorrect coarse part to display\n");
    return;
  }

  _part_display(coarse_part->part, n_part, n_total_part, n_face_group);

  if (coarse_part->coarse_cell_cell_idx != NULL) {
    PDM_printf("\nContent of coarse_cell_cell_idx\n");
    for(int i = 0; i < coarse_part->part->n_cell + 1; i++) {
      PDM_printf(" %d ", coarse_part->coarse_cell_cell_idx[i]);
    }
    PDM_printf("\n");
  }

  if (coarse_part->coarse_cell_cell != NULL) {
    PDM_printf("\nContent of coarse_cell_cell\n");
    for(int i = 0; i < coarse_part->coarse_cell_cell_idx[coarse_part->part->n_cell]; i++) {
      PDM_printf(" %d ", coarse_part->coarse_cell_cell[i]);
    }
    PDM_printf("\n");
  }

  if (coarse_part->coarse_face_group_to_fine_face_group != NULL) {
    PDM_printf("\nContent of coarse_face_group_to_fine_face_group\n");
    for(int i = 0; i < coarse_part->part->face_group_idx[n_face_group]; i++) {
      PDM_printf(" %d ", coarse_part->coarse_face_group_to_fine_face_group[i]);
    }
    PDM_printf("\n");
  }

  if (coarse_part->coarse_face_to_fine_face != NULL) {
    PDM_printf("\nContent of coarse_face_to_fine_face\n");
    for(int i = 0; i < coarse_part->part->n_face; i++) {
      PDM_printf(" %d ", coarse_part->coarse_face_to_fine_face[i]);
    }
    PDM_printf("\n");
  }

  if (coarse_part->coarse_vtx_to_fine_vtx != NULL) {
    PDM_printf("\nContent of coarse_vtx_to_fine_vtx\n");
    for(int i = 0; i < coarse_part->part->n_vtx; i++) {
      PDM_printf(" %d ", coarse_part->coarse_vtx_to_fine_vtx[i]);
    }
    PDM_printf("\n");
  }
}

/**
 *
 * \brief Free partition
 *
 * \param [in]   part      partition
 *
 */

static void
_part_free
(
 _part_t *part
)
{
  if (part == NULL) {
    PDM_printf("Incorrect part to free");
    return;
  }

  if (part->cell_face_idx != NULL)
    free(part->cell_face_idx);
  part->cell_face_idx = NULL;

  if (part->gcell_face != NULL)
    free(part->gcell_face);
  part->gcell_face = NULL;

  if (part->cell_face != NULL)
    free(part->cell_face);
  part->cell_face = NULL;

  if (part->cell_ln_to_gn != NULL)
    free(part->cell_ln_to_gn);
  part->cell_ln_to_gn = NULL;

  if (part->cell_tag != NULL)
    free(part->cell_tag);
  part->cell_tag = NULL;

  if (part->face_cell != NULL)
    free(part->face_cell);
  part->face_cell = NULL;

  if (part->face_vtx_idx != NULL)
    free(part->face_vtx_idx);
  part->face_vtx_idx = NULL;

  if (part->gface_vtx != NULL)
    free(part->gface_vtx);
  part->gface_vtx = NULL;

  if (part->face_vtx != NULL)
    free(part->face_vtx);
  part->face_vtx = NULL;

  if (part->face_ln_to_gn != NULL)
    free(part->face_ln_to_gn);
  part->face_ln_to_gn = NULL;

  if (part->face_tag != NULL)
    free(part->face_tag);
  part->face_tag = NULL;

  if (part->face_part_bound_proc_idx != NULL)
    free(part->face_part_bound_proc_idx);
  part->face_part_bound_proc_idx = NULL;

  if (part->face_part_bound_part_idx != NULL)
    free(part->face_part_bound_part_idx);
  part->face_part_bound_part_idx = NULL;

  if (part->face_part_bound != NULL)
    free(part->face_part_bound);
  part->face_part_bound = NULL;

  if (part->face_group_idx != NULL)
    free(part->face_group_idx);
  part->face_group_idx = NULL;

  if (part->face_group != NULL)
    free(part->face_group);
  part->face_group = NULL;

  if (part->face_group_ln_to_gn != NULL)
    free(part->face_group_ln_to_gn);
  part->face_group_ln_to_gn = NULL;

  if (part->vtx != NULL)
    free(part->vtx);
  part->vtx = NULL;

  if (part->vtx_ln_to_gn != NULL)
    free(part->vtx_ln_to_gn);
  part->vtx_ln_to_gn = NULL;

  if (part->vtx_tag != NULL)
    free(part->vtx_tag);
  part->vtx_tag = NULL;

  free(part);
}

/**
 *
 * \brief Free coarse partition
 *
 * \param [in]   coarse_part     coarse partition
 *
 */

static void
_coarse_part_free
(
 _coarse_part_t *coarse_part
)
{
  _part_free(coarse_part->part);

  if (coarse_part->specific_data != NULL) {
    free (coarse_part->specific_data);
  }

  if (coarse_part->coarse_cell_cell != NULL)
    free(coarse_part->coarse_cell_cell);
  coarse_part->coarse_cell_cell = NULL;

  if (coarse_part->coarse_cell_cell_idx != NULL)
    free(coarse_part->coarse_cell_cell_idx);
  coarse_part->coarse_cell_cell_idx = NULL;

  if (coarse_part->coarse_face_group_to_fine_face_group != NULL)
    free(coarse_part->coarse_face_group_to_fine_face_group);
  coarse_part->coarse_face_group_to_fine_face_group = NULL;

  if (coarse_part->coarse_face_to_fine_face != NULL)
    free(coarse_part->coarse_face_to_fine_face);
  coarse_part->coarse_face_to_fine_face = NULL;

  if (coarse_part->coarse_vtx_to_fine_vtx != NULL)
    free(coarse_part->coarse_vtx_to_fine_vtx);
  coarse_part->coarse_vtx_to_fine_vtx = NULL;

  free(coarse_part);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Return an initialized coarse mesh object
 *
 * \param [out]  cmId              Coarse mesh identifier
 *
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   n_part             Number of partitions
 * \param [in]   n_total_part            Total number of partitions
 * \param [in]   n_face_group        Total number of groups
 * \param [in]   have_cell_tag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_face_tag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtx_tag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cell_weight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_face_weight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_face_group    Presence des tableaux de groupes de faces
 */

void
PDM_part_coarse_mesh_create
(
 int                *cmId,
 PDM_MPI_Comm        comm,
 const char*         method,
 const char         *renum_cell_method,
 const char         *renum_face_method,
 const int           n_property_cell,
 const int          *renum_properties_cell,
 const int           n_property_face,
 const int          *renum_properties_face,
 const int           n_part,
 const int           n_total_part,
 const int           n_face_group,
 const int           have_cell_tag,
 const int           have_face_tag,
 const int           have_vtx_tag,
 const int           have_cell_weight,
 const int           have_face_weight,
 const int           have_face_group
)
{
  if (_cm == NULL) {
    _cm = PDM_Handles_create (4);
  }

  _coarse_mesh_t *cm  = _coarse_mesh_create (comm,
                                             method,
                                             renum_cell_method,
                                             renum_face_method,
                                             n_property_cell,
                                             renum_properties_cell,
                                             n_property_face,
                                             renum_properties_face,
                                             n_part,
                                             n_total_part,
                                             n_face_group,
                                             have_cell_tag,
                                             have_face_tag,
                                             have_vtx_tag,
                                             have_cell_weight,
                                             have_face_weight,
                                             have_face_group);

  *cmId = PDM_Handles_store (_cm, cm);
}

void
PROCF (pdm_part_coarse_mesh_create_cf, PDM_PART_COARSE_MESH_CREATE_CF)
(
 int                *cmId,
 PDM_MPI_Fint       *fcomm,
 const char         *method,
 const int          *l_method,
 const char         *renum_cell_method,
 const int          *l_renum_cell_method,
 const char         *renum_face_method,
 const int          *l_renum_face_method,
 const int          *n_property_cell,
 const int          *renum_properties_cell,
 const int          *n_property_face,
 const int          *renum_properties_face,
 const int          *n_part,
 const int          *n_total_part,
 const int          *n_face_group,
 const int          *have_cell_tag,
 const int          *have_face_tag,
 const int          *have_vtx_tag,
 const int          *have_cell_weight,
 const int          *have_face_weight,
 const int          *have_face_group
)
{

  PDM_MPI_Comm comm = PDM_MPI_Comm_f2c (*fcomm);

  char *_method = PDM_fortran_to_c_string (method, *l_method);

  char *_renum_cell_method =
      PDM_fortran_to_c_string(renum_cell_method, *l_renum_cell_method);

  char *_renum_face_method =
      PDM_fortran_to_c_string(renum_face_method, *l_renum_face_method);

  PDM_part_coarse_mesh_create (cmId,
                               comm,
                               _method,
                               _renum_cell_method,
                               _renum_face_method,
                               *n_property_cell,
                               renum_properties_cell,
                               *n_property_face,
                               renum_properties_face,
                               *n_part,
                               *n_total_part,
                               *n_face_group,
                               *have_cell_tag,
                               *have_face_tag,
                               *have_vtx_tag,
                               *have_cell_weight,
                               *have_face_weight,
                               *have_face_group);

  free (_method);

}

/**
 *
 * \brief Build a coarse mesh
 *
 * \param [in]  cmId               Coarse mesh identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  nCoarseCell        Number of cells in the coarse grid
 * \param [in]  n_cell              Number of cells
 * \param [in]  n_face              Number of faces
 * \param [in]  n_face_part_bound     Number of partitioning boundary faces
 * \param [in]  n_vtx               Number of vertices
 * \param [in]  n_face_group         Number of face groups
 * \param [in]  cell_face_idx        Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [in]  cell_face           Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                             numbering : 1 to n)
 * \param [in]  cell_tag            Cell tag (size = n_cell)
 * \param [in]  cell_ln_to_gn         Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [in]  cellWeight         Cell weight (size = n_cell)
 * \param [in]  faceWeight         Face weight (size = n_face)
 * \param [in]  face_cell           Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [in]  face_vtx_idx         Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [in]  face_vtx            Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [in]  face_tag            Face tag (size = n_face)
 * \param [in]  face_ln_to_gn         Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [in]  vtxCoord           Vertex coordinates (size = 3 * nVertex)
 * \param [in]  vtx_tag             Vertex tag (size = nVertex)
 * \param [in]  vtx_ln_to_gn          Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [in]  face_group_idx       Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [in]  face_group          faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 */

void
PDM_part_coarse_mesh_input
(
 int                 cmId,
 int                 i_part,
 const int           nCoarseCellWanted,
 const int           n_cell,
 const int           n_face,
 const int           n_vtx,
 const int           n_face_group, //FIXME: Argument a eliminer : Information deja donnee
 const int           n_face_part_bound,
 const int          *cell_face_idx,
 const int          *cell_face,
 const int          *cell_tag,
 const int          *cellWeight,
 const int          *faceWeight,
 const PDM_g_num_t  *cell_ln_to_gn,
 const int          *face_cell,
 const int          *face_vtx_idx,
 const int          *face_vtx,
 const int          *face_tag,
 const PDM_g_num_t  *face_ln_to_gn,
 const double       *vtxCoord,
 const int          *vtx_tag,
 const PDM_g_num_t  *vtx_ln_to_gn,
 const int          *face_group_idx,
 const int          *face_group,
 const PDM_g_num_t  *face_group_ln_to_gn,
 const int          *face_part_bound_proc_idx,
 const int          *face_part_bound_part_idx,
 const int          *face_part_bound
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  _coarse_grid_mesh_input (cm,
                           i_part,
                           nCoarseCellWanted,
                           n_cell,
                           n_face,
                           n_vtx,
                           n_face_group,
                           n_face_part_bound,
                           cell_face_idx,
                           cell_face,
                           cell_tag,
                           cellWeight,
                           faceWeight,
                           cell_ln_to_gn,
                           face_cell,
                           face_vtx_idx,
                           face_vtx,
                           face_tag,
                           face_ln_to_gn,
                           vtxCoord,
                           vtx_tag,
                           vtx_ln_to_gn,
                           face_group_idx,
                           face_group,
                           face_group_ln_to_gn,
                           face_part_bound_proc_idx,
                           face_part_bound_part_idx,
                           face_part_bound);
}

void
PROCF (pdm_part_coarse_mesh_input, PDM_PART_COARSE_MESH_INPUT)
(
 int                *cmId,
 int                *i_part,
 const int          *nCoarseCellWanted,
 const int          *n_cell,
 const int          *n_face,
 const int          *n_vtx,
 const int          *n_face_group,
 const int          *n_face_part_bound,
 const int          *cell_face_idx,
 const int          *cell_face,
 const int          *have_cell_tag,
 const int          *cell_tag,
 const int          *have_cell_weight,
 const int          *cellWeight,
 const int          *have_face_weight,
 const int          *faceWeight,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          *face_cell,
 const int          *face_vtx_idx,
 const int          *face_vtx,
 const int          *have_face_tag,
 const int          *face_tag,
 const PDM_g_num_t *face_ln_to_gn,
 const double       *vtxCoord,
 const int          *have_vtx_tag,
 const int          *vtx_tag,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int          *have_face_group,
 const int          *face_group_idx,
 const int          *face_group,
 const PDM_g_num_t *face_group_ln_to_gn,
 const int          *face_part_bound_proc_idx,
 const int          *face_part_bound_part_idx,
 const int          *face_part_bound
)
{

  int *_cell_tag = (int *) cell_tag;
  int *_face_tag = (int *) face_tag;
  int *_vtx_tag = (int *) vtx_tag;
  int *_cellWeight = (int *) cellWeight;
  int *_faceWeight = (int *) faceWeight;
  int *_faceGroupIdx = (int *) face_group_idx;
  int *_faceGroup = (int *) face_group;
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2312)
#endif
  PDM_g_num_t *_faceGroupLNToGN = (PDM_g_num_t *) face_group_ln_to_gn;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

  if (*have_cell_tag == 0) {
    _cell_tag = NULL;
  }

  if (*have_face_tag == 0) {
    _face_tag = NULL;
  }

  if (*have_vtx_tag == 0) {
    _vtx_tag = NULL;
  }

  if (*have_cell_weight == 0) {
    _cellWeight = NULL;
  }

  if (*have_face_weight == 0) {
    _faceWeight = NULL;
  }

  if (*have_face_group == 0) {
    _faceGroupIdx = NULL;
    _faceGroup = NULL;
    _faceGroupLNToGN = NULL;
  }

  PDM_part_coarse_mesh_input(*cmId,
                             *i_part,
                             *nCoarseCellWanted,
                             *n_cell,
                             *n_face,
                             *n_vtx,
                             *n_face_group,
                             *n_face_part_bound,
                             cell_face_idx,
                             cell_face,
                             _cell_tag,
                             _cellWeight,
                             _faceWeight,
                             cell_ln_to_gn,
                             face_cell,
                             face_vtx_idx,
                             face_vtx,
                             _face_tag,
                             face_ln_to_gn,
                             vtxCoord,
                             _vtx_tag,
                             vtx_ln_to_gn,
                             _faceGroupIdx,
                             _faceGroup,
                             _faceGroupLNToGN,
                             face_part_bound_proc_idx,
                             face_part_bound_part_idx,
                             face_part_bound);
}

/**
 *
 * \brief Updates all the arrays dealing with MPI exchanges
 *
 * \param [in] cmId               Coarse mesh identifier
 */

void
PDM_part_coarse_mesh_compute
(
                                int cmId
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  /* First step : Manage independently coarse grid generation */

  for (int i_part = 0; i_part < cm->n_part; i_part++) {
    _coarse_grid_compute(cm, i_part);
  }

  /* Second step : Manage MPI */
  int itime = 13;

  //    PDM_part_coarse_mesh_display(cmId);
  PDM_timer_resume(cm->timer);

  _build_coarsecell_ln_to_gn(cm);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);

  _build_face_ln_to_gn(cm);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);

  _build_vtx_ln_to_gn(cm);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);

  _build_faceGroupLNToGN(cm);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);
  itime += 1;

  PDM_timer_resume(cm->timer);

  // PDM_printf(" ------------------------------------------- \n");
  // for (int i = 0; i < part->face_group_idx[part->n_face_group]; i++) {
  //   PDM_printf("part->face_group[%i] = %i  \n", i, part->face_group[i]);
  //   int iFace = part->face_group[i];
  //   PDM_printf("part->face_cell = %i/%i  \n", i, part->face_cell[2*iFace],part->face_cell[2*iFace+1] );
  // }
  /* Renumbering */

  // Demander a Eric :
  // Il va manquer le reordering des tableau specific au multigrille ?
  // coarse_cell_cell, coarse_cell_cell_idx
  // coarse_cell_cell, coarse_vtx_to_fine_vtx
  printf("Renumbering Coarse mesh \n");

  for (int i_part = 0; i_part < cm->n_part; i_part++) {

    /* Cell renumbering */
    PDM_part_renum_cell (        &cm->part_res[i_part]->part,
                                 1,
                                 cm->renum_cell_method,
                         (void*) cm->renum_properties_cell);

    printf("PDM_part_renum_connectivities \n");
    if(cm->part_res[i_part]->part->new_to_old_order_cell != NULL)
    {
      /* Verbose */
      if( 0 == 1){
        printf("PDM_part_renum_connectivities end %i \n", cm->part_res[i_part]->part->n_cell);
        for (int i = 0; i < cm->part_res[i_part]->part->n_cell; i++){
          printf("part->new_to_old_order_cell[%i] = %i\n ", i,cm->part_res[i_part]->part->new_to_old_order_cell[i] );
        }
      }

      /* Renum CoarseGrid connectivity */
      PDM_part_renum_connectivities(cm->part_res[i_part]->part->n_cell,
                                    cm->part_res[i_part]->part->new_to_old_order_cell,
                                    cm->part_res[i_part]->coarse_cell_cell_idx,
                                    cm->part_res[i_part]->coarse_cell_cell);

      /* Si agglomeration method = 3 il faut changer les tableaux */
      // _coarse_part_t *part_res = cm->part_res[i_part];
      // _part_aniso_agglo_data_t *part_aniso_agglo_data  = (_part_aniso_agglo_data_t *) part_res->specific_data;



    }

    /* Face renumbering */
    PDM_part_renum_face (        &cm->part_res[i_part]->part,
                                 1,
                                 cm->renum_face_method,
                         (void*) cm->renum_properties_face);

    if(cm->part_res[i_part]->part->new_to_old_order_face != NULL)
    {
      /* Verbose */
      if( 0 == 1){
        printf("PDM_part_renum_connectivities end %i \n", cm->part_res[i_part]->part->n_face);
        for (int i = 0; i < cm->part_res[i_part]->part->n_face; i++){
          printf("part->new_to_old_order_face[%i] = %i\n ", i,cm->part_res[i_part]->part->new_to_old_order_face[i] );
        }
      }

      // printf("PDM_order_array \n");

      /* Renum CoarseGrid connectivity */
      PDM_order_array(cm->part_res[i_part]->part->n_face,
                      sizeof(PDM_l_num_t),
                      cm->part_res[i_part]->part->new_to_old_order_face,
                      cm->part_res[i_part]->coarse_face_to_fine_face);
      // printf("PDM_order_array end\n");

    }

    /* Compute renumbering to specific agglomeration method */
    if(cm->specific_func != NULL){
      (*cm->specific_func)(cm->part_res[i_part]);
    }


  }

  printf("Renumbering Coarse mesh end ... \n");

  _build_facePartBound(cm);

  PDM_timer_hang_on(cm->timer);
  cm->times_elapsed[itime] = PDM_timer_elapsed(cm->timer);
  cm->times_cpu[itime]     = PDM_timer_cpu(cm->timer);
  cm->times_cpu_u[itime]   = PDM_timer_cpu_user(cm->timer);
  cm->times_cpu_s[itime]   = PDM_timer_cpu_sys(cm->timer);

  cm->times_elapsed[0]     = cm->times_elapsed[itime];
  cm->times_cpu[0]         = cm->times_cpu[itime];
  cm->times_cpu_u[0]       = cm->times_cpu_u[itime];
  cm->times_cpu_s[0]       = cm->times_cpu_s[itime];

  for (int i = itime; i > 1; i--) {
    cm->times_elapsed[i] -= cm->times_elapsed[i-1];
    cm->times_cpu[i]     -= cm->times_cpu[i-1];
    cm->times_cpu_u[i]   -= cm->times_cpu_u[i-1];
    cm->times_cpu_s[i]   -= cm->times_cpu_s[i-1];
  }
}

void
PROCF (pdm_part_coarse_mesh_compute, PDM_PART_COARSE_MESH_COMPUTE)
(
 int *cmId
)
{
  PDM_part_coarse_mesh_compute(*cmId);
}

/**
 *
 * \brief Return a coarse mesh partition dimensions
 *
 * \param [in]   ppart_id            Coarse mesh identifier
 * \param [in]   i_part              Current partition
 *
 * \param [out]  n_cell              Number of cells
 * \param [out]  n_face              Number of faces
 * \param [out]  n_face_part_bound     Number of partitioning boundary faces
 * \param [out]  n_vtx               Number of vertices
 * \param [out]  n_proc              Number of processus
 * \param [out]  n_total_part             Number of partitions
 * \param [out]  scell_face          Size of cell-face connectivity
 * \param [out]  sface_vtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound     Size of face_part_bound array
 * \param [out]  sface_group         Size of face_group array
 * \param [out]  sCoarseCellToFineCell  Size of coarseCellToFineCell array
 *
 */

void
PDM_part_coarse_mesh_part_dim_get
(
 const int      cmId,
 int            i_part,
 int           *n_cell,
 int           *n_face,
 int           *n_face_part_bound,
 int           *n_vtx,
 int           *n_proc,
 int           *n_total_part,
 int           *n_face_group,
 int           *scell_face,
 int           *sface_vtx,
 int           *sface_group,
 int           *sCoarseCellToFineCell
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  _coarse_part_t *part_res = NULL;

  int numProcs;
  PDM_MPI_Comm_size(cm->comm, &numProcs);

  if (i_part < cm->n_part) {
    part_res = cm->part_res[i_part];
  }

  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_dim_get error : unknown partition\n");
    exit(1);
  }

  *n_face_group = cm->n_face_group;

  *n_cell = part_res->part->n_cell;
  *n_face = part_res->part->n_face;
  *n_vtx  = part_res->part->n_vtx;

  *n_face_part_bound  = part_res->part->n_face_part_bound;
  *n_proc           = numProcs;
  *n_total_part          = cm->n_total_part;
  *scell_face       = part_res->part->cell_face_idx[*n_cell];
  *sface_vtx        = part_res->part->face_vtx_idx[*n_face];
  *sCoarseCellToFineCell = part_res->coarse_cell_cell_idx[*n_cell];
  *sface_group      = 0;
  if (cm->n_face_group > 0 && part_res->part->face_group_idx != NULL) {
    *sface_group    = part_res->part->face_group_idx[cm->n_face_group];
  }
}

void
PROCF (pdm_part_coarse_mesh_part_dim_get, PDM_PART_COARSE_MESH_PART_DIM_GET)
(
 const int     *cmId,
 int           *i_part,
 int           *n_cell,
 int           *n_face,
 int           *n_face_part_bound,
 int           *n_vtx,
 int           *n_proc,
 int           *n_total_part,
 int           *n_face_group,
 int           *scell_face,
 int           *sface_vtx,
 int           *sface_group,
 int           *sCoarseCellToFineCell
)
{
  PDM_part_coarse_mesh_part_dim_get(*cmId,
                                    *i_part,
                                    n_cell,
                                    n_face,
                                    n_face_part_bound,
                                    n_vtx,
                                    n_proc,
                                    n_total_part,
                                    n_face_group,
                                    scell_face,
                                    sface_vtx,
                                    sface_group,
                                    sCoarseCellToFineCell);
}

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   cmId               Coarse mesh identifier
 * \param [in]   i_part              Current partition
 *
 * \param [out]  cell_face_idx        Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [out]  cell_face           Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                             numbering : 1 to n)
 * \param [out]  cell_tag            Cell tag (size = n_cell)
 * \param [out]  cell_ln_to_gn         Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [out]  cellInitCellIdx    Array of indexes of the connected partitions (size : nCoarseCell + 1)
 * \param [out]  cellInitCell       Partitioning array (size : cellInitCellIdx[nCoarseCell])
 *
 * \param [out]  face_cell           Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [out]  face_vtx_idx         Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [out]  face_vtx            Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [out]  face_tag            Face tag (size = n_face)
 * \param [out]  face_ln_to_gn         Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [out]  faceInitFace       Coarse face - fine face connectivity (size = nCoarseFace)
 *
 * \param [out]  vtxCoord           Vertex coordinates (size = 3 * n_vtx)
 * \param [out]  vtx_tag             Vertex tag (size = n_vtx)
 * \param [out]  vtx_ln_to_gn          Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [out]  vtxInitVtx         Coarse vertex - fine vertex connectivity (size = nCoarseVtx)
 *
 * \param [out]  face_group_idx       Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [out]  face_group          faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [out]  face_group_ln_to_gn    Faces global numbering for each group
 *                                  (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 *
 * \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 * \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 * \param [out]  face_part_bound      Partitioning boundary faces (size = 4 * n_face_part_bound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number
 *                                          in the connected partition (numbering :1 to n)
 *
 */

void
PDM_part_coarse_mesh_part_get
(
 const int    cmId,
 const int    i_part,
 int          **cell_face_idx,
 int          **cell_face,
 int          **cell_tag,
 PDM_g_num_t **cell_ln_to_gn,
 int          **cellInitCellIdx,
 int          **cellInitCell,
 int          **face_cell,
 int          **face_vtx_idx,
 int          **face_vtx,
 int          **face_tag,
 PDM_g_num_t **face_ln_to_gn,
 int          **faceGroupInitFaceGroup,
 int          **faceInitFace,
 double       **vtxCoord,
 int          **vtx_tag,
 PDM_g_num_t **vtx_ln_to_gn,
 int          **vtxInitVtx,
 int          **face_group_idx,
 int          **face_group,
 PDM_g_num_t **face_group_ln_to_gn,
 int          **face_part_bound_proc_idx,
 int          **face_part_bound_part_idx,
 int          **face_part_bound
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  _coarse_part_t *part_res = NULL;

  if (i_part < cm->n_part) {
    part_res = cm->part_res[i_part];
  }

  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_get error : unknown partition\n");
    exit(1);
  }

  *cellInitCellIdx        = part_res->coarse_cell_cell_idx;
  *cellInitCell           = part_res->coarse_cell_cell;
  *faceGroupInitFaceGroup = part_res->coarse_face_group_to_fine_face_group;
  *faceInitFace           = part_res->coarse_face_to_fine_face;
  *vtxInitVtx             = part_res->coarse_vtx_to_fine_vtx;

  *cell_face_idx          = part_res->part->cell_face_idx;
  *cell_face             = part_res->part->cell_face;
  *cell_tag              = part_res->part->cell_tag;
  *cell_ln_to_gn           = part_res->part->cell_ln_to_gn;
  *face_cell             = part_res->part->face_cell;
  *face_vtx_idx           = part_res->part->face_vtx_idx;
  *face_vtx              = part_res->part->face_vtx;
  *face_tag              = part_res->part->face_tag;
  *face_ln_to_gn           = part_res->part->face_ln_to_gn;
  *vtxCoord             = part_res->part->vtx;
  *vtx_tag               = part_res->part->vtx_tag;
  *vtx_ln_to_gn            = part_res->part->vtx_ln_to_gn;
  *face_group_idx         = part_res->part->face_group_idx;
  *face_group            = part_res->part->face_group;
  *face_group_ln_to_gn      = part_res->part->face_group_ln_to_gn;
  *face_part_bound_proc_idx = part_res->part->face_part_bound_proc_idx;
  *face_part_bound_part_idx = part_res->part->face_part_bound_part_idx;
  *face_part_bound        = part_res->part->face_part_bound;

}


void
PROCF (pdm_part_coarse_mesh_part_get, PDM_PART_COARSE_MESH_PART_GET)
(
 int          *cmId,
 int          *i_part,
 int          *cell_face_idx,
 int          *cell_face,
 int          *cell_tag,
 PDM_g_num_t *cell_ln_to_gn,
 int          *cellInitCellIdx,
 int          *cellInitCell,
 int          *face_cell,
 int          *face_vtx_idx,
 int          *face_vtx,
 int          *face_tag,
 PDM_g_num_t *face_ln_to_gn,
 int          *faceGroupInitFaceGroup,
 int          *faceInitFace,
 double       *vtxCoord,
 int          *vtx_tag,
 PDM_g_num_t *vtx_ln_to_gn,
 int          *vtxInitVtx,
 int          *face_group_idx,
 int          *face_group,
 PDM_g_num_t *face_group_ln_to_gn,
 int          *face_part_bound_proc_idx,
 int          *face_part_bound_part_idx,
 int          *face_part_bound
)
{
  _coarse_mesh_t * cm = _get_from_id (*cmId);

  int numProcs;
  PDM_MPI_Comm_size(cm->comm, &numProcs);

  _coarse_part_t *part_res = NULL;

  if (*i_part < cm->n_part) {
    part_res = cm->part_res[*i_part];
  }

  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_get error : unknown partition\n");
    exit(1);
  }

  for (int i = 0; i < part_res->part->n_cell + 1; i++)
    cell_face_idx[i] = part_res->part->cell_face_idx[i];

  for (int i = 0; i < part_res->part->cell_face_idx[part_res->part->n_cell]; i++)
    cell_face[i] = part_res->part->cell_face[i];

  for (int i = 0; i < part_res->part->n_cell; i++){
    cell_ln_to_gn[i] = part_res->part->cell_ln_to_gn[i];

    if (part_res->part->cell_tag != NULL)
        cell_tag[i] = part_res->part->cell_tag[i];
  }

  for (int i = 0; i < part_res->part->n_cell + 1; i++)
    cellInitCellIdx[i] = part_res->coarse_cell_cell_idx[i];

  for (int i = 0; i < part_res->coarse_cell_cell_idx[part_res->part->n_cell]; i++)
    cellInitCell[i] = part_res->coarse_cell_cell[i];

  for (int i = 0; i < 2 * part_res->part->n_face; i++){
    face_cell[i] = part_res->part->face_cell[i];
  }

  for (int i = 0; i < part_res->part->n_face + 1; i++)
    face_vtx_idx[i] = part_res->part->face_vtx_idx[i];

  for (int i = 0; i < part_res->part->face_vtx_idx[part_res->part->n_face]; i++)
    face_vtx[i] = part_res->part->face_vtx[i];

  for (int i = 0; i < part_res->part->n_face; i++){
    face_ln_to_gn[i] = part_res->part->face_ln_to_gn[i];

    if (part_res->part->face_tag != NULL)
      face_tag[i] = part_res->part->face_tag[i];
  }

  if (part_res->part->face_group_idx != NULL) {
    for (int i = 0; i < part_res->part->face_group_idx[cm->n_face_group]; i++) {
      faceGroupInitFaceGroup[i] = part_res->coarse_face_group_to_fine_face_group[i];
    }
  }

  for (int i = 0; i < part_res->part->n_face; i++)
    faceInitFace[i] = part_res->coarse_face_to_fine_face[i];

  for (int i = 0; i < 3 * part_res->part->n_vtx; i++){
    vtxCoord[i] = part_res->part->vtx[i];
  }

  for (int i = 0; i < part_res->part->n_vtx; i++){
    vtx_ln_to_gn[i] = part_res->part->vtx_ln_to_gn[i];

    if (part_res->part->vtx_tag != NULL)
      vtx_tag[i] = part_res->part->vtx_tag[i];

  }

  for (int i = 0; i < part_res->part->n_vtx; i++)
    vtxInitVtx[i] = part_res->coarse_vtx_to_fine_vtx[i];

  if (part_res->part->face_group_idx != NULL) {
    for (int i = 0; i < cm->n_face_group + 1; i++) {
      face_group_idx[i] = part_res->part->face_group_idx[i];
    }

    for (int i = 0; i < part_res->part->face_group_idx[cm->n_face_group]; i++) {
      face_group[i]       = part_res->part->face_group[i];
      face_group_ln_to_gn[i] = part_res->part->face_group_ln_to_gn[i];
    }
  }

  for (int i = 0; i < 4 * part_res->part->n_face_part_bound; i++)
    face_part_bound[i] = part_res->part->face_part_bound[i];

  for (int i = 0; i < numProcs + 1; i++)
    face_part_bound_proc_idx[i] = part_res->part->face_part_bound_proc_idx[i];

  for (int i = 0; i < cm->n_total_part + 1; i++)
    face_part_bound_part_idx[i] = part_res->part->face_part_bound_part_idx[i];

}

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   ppart_id            ppart identifier
 * \param [in]   i_part              Current partition
 * \param [out]  cellColor          Cell tag (size = n_cell)
 * \param [out]  faceColor          Face tag (size = n_face)

 */

void PDM_part_coarse_color_get
(
 const int   cmId,
 const int   i_part,
       int **cellColor,
       int **faceColor,
       int **threadColor,
       int **hyperPlaneColor
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  _coarse_part_t *part_res = NULL;

  if (i_part < cm->n_part) {
    part_res = cm->part_res[i_part];
  }

  if (part_res == NULL) {
    PDM_printf("PDM_part_coarse_mesh_part_get error : unknown partition\n");
    exit(1);
  }

  *cellColor       = part_res->part->cellColor;
  *faceColor       = part_res->part->faceColor;
  *threadColor     = part_res->part->threadColor;
  *hyperPlaneColor = part_res->part->hyperPlaneColor;
}

/**
 *
 * \brief Free coarse mesh
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

void
PDM_part_coarse_mesh_free
(
 const int    cmId
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  for (int i = 0; i < cm->n_part; i++) {
    free(cm->part_ini[i]);
    _coarse_part_free(cm->part_res[i]);
    cm->part_ini[i] = NULL;
    cm->part_res[i] = NULL;
  }

  if (cm->specific_data != NULL) {
    free (cm->specific_data);
  }

  free(cm->part_ini);
  free(cm->part_res);

  cm->part_ini = NULL;
  cm->part_res = NULL;

  PDM_timer_free(cm->timer);
  cm->timer = NULL;

  free(cm);

  PDM_Handles_handle_free (_cm, cmId, PDM_FALSE);

  const int n_cm = PDM_Handles_n_get (_cm);

  if (n_cm == 0) {
    _cm = PDM_Handles_free (_cm);
  }

}

void
PROCF (pdm_part_coarse_mesh_free, PDM_PART_COARSE_MESH_FREE)
(
  int *cmId
)
{
  PDM_part_coarse_mesh_free(*cmId);
}

/**
 *
 * \brief Return times
 *
 * \param [in]   cmId        coarse mesh identifier
 * \param [out]  elapsed     elapsed times (size = 18)
 * \param [out]  cpu         cpu times (size = 18)
 * \param [out]  cpu_user    user cpu times (size = 18)
 * \param [out]  cpu_sys     system cpu times (size = 18)
 *
 */

void PDM_part_coarse_mesh_time_get
(
 int       cmId,
 double  **elapsed,
 double  **cpu,
 double  **cpu_user,
 double  **cpu_sys
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  *elapsed  = cm->times_elapsed;
  *cpu      = cm->times_cpu;
  *cpu_user = cm->times_cpu_u;
  *cpu_sys  = cm->times_cpu_s;
}

void
PROCF (pdm_part_coarse_mesh_time_get, PDM_PART_COARSE_MESH_TIME_GET)
(
 int      *cmId,
 double   *elapsed,
 double   *cpu,
 double   *cpu_user,
 double   *cpu_sys
 )
{
  _coarse_mesh_t * cm = _get_from_id (*cmId);

  for (int i = 0; i < 18; i++) {
    elapsed[i]  = cm->times_elapsed[i];
    cpu[i]      = cm->times_cpu[i];
    cpu_user[i] = cm->times_cpu_u[i];
    cpu_sys[i]  = cm->times_cpu_s[i];
  }
}


/**
 *
 * \brief Displays all the arrays of a coarse mesh
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

void
PDM_part_coarse_mesh_display
(
 const int    cmId
)
{
  _coarse_mesh_t * cm = _get_from_id (cmId);

  //Display all the elements of the structure (not part of the other structures)

  PDM_printf("\n");
  PDM_printf("Value of n_part : %d \n", cm->n_part);
  PDM_printf("Value of comm : %d \n", cm->comm);
  PDM_printf("Value of method : %d \n", cm->method);
  PDM_printf("Value of n_total_part : %d \n", cm->n_total_part);
  PDM_printf("Value of n_face_group : %d \n", cm->n_face_group);

  for (int i = 0; i < cm->n_part; i++) {
    PDM_printf("\n=============================================\n");
    PDM_printf("Valeur de i : %d \n", i);
    PDM_printf("\n=============================================\n");

    PDM_printf("\n----------Affichage de part_ini-----------------\n");
    _part_display(cm->part_ini[i],cm->n_part,cm->n_total_part,cm->n_face_group);

    PDM_printf("\n----------Affichage de part_res-----------------\n");
    _coarse_part_display(cm->part_res[i], cm->n_part,cm->n_total_part,cm->n_face_group);
    PDM_printf("\n-----------------------------------------\n");
  }
}

void
PROCF (pdm_part_coarse_mesh_display, PDM_PART_COARSE_MESH_DISPLAY)
(
  int *cmId
)
{
  PDM_part_coarse_mesh_display (*cmId);
}


/**
 *
 * \brief Add a new coarse mesh method
 *
 * \param [in]      name          Mesh entity to renumber
 * \param [in]      fct           Function
 *
 */

int
PDM_coarse_mesh_method_add
(
 const char                 *name,     /*!< Name          */
 PDM_coarse_mesh_fct_t       fct       /*!< Function      */
)
{
  if (_coarse_mesh_methods == NULL) {
      PDM_coarse_mesh_method_load_local();
  }

  _coarse_mesh_method_t *method_ptr = malloc (sizeof(_coarse_mesh_method_t));

  int idx = PDM_Handles_store  (_coarse_mesh_methods, method_ptr);

  method_ptr->name = malloc (sizeof(char) * (strlen(name) + 1));
  strcpy (method_ptr->name, name);

  method_ptr->fct = fct;

  return idx;
}


/**
 *
 * \brief Get index of a coarse mesh method from it's name
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PROCF (pdm_coarse_mesh_method_idx_get_cf, PDM_COARSE_MESH_METHOD_IDX_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  char *_name = PDM_fortran_to_c_string (name, *l_name);

  *idx = PDM_coarse_mesh_method_idx_get (_name);

  free (_name);

}

int
PDM_coarse_mesh_method_idx_get
(
const char *name
)
{
  if (_coarse_mesh_methods == NULL) {
    PDM_coarse_mesh_method_load_local();
  }
  int idx = -1;

  if (_coarse_mesh_methods != NULL) {
    int n_methods = PDM_Handles_n_get (_coarse_mesh_methods);
    const int *index =  PDM_Handles_idx_get (_coarse_mesh_methods);

    for (int i = 0; i < n_methods; i++) {
      _coarse_mesh_method_t *method_ptr =
              (_coarse_mesh_method_t *) PDM_Handles_get (_coarse_mesh_methods, index[i]);
      if (!strcmp(method_ptr->name, name)) {
        idx = index[i];
        break;
      }
    }
  }
  return idx;
}


/**
 *
 * \brief Get name of a coarse mesh method from it's index
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PROCF (pdm_coarse_mesh_method_name_get_cf, PDM_COARSE_MESH_METHOD_NAME_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  const char *_name = PDM_coarse_mesh_method_name_get (*idx);

  const int _l_name = strlen(_name);

  *l_name = PDM_MAX (_l_name, PDM_MAX_CHAR_LENGTH);

  strncpy (name, _name, *l_name);
}

char *
PDM_coarse_mesh_method_name_get
(
const int id
)
{
  if (_coarse_mesh_methods == NULL) {
    PDM_coarse_mesh_method_load_local();
  }

  int n_methods = PDM_Handles_n_get (_coarse_mesh_methods);

  if (id >= n_methods) {
    return NULL;
  }

  const int *index =  PDM_Handles_idx_get (_coarse_mesh_methods);

  _coarse_mesh_method_t *method_ptr =
            (_coarse_mesh_method_t *) PDM_Handles_get (_coarse_mesh_methods, index[id]);

  return method_ptr->name;
}


/**
 *
 * \brief Get the number of coarse mesh method
 *
 * \return Number of methods
 *
 */

void
PROCF (pdm_coarse_mesh_method_n_get, PDM_COARSE_MESH_METHOD_N_GET)
(
 int  *n_method
 )
{
  *n_method = PDM_coarse_mesh_method_n_get ();
}

int
PDM_coarse_mesh_method_n_get
(
void
)
{
  if (_coarse_mesh_methods == NULL) {
    PDM_coarse_mesh_method_load_local();
  }

  return PDM_Handles_n_get (_coarse_mesh_methods);

}

/**
 *
 * \brief Purge coarse mesh methods catalog
 *
 */

void
PDM_coarse_mesh_method_purge
(
void
)
{
  if (_coarse_mesh_methods != NULL) {

    const int *index =  PDM_Handles_idx_get (_coarse_mesh_methods);
    int n_methods = PDM_Handles_n_get (_coarse_mesh_methods);

    while (n_methods > 0) {
      int idx = index[0];
      _coarse_mesh_method_t *method_ptr =
              (_coarse_mesh_method_t *) PDM_Handles_get (_coarse_mesh_methods, idx);
      free (method_ptr->name);
      PDM_Handles_handle_free (_coarse_mesh_methods, idx, PDM_TRUE);
      n_methods = PDM_Handles_n_get (_coarse_mesh_methods);
    }

    _coarse_mesh_methods = PDM_Handles_free (_coarse_mesh_methods);

  }
}

/**
 *
 * \brief Load local coarse mesh methods
 *
 */

void
PDM_coarse_mesh_method_load_local
(
void
)
{
  if (_coarse_mesh_methods == NULL)  {

    const int n_default_methods = 2;
    _coarse_mesh_methods = PDM_Handles_create (n_default_methods);

    PDM_coarse_mesh_method_add ("PDM_COARSE_MESH_SCOTCH",
                             _coarse_from_scotch);
    PDM_coarse_mesh_method_add ("PDM_COARSE_MESH_METIS",
                             _coarse_from_metis);

  }

}



/**
 *
 * \brief Return coarse mesh object from its identifier
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

_coarse_mesh_t *
PDM_part_coarse_mesh_get_from_id
(
 int  cmId
 )
{
  return _get_from_id (cmId);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

