/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_distant_neighbor.h"
#include "pdm_logging.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_order.h"
#include "pdm_part_extension.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static inline
int
_is_same_triplet
(
int iproc1, int ipart1, int ielt1,
int iproc2, int ipart2, int ielt2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        return 1;
      }
    }
  }
  return 0;
}

static
void
_create_cell_cell_graph
(
  PDM_part_extension_t *part_ext,
  PDM_extend_type_t     extend_type
)
{
  assert(extend_type == PDM_EXTEND_FROM_FACE);

  part_ext->cell_cell_idx = malloc( part_ext->n_domain * sizeof(int ***));
  part_ext->cell_cell     = malloc( part_ext->n_domain * sizeof(int ***));

  if(extend_type == PDM_EXTEND_FROM_FACE ){

    /* In order to call generic fonction we tranform face_cell with idx */
    int ***face_cell_idx = malloc( part_ext->n_domain * sizeof(int ***));
    int ***face_cell     = malloc( part_ext->n_domain * sizeof(int ***));

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      int  *pn_face        = (int * ) malloc( part_ext->n_part[i_domain] * sizeof(int  ));
      int  *pn_cell        = (int * ) malloc( part_ext->n_part[i_domain] * sizeof(int  ));
      int **pface_cell     = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      int **pcell_face     = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      int **pcell_face_idx = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_face       [i_part] = part_ext->parts[i_domain][i_part].n_face;
        pn_cell       [i_part] = part_ext->parts[i_domain][i_part].n_cell;
        pface_cell    [i_part] = part_ext->parts[i_domain][i_part].face_cell;
        pcell_face_idx[i_part] = part_ext->parts[i_domain][i_part].cell_face_idx;
        pcell_face    [i_part] = part_ext->parts[i_domain][i_part].cell_face;
      }

      PDM_part_connectivity_to_connectity_idx(part_ext->n_part[i_domain],
                                              pn_face,
                                              pface_cell,
                                              &face_cell_idx[i_domain],
                                              &face_cell[i_domain]);

      part_ext->cell_cell_idx[i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      part_ext->cell_cell    [i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        // PDM_log_trace_array_int(face_cell_idx[i_domain][i_part], pn_face[i_part], "face_cell_idx::");
        // PDM_log_trace_array_int(face_cell    [i_domain][i_part], face_cell_idx[i_domain][i_part][pn_face[i_part]], "face_cell::");
        // printf("[%i] face_cell_idx -> \n ", i_part);
        // for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
        //   printf(" (%i) = ", i_face);
        //   for(int idx_cell = face_cell_idx[i_domain][i_part][i_face]; idx_cell < face_cell_idx[i_domain][i_part][i_face+1]; ++idx_cell) {
        //     printf(" %i", face_cell[i_domain][i_part][idx_cell]-1);
        //   }
        //   printf("\n");
        // }

        // On fait égalemnt le cell_cell
        // printf("PDM_combine_connectivity[%i] -> %i %i \n", i_part, pn_cell[i_part], pn_face[i_part]);
        PDM_combine_connectivity(pn_cell[i_part],
                                 pcell_face_idx[i_part],
                                 pcell_face[i_part],
                                 face_cell_idx[i_domain][i_part],
                                 face_cell[i_domain][i_part],
                                 &part_ext->cell_cell_idx[i_domain][i_part],
                                 &part_ext->cell_cell[i_domain][i_part]);

        /*
         * Setup shortcut and free useless memory
         */
        part_ext->entity_cell_idx[i_part+shift_part] = face_cell_idx[i_domain][i_part];
        part_ext->entity_cell    [i_part+shift_part] = face_cell    [i_domain][i_part];
        part_ext->entity_cell_n  [i_part+shift_part] = (int * ) malloc( pn_face[i_part] * sizeof(int));
        for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
          int n_connected = part_ext->entity_cell_idx[i_part+shift_part][i_face+1] - part_ext->entity_cell_idx[i_part+shift_part][i_face];
          part_ext->entity_cell_n[i_part+shift_part][i_face] = n_connected;
        }
      }

      free(pn_face);
      free(pn_cell);
      free(pface_cell);
      free(pcell_face);
      free(pcell_face_idx);

      free(face_cell_idx[i_domain]);
      free(face_cell[i_domain]);

      shift_part += part_ext->n_part[i_domain];
    }

    free(face_cell_idx);
    free(face_cell);
  } else {
    abort();
  }
}

static
void
_create_cell_graph_comm
(
  PDM_part_extension_t *part_ext
)
{
  /*
   * The idea is to use the first communicator graph to exchange directly the entity_cell
   *     -> We deduced first from this graph the data necessary for distant_neigbor
   *     -> We apply exchange with the connectivity entitiy_cell
   *        The information in graph comm is usefull also for directly adressing the results
   */
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  assert(part_ext->n_domain == 1);
  /* Si multidomain on fait un shift et tt roule */

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }
  // printf(" n_tot_all_domain      = %i \n", n_tot_all_domain);
  // printf(" n_part_loc_all_domain = %i \n", n_part_loc_all_domain);

  /* We flat all partition */
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->n_entity_bound == NULL);

  part_ext->neighbor_idx   = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->neighbor_desc  = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->n_entity_bound = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  part_ext->entity_cell_opp_idx = (int  **) malloc( n_part_loc_all_domain * sizeof(int  *) );

  part_ext->dist_neighbor_cell_n    = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_idx  = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_desc = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );

  // Begin with exchange by the connectivity the cell opposite number
  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    int n_part_total = part_ext->n_tot_part_by_domain[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_part_entity_bound_tot = part_ext->parts[i_domain][i_part].face_part_bound_part_idx[n_part_total];

      // printf(" i_part+shift_part = %i \n", i_part+shift_part);
      // printf(" n_part_entity_bound_tot = %i \n", n_part_entity_bound_tot);
      int* _entity_part_bound = NULL;
      if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        // part_ext->n_entity_bound[i_part+shift_part] = n_part_entity_bound_tot;
        part_ext->n_entity_bound[i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
        _entity_part_bound = part_ext->parts[i_domain][i_part].face_part_bound;
      } else {
        abort();
      }

      part_ext->neighbor_idx       [i_part+shift_part] = (int *) malloc(    (part_ext->n_entity_bound[i_part+shift_part]+1) * sizeof(int) );

      int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      part_ext->n_cell[i_part+shift_part] = n_cell;

      int* _neighbor_idx  = part_ext->neighbor_idx [i_part+shift_part];

      int* _neighbor_n = (int * ) malloc( part_ext->n_entity_bound[i_part+shift_part] * sizeof(int));
      for(int i_entity = 0; i_entity < part_ext->n_entity_bound[i_part+shift_part]; ++i_entity) {
        _neighbor_n[i_entity] = 0;
      }

      /* Just copy to remove the 4 indices to 3 indices */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        _neighbor_n[i_entity] += 1;
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      for(int i_entity = 0; i_entity < part_ext->n_entity_bound[i_part+shift_part]; ++i_entity) {
        _neighbor_idx[i_entity+1] = _neighbor_idx[i_entity] + _neighbor_n[i_entity];
        _neighbor_n[i_entity] = 0;
      }

      /* Ici il faut faire les raccords entre domaine ---> Count */
      part_ext->neighbor_desc[i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]] * sizeof(int) );
      int* _neighbor_desc = part_ext->neighbor_desc[i_part+shift_part];

      /* Second loop for fill */
      /* Just copy to remove the 4 indices to 3 indices */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
        _neighbor_desc[3*idx_write  ] = _entity_part_bound[4*idx_entity+1];   // i_proc_opp;
        _neighbor_desc[3*idx_write+1] = _entity_part_bound[4*idx_entity+2]-1; // i_part_opp
        _neighbor_desc[3*idx_write+2] = _entity_part_bound[4*idx_entity+3]-1; // i_entity_opp
      }

      free(_neighbor_n);

      /* Ici il faut faire les raccords entre domaine ---> Fill - Count _distant_neighbor_cell also */


    }

    shift_part += part_ext->n_part[i_domain];
  }


  /*
   * All partition of all domain is treated in the same time
   *    -> All partition are flatten along i_domain
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_entity_bound,
                                                           part_ext->neighbor_idx,
                                                           part_ext->neighbor_desc);


  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR,
                            -1,
                            part_ext->entity_cell_n,
                  (void **) part_ext->entity_cell,
                           &part_ext->entity_cell_opp_n,
                 (void ***)&part_ext->entity_cell_opp);
  // PDM_distant_neighbor_exch_int(dn,
  //                           sizeof(int),
  //                           PDM_STRIDE_VAR,
  //                           -1,
  //                           part_ext->entity_cell_n,
  //                 (void **) part_ext->entity_cell,
  //                          &part_ext->entity_cell_opp_n,
  //                (void ***)&part_ext->entity_cell_opp);

  /* Compute idx */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    int n_part_total = part_ext->n_tot_part_by_domain[i_domain];
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_part_entity_bound_tot = part_ext->parts[i_domain][i_part].face_part_bound_part_idx[n_part_total];

      int* _entity_part_bound = NULL;
      if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        _entity_part_bound = part_ext->parts[i_domain][i_part].face_part_bound;
      } else {
        abort();
      }

      int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      int* _neighbor_idx        = part_ext->neighbor_idx       [i_part+shift_part];
      int* _entity_cell_opp_n   = part_ext->entity_cell_opp_n  [i_part+shift_part];
      int* _entity_cell_opp     = part_ext->entity_cell_opp    [i_part+shift_part];

      int* _entity_cell_idx = part_ext->entity_cell_idx[i_part+shift_part];
      int* _entity_cell     = part_ext->entity_cell    [i_part+shift_part];

      int n_entity = part_ext->n_entity_bound[i_part+shift_part]; // n_face / n_edge / n_vtx
      int n_elmt   = _neighbor_idx[n_entity];

      part_ext->entity_cell_opp_idx   [i_part+shift_part] = (int *) malloc((n_elmt+1) * sizeof(int) );

      // int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      part_ext->dist_neighbor_cell_n  [i_part+shift_part] = (int *) malloc(  n_cell    * sizeof(int) );
      part_ext->dist_neighbor_cell_idx[i_part+shift_part] = (int *) malloc( (n_cell+1) * sizeof(int) );

      int* _entity_cell_opp_idx    = part_ext->entity_cell_opp_idx   [i_part+shift_part];
      int* _dist_neighbor_cell_n   = part_ext->dist_neighbor_cell_n  [i_part+shift_part];
      int* _dist_neighbor_cell_idx = part_ext->dist_neighbor_cell_idx[i_part+shift_part];

      // printf(" n_elmt = %i \n", n_elmt);
      _entity_cell_opp_idx[0] = 0;
      for(int i_elmt = 0; i_elmt < n_elmt; ++i_elmt) {
        _entity_cell_opp_idx[i_elmt+1] = _entity_cell_opp_idx[i_elmt] + _entity_cell_opp_n[i_elmt];
      }

      /* Init in order to count after - two pass to count = interior then border */
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _dist_neighbor_cell_n[i_cell] = 0;
      }

      part_ext->border_cell_list[i_part+shift_part] = malloc( n_elmt * sizeof(int));
      // for(int i = 0; i < n_elmt; ++i) {
      //   part_ext->border_cell_list[i_part+shift_part][i] = -1000;
      // }

      /* For each border count number of opposite */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        for(int idx_cell = _entity_cell_idx[i_entity]; idx_cell < _entity_cell_idx[i_entity+1]; ++idx_cell) {
          int i_cell = PDM_ABS(_entity_cell[idx_cell])-1;
          _dist_neighbor_cell_n[i_cell] += _entity_cell_opp_n[idx_entity];
        }
      }

      // PDM_log_trace_array_int(_dist_neighbor_cell_n , n_cell, "_dist_neighbor_cell_n::");
      // PDM_log_trace_array_int(_entity_cell_opp_n , n_part_entity_bound_tot, "_entity_cell_opp_n::");
      // PDM_log_trace_array_int(_entity_cell_opp , n_part_entity_bound_tot, "_entity_cell_opp::");

      /* Ici il faut faire les raccords entre domaine ---> Count */

      /* Compute idx */
      int idx_indic = 0;
      _dist_neighbor_cell_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _dist_neighbor_cell_idx[i_cell+1] = _dist_neighbor_cell_idx[i_cell] + _dist_neighbor_cell_n[i_cell];
        if(_dist_neighbor_cell_n[i_cell] > 0){
          part_ext->border_cell_list[i_part+shift_part][idx_indic++] = i_cell;
        }
      }
      /* Because cell can be associated twice */
      part_ext->n_cell_border[i_part+shift_part] = idx_indic;

      // printf(" idx_indic = %i\n", idx_indic);
      // PDM_log_trace_array_int(part_ext->border_cell_list[i_part+shift_part]  , idx_indic  , "border_cell_list::");

      /* Reset */
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _dist_neighbor_cell_n[i_cell] = 0;
      }

      /* Allocate */
      part_ext->dist_neighbor_cell_desc[i_part+shift_part] = (int * ) malloc( 3 * _dist_neighbor_cell_idx[n_cell] * sizeof(int));
      int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

      /* For each border count number of opposite */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity   = _entity_part_bound[4*idx_entity  ]-1;
        int i_proc_opp = _entity_part_bound[4*idx_entity+1];
        int i_part_opp = _entity_part_bound[4*idx_entity+2]-1;
        // int i_entity_opp = _entity_part_bound[4*idx_entity+3]-1;
        for(int idx_cell = _entity_cell_idx[i_entity]; idx_cell < _entity_cell_idx[i_entity+1]; ++idx_cell) {
          int i_cell = PDM_ABS(_entity_cell[idx_cell])-1;

          for(int idx_cell_opp = _entity_cell_opp_idx[idx_entity]; idx_cell_opp < _entity_cell_opp_idx[idx_entity+1]; ++idx_cell_opp) {
            int idx_write = _dist_neighbor_cell_idx[i_cell] + _dist_neighbor_cell_n[i_cell]++;
            int i_opp_cell = _entity_cell_opp[idx_cell_opp]-1; // Commence à zero
            // printf("[%i] - _entity_cell_opp[%i] = %i \n", i_part, idx_cell_opp,  _entity_cell_opp[idx_cell_opp]);
            _dist_neighbor_cell_desc[3*idx_write  ] = i_proc_opp;
            _dist_neighbor_cell_desc[3*idx_write+1] = i_part_opp;
            _dist_neighbor_cell_desc[3*idx_write+2] = i_opp_cell;
            // printf("[%i][%i] _dist_neighbor_cell[%i] = %i %i %i - i_entity = %i | i_entity_opp = %i \n", i_part, i_cell, idx_write, i_proc_opp, i_part_opp, i_opp_cell, i_entity, i_entity_opp);
          }
        }
      }
    }
    shift_part += part_ext->n_part[i_domain];
  }

  PDM_distant_neighbor_free(dn);

  /*
   * Panic verbose
   */
  if(0 == 1) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        int* _neighbor_idx = part_ext->neighbor_idx                [i_part+shift_part];
        int n_elmt         = _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]];
        int n_data         = part_ext->entity_cell_opp_idx[i_part+shift_part][n_elmt];

        PDM_log_trace_array_int(part_ext->entity_cell_opp_n  [i_part+shift_part], n_elmt  , "entity_cell_opp_n::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp_idx[i_part+shift_part], n_elmt+1, "entity_cell_opp_idx::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp    [i_part+shift_part], n_data  , "entity_cell_opp::");

        int n_cell = part_ext->parts[i_domain][i_part].n_cell;
        int* _dist_neighbor_cell_idx  = part_ext->dist_neighbor_cell_idx [i_part+shift_part];
        int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

        printf(" _dist_neighbor_cell_idx ---------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          for(int idx = _dist_neighbor_cell_idx[i_cell]; idx < _dist_neighbor_cell_idx[i_cell+1]; ++idx) {
            printf("(%i, %i) ", _dist_neighbor_cell_desc[3*idx+1], _dist_neighbor_cell_desc[3*idx+2]);
          }
          printf("\n");
        }
        printf(" _dist_neighbor_cell_idx ---------------- END \n");

        PDM_log_trace_array_int(_dist_neighbor_cell_idx , n_cell+1                       , "_dist_neighbor_cell_idx::");
        PDM_log_trace_array_int(_dist_neighbor_cell_desc, 3 * _dist_neighbor_cell_idx[n_cell], "_dist_neighbor_cell_desc::");
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }
  // exit(1);

  /*
   *   Pour la reconstruction des ghost cell
   *     -> On peut à mon avis faire une astuce pour faire 1 distant_neighbor avec
   *         "1 voisins = 1 rangées "
   *     -> Assez pratique car par partitione on va récupérer une liste de rangés également -_-
   *     -> Pb possible quand une cellule apparait dans une partition de rang 2 et dans une autre en rang 1
   */

}


static
void
_compute_dual_graph
(
  PDM_part_extension_t   *part_ext,
  PDM_distant_neighbor_t *dn,
  int                     i_depth
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  assert(i_depth > 0);
  int** prev_cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth-1];
  int** prev_cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth-1];
  int** prev_cell_cell_extended     = part_ext->cell_cell_extended    [i_depth-1];

  // int** next_cell_cell_extended_idx = NULL;
  int** next_cell_cell_extended_n   = NULL;
  int** next_cell_cell_extended     = NULL;

  /*
   * Exchange of the previous rank
   */
  PDM_distant_neighbor_exch(dn,
                            3 * sizeof(int),
                            PDM_STRIDE_VAR,
                            -1,
                            prev_cell_cell_extended_n,
                  (void **) prev_cell_cell_extended,
                           &next_cell_cell_extended_n,
                 (void ***)&next_cell_cell_extended);

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell        = part_ext->n_cell       [i_part+shift_part];
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int* _border_cell_cell_extended_n   = next_cell_cell_extended_n[i_part+shift_part];
      int* _border_cell_cell_extended     = next_cell_cell_extended  [i_part+shift_part];
      int* _border_cell_cell_extended_idx = (int * ) malloc( (n_cell_border+1) * sizeof(int) );


      if( 0 == 1) {
        printf("prev_cell_cell_extended :: --------------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          for(int idx = prev_cell_cell_extended_idx[i_part+shift_part][i_cell]; idx < prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1]; ++idx) {
            printf("(%i, %i) ", prev_cell_cell_extended[i_part+shift_part][3*idx+1], prev_cell_cell_extended[i_part+shift_part][3*idx+2]);
          }
          printf("\n");
        }
        printf("prev_cell_cell_extended :: --------------------- END \n");
      }

      _border_cell_cell_extended_idx[0] = 0;
      for(int i = 0; i < n_cell_border; ++i) {
        _border_cell_cell_extended_idx[i+1] = _border_cell_cell_extended_idx[i] + _border_cell_cell_extended_n[i];
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_border_cell_cell_extended_idx, n_cell_border+1, "next_border_cell_cell_extended_idx::");
        PDM_log_trace_array_int(_border_cell_cell_extended_n  , n_cell_border,   "next_border_cell_cell_extended_n::");
        PDM_log_trace_array_int(_border_cell_cell_extended    , 3 * _border_cell_cell_extended_idx[n_cell_border], "next_border_cell_cell_extended::");

        printf("_border_cell_cell_extended :: --------------------- \n");
        for(int i = 0; i < n_cell_border; ++i) {
          int i_cell = part_ext->border_cell_list[i_part+shift_part][i];
          printf("i_cell -> %i -->  ", i_cell);
          for(int idx = _border_cell_cell_extended_idx[i]; idx < _border_cell_cell_extended_idx[i+1]; ++idx) {
            printf("(%i, %i) ", _border_cell_cell_extended[3*idx+1], _border_cell_cell_extended[3*idx+2]);
          }
          printf("\n");
        }
        printf("_border_cell_cell_extended :: --------------------- END \n");

      }

      /* Now we have the opposite border_cell_cell_extended
       *   -> We need to reput all together
       *   -> The current cell has to be also extendted but caution !!!
       *      the extension depend also of border
       */

      /*
       * Generate tag to know is a local cell is a border or not
       */
      int* idx_border_cell = (int *) malloc( n_cell * sizeof(int));
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        idx_border_cell[i_cell] = -1;
      }

      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        idx_border_cell[i_cell] = idx_cell; /* Keep idx_cell for imediate loop */
      }

      /* Allocate current depth and shorcut */
      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = (int *) malloc( (n_cell + 1 ) * sizeof(int));
      part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part] = (int *) malloc( (n_cell     ) * sizeof(int));

      int* _pcell_cell_extended_n   = prev_cell_cell_extended_n  [i_part+shift_part];
      int* _pcell_cell_extended_idx = prev_cell_cell_extended_idx[i_part+shift_part];
      int* _pcell_cell_extended     = prev_cell_cell_extended    [i_part+shift_part];

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];
      int* _cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part];

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell] = 0;
        _cell_cell_extended_n  [i_cell] = 0;
      }
      _cell_cell_extended_idx[n_cell] = 0;


      /* Let's go - First pass to count */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];

        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior - we add the previous rank */
        _cell_cell_extended_n[i_cell] = _pcell_cell_extended_n[i_cell];

        /* From border */
        // printf("_border_cell_cell_extended_n[%i][%i] = %i\n", i_part, idx_cell, _border_cell_cell_extended_n[idx_cell]);
        assert(_border_cell_cell_extended_n[idx_cell] > 0);
        _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_cell];

        /* Now we have to extend the interior */
        for(int idx_neight = _pcell_cell_extended_idx[i_cell]; idx_neight < _pcell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _pcell_cell_extended[3*idx_neight  ];
          int i_part_neight = _pcell_cell_extended[3*idx_neight+1];
          /* We add stencil only if it's local */
          if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {
            int i_cell_neight = _pcell_cell_extended[3*idx_neight+2];

            /* From interior */
            _cell_cell_extended_n[i_cell] += _pcell_cell_extended_n[i_cell_neight];

            int idx_border_neight = idx_border_cell[i_cell_neight];
            if(idx_border_neight != -1) {
              // Il faut rajouter les voisins aussi
              _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_border_neight];
            }
          } /* End if same part and same proc */
        } /* End loop neighbor */
      } /* End loop border */

      /* Setup idx and reset */
      int max_neight = 0;
      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
        max_neight = PDM_MAX(max_neight, _cell_cell_extended_n[i_cell]);
        _cell_cell_extended_n[i_cell] = 0;
      }

      /* Allocate */
      part_ext->cell_cell_extended[i_depth][i_part+shift_part] = (int *) malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth][i_part+shift_part];

      /* Let's go - Second pass to fill */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];

        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior - we add the previous rank */
        for(int idx_neight = _pcell_cell_extended_idx[i_cell]; idx_neight < _pcell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = _pcell_cell_extended[3*idx_neight  ];
          _cell_cell_extended[3*idx_write+1] = _pcell_cell_extended[3*idx_neight+1];
          _cell_cell_extended[3*idx_write+2] = _pcell_cell_extended[3*idx_neight+2];
        }

        /* From border */
        assert(_border_cell_cell_extended_n[idx_cell] > 0);
        for(int idx_neight = _border_cell_cell_extended_idx[idx_cell]; idx_neight < _border_cell_cell_extended_idx[idx_cell+1]; ++idx_neight) {
          int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = _border_cell_cell_extended[3*idx_neight  ];
          _cell_cell_extended[3*idx_write+1] = _border_cell_cell_extended[3*idx_neight+1];
          _cell_cell_extended[3*idx_write+2] = _border_cell_cell_extended[3*idx_neight+2];
        }

        /* Now we have to extend the interior */
        for(int idx_neight = _pcell_cell_extended_idx[i_cell]; idx_neight < _pcell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _pcell_cell_extended[3*idx_neight  ];
          int i_part_neight = _pcell_cell_extended[3*idx_neight+1];
          /* We add stencil only if it's local */
          if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {
            int i_cell_neight = _pcell_cell_extended[3*idx_neight+2];

            /* From interior */
            for(int idx_neight2 = _pcell_cell_extended_idx[i_cell_neight]; idx_neight2 < _pcell_cell_extended_idx[i_cell_neight+1]; ++idx_neight2) {
              int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
              _cell_cell_extended[3*idx_write  ] = _pcell_cell_extended[3*idx_neight2  ];
              _cell_cell_extended[3*idx_write+1] = _pcell_cell_extended[3*idx_neight2+1];
              _cell_cell_extended[3*idx_write+2] = _pcell_cell_extended[3*idx_neight2+2];
            }

            int idx_border_neight = idx_border_cell[i_cell_neight];
            if(idx_border_neight != -1) {
              // Il faut rajouter les voisins aussi
              for(int idx_neight2 = _border_cell_cell_extended_idx[idx_border_neight]; idx_neight2 < _border_cell_cell_extended_idx[idx_border_neight+1]; ++idx_neight2) {
                int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                _cell_cell_extended[3*idx_write  ] = _border_cell_cell_extended[3*idx_neight2  ];
                _cell_cell_extended[3*idx_write+1] = _border_cell_cell_extended[3*idx_neight2+1];
                _cell_cell_extended[3*idx_write+2] = _border_cell_cell_extended[3*idx_neight2+2];
              }
            }
          } /* End if same part and same proc */

        } /* End loop neighbor */
      } /* End loop border */


      /* The _cell_cell_extended need to be sorted because many entry is duplicated */
      int* order = malloc( max_neight * sizeof(int));
      int* _ncell_cell_extended_n   = malloc( (n_cell    ) * sizeof(int));
      int* _ncell_cell_extended_idx = malloc( (n_cell + 1) * sizeof(int));
      int* _ncell_cell_extended     = malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      _ncell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

        int beg = _cell_cell_extended_idx[i_cell];
        int n_connect = _cell_cell_extended_idx[i_cell+1] - beg;
        assert(n_connect == _cell_cell_extended_n[i_cell]);

        PDM_order_lnum_s(&_cell_cell_extended[3*beg], 3, order, n_connect);

        // printf(" order[%i] = ", i_cell);
        // for(int i = 0; i < n_connect; ++i) {
        //   printf(" %i", order[i]);
        // }
        // printf("\n");

        _ncell_cell_extended_n  [i_cell  ] = 0;
        _ncell_cell_extended_idx[i_cell+1] = _ncell_cell_extended_idx[i_cell];

        // int idx_unique = -1;
        int last_proc  = -1;
        int last_part  = -1;
        int last_elmt  = -1;
        for(int i = 0; i < n_connect; ++i) {
          int old_order = order[i];
          int curr_proc = _cell_cell_extended[3*(beg+old_order)  ];
          int curr_part = _cell_cell_extended[3*(beg+old_order)+1];
          int curr_cell = _cell_cell_extended[3*(beg+old_order)+2];
          int is_same  = _is_same_triplet(last_proc, last_part, last_elmt,
                                          curr_proc, curr_part, curr_cell);

          if(is_same == 0){ // N'est pas le meme
            // idx_unique++;
            last_proc = curr_proc;
            last_part = curr_part;
            last_elmt = curr_cell;

            int beg_write = 3 * _ncell_cell_extended_idx[i_cell+1];
            // printf(" write in = %i  | beg_write = %i | idx_unique = %i\n", beg_write + idx_unique, beg_write, idx_unique);
            _ncell_cell_extended[beg_write  ] = curr_proc;
            _ncell_cell_extended[beg_write+1] = curr_part;
            _ncell_cell_extended[beg_write+2] = curr_cell;

            /* Increment the new counter */
            _ncell_cell_extended_idx[i_cell+1]++;
            _ncell_cell_extended_n  [i_cell  ]++;
          }
        }
      }

      /* Free old ptr and assign the sort one */
      free(part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part]);
      free(part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part]);
      free(part_ext->cell_cell_extended    [i_depth][i_part+shift_part]);

      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = _ncell_cell_extended_idx;
      part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part] = _ncell_cell_extended_n;
      part_ext->cell_cell_extended    [i_depth][i_part+shift_part] = _ncell_cell_extended;

      if(0 == 1) {
        PDM_log_trace_array_int(_ncell_cell_extended_idx, n_cell+1, "_ncell_cell_extended_idx:: ");
        PDM_log_trace_array_int(_ncell_cell_extended_n  , n_cell  , "_ncell_cell_extended_n:: ");
        PDM_log_trace_array_int(_ncell_cell_extended    , 3 * _ncell_cell_extended_idx[n_cell], "_ncell_cell_extended:: ");
      }

      /* Free */
      free(idx_border_cell);
      free(_border_cell_cell_extended_idx);
      free(order);

      /* Free allocated memory in distant neigbor exhange */
      free(next_cell_cell_extended_n[i_part+shift_part]);
      free(next_cell_cell_extended  [i_part+shift_part]);

    }
    shift_part += part_ext->n_part[i_domain];
  }

  free(next_cell_cell_extended_n);
  free(next_cell_cell_extended);

}



static
void
_compute_first_extended_cell_graph
(
 PDM_part_extension_t *part_ext
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int i_depth_cur = 0;

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell = part_ext->n_cell[i_part+shift_part];
      part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part] = (int *) malloc( (n_cell + 1 ) * sizeof(int));
      part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part] = (int *) malloc( (n_cell     ) * sizeof(int));

      int* _dist_neighbor_cell_n    = part_ext->dist_neighbor_cell_n   [i_part+shift_part];
      int* _dist_neighbor_cell_idx  = part_ext->dist_neighbor_cell_idx [i_part+shift_part];
      int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

      /* Uniquement besoin du graph sur les cellules de bords */
      // int n_cell_border = _dist_neighbor_cell_idx[n_cell];
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part];
      int* _cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part];

      int* _cell_cell_idx = part_ext->cell_cell_idx[i_depth_cur][i_part+shift_part];
      int* _cell_cell     = part_ext->cell_cell    [i_depth_cur][i_part+shift_part];

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell] = 0;
        _cell_cell_extended_n  [i_cell] = 0;
      }
      _cell_cell_extended_idx[n_cell] = 0;

      // printf(" n_cell_border = %i \n", n_cell_border);

      /* First pass to count */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior */
        _cell_cell_extended_n[i_cell] = _cell_cell_idx[i_cell+1] - _cell_cell_idx[i_cell];

        /* From border */
        assert(_dist_neighbor_cell_n[i_cell] > 0);
        _cell_cell_extended_n[i_cell] += _dist_neighbor_cell_n[i_cell];
      }

      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
        _cell_cell_extended_n[i_cell] = 0;
      }

      part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part] = (int *) malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part];

      /* Second pass to fill */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];

        /* From interior */
        for(int idx_neight = _cell_cell_idx[i_cell]; idx_neight < _cell_cell_idx[i_cell+1]; ++idx_neight ) {
          int i_cell_neight = _cell_cell[idx_neight];
          int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = i_rank;
          _cell_cell_extended[3*idx_write+1] = i_part;
          _cell_cell_extended[3*idx_write+2] = i_cell_neight-1;
        }

        /* From border */
        for(int idx_neight = _dist_neighbor_cell_idx[i_cell]; idx_neight < _dist_neighbor_cell_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _dist_neighbor_cell_desc[3*idx_neight  ];
          int i_part_neight = _dist_neighbor_cell_desc[3*idx_neight+1];
          int i_cell_neight = _dist_neighbor_cell_desc[3*idx_neight+2];
          int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = i_rank_neight;
          _cell_cell_extended[3*idx_write+1] = i_part_neight;
          _cell_cell_extended[3*idx_write+2] = i_cell_neight;
        }

        // printf("[%i] _cell_cell_extended_n[%i] = %i\n", i_part, i_cell, _cell_cell_extended_n[i_cell]);
      }

      if(1 == 1) {
        PDM_log_trace_array_int(_cell_cell_extended_n  , n_cell  , "_cell_cell_extended_n::");
        PDM_log_trace_array_int(_cell_cell_extended_idx, n_cell+1, "_cell_cell_extended_idx::");
        PDM_log_trace_array_int(_cell_cell_extended, 3 * _cell_cell_extended_idx[n_cell]  , "_cell_cell_extended::");
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }
}

static
void
_rebuild_faces
(
  PDM_part_extension_t *part_ext,
  int i_depth
)
{

  printf("_rebuild_faces \n");

  // Besoin
  // _cell_cell_extended + _cell_cell_extended_idx --> distant_neighbor
  // int* cell_face_idx = part_ext->cell_face_idx;
  // int* cell_face     = part_ext->cell_face;

  /* In order to get the face connectivity we need to prepare array of cell_face AND face_group */
  // _rebuild_part_extension();

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  if(0 == 1) {
    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        printf("cell_cell_extended[i_depth=%i, %i] :: --------------------- \n", i_depth, i_part+shift_part);
        int n_cell_border = part_ext->n_cell_border[i_part+shift_part];
        int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];
        int* _cell_cell_extended     = part_ext->cell_cell_extended    [i_depth][i_part+shift_part];
        for(int i = 0; i < n_cell_border; ++i) {
          int i_cell = part_ext->border_cell_list[i_part+shift_part][i];
          printf("i_cell -> %i -->  ", i_cell);
          for(int idx = _cell_cell_extended_idx[i]; idx < _cell_cell_extended_idx[i+1]; ++idx) {
            printf("(%i, %i) ", _cell_cell_extended[3*idx+1], _cell_cell_extended[3*idx+2]);
          }
          printf("\n");
        }
        printf("cell_cell_extended :: --------------------- END \n");
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_cell,
                                                           part_ext->cell_cell_extended_idx[i_depth],
                                                           part_ext->cell_cell_extended    [i_depth]);

  /* Prepare */
  int         **cell_face_n = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **gcell_face  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **cell_flags    = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int* cell_face_idx =  part_ext->parts[i_domain][i_part].cell_face_idx;
      int* cell_face     =  part_ext->parts[i_domain][i_part].cell_face;

      int n_cell      = part_ext->parts[i_domain][i_part].n_cell;
      // int n_face      = part_ext->parts[i_domain][i_part].n_face;
      int s_cell_face = cell_face_idx[n_cell];

      cell_face_n[i_part+shift_part] = (int         *) malloc( n_cell      * sizeof(int        ));
      gcell_face [i_part+shift_part] = (PDM_g_num_t *) malloc( s_cell_face * sizeof(PDM_g_num_t));

      PDM_g_num_t* face_ln_to_gn = part_ext->parts[i_domain][i_part].face_ln_to_gn;

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        cell_face_n[i_part+shift_part][i_cell] = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];
        for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {
          int sgn    = PDM_SIGN(cell_face[idx_face]);
          int i_face = PDM_ABS (cell_face[idx_face])-1;
          // gcell_face[i_part+shift_part][idx_face] = sgn * face_ln_to_gn[i_face];
          if(i_part == 0) {
            gcell_face[i_part+shift_part][idx_face] = (1+idx_face);
          } else {
            gcell_face[i_part+shift_part][idx_face] = -(1+idx_face);
          }
          gcell_face[i_part+shift_part][idx_face] = sgn * face_ln_to_gn[i_face];
          printf("gcell_face[%i][%i] = %i \n", i_part+shift_part, idx_face, i_part);
          // gcell_face[i_part+shift_part][idx_face] = i_part;
          // gcell_face[i_part+shift_part][idx_face] = i_cell;
        }
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Exchange */
  int         **border_gcell_face_n;
  PDM_g_num_t **border_gcell_face;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            -1,
                            cell_face_n,
                 (void **)  gcell_face,
                           &border_gcell_face_n,
                (void ***) &border_gcell_face);

  // PDM_distant_neighbor_exch_int(dn,
  //                           sizeof(PDM_g_num_t),
  //                           PDM_STRIDE_VAR,
  //                           -1,
  //                           cell_face_n,
  //                (void **)  gcell_face,
  //                          &border_gcell_face_n,
  //               (void ***) &border_gcell_face);


  /* Post treatment */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell        = part_ext->n_cell       [i_part+shift_part];
      int n_face        = part_ext->parts[i_domain][i_part].n_face;
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int         *_cell_cell_extended     = part_ext->cell_cell_extended    [i_depth][i_part+shift_part];
      int         *_cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];

      int         *_border_gcell_face_n   = border_gcell_face_n[i_part+shift_part];
      PDM_g_num_t *_border_gcell_face     = border_gcell_face[i_part+shift_part];
      // int         *_border_gcell_face_idx = (int * ) malloc( (n_cell_border+1) * sizeof(int) );
      int         *_border_gcell_face_idx = (int * ) malloc( (_cell_cell_extended_idx[n_cell]+1) * sizeof(int) );

      _border_gcell_face_idx[0] = 0;
      // for(int i = 0; i < n_cell_border; ++i) {
      // --> PAS BON
      // for(int i = 0; i < _cell_cell_extended_idx[n_cell_border]; ++i) {
      //   _border_gcell_face_idx[i+1] = _border_gcell_face_idx[i] + _border_gcell_face_n[i];
      // }
      int n_neight_tot = _cell_cell_extended_idx[n_cell];
      int s_tot = 0;
      for(int i = 0; i < n_cell_border; ++i) {
        printf("_border_gcell_face_idx[%i] -> \n", i);
        for(int idx_neight = _cell_cell_extended_idx[i]; idx_neight < _cell_cell_extended_idx[i+1]; ++idx_neight) {
          printf("(%i,%i) -> %i ", _cell_cell_extended[3*idx_neight+1],  _cell_cell_extended[3*idx_neight+2], _border_gcell_face_n[idx_neight]);
          s_tot += _border_gcell_face_n[idx_neight];
          _border_gcell_face_idx[idx_neight+1] =  _border_gcell_face_idx[idx_neight] + _border_gcell_face_n[idx_neight];
        }
        printf("\n");
      }
      printf(" s_tot = %i \n", s_tot);
      printf(" _cell_cell_extended_idx[n_cell] = %i \n", _cell_cell_extended_idx[n_cell]);

      if(1 == 1) {
        printf("n_cell_border = %i \n", n_cell_border);
        // PDM_log_trace_array_int (_border_gcell_face_idx, n_cell_border+1, "_border_gcell_face_idx::");
        // PDM_log_trace_array_int (_border_gcell_face_n  , n_cell_border,   "_border_gcell_face_n::");
        // Pas la bonne taille sur
        PDM_log_trace_array_int (_border_gcell_face_idx, n_neight_tot+1, "_border_gcell_face_idx::");
        PDM_log_trace_array_int (_border_gcell_face_n  , n_neight_tot  , "_border_gcell_face_n::");
        // PDM_log_trace_array_long(_border_gcell_face    , _border_gcell_face_idx[_cell_cell_extended_idx[n_cell_border]]  , "_border_gcell_face::");
        PDM_log_trace_array_long(_border_gcell_face, s_tot, "_border_gcell_face::");
      }

      /* Preparation du tableau des faces venant de l'exterieur */
      // PDM_g_num_t* _border_face_ln_to_gn = (PDM_g_num_t * ) malloc( _border_gcell_face_idx[n_cell_border] * sizeof(PDM_g_num_t));
      PDM_g_num_t* _border_face_ln_to_gn = (PDM_g_num_t * ) malloc( s_tot * sizeof(PDM_g_num_t));
      PDM_g_num_t* face_ln_to_gn = part_ext->parts[i_domain][i_part].face_ln_to_gn;

      PDM_g_num_t* _sorted_face_ln_to_gn = (PDM_g_num_t * ) malloc( n_face * sizeof(PDM_g_num_t));
      for(int i_face = 0; i_face < n_face; ++i_face ) {
        _sorted_face_ln_to_gn[i_face] = face_ln_to_gn[i_face];
      }
      PDM_sort_long(_sorted_face_ln_to_gn, 0, n_face-1);

      // for(int i = 0; i < _border_gcell_face_idx[n_cell_border]; ++i) {
      for(int i = 0; i < s_tot; ++i) {
        _border_face_ln_to_gn[i] = _border_gcell_face[i];
      }
      // int n_face_extended = PDM_inplace_unique_long(_border_face_ln_to_gn, 0, _border_gcell_face_idx[n_cell_border]-1);
      // int n_face_extended = PDM_inplace_unique_long(_border_face_ln_to_gn, 0, _border_gcell_face_idx[_cell_cell_extended_idx[n_cell_border]]-1);
      int n_face_extended = PDM_inplace_unique_long(_border_face_ln_to_gn, 0, s_tot-1);

      if(1 == 1) {
        PDM_log_trace_array_long(_border_face_ln_to_gn, n_face_extended, "_border_face_ln_to_gn::");
      }

      /* Pour chaque elements on chercher si il est dans les face_ln_to_gn */
      // for(int i_face = 0; i_face < n_face_extended; ++i_face) {
      for(int i_face = 0; i_face < s_tot; ++i_face) {
        PDM_g_num_t g_face = _border_face_ln_to_gn[i_face];
        int pos = PDM_binary_search_long(g_face, _sorted_face_ln_to_gn, n_face);
        printf(" [%i] found [%i] = %i\n", i_part+shift_part, i_face, pos);
      }

      free(_border_gcell_face_idx);
      free(_border_face_ln_to_gn);
      free(_sorted_face_ln_to_gn);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Pour les faces group on peut faire aussi le gnum location --> Marche pas en multidomain (ou il faut shifter )*/
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {


      free(cell_face_n[i_part+shift_part]);
      free(gcell_face[i_part+shift_part]);
      free(border_gcell_face_n[i_part+shift_part]);
      free(border_gcell_face[i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  PDM_distant_neighbor_free(dn);
  free(cell_face_n);
  free(gcell_face);
  free(border_gcell_face_n);
  free(border_gcell_face);
  printf("_rebuild_faces end \n");
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a part extension structure
 *
 * \param [in]   comm           Communicator
 *
 */

PDM_part_extension_t*
PDM_part_extension_create
(
 const int                n_domain,
 const int               *n_part,
       PDM_extend_type_t  extend_type,
       int                depth,
 const PDM_MPI_Comm       comm,
 const PDM_ownership_t    owner
)
{
  PDM_part_extension_t *part_ext = (PDM_part_extension_t *) malloc(sizeof(PDM_part_extension_t));

  part_ext->n_domain    = n_domain;
  part_ext->n_part      = n_part;
  part_ext->comm        = comm;
  part_ext->owner       = owner;
  part_ext->extend_type = extend_type;
  part_ext->depth       = depth;

  part_ext->parts = malloc(n_domain * sizeof(_part_t *));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    part_ext->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
  }

  part_ext->neighbor_idx     = NULL;
  part_ext->neighbor_desc    = NULL;
  part_ext->n_entity_bound   = NULL;
  part_ext->border_cell_list = NULL;

  part_ext->dist_neighbor_cell_n     = NULL;
  part_ext->dist_neighbor_cell_idx   = NULL;
  part_ext->dist_neighbor_cell_desc  = NULL;

  part_ext->n_tot_part_by_domain = NULL;

  part_ext->entity_cell_idx     = NULL;
  part_ext->entity_cell         = NULL;
  part_ext->entity_cell_opp_idx = NULL;
  part_ext->entity_cell_opp     = NULL;

  return part_ext;
}

void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                   n_cell,
  int                   n_face,
  int                   n_face_part_bound,
  int                   n_vtx,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn
)
{

  part_ext->parts[i_domain][i_part].n_cell            = n_cell;
  part_ext->parts[i_domain][i_part].n_face            = n_face;
  part_ext->parts[i_domain][i_part].n_face_part_bound = n_face_part_bound;
  part_ext->parts[i_domain][i_part].n_vtx             = n_vtx;

  part_ext->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  part_ext->parts[i_domain][i_part].cell_face     = cell_face;
  part_ext->parts[i_domain][i_part].face_cell     = face_cell;

  part_ext->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  part_ext->parts[i_domain][i_part].face_edge     = face_edge;

  part_ext->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  part_ext->parts[i_domain][i_part].face_vtx      = face_vtx;

  part_ext->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  part_ext->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  part_ext->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  part_ext->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  part_ext->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = face_part_bound_proc_idx;
  part_ext->parts[i_domain][i_part].face_part_bound_part_idx = face_part_bound_part_idx;
  part_ext->parts[i_domain][i_part].face_part_bound          = face_part_bound;

  part_ext->parts[i_domain][i_part].face_bound_idx      = face_bound_idx;
  part_ext->parts[i_domain][i_part].face_bound          = face_bound;
  part_ext->parts[i_domain][i_part].face_join_idx       = face_join_idx;
  part_ext->parts[i_domain][i_part].face_join           = face_join;
  part_ext->parts[i_domain][i_part].face_bound_ln_to_gn = face_group_ln_to_gn;
}


/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext
)
{
  /*
   *  A prevoir : reconstruire les "ghost" sans toute la topologie
   *              par exemple pour l'algébre linéaire
   */
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int depth = part_ext->depth;

  part_ext->n_tot_part_by_domain = (int *) malloc( part_ext->n_domain * sizeof(int));
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    part_ext->n_tot_part_by_domain[i_domain] = -1;
    int n_part_loc = part_ext->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &part_ext->n_tot_part_by_domain[i_domain], 1, PDM_MPI_INT, PDM_MPI_SUM, part_ext->comm);

    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  part_ext->entity_cell_n    = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell_idx  = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell      = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->n_cell           = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  part_ext->n_cell_border    = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  part_ext->border_cell_list = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

  // int shift_part = 0;
  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

  //   /* Donc face_cell ou vtx_cell */
  //   part_ext->entity_cell_idx[i_part+shift_part] = ;
  //   part_ext->entity_cell    [i_part+shift_part] = ;

  //   shift_part += part_ext->n_part[i_domain];
  // }


  assert(part_ext != NULL);
  printf(" PDM_part_extension_compute : depth = %i | extend_type = %i \n", depth, part_ext->extend_type);

  assert(part_ext->extend_type == PDM_EXTEND_FROM_FACE);

  _create_cell_cell_graph(part_ext, part_ext->extend_type);

  part_ext->cell_cell_extended_idx  = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended_n    = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended      = (int *** ) malloc( (depth + 1) * sizeof(int  **));

  for(int i_depth = 0; i_depth < depth+1; ++i_depth) {
    part_ext->cell_cell_extended_idx [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_extended_n   [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_extended     [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  }

  // TODO : vtx_cell
  _create_cell_graph_comm(part_ext);

  /*
   * Create for first level the proper graph
   */
  _compute_first_extended_cell_graph(part_ext);


  /*
   *   Step 3 : Compute the graph cell with triplet
   *      -> Now we have all the things ok to exchange direclty on cells
   *      -> In order to build the next depth of ghost cell
   *         we need to prepare an extended graph with triplet
   *      -> With the same distant neigbor we exchange for each depth this graph (containing the information of surrounding cells)
   */



  /*
   * Init the distant_neightbor from cell
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_cell,
                                                           part_ext->dist_neighbor_cell_idx,
                                                           part_ext->dist_neighbor_cell_desc);

  for(int i_depth = 1; i_depth < depth+1; ++i_depth) {

    /* Graph compute = local + distant */
    _compute_dual_graph(part_ext, dn, i_depth);

    /* Exchange for next rank */

    /* Post */

  }

  PDM_distant_neighbor_free(dn);

  /*
   * At this step we have for each level the opposite connected cell
   *   In order to setup all ghost cell we need to deduced all descending connectivity
   *       cell -> face -> edge -> vtx
   * Another step is to compute the ln_to_gn in mulitpart context -> (i_domain, ln_to_gn)
   */

  /*
   *  Implem du binary search bi-niveaux -> i_domain + ln_to_gn
   */
  printf(" ATTENTION CHANGEMENT DEPTH DEBUG \n");
  _rebuild_faces(part_ext, depth-1);
  printf(" ATTENTION CHANGEMENT DEPTH DEBUG \n");



  printf(" PDM_part_extension_compute end \n");
}


/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
)
{

  if(part_ext->n_tot_part_by_domain != NULL) {
    free(part_ext->n_tot_part_by_domain);
    part_ext->n_tot_part_by_domain = NULL;
  }

  if(part_ext->neighbor_idx != NULL) {
    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->neighbor_idx       [i_part+shift_part]);
        free(part_ext->neighbor_desc      [i_part+shift_part]);

        free(part_ext->dist_neighbor_cell_n   [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_idx [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_desc[i_part+shift_part]);

        free(part_ext->entity_cell_opp_idx[i_part+shift_part]);
        free(part_ext->entity_cell_opp_n  [i_part+shift_part]);
        free(part_ext->entity_cell_opp    [i_part+shift_part]);

        free(part_ext->border_cell_list    [i_part+shift_part]);

        if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
          free(part_ext->entity_cell_idx[i_part+shift_part]);
          free(part_ext->entity_cell_n  [i_part+shift_part]);
          free(part_ext->entity_cell    [i_part+shift_part]);
        }

        free(part_ext->cell_cell_idx[i_domain][i_part]);
        free(part_ext->cell_cell    [i_domain][i_part]);

      }
      free(part_ext->cell_cell_idx[i_domain]);
      free(part_ext->cell_cell    [i_domain]);

      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->n_entity_bound);
  }
  free(part_ext->neighbor_idx);
  free(part_ext->neighbor_desc);
  free(part_ext->n_cell);
  free(part_ext->n_cell_border);
  free(part_ext->border_cell_list);

  free(part_ext->cell_cell_idx);
  free(part_ext->cell_cell);


  for(int i_depth = 0; i_depth < part_ext->depth+1; ++i_depth) {

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_extended    [i_depth][i_part+shift_part]);
      }
      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->cell_cell_extended_idx[i_depth]);
    free(part_ext->cell_cell_extended_n  [i_depth]);
    free(part_ext->cell_cell_extended    [i_depth]);
  }

  free(part_ext->cell_cell_extended_idx);
  free(part_ext->cell_cell_extended_n);
  free(part_ext->cell_cell_extended);

  /* Only shortcut of user data */
  free(part_ext->entity_cell_idx    );
  free(part_ext->entity_cell        );
  free(part_ext->entity_cell_n);

  free(part_ext->dist_neighbor_cell_n   );
  free(part_ext->dist_neighbor_cell_idx );
  free(part_ext->dist_neighbor_cell_desc);

  /* Allocated by distant neightbor */
  free(part_ext->entity_cell_opp_idx);
  free(part_ext->entity_cell_opp_n  );
  free(part_ext->entity_cell_opp    );

  part_ext->n_part = NULL;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(part_ext->parts[i_domain]);
  }
  free(part_ext->parts);


  free(part_ext);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

