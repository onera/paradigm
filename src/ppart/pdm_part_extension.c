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

// static
// int
// _setup_neighbor
// (
//  int           i_proc_cur,
//  int           i_part_cur,
//  int           ntot_part,
//  int          *entity_part_bound_part_idx,
//  int          *entity_part_bound,
//  int          *entity_cell_idx,
//  int          *entity_cell,
//  int          *cell_list,
//  int          *graph_comm_cell,
//  PDM_g_num_t **cell_border_ln_to_gn,
//  int         **neighbor_idx,
//  int         **neighbor_desc
// )
// {
//   abort();
//   *neighbor_idx         = malloc(     ( entity_part_bound_part_idx[ntot_part] + 1 ) * sizeof(int        ));
//   *neighbor_desc        = malloc(  3 *  entity_part_bound_part_idx[ntot_part]       * sizeof(int        )); // i_proc, i_part, i_entity

//   int* _neighbor_idx  = *neighbor_idx;
//   int* _neighbor_desc = *neighbor_desc;

//   PDM_g_num_t* _cell_border_ln_to_gn = *cell_border_ln_to_gn;

//   int n_entity_bound = 0;
//   for(int i_part = 0; i_part < ntot_part; ++i_part ){
//     _neighbor_idx[0] = 0;
//     int idx_opp = 0;
//     for(int idx_entity = entity_part_bound_part_idx[i_part]; idx_entity < entity_part_bound_part_idx[i_part+1]; ++idx_entity) {

//       int i_entity     = entity_part_bound[4*idx_entity  ]-1;
//       int i_proc_opp   = entity_part_bound[4*idx_entity+1];
//       int i_part_opp   = entity_part_bound[4*idx_entity+2]-1;
//       // int i_entity_opp = entity_part_bound[4*idx_entity+3];

//       _neighbor_idx [idx_entity+1  ] = _neighbor_idx[idx_entity] + 1;
//       _neighbor_desc[3*idx_entity  ] = i_proc_opp;
//       _neighbor_desc[3*idx_entity+1] = i_part_opp;
//       _neighbor_desc[3*idx_entity+2] = idx_opp; // Car symétrique

//       /* Each bound can be connected to multiple cell */
//       int icell1 = entity_cell[2*i_entity]-1;

//       cell_list[idx_entity] = icell1;
//       graph_comm_cell[3*idx_entity  ] = i_proc_cur;
//       graph_comm_cell[3*idx_entity+1] = i_part_cur;
//       graph_comm_cell[3*idx_entity+2] = icell1;

//       idx_opp++;
//       n_entity_bound++;
//     }
//   }

//   return n_entity_bound;
// }


// static
// int
// _setup_neighbor_gen
// (
//  int           i_proc_cur,
//  int           i_part_cur,
//  int           ntot_part,
//  int          *entity_part_bound_part_idx,
//  int          *entity_part_bound,
//  int          *entity_cell_idx,
//  int          *entity_cell,
//  PDM_g_num_t  *cell_ln_to_gn,
//  int          *cell_list,
//  int         **graph_comm_cell_idx,
//  int         **graph_comm_cell,
//  int         **neighbor_idx,
//  int         **neighbor_desc
// )
// {

// }



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
  } else {
    abort();
  }
}


// static
// void
// _create_cell_graph_comm
// (
//   PDM_part_extension_t *part_ext
// )
// {

//   int i_rank;
//   int n_rank;
//   PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
//   PDM_MPI_Comm_size(part_ext->comm, &n_rank);

//   int**  n_face_part_bound = (int  **) malloc( part_ext->n_domain * sizeof(int  * ) );

//   // Begin with exchange by the connectivity the cell opposite number
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

//     int n_part_total = -1;
//     int n_part_loc = part_ext->n_part[i_domain];
//     PDM_MPI_Allreduce(&n_part_loc, &n_part_total, 1, PDM_MPI_INT, PDM_MPI_SUM, part_ext->comm);

//     part_ext->neighbor_idx [0][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *) );
//     part_ext->neighbor_desc[0][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *) );
//     n_face_part_bound[i_domain] = (int  *) malloc( part_ext->n_part[i_domain] * sizeof(int  ) );

//     part_ext->n_cell_per_bound[0][i_domain] = (int  *) malloc( part_ext->n_part[i_domain] * sizeof(int  ) );

//     printf(" n_part_total = %i \n", n_part_total);
//     PDM_g_num_t** cell_border_ln_to_gn = malloc( part_ext->n_part [i_domain] * sizeof(PDM_g_num_t * ));

//     part_ext->cell_list      [0][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
//     part_ext->graph_comm_cell[0][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

//       int n_part_face_bound_tot = part_ext->parts[i_domain][i_part].face_part_bound_part_idx[n_part_total];
//       part_ext->cell_list      [0][i_domain][i_part] = (int *) malloc(     n_part_face_bound_tot * sizeof(int));
//       part_ext->graph_comm_cell[0][i_domain][i_part] = (int *) malloc( 3 * n_part_face_bound_tot * sizeof(int));

//       // Pour chaque partition on hook le graph de comm
//       // n_face_part_bound[i_domain][i_part] = _setup_neighbor(i_rank,
//       //                                                       i_part,
//       //                                                       n_part_total,
//       //                                                       part_ext->parts[i_domain][i_part].face_part_bound_part_idx,
//       //                                                       part_ext->parts[i_domain][i_part].face_part_bound,
//       //                                                       part_ext->parts[i_domain][i_part].face_cell,
//       //                                                       part_ext->parts[i_domain][i_part].cell_ln_to_gn,
//       //                                                       part_ext->cell_list[0][i_domain][i_part],
//       //                                                       part_ext->graph_comm_cell[0][i_domain][i_part],
//       //                                                       &cell_border_ln_to_gn[i_part],
//       //                                                       &part_ext->neighbor_idx [0][i_domain][i_part],
//       //                                                       &part_ext->neighbor_desc[0][i_domain][i_part]);

//       // For the first depth it's the same
//       part_ext->n_cell_per_bound[0][i_domain][i_part] = n_face_part_bound[i_domain][i_part];
//       printf(" n_face_part_bound[%i][%i] = %i \n", i_domain, i_part, n_face_part_bound[i_domain][i_part]);
//     }

//     /*
//      * Create distant neighbor
//      */
//     PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
//                                                              part_ext->n_part [i_domain],
//                                                              n_face_part_bound[i_domain],
//                                                              part_ext->neighbor_idx [0][i_domain],
//                                                              part_ext->neighbor_desc[0][i_domain]);

//     /*
//      * Echange du numero :
//      *     - i_domain (useless now)
//      *     - cell_ln_to_gn
//      *     - cell_face (in gn )
//      *     - face_ln_to_gn
//      */
//     PDM_g_num_t **cell_ln_to_gn_opp;
//     PDM_distant_neighbor_exch(dn,
//                               sizeof(PDM_g_num_t),
//                               PDM_STRIDE_CST,
//                               1,
//                               NULL,
//                     (void **) cell_border_ln_to_gn,
//                               NULL,
//                    (void ***)&cell_ln_to_gn_opp);

//     PDM_distant_neighbor_exch(dn,
//                               sizeof(int),
//                               PDM_STRIDE_CST,
//                               3,
//                               NULL,
//                     (void **) part_ext->graph_comm_cell[0][i_domain],
//                               NULL,
//                    (void ***)&part_ext->graph_comm_cell_opp[0][i_domain]);

//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
//       PDM_log_trace_array_long(cell_border_ln_to_gn[i_part], n_face_part_bound[i_domain][i_part], "cell_border_ln_to_gn = ");
//       PDM_log_trace_array_long(cell_ln_to_gn_opp[i_part], n_face_part_bound[i_domain][i_part], "cell_ln_to_gn_opp = ");
//       PDM_log_trace_array_int(part_ext->graph_comm_cell_opp[0][i_domain][i_part], 3 * n_face_part_bound[i_domain][i_part], "graph_comm_cell_opp = ");
//       free(cell_border_ln_to_gn[i_part]);
//       free(cell_ln_to_gn_opp[i_part]);
//     }
//     free(cell_border_ln_to_gn);
//     free(cell_ln_to_gn_opp);

//     /* Setup for the next distant neighbor */
//     part_ext->graph_comm_cell_opp_idx[0][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
//       part_ext->graph_comm_cell_opp_idx[0][i_domain][i_part] = (int *) malloc( (part_ext->n_cell_per_bound[0][i_domain][i_part]+1) * sizeof(int));
//       part_ext->graph_comm_cell_opp_idx[0][i_domain][i_part][0] = 0;
//       for(int i = 0; i < part_ext->n_cell_per_bound[0][i_domain][i_part]; ++i) {
//         part_ext->graph_comm_cell_opp_idx[0][i_domain][i_part][i+1] = part_ext->graph_comm_cell_opp_idx[0][i_domain][i_part][i] + 1;
//       }
//     }
//     PDM_distant_neighbor_free(dn);
//   }
// }

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
  printf(" n_tot_all_domain      = %i \n", n_tot_all_domain);
  printf(" n_part_loc_all_domain = %i \n", n_part_loc_all_domain);

  /* We flat all partition */
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->n_entity_bound == NULL);

  part_ext->neighbor_idx   = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->neighbor_desc  = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->n_entity_bound = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  part_ext->entity_cell_opp_idx = (int  **) malloc( n_part_loc_all_domain * sizeof(int  *) );

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

      /* Ici il faut faire les raccords entre domaine ---> Fill */


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

  /* Compute idx */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      int* _neighbor_idx        = part_ext->neighbor_idx       [i_part+shift_part];
      int* _entity_cell_opp_n   = part_ext->entity_cell_opp_n  [i_part+shift_part];

      int n_elmt = _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]];

      part_ext->entity_cell_opp_idx[i_part+shift_part] = (int *) malloc((n_elmt+1) * sizeof(int) );
      int* _entity_cell_opp_idx = part_ext->entity_cell_opp_idx[i_part+shift_part];

      printf(" n_elmt = %i \n", n_elmt);
      _entity_cell_opp_idx[0] = 0;
      for(int i_elmt = 0; i_elmt < n_elmt; ++i_elmt) {
        _entity_cell_opp_idx[i_elmt+1] = _entity_cell_opp_idx[i_elmt] + _entity_cell_opp_n[i_elmt];
      }
    }
    shift_part += part_ext->n_part[i_domain];
  }


  if(1 == 1) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        int* _neighbor_idx = part_ext->neighbor_idx                [i_part+shift_part];
        int n_elmt         = _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]];
        int n_data         = part_ext->entity_cell_opp_idx[i_part+shift_part][n_elmt];
        PDM_log_trace_array_int(part_ext->entity_cell_opp_n  [i_part+shift_part], n_elmt  , "entity_cell_opp_n::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp_idx[i_part+shift_part], n_elmt+1, "entity_cell_opp_idx::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp    [i_part+shift_part], n_data  , "entity_cell_opp::");
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }


  PDM_distant_neighbor_free(dn);

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
_extend_cell_graph
(
  PDM_part_extension_t *part_ext,
  int                   i_depth
)
{
  /*
   *  Create the cell_cell and cell_cell_idx with i_proc, i_part, i_entitiy
   *   Only for the border
   */
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell_bound = part_ext->n_cell_per_bound[i_depth-1][i_domain][i_part];
      int n_cell       = part_ext->parts[i_domain][i_part].n_cell;

      /*
       *  Compute size of the associated graph
       */
      int* cell_cell_tmp_idx = malloc( (n_cell + 1) * sizeof(int));
      int* _cell_cell_idx          = part_ext->cell_cell_idx[i_domain][i_part];
      // int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_domain][i_part];

      for(int i_cell = 0; i_cell < n_cell+1; ++i_cell) {
        cell_cell_tmp_idx[i_cell] = 0;
      }

      for(int idx_cell = 0; idx_cell < n_cell_bound; ++idx_cell) {
        int i_cell = part_ext->cell_list[i_depth-1][i_domain][i_part][idx_cell];

        /* From interior */
        int n_cell_connected = _cell_cell_idx[i_cell+1] - _cell_cell_idx[i_cell];
        cell_cell_tmp_idx[i_cell+1] += n_cell_connected;

        /* From border */


      }

      /* */
      for(int idx_cell = 0; idx_cell < n_cell_bound; ++idx_cell) {
        printf(" cell_list [%i] = %i \n", idx_cell, part_ext->cell_list[i_depth-1][i_domain][i_part][idx_cell]);
      }

      if(part_ext->cell_cell_extended[i_domain][i_part] == NULL) {
        part_ext->cell_cell_extended[i_depth][i_domain][i_part] = malloc( (part_ext->cell_cell_extended[i_depth][i_domain][i_part][n_cell] + 1) * sizeof(int));
      } else {
        part_ext->cell_cell_extended[i_depth][i_domain][i_part] = realloc(part_ext->cell_cell_extended[i_depth][i_domain][i_part], (part_ext->cell_cell_extended[i_depth][i_domain][i_part][n_cell] + 1) * sizeof(int));
      }


    }
  }


  /*
   *  On fait 2 échanges :
   *     - Le premier permet d'avoir le graph de comm pour le prochain échange
   *     - Le deuxième permet de connaitre le graphe associé
   */


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

  part_ext->parts = malloc(n_domain * sizeof(_part_t *));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    part_ext->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
  }

  part_ext->neighbor_idx   = NULL;
  part_ext->neighbor_idx   = NULL;
  part_ext->n_entity_bound = NULL;

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
  PDM_part_extension_t *part_ext,
  int                   depth
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

  part_ext->n_tot_part_by_domain = (int *) malloc( part_ext->n_domain * sizeof(int));
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    part_ext->n_tot_part_by_domain[i_domain] = -1;
    int n_part_loc = part_ext->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &part_ext->n_tot_part_by_domain[i_domain], 1, PDM_MPI_INT, PDM_MPI_SUM, part_ext->comm);

    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  part_ext->entity_cell_n   = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell_idx = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell     = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

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

  // part_ext->neighbor_idx    = (int ****) malloc( (depth + 1) * sizeof(int ***));
  // part_ext->neighbor_desc   = (int ****) malloc( (depth + 1) * sizeof(int ***));
  part_ext->graph_comm_cell = (int ****) malloc( (depth + 1) * sizeof(int ***));

  part_ext->cell_list               = (int ****) malloc( (depth + 1) * sizeof(int ***));
  part_ext->graph_comm_cell_opp_idx = (int ****) malloc( (depth + 1) * sizeof(int ***));
  part_ext->graph_comm_cell_opp     = (int ****) malloc( (depth + 1) * sizeof(int ***));
  part_ext->n_cell_per_bound        = (int  ***) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended_idx  = (int ****) malloc( (depth + 1) * sizeof(int ***));
  part_ext->cell_cell_extended      = (int ****) malloc( (depth + 1) * sizeof(int ***));

  for(int i_depth = 0; i_depth < depth+1; ++i_depth) {
    // part_ext->neighbor_idx           [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    // part_ext->neighbor_desc          [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    part_ext->graph_comm_cell        [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    part_ext->cell_list              [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    part_ext->graph_comm_cell_opp_idx[i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    part_ext->graph_comm_cell_opp    [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    part_ext->n_cell_per_bound       [i_depth] = (int  **) malloc( part_ext->n_domain * sizeof(int  *));
    part_ext->cell_cell_extended_idx [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
    part_ext->cell_cell_extended     [i_depth] = (int ***) malloc( part_ext->n_domain * sizeof(int **));
  }

  // TODO : vtx_cell
  _create_cell_graph_comm(part_ext);

  /*
   * Creation du premier graph donc on utilise le cell_cell ET ce qui vient du graph de comm
   * Cette étape initialize une sorte de recurence
   */
  int i_depth_cur = 1;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    part_ext->cell_cell_extended_idx[i_depth_cur][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
    part_ext->cell_cell_extended    [i_depth_cur][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
    part_ext->n_cell_per_bound      [i_depth_cur][i_domain] = (int  *) malloc( part_ext->n_part[i_domain] * sizeof(int  ));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell_bound = part_ext->n_cell_per_bound[i_depth_cur-1][i_domain][i_part];
      int n_cell       = part_ext->parts[i_domain][i_part].n_cell;

      part_ext->cell_cell_extended_idx[i_depth_cur][i_domain][i_part] = malloc( (n_cell + 1) * sizeof(int));
      part_ext->cell_cell_extended    [i_depth_cur][i_domain][i_part] = NULL;

      int* _cell_cell_idx          = part_ext->cell_cell_idx[i_domain][i_part];
      int* _cell_cell              = part_ext->cell_cell    [i_domain][i_part];
      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth_cur][i_domain][i_part];
      int* _cell_cell_extended_n   = malloc( n_cell * sizeof(int));
      int* _graph_comm_cell_opp    = part_ext->graph_comm_cell_opp[i_depth_cur-1][i_domain][i_part];

      int* cell_flag = (int *) malloc( n_cell * sizeof(int));
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        cell_flag[i_cell] = 0;
        _cell_cell_extended_n[i_cell] = 0;
      }

      /*
       *  Compute size of the associated graph
       */
      part_ext->n_cell_per_bound[i_depth_cur][i_domain][i_part] = 0;
      for(int idx_cell = 0; idx_cell < n_cell_bound; ++idx_cell) {
        int i_cell = part_ext->cell_list[i_depth_cur-1][i_domain][i_part][idx_cell];

        /* Si la cellule n'est pas flaggé on rajoute son voisinage */
        if(cell_flag[i_cell] == 0) {
          int n_cell_connected = _cell_cell_idx[i_cell+1] - _cell_cell_idx[i_cell];
          cell_flag[i_cell] = 1;
          part_ext->n_cell_per_bound[i_depth_cur][i_domain][i_part] += 1;
          _cell_cell_extended_n[i_cell] += n_cell_connected;
          _cell_cell_extended_n[i_cell] += 1;
        } else {
          // Count only the border (for the first pass one face )
          // If graph is deduced by vertex connectivity we pass multiple time here (also true for cell in corner )
          _cell_cell_extended_n[i_cell] += 1;
        }
      }

      //
      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
      }

      /* Reset */
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        cell_flag[i_cell] = 0;
        _cell_cell_extended_n[i_cell] = 0;
      }

      part_ext->cell_cell_extended[i_depth_cur][i_domain][i_part] = (int*) malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int) );
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth_cur][i_domain][i_part];


      for(int idx_cell = 0; idx_cell < n_cell_bound; ++idx_cell) {
        int i_cell = part_ext->cell_list[i_depth_cur-1][i_domain][i_part][idx_cell];

        if(cell_flag[i_cell] == 0) {
          cell_flag[i_cell] = 1;

          for(int idx_connect_cell = _cell_cell_idx[i_cell]; idx_connect_cell < _cell_cell_idx[i_cell+1]; ++idx_connect_cell) {
            int i_cell_connected = PDM_ABS(_cell_cell[idx_connect_cell])-1;
            int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];

            _cell_cell_extended[3*idx_write  ] = i_rank;
            _cell_cell_extended[3*idx_write+1] = i_part;
            _cell_cell_extended[3*idx_write+2] = i_cell_connected;

            _cell_cell_extended_n[i_cell] += 1;
          }

          int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
          _cell_cell_extended[3*idx_write  ] = _graph_comm_cell_opp[3*idx_cell  ];
          _cell_cell_extended[3*idx_write+1] = _graph_comm_cell_opp[3*idx_cell+1];
          _cell_cell_extended[3*idx_write+2] = _graph_comm_cell_opp[3*idx_cell+2];

          _cell_cell_extended_n[i_cell] += 1;

        } else {

          int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
          _cell_cell_extended[3*idx_write  ] = _graph_comm_cell_opp[3*idx_cell  ];
          _cell_cell_extended[3*idx_write+1] = _graph_comm_cell_opp[3*idx_cell+1];
          _cell_cell_extended[3*idx_write+2] = _graph_comm_cell_opp[3*idx_cell+2];

          // Count only the border
          _cell_cell_extended_n[i_cell] += 1;
        }

      }

      printf("part_ext->n_cell_per_bound[%i][%i][%i] = %i\n", i_depth_cur, i_domain, i_part, part_ext->n_cell_per_bound[i_depth_cur][i_domain][i_part]);

      PDM_log_trace_array_int(_cell_cell_extended_idx, n_cell+1                           , "_cell_cell_extended_idx::");
      PDM_log_trace_array_int(_cell_cell_extended_n  , n_cell                             , "_cell_cell_extended_n::"  );
      PDM_log_trace_array_int(_cell_cell_extended    , 3 * _cell_cell_extended_idx[n_cell], "cell_cell_extended::");

      printf(" ------------------------ \n");
      for(int i = 0; i < n_cell; ++i) {
        printf(" cell_cell_extended[%i][%i] = ", i_part, i);
        for(int idx = _cell_cell_extended_idx[i]; idx < _cell_cell_extended_idx[i+1]; ++idx) {
          printf(" ( i_part_opp = %i | i_entity_opp = %i ) ", _cell_cell_extended[3*idx+1], _cell_cell_extended[3*idx+2]);
        }
        printf("\n");
      }

      free(_cell_cell_extended_n);
      free(cell_flag);
    }
  }

  for(int i_depth = 1; i_depth < depth+1; ++i_depth) {
    _extend_cell_graph(part_ext, i_depth);
  }

  for(int i_depth = 0; i_depth < depth+1; ++i_depth) {
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      free(part_ext->n_cell_per_bound[i_depth][i_domain]);
    }
    free(part_ext->n_cell_per_bound[i_depth]);
  }
  free(part_ext->n_cell_per_bound);

  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //   for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //     free(face_cell[i_domain][i_part]);
  //     free(face_cell_idx[i_domain][i_part]);
  //     free(cell_cell_idx[i_domain][i_part]);
  //     free(cell_cell[i_domain][i_part]);
  //   }
  //   free(face_cell[i_domain]);
  //   free(face_cell_idx[i_domain]);
  //   free(cell_cell_idx[i_domain]);
  //   free(cell_cell[i_domain]);
  // }
  // free(face_cell);
  // free(face_cell_idx);
  // free(cell_cell_idx);
  // free(cell_cell);

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
        free(part_ext->entity_cell_opp_idx[i_part+shift_part]);
        free(part_ext->entity_cell_opp_n  [i_part+shift_part]);
        free(part_ext->entity_cell_opp    [i_part+shift_part]);

        if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
          free(part_ext->entity_cell_idx[i_part+shift_part]);
          free(part_ext->entity_cell_n  [i_part+shift_part]);
          free(part_ext->entity_cell    [i_part+shift_part]);
        }

      }
      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->n_entity_bound);
  }

  /* Only shortcut of user data */
  free(part_ext->entity_cell_idx    );
  free(part_ext->entity_cell        );

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





  /*
   *  So we have the comm graph between cell directly and the cell_list
   *    we need to prepare the newt send
   *    for the first rank this step is enought depth = 0
   */
  // int*** graph_comm_cell_opp = (int ***) malloc( depth * sizeof(int **));
  // for(int i_depth = 1; i_depth < depth+1; ++i_depth) {
  //   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

  //     n_cell_per_bound[i_depth][i_domain] = (int  *) malloc( part_ext->n_part[i_domain] * sizeof(int  ) );
  //     cell_list       [i_depth][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

  //     /*
  //      *  Collect data connexion AND next graph_comm
  //      */
  //     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //       n_cell_per_bound[i_depth][i_domain][i_part] = 0;
  //       for(int i = 0; i < n_cell_per_bound[i_depth-1][i_domain][i_part]; ++i) {
  //         int i_cell = cell_list[i_depth-1][i_domain][i_part][i];
  //         int beg = part_ext->parts[i_domain][i_part].cell_face_idx[i_cell  ];
  //         int end = part_ext->parts[i_domain][i_part].cell_face_idx[i_cell+1];
  //         for(int idx_face = beg; idx_face < end; ++idx_face) {
  //           int i_face = PDM_ABS(part_ext->parts[i_domain][i_part].cell_face[idx_face])-1;
  //           int i_cell2 = part_ext->parts[i_domain][i_part].face_cell[2*i_face+1]-1;
  //           n_cell_per_bound[i_depth][i_domain][i_part] += 1;
  //           if(i_cell2 > -1) {
  //             n_cell_per_bound[i_depth][i_domain][i_part] += 1;
  //           }
  //         }
  //       }
  //       printf(" n_cell_per_bound[%i][%i][%i] = %i \n", i_depth, i_domain, i_part, n_cell_per_bound[i_depth][i_domain][i_part]);
  //     }

  //     graph_comm_cell_opp    [i_depth][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
  //     graph_comm_cell_opp_idx[i_depth][i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

  //     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //       cell_list              [i_depth][i_domain][i_part] = (int *) malloc(     (n_cell_per_bound[i_depth][i_domain][i_part]  ) * sizeof(int));
  //       graph_comm_cell_opp    [i_depth][i_domain][i_part] = (int *) malloc( 3 * (n_cell_per_bound[i_depth][i_domain][i_part]  ) * sizeof(int));
  //       graph_comm_cell_opp_idx[i_depth][i_domain][i_part] = (int *) malloc(     (n_cell_per_bound[i_depth][i_domain][i_part]+1) * sizeof(int));
  //       graph_comm_cell_opp_idx[i_depth][i_domain][i_part][0] = 0;

  //       int idx = 0;
  //       for(int i = 0; i < n_cell_per_bound[i_depth-1][i_domain][i_part]; ++i) {
  //         int i_cell = cell_list[i_depth-1][i_domain][i_part][i];
  //         int beg = part_ext->parts[i_domain][i_part].cell_face_idx[i_cell  ];
  //         int end = part_ext->parts[i_domain][i_part].cell_face_idx[i_cell+1];
  //         for(int idx_face = beg; idx_face < end; ++idx_face) {
  //           int i_face = PDM_ABS(part_ext->parts[i_domain][i_part].cell_face[idx_face])-1;
  //           int i_cell1 = part_ext->parts[i_domain][i_part].face_cell[2*i_face  ]-1;
  //           int i_cell2 = part_ext->parts[i_domain][i_part].face_cell[2*i_face+1]-1;
  //           n_cell_per_bound[i_depth][i_domain][i_part] += 1;

  //           graph_comm_cell_opp[i_depth][i_domain][i_part][3*idx  ] = idx;
  //           graph_comm_cell_opp[i_depth][i_domain][i_part][3*idx+1] = i_rank;
  //           graph_comm_cell_opp[i_depth][i_domain][i_part][3*idx+2] = i_part;

  //           cell_list[i_depth][i_domain][i_part][idx++] = i_cell1;

  //           if(i_cell2 > -1) {
  //             n_cell_per_bound[i_depth][i_domain][i_part] += 1;
  //             graph_comm_cell_opp[i_depth][i_domain][i_part][3*idx  ] = idx;
  //             graph_comm_cell_opp[i_depth][i_domain][i_part][3*idx+1] = i_rank;
  //             graph_comm_cell_opp[i_depth][i_domain][i_part][3*idx+2] = i_part;

  //             cell_list[i_depth][i_domain][i_part][idx++] = i_cell2;

  //           }
  //         }
  //       }

  //     }


  //     printf(" manage depth = %i \n", i_depth);

  //     PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
  //                                                              part_ext->n_part[i_domain],
  //                                                              n_cell_per_bound       [i_depth-1][i_domain],
  //                                                              graph_comm_cell_opp_idx[i_depth-1][i_domain],
  //                                                              graph_comm_cell_opp    [i_depth-1][i_domain]);

  //     /*
  //      *  Prepare exchange with cell_list of previous depth
  //      */

  //     PDM_distant_neighbor_free(dn);

  //   }
  // }



  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //   for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //     free(neighbor_idx [i_domain][i_part]);
  //     free(neighbor_desc[i_domain][i_part]);
  //     free(graph_comm_cell[i_domain][i_part]);
  //   }
  //   free(neighbor_idx[i_domain]);
  //   free(neighbor_desc[i_domain]);
  //   free(n_face_part_bound[i_domain]);
  //   free(graph_comm_cell[i_domain]);
  // }
  // free(neighbor_idx);
  // free(neighbor_desc);
  // free(n_face_part_bound);
  // free(graph_comm_cell);


  // for(int i_depth = 0; i_depth < depth+1; ++i_depth) {
  // // for(int i_depth = 0; i_depth < 1; ++i_depth) {
  //   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //       free(graph_comm_cell_opp    [i_depth][i_domain][i_part]);
  //       free(graph_comm_cell_opp_idx[i_depth][i_domain][i_part]);
  //       free(cell_list[i_depth][i_domain][i_part]);
  //     }
  //     free(graph_comm_cell_opp    [i_depth][i_domain]);
  //     free(graph_comm_cell_opp_idx[i_depth][i_domain]);
  //     free(cell_list[i_depth][i_domain]);
  //   }
  //   free(graph_comm_cell_opp    [i_depth]);
  //   free(graph_comm_cell_opp_idx[i_depth]);
  //   free(cell_list[i_depth]);
  // }
  // free(graph_comm_cell_opp);
  // free(graph_comm_cell_opp_idx);
  // free(cell_list);
