#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_partitioning.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_partitioning_method_t method  = PDM_PARTITIONING_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_partitioning_method_t method  = PDM_PARTITIONING_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int          id;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_gen_init(&id,
                      comm,
                      n_vtx_seg,
                      length,
            		      0.,
		                  0.,
		                  0.);

  PDM_dcube_gen_dim_get(id,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtxL,
                         &dFaceGroupL);

  PDM_dcube_gen_data_get(id,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  if (0 == 1) {

    PDM_printf("[%i] n_face_group    : %i\n", i_rank, n_face_group);
    PDM_printf("[%i] dn_cell        : %i\n", i_rank, dn_cell);
    PDM_printf("[%i] dn_face        : %i\n", i_rank, dn_face);
    PDM_printf("[%i] dn_vtx         : %i\n", i_rank, dn_vtx);

    PDM_printf("[%i] dface_cell     : ", i_rank);
    for (int i = 0; i < 2 * dn_face; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_cell[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx_idx   : ", i_rank);
    for (int i = 0; i < dn_face + 1; i++)
      PDM_printf(" %i", dface_vtx_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx      : ", i_rank);
    for (int i = 0; i < dface_vtx_idx[dn_face]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_vtx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dvtx_coord     : ", i_rank);
    for (int i = 0; i < 3*dn_vtx; i++)
      PDM_printf(" %12.5e", dvtx_coord[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group_idx : ", i_rank);
    for (int i = 0; i < n_face_group + 1; i++)
      PDM_printf(" %i", dface_group_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group    : ", i_rank);
    for (int i = 0; i < dface_group_idx[n_face_group]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_group[i]);
    PDM_printf("\n");

  }

  // printf(" PDM_PART_NONE::%d\n"     , PDM_PART_NONE);
  // printf(" PDM_PART_NULL::%d\n"     , PDM_PART_NULL);
  // printf(" PDM_PART_FACE_CELL::%d\n", PDM_PART_FACE_CELL);
  // printf(" PDM_PART_CELL_FACE::%d\n", PDM_PART_CELL_FACE);

  // int flags = PDM_PART_FACE_CELL|PDM_PART_CELL_FACE;
  // printf("PDM_HASFLAG(flags, PDM_PART_FACE_CELL) :: %d\n", PDM_HASFLAG(flags, PDM_PART_FACE_CELL) );
  // printf("PDM_HASFLAG(flags, PDM_PART_CELL_FACE) :: %d\n", PDM_HASFLAG(flags, PDM_PART_CELL_FACE) );
  // printf("PDM_HASFLAG(flags, PDM_PART_FACE_VTX) :: %d\n", PDM_HASFLAG(flags, PDM_PART_FACE_VTX) );
  // printf("x::PDM_HASFLAG(flags, PDM_PART_FACE_VTX) :: %x\n", PDM_PART_FACE_VTX);


  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_FACE_CELL, NULL);
  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_CELL_FACE, NULL);
  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_FACE_VTX, NULL);
  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_FACE_VTX|PDM_PART_FACE_CELL, NULL);

  // PDM_dmesh_partitioning_get(1, PDM_PART_FACE_CELL, NULL);
  // PDM_dmesh_partitioning_get(1, PDM_PART_CELL_FACE, NULL);
  // PDM_dmesh_partitioning_get(1, PDM_PART_FACE_VTX, NULL);
  // PDM_dmesh_partitioning_get(1, PDM_PART_FACE_VTX|PDM_PART_FACE_CELL, NULL);

  // int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */
  PDM_g_num_t* cell_distribution = PDM_compute_entity_distribution(comm, dn_cell);
  PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t* part_distribution = PDM_compute_entity_distribution(comm, n_part );

  // printf("part_distribution::\n");
  // for(int i_part = 0; i_part < n_rank+1; ++i_part){
  //   printf("%d ", part_distribution[i_part]);
  // }
  // printf("\n");

  /*
   * Compute dual graph
   */
  int* dual_graph_idx;
  PDM_g_num_t* dual_graph;
  int* dcell_face_idx;
  PDM_g_num_t* dcell_face;
  PDM_para_graph_dual_from_face_cell(comm,
                                     cell_distribution,
                                     face_distribution,
                                     dface_cell,
                     (int        **) &dual_graph_idx,
                     (PDM_g_num_t**) &dual_graph,
                                     1,
                     (int        **) &dcell_face_idx,
                     (PDM_g_num_t**) &dcell_face);

  /*
   * Split it !!! CAUTION dn_cell can be different of the size of dual graph !!!
   */
  // printf("PDM_split_graph\n");
  int* cell_part    = (int *) malloc( sizeof(int) * dn_cell );
  int* dcell_weight = (int *) malloc( sizeof(int) * dn_cell );
  for(int i = 0; i < dn_cell; ++i){
    dcell_weight[i] = dual_graph_idx[i+1] - dual_graph_idx[i];
  }

  if( 0 == 1 ){
    printf("n_cell_block:: %d \n", dn_cell);
    for(int i = 0; i < dn_cell; ++i){
      printf(" dual_graph_idx = %d ---> \n", dual_graph_idx[i]);
      for(int i_data = dual_graph_idx[i]; i_data < dual_graph_idx[i+1]; ++i_data){
        // printf("%d ", dual_graph[i_data]);
        printf("\t dual_graph[%d] = %d \n", i_data, dual_graph[i_data]);
      }
      printf("\n");
    }
  }

  // PDM_split_graph(comm, dual_graph_idx, dual_graph, NULL, cell_part, dn_cell, n_rank);
  // PDM_split_graph(comm, dual_graph_idx, dual_graph, NULL, cell_part, dn_cell, n_part);
  PDM_split_graph(comm, dual_graph_idx, dual_graph, NULL, cell_part, dn_cell, part_distribution[n_rank]-1);
  // PDM_split_graph(comm, dual_graph_idx, dual_graph, NULL, cell_part, dn_cell, n_rank);


  // printf("cell_part[%d]::", dn_cell);
  // for(int i = 0; i < dn_cell; ++i){
  //   printf("%d ", cell_part[i]);
  // }
  // printf("\n");

  /*
   * On dispose pour chaque cellule de la partition associé : il faut retrouver le
   *  cell_ln_to_gn
   *    Idée : faire un part_to_block avec ln_to_gn == cell_part et distribution imposé = dpart_proc
   *    On devrait recupérer n_part block de données avec le ln_to_gn !
   */
  int** pcell_ln_to_gn;
  int*  pn_cell;

  int n_res_part = PDM_generate_part_cell_ln_to_gn(comm,
                                                   part_distribution,
                                                   cell_distribution,
                                                   cell_part,
                                       (int ** )  &pn_cell,
                                       (int ***)  &pcell_ln_to_gn);

  /*
   *  At this stage we have the cell_ln_to_gn :
   *      --> We need to deduce the other if needed
   *           --> face_ln_to_gn
   *           --> vtx_ln_to_gn
   *  We can do by the relationship of each varibles (connectivity)
   *  For example face_ln_to_gn can be deduce with cell_face connectivity if we have cell_ln_to_gn
   */
  int** pcell_face_idx;
  int** pcell_face;
  int*  pn_faces;
  PDM_g_num_t** pface_ln_to_gn = NULL;

  PDM_generate_part_entity_ln_to_gn(comm,
                                    part_distribution,
                                    cell_distribution,
                                    dcell_face_idx,
                                    dcell_face,
                                    n_res_part,
                                    pn_cell,
                                    pcell_ln_to_gn,
                (int         ** )  &pn_faces,
                (PDM_g_num_t ***)  &pface_ln_to_gn,
                (int         ***)  &pcell_face_idx,
                (int         ***)  &pcell_face);

  /*
   * Generate vtx
   */

  int** pface_vtx_idx;
  int** pface_vtx;
  int*  pn_vtx;
  PDM_g_num_t** pvtx_ln_to_gn = NULL;

  PDM_generate_part_entity_ln_to_gn(comm,
                                    part_distribution,
                                    face_distribution,
                                    dface_vtx_idx,
                                    dface_vtx,
                                    n_res_part,
                                    pn_cell,
                                    pcell_ln_to_gn,
                (int         ** )  &pn_vtx,
                (PDM_g_num_t ***)  &pvtx_ln_to_gn,
                (int         ***)  &pface_vtx_idx,
                (int         ***)  &pface_vtx);


  /*
   * On doit calculer le dcell_face car avec lui et le cell_ln_to_gn on retrouve facilement
   *   le dcell_face sur la partition donc on n'a plus qu'a trié pour avoir le face_ln_to_gn
   *      Une fois le face_ln_to_gn trouvé on doit avoir facilement
   *      le dface_cell mais en numéro absolu (attention au signe ) et la on a juste à chercher
   *      dans le cell_ln_to_gn l'index correspondant pour transformé le numero global en numéro local
   */

  /*
   *  L'astuce de tri pour le face_cell et le meme que pour le face_vtx -->
   *    permet d'avoir le vtx_ln_to_gn
   */

  /*
   *  Boundary condition (face group )
   */
  PDM_g_num_t** pface_group_ln_to_gn;
  int** pface_group;
  int** pface_group_idx;
  PDM_generate_part_face_group_ln_to_gn(comm,
                                        face_distribution,
                                        dface_group_idx,
                                        dface_group,
                                        n_res_part,
                                        n_face_group,
                                        pn_faces,
                                        pface_ln_to_gn,
                     (PDM_g_num_t ***) &pface_group_ln_to_gn,
                     (int         ***) &pface_group,
                     (int         ***) &pface_group_idx);

  /*
   *  Edge
   */


  /*
   * Graph communication build
   */
  int** pproc_face_bound_idx;
  int** ppart_face_bound_idx;
  int** pface_bound_idx;
  PDM_generate_entity_graph_comm(comm,
                                 part_distribution,
                                 face_distribution,
                                 n_part,
                                 pn_faces,
                                 pface_ln_to_gn,
                      (int ***) &pproc_face_bound_idx,
                      (int ***) &ppart_face_bound_idx,
                      (int ***) &pface_bound_idx);





  // Attention on veut garder l'orientation donc il y a un signe dans le face_cell / cell_face
  // Reflechir sur les connectivité d'edge également ...

  /*
   * Free
   */
  free(dual_graph_idx);
  free(dual_graph);
  free(cell_part);
  free(dcell_face);
  free(dcell_face_idx);
  free(dcell_weight);
  free(cell_distribution);
  free(face_distribution);
  free(part_distribution);
  for(int i_part = 0; i_part < n_res_part; ++i_part){
    free(pface_ln_to_gn[i_part]);
    free(pcell_ln_to_gn[i_part]);
    free(pvtx_ln_to_gn[i_part]);
    free(pcell_face[i_part]);
    free(pcell_face_idx[i_part]);
    free(pface_vtx_idx[i_part]);
    free(pface_vtx[i_part]);
    free(pface_group_ln_to_gn[i_part]);
    free(pface_group[i_part]);
    free(pface_group_idx[i_part]);
    free(pproc_face_bound_idx[i_part]);
    free(ppart_face_bound_idx[i_part]);
    free(pface_bound_idx[i_part]);
  }
  free(pcell_face);
  free(pcell_face_idx);
  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_ln_to_gn);
  free(pproc_face_bound_idx);
  free(ppart_face_bound_idx);
  free(pface_bound_idx);
  free(pface_ln_to_gn);
  free(pn_cell);
  free(pn_faces);
  free(pn_vtx);
  free(pface_group_ln_to_gn);
  free(pface_group);
  free(pface_group_idx);



  PDM_dcube_gen_free(id);

  PDM_MPI_Finalize();

  return 0;
}


  // Unit test ...
  // int test_unique[10] = {10, 2, 1, 12, 31, 2, 31, 4, 5, 2};
  // // int ns = PDM_inplace_unique_long(test_unique, 4, 10);
  // // PDM_quick_sort_long(test_unique, 0, 9);
  // int ns = PDM_inplace_unique_long(test_unique, 0, 9);
  // printf("ns::%d\n", ns);
  // printf("test_unique::");
  // for(int i = 0; i < ns; ++i){
  //   printf("%d ", test_unique[i]);
  // }
  // printf("\n");

  // Unit test ...
  // int cell_cell_n[4] = {4, 5, 3, 2};
  // PDM_g_num_t dual_graph_test[14] = {31, 10, 5, 6,
  //                                    8, 8, 2, 1, 4,
  //                                    2, 3, 2,
  //                                    91, 92 };
  // int dual_graph_test_idx[5];

  // PDM_compress_connectivity(dual_graph_test, dual_graph_test_idx, cell_cell_n, 4);


  // if( 1 == 1 ){
  //   // printf("n_cell_block:: %d \n", n_cell_block);
  //   for(int i = 0; i < 4; ++i){
  //     printf(" dual_graph_test_idx = %d ---> \n", dual_graph_test_idx[i]);
  //     for(int i_data = dual_graph_test_idx[i]; i_data < dual_graph_test_idx[i+1]; ++i_data){
  //       // printf("%d ", _dual_graph[i_data]);
  //       printf("\t _dual_graph[%d] = %d \n", i_data, dual_graph_test[i_data]);
  //     }
  //     printf("\n");
  //   }
  // }
