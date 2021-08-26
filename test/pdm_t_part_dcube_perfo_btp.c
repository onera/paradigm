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
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
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
  double             length  = 1.;
  int                n_part   = 1;
  int                post    = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
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

  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtxL,
                         &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
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
  int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_create(&ppart_id,
                  comm,
                  method,
                  "PDM_PART_RENUM_CELL_NONE",
                  "PDM_PART_RENUM_FACE_NONE",
                  n_property_cell,
                  renum_properties_cell,
                  n_property_face,
                  renum_properties_face,
                  n_part,
                  dn_cell,
                  dn_face,
                  dn_vtx,
                  n_face_group,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  have_dcell_part,
                  dcell_part,
                  dface_cell,
                  dface_vtx_idx,
                  dface_vtx,
                  NULL,
                  dvtx_coord,
                  NULL,
                  dface_group_idx,
                  dface_group);

  double  *elapsed  = NULL;
  double  *cpu      = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys  = NULL;

  PDM_part_time_get(ppart_id,
                    &elapsed,
                    &cpu,
                    &cpu_user,
                    &cpu_sys);

  PDM_printf("[%i]   - elapsed total                    : %12.5e\n", i_rank, elapsed[0]);
  PDM_printf("[%i]   - elapsed building graph           : %12.5e\n", i_rank, elapsed[1]);
  PDM_printf("[%i]   - elapsed splitting graph          : %12.5e\n", i_rank, elapsed[2]);
  PDM_printf("[%i]   - elapsed building mesh partitions : %12.5e\n", i_rank, elapsed[3]);

  PDM_printf("[%i]   - cpu total                        : %12.5e\n", i_rank, cpu[0]);
  PDM_printf("[%i]   - cpu building graph               : %12.5e\n", i_rank, cpu[1]);
  PDM_printf("[%i]   - cpu splitting graph              : %12.5e\n", i_rank, cpu[2]);
  PDM_printf("[%i]   - cpu building mesh partitions     : %12.5e\n", i_rank, cpu[3]);

  PDM_printf("[%i]   - cpu_user total                   : %12.5e\n", i_rank, cpu_user[0]);
  PDM_printf("[%i]   - cpu_user building graph          : %12.5e\n", i_rank, cpu_user[1]);
  PDM_printf("[%i]   - cpu_user splitting graph         : %12.5e\n", i_rank, cpu_user[2]);
  PDM_printf("[%i]   - cpu_user building mesh partitions: %12.5e\n", i_rank, cpu_user[3]);

  PDM_printf("[%i]   - cpu_sys total                    : %12.5e\n", i_rank, cpu_sys[0]);
  PDM_printf("[%i]   - cpu_sys building graph           : %12.5e\n", i_rank, cpu_sys[1]);
  PDM_printf("[%i]   - cpu_sys splitting graph          : %12.5e\n", i_rank, cpu_sys[2]);
  PDM_printf("[%i]   - cpu_sys building mesh partitions : %12.5e\n", i_rank, cpu_sys[3]);

  struct timeval t_elaps_fin;
  gettimeofday(&t_elaps_fin, NULL);

  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
                         (t_elaps_debut.tv_usec + 1000000 *
                          t_elaps_debut.tv_sec);
  long tranche_elapsed_max = tranche_elapsed;
  double t_elapsed = (double) tranche_elapsed_max/1000000.;

  PDM_printf("[%i]   - TEMPS DANS PART_CUBE  : %12.5e\n", i_rank,  t_elapsed);

  int *n_cell            = (int *) malloc(n_part * sizeof(int));
  int *n_face            = (int *) malloc(n_part * sizeof(int));
  int *n_face_part_bound = (int *) malloc(n_part * sizeof(int));
  int *n_vtx             = (int *) malloc(n_part * sizeof(int));
  int *n_proc            = (int *) malloc(n_part * sizeof(int));
  int *n_total_part      = (int *) malloc(n_part * sizeof(int));
  int *scell_face        = (int *) malloc(n_part * sizeof(int));
  int *sface_vtx         = (int *) malloc(n_part * sizeof(int));
  int *sface_group       = (int *) malloc(n_part * sizeof(int));
  int *n_face_group2     = (int *) malloc(n_part * sizeof(int));

  int          **cell_tag                 = (int          **) malloc( n_part * sizeof(int          *));
  int          **cell_face_idx            = (int          **) malloc( n_part * sizeof(int          *));
  int          **cell_face                = (int          **) malloc( n_part * sizeof(int          *));
  PDM_g_num_t  **cell_ln_to_gn            = (PDM_g_num_t  **) malloc( n_part * sizeof(PDM_g_num_t  *));
  int          **face_tag                 = (int          **) malloc( n_part * sizeof(int          *));
  int          **face_cell                = (int          **) malloc( n_part * sizeof(int          *));
  int          **face_vtx_idx             = (int          **) malloc( n_part * sizeof(int          *));
  int          **face_vtx                 = (int          **) malloc( n_part * sizeof(int          *));
  PDM_g_num_t  **face_ln_to_gn            = (PDM_g_num_t  **) malloc( n_part * sizeof(PDM_g_num_t  *));
  int          **face_part_bound_proc_idx = (int          **) malloc( n_part * sizeof(int          *));
  int          **face_part_bound_part_idx = (int          **) malloc( n_part * sizeof(int          *));
  int          **face_part_bound          = (int          **) malloc( n_part * sizeof(int          *));
  int          **vtx_tag                  = (int          **) malloc( n_part * sizeof(int          *));
  double       **vtx                      = (double       **) malloc( n_part * sizeof(double       *));
  PDM_g_num_t  **vtx_ln_to_gn             = (PDM_g_num_t  **) malloc( n_part * sizeof(PDM_g_num_t  *));
  int          **face_group_idx           = (int          **) malloc( n_part * sizeof(int          *));
  int          **face_group               = (int          **) malloc( n_part * sizeof(int          *));
  PDM_g_num_t  **face_group_ln_to_gn      = (PDM_g_num_t  **) malloc( n_part * sizeof(PDM_g_num_t  *));

  for (int i_part = 0; i_part < n_part; i_part++) {


    PDM_part_part_dim_get(ppart_id,
                          i_part,
                          &n_cell[i_part],
                          &n_face[i_part],
                          &n_face_part_bound[i_part],
                          &n_vtx[i_part],
                          &n_proc[i_part],
                          &n_total_part[i_part],
                          &scell_face[i_part],
                          &sface_vtx[i_part],
                          &sface_group[i_part],
                          &n_face_group2[i_part]);


    PDM_part_part_val_get(ppart_id,
                          i_part,
                          &cell_tag[i_part],
                          &cell_face_idx[i_part],
                          &cell_face[i_part],
                          &cell_ln_to_gn[i_part],
                          &face_tag[i_part],
                          &face_cell[i_part],
                          &face_vtx_idx[i_part],
                          &face_vtx[i_part],
                          &face_ln_to_gn[i_part],
                          &face_part_bound_proc_idx[i_part],
                          &face_part_bound_part_idx[i_part],
                          &face_part_bound[i_part],
                          &vtx_tag[i_part],
                          &vtx[i_part],
                          &vtx_ln_to_gn[i_part],
                          &face_group_idx[i_part],
                          &face_group[i_part],
                          &face_group_ln_to_gn[i_part]);

  }

  /* Calculs statistiques */

  int    cells_average;
  int    cells_median;
  double cells_std_deviation;
  int    cells_min;
  int    cells_max;
  int    bound_part_faces_average;
  int    bound_part_faces_median;
  double bound_part_faces_std_deviation;
  int    bound_part_faces_min;
  int    bound_part_faces_max;
  int    bound_part_faces_sum;

  PDM_part_stat_get(ppart_id,
                    &cells_average,
                    &cells_median,
                    &cells_std_deviation,
                    &cells_min,
                    &cells_max,
                    &bound_part_faces_average,
                    &bound_part_faces_median,
                    &bound_part_faces_std_deviation,
                    &bound_part_faces_min,
                    &bound_part_faces_max,
                    &bound_part_faces_sum);

  if (i_rank == 0) {
    PDM_printf("Statistics :\n");
    PDM_printf("  - Number of cells :\n");
    PDM_printf("       * average            : %i\n", cells_average);
    PDM_printf("       * median             : %i\n", cells_median);
    PDM_printf("       * standard deviation : %12.5e\n", cells_std_deviation);
    PDM_printf("       * min                : %i\n", cells_min);
    PDM_printf("       * max                : %i\n", cells_max);
    PDM_printf("  - Number of faces exchanging with another partition :\n");
    PDM_printf("       * average            : %i\n", bound_part_faces_average);
    PDM_printf("       * median             : %i\n", bound_part_faces_median);
    PDM_printf("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
    PDM_printf("       * min                : %i\n", bound_part_faces_min);
    PDM_printf("       * max                : %i\n", bound_part_faces_max);
    PDM_printf("       * total              : %i\n", bound_part_faces_sum);
  }

  /*
   * btp to bench
   */
  int n_field = 5;
  PDM_g_num_t* dcell_distrib = PDM_compute_entity_distribution(comm, dn_cell);
  PDM_g_num_t* dface_distrib = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t* dvtx_distrib  = PDM_compute_entity_distribution(comm, dn_vtx );

  PDM_block_to_part_t *btp_cell = PDM_block_to_part_create (dcell_distrib,
                                     (const PDM_g_num_t **) cell_ln_to_gn,
                                                            n_cell,
                                                            n_part,
                                                            comm);

  double* dcell_field = malloc( n_field * dn_cell * sizeof(double));
  for(int i = 0; i < dn_cell; ++i ) {
    for(int i_field = 0; i_field < n_field; ++i_field) {
      dcell_field[i_field*n_field+i] = 100000000*i_field + i;
    }
  }

  double **pcell_field;
  int stride_one = n_field;
  PDM_block_to_part_exch2 (btp_cell,
                           sizeof(double),
                           PDM_STRIDE_CST,
                           &stride_one,
                  (void *) dcell_field,
                           NULL,
                (void ***) &pcell_field);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pcell_field[i_part]);
  }
  free(pcell_field);


  free(dcell_field);

  PDM_block_to_part_free(btp_cell);

  PDM_block_to_part_t *btp_face = PDM_block_to_part_create (dface_distrib,
                                     (const PDM_g_num_t **) face_ln_to_gn,
                                                            n_face,
                                                            n_part,
                                                            comm);

  PDM_block_to_part_free(btp_face);

  PDM_block_to_part_t *btp_vtx = PDM_block_to_part_create (dvtx_distrib,
                                     (const PDM_g_num_t **) vtx_ln_to_gn,
                                                            n_vtx,
                                                            n_part,
                                                            comm);

  PDM_block_to_part_free(btp_vtx);

  free(dcell_distrib);
  free(dface_distrib);
  free(dvtx_distrib );


  free(n_cell           );
  free(n_face           );
  free(n_face_part_bound);
  free(n_vtx            );
  free(n_proc           );
  free(n_total_part     );
  free(scell_face       );
  free(sface_vtx        );
  free(sface_group      );
  free(n_face_group2    );

  free(cell_tag                );
  free(cell_face_idx           );
  free(cell_face               );
  free(cell_ln_to_gn           );
  free(face_tag                );
  free(face_cell               );
  free(face_vtx_idx            );
  free(face_vtx                );
  free(face_ln_to_gn           );
  free(face_part_bound_proc_idx);
  free(face_part_bound_part_idx);
  free(face_part_bound         );
  free(vtx_tag                 );
  free(vtx                     );
  free(vtx_ln_to_gn            );
  free(face_group_idx          );
  free(face_group              );
  free(face_group_ln_to_gn     );
  free(dcell_part);
  PDM_part_free(ppart_id);

  PDM_dcube_gen_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
