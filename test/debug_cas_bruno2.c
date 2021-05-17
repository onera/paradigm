#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void
_read_data
(
 PDM_MPI_Comm  comm,
 PDM_g_num_t  *surf_ng_face,
 PDM_g_num_t  *surf_ng_vtx,
 int          *surf_n_face,
 PDM_g_num_t **surf_face_g_num,
 int         **surf_face_vtx_idx,
 int         **surf_face_vtx,
 int          *surf_n_vtx,
 PDM_g_num_t **surf_vtx_g_num,
 double      **surf_vtx_coord,
 int          *n_cell,
 PDM_g_num_t **cell_g_num,
 int         **cell_face_idx,
 int         **cell_face,
 double      **cell_center,
 int          *n_face,
 PDM_g_num_t **face_g_num,
 int         **face_vtx_idx,
 int         **face_vtx,
 int          *n_vtx,
 PDM_g_num_t **vtx_g_num,
 double      **vtx_coord
 )
{
  int i_rank, n_rank;

  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  FILE *f = NULL;

  char directory[999];
  sprintf(directory, "/home/bastien/workspace/debug/dist_cloud_surf/BMaugars/zz_DEBUG_%drank", n_rank);

  char filename[9999];


  /*
   *  Read surface mesh
   */
  *surf_ng_face = 0;
  *surf_ng_vtx  = 0;
  sprintf(filename, "%s/surf_size_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  fscanf(f, "n_g_face = "PDM_FMT_G_NUM"\n", surf_ng_face);
  fscanf(f, "n_g_vtx = "PDM_FMT_G_NUM"\n",  surf_ng_vtx);

  fclose(f);




  // Read surf_face_g_num
  *surf_n_face = 0;
  *surf_face_g_num = malloc (sizeof (PDM_g_num_t) * (*surf_ng_face));
  sprintf(filename, "%s/surf_face_ln_to_gn_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (PDM_g_num_t i = 0; i < (*surf_ng_face); i++) {
    int stat = fscanf(f, PDM_FMT_G_NUM, *surf_face_g_num + i);
    if (stat == EOF) break;
    (*surf_n_face)++;
  }
  fclose(f);

  *surf_face_g_num = realloc (*surf_face_g_num, sizeof(PDM_g_num_t) * (*surf_n_face));


  // Read surf_face_vtx_idx
  *surf_face_vtx_idx = malloc (sizeof(int) * (*surf_n_face + 1));
  sprintf(filename, "%s/surf_face_vtx_idx_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i <= (*surf_n_face); i++) {
    int stat = fscanf(f, "%d", *surf_face_vtx_idx + i);
    assert (stat != EOF);
  }
  fclose(f);


  // Read surf face_vtx
  *surf_face_vtx = malloc (sizeof(int) * (*surf_face_vtx_idx)[*surf_n_face]);
  sprintf(filename, "%s/surf_face_vtx_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i < (*surf_face_vtx_idx)[*surf_n_face]; i++) {
    int stat = fscanf(f, "%d", *surf_face_vtx + i);
    assert (stat != EOF);
  }
  fclose(f);



  // Read surf_vtx_g_num
  *surf_n_vtx = 0;
  *surf_vtx_g_num = malloc (sizeof(PDM_g_num_t) * (*surf_ng_vtx));
  sprintf(filename, "%s/surf_vtx_ln_to_gn_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (PDM_g_num_t i = 0; i < (*surf_ng_vtx); i++) {
    int stat = fscanf(f, PDM_FMT_G_NUM, *surf_vtx_g_num + i);
    if (stat == EOF) break;
    (*surf_n_vtx)++;
  }

  *surf_vtx_g_num = realloc (*surf_vtx_g_num, sizeof(PDM_g_num_t) * (*surf_n_vtx));

  // Read surf vtx_coord
  *surf_vtx_coord = malloc (sizeof(double) * (*surf_n_vtx) * 3);
  sprintf(filename, "%s/surf_coords_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i < (*surf_n_vtx); i++) {
    int stat = fscanf(f, "%lf %lf %lf", *surf_vtx_coord + 3*i, *surf_vtx_coord + 3*i+1, *surf_vtx_coord + 3*i+2);
    if (stat == EOF) break;
  }
  fclose(f);



  /*
   *  Read volume mesh
   */
  // Read cell_g_num
  *n_cell = 0;
  int n_tmp = 1000;

  *cell_g_num = malloc (sizeof(PDM_g_num_t) * n_tmp);
  sprintf(filename, "%s/cell_ln_to_gn_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  while (1) {
    if (*n_cell >= n_tmp) {
      n_tmp = PDM_MAX (2*n_tmp, *n_cell);
      *cell_g_num = realloc (*cell_g_num, sizeof(PDM_g_num_t) * n_tmp);
    }
    int stat = fscanf(f, PDM_FMT_G_NUM, *cell_g_num + (*n_cell));
    if (stat == EOF) break;
    (*n_cell)++;
  }

  fclose(f);

  *cell_g_num = realloc (*cell_g_num, sizeof(PDM_g_num_t) * (*n_cell));


  // Read cell_face_idx
  *cell_face_idx = malloc (sizeof(int) * (*n_cell + 1));
  sprintf(filename, "%s/cell_face_idx_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i <= (*n_cell); i++) {
    int stat = fscanf(f, "%d", *cell_face_idx + i);
    assert (stat != EOF);
  }
  fclose(f);


  // Read cell_face
  *cell_face = malloc (sizeof(int) * (*cell_face_idx)[*n_cell]);
  sprintf(filename, "%s/cell_face_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i < (*cell_face_idx)[*n_cell]; i++) {
    int stat = fscanf(f, "%d", *cell_face + i);
    assert (stat != EOF);
  }
  fclose(f);


  // Read cell_center
  *cell_center = malloc (sizeof(double) * (*n_cell) * 3);
  sprintf(filename, "%s/cell_center_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i < (*n_cell); i++) {
    int stat = fscanf(f, "%lf %lf %lf", *cell_center + 3*i, *cell_center + 3*i+1, *cell_center + 3*i+2);
    assert (stat != EOF);
  }
  fclose(f);


  // Read face_g_num
  *n_face = 0;
  n_tmp = 1000;

  *face_g_num = malloc (sizeof(PDM_g_num_t) * n_tmp);
  sprintf(filename, "%s/face_ln_to_gn_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  while (1) {
    if (*n_face >= n_tmp) {
      n_tmp = PDM_MAX (2*n_tmp, *n_face);
      *face_g_num = realloc (*face_g_num, sizeof(PDM_g_num_t) * n_tmp);
    }
    int stat = fscanf(f, PDM_FMT_G_NUM, *face_g_num + (*n_face));
    if (stat == EOF) break;
    (*n_face)++;
  }

  fclose(f);

  *face_g_num = realloc (*face_g_num, sizeof(PDM_g_num_t) * (*n_face));


  // Read face_vtx_idx
  *face_vtx_idx = malloc (sizeof(int) * (*n_face + 1));
  sprintf(filename, "%s/face_vtx_idx_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i <= (*n_face); i++) {
    int stat = fscanf(f, "%d", *face_vtx_idx + i);
    assert (stat != EOF);
  }
  fclose(f);


  // Read face_vtx
  *face_vtx = malloc (sizeof(int) * (*face_vtx_idx)[*n_face]);
  sprintf(filename, "%s/face_vtx_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i < (*face_vtx_idx)[*n_face]; i++) {
    int stat = fscanf(f, "%d", *face_vtx + i);
    assert (stat != EOF);
  }
  fclose(f);


  // Read vtx_g_num
  *n_vtx = 0;
  n_tmp = 1000;

  *vtx_g_num = malloc (sizeof(PDM_g_num_t) * n_tmp);
  sprintf(filename, "%s/vtx_ln_to_gn_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  while (1) {
    if (*n_vtx >= n_tmp) {
      n_tmp = PDM_MAX (2*n_tmp, *n_vtx);
      *vtx_g_num = realloc (*vtx_g_num, sizeof(PDM_g_num_t) * n_tmp);
    }
    int stat = fscanf(f, PDM_FMT_G_NUM, *vtx_g_num + (*n_vtx));
    if (stat == EOF) break;
    (*n_vtx)++;
  }

  fclose(f);

  *vtx_g_num = realloc (*vtx_g_num, sizeof(PDM_g_num_t) * (*n_vtx));


  // Read vtx_coord
  *vtx_coord = malloc (sizeof(double) * (*n_vtx) * 3);
  sprintf(filename, "%s/coords_%d.txt", directory, i_rank);
  f = fopen(filename, "r");
  assert (f != NULL);

  for (int i = 0; i < (*n_vtx); i++) {
    int stat = fscanf(f, "%lf %lf %lf", *vtx_coord + 3*i, *vtx_coord + 3*i+1, *vtx_coord + 3*i+2);
    assert (stat != EOF);
  }
  fclose(f);
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int dist_function = 0;
  if (argc > 1) {
    dist_function = atoi(argv[1]);
  }

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  assert (n_rank <= 2);


  /*
   *  Read data
   */
  if (i_rank == 0) {
    printf("-- Read data\n");
    fflush(stdout);
  }
  PDM_g_num_t  surf_ng_face;
  PDM_g_num_t  surf_ng_vtx;
  int          surf_n_face;
  PDM_g_num_t *surf_face_g_num = NULL;
  int         *surf_face_vtx_idx = NULL;
  int         *surf_face_vtx = NULL;
  int          surf_n_vtx;
  PDM_g_num_t *surf_vtx_g_num = NULL;
  double      *surf_vtx_coord = NULL;
  int          n_cell;
  PDM_g_num_t *cell_g_num = NULL;
  int         *cell_face_idx = NULL;
  int         *cell_face = NULL;
  double      *cell_center = NULL;
  int          n_face;
  PDM_g_num_t *face_g_num = NULL;
  int         *face_vtx_idx = NULL;
  int         *face_vtx = NULL;
  int          n_vtx;
  PDM_g_num_t *vtx_g_num = NULL;
  double      *vtx_coord = NULL;

  _read_data (PDM_MPI_COMM_WORLD,
              &surf_ng_face,
              &surf_ng_vtx,
              &surf_n_face,
              &surf_face_g_num,
              &surf_face_vtx_idx,
              &surf_face_vtx,
              &surf_n_vtx,
              &surf_vtx_g_num,
              &surf_vtx_coord,
              &n_cell,
              &cell_g_num,
              &cell_face_idx,
              &cell_face,
              &cell_center,
              &n_face,
              &face_g_num,
              &face_vtx_idx,
              &face_vtx,
              &n_vtx,
              &vtx_g_num,
              &vtx_coord);

  int n_part = 1;

  PDM_g_num_t n_g_vol_cell = -1;
  PDM_g_num_t n_g_vol_face = -1;
  PDM_g_num_t n_g_vol_vtx = -1;
  for (int i = 0; i < n_face; i++) {
    n_g_vol_face = PDM_MAX(face_g_num[i], n_g_vol_face);
  }

  for (int i = 0; i < n_vtx; i++) {
    n_g_vol_vtx = PDM_MAX(vtx_g_num[i], n_g_vol_vtx);
  }

  for (int i = 0; i < n_cell; i++) {
    n_g_vol_cell = PDM_MAX(cell_g_num[i], n_g_vol_cell);
  }


  /*
   *  Create dist structure
   */
  int n_point_cloud = 1;
  int id_dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                            n_point_cloud,
                                            PDM_MPI_COMM_WORLD,
                                            PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (id_dist,
                                                 surf_ng_face,
                                                 surf_ng_vtx,
                                                 n_part);

  PDM_dist_cloud_surf_n_part_cloud_set (id_dist,
                                        0,
                                        n_part);



  PDM_dist_cloud_surf_surf_mesh_part_set (id_dist,
                                          0,
                                          surf_n_face,
                                          surf_face_vtx_idx,
                                          surf_face_vtx,
                                          surf_face_g_num,
                                          surf_n_vtx,
                                          surf_vtx_coord,
                                          surf_vtx_g_num);

  PDM_dist_cloud_surf_cloud_set (id_dist,
                                 0,
                                 0,
                                 n_cell,
                                 cell_center,
                                 cell_g_num);

  /*
   *  Compute distance
   */
  if (i_rank == 0) {
    printf("-- Dist compute (function #%d)\n", dist_function);
    fflush(stdout);
  }

  if (dist_function == 2) {
    PDM_dist_cloud_surf_compute2 (id_dist);
  }
  else if (dist_function == 2) {
    PDM_dist_cloud_surf_compute3 (id_dist);
  }
  else {
    PDM_dist_cloud_surf_compute (id_dist);
  }
  PDM_dist_cloud_surf_dump_times (id_dist);

  /*
   *  Free memory
   */
  PDM_dist_cloud_surf_free (id_dist);

  free (surf_face_g_num);
  free (surf_face_vtx_idx);
  free (surf_face_vtx);
  free (surf_vtx_g_num);
  free (surf_vtx_coord);
  free (cell_g_num);
  free (cell_face_idx);
  free (cell_face);
  free (cell_center);
  free (face_g_num);
  free (face_vtx_idx);
  free (face_vtx);
  free (vtx_g_num);
  free (vtx_coord);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
