#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_dmesh.h"
#include "pdm_writer.h"
#include "pdm_distrib.h"

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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -hilbert         Call PT-Hilbert.\n\n"
     "  -post            Write output in Ensight format.\n\n"
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
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch, 3 Hilbert)
 * \param [inout]   post     Write output or not
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           int           *n_part,
           int           *post,
           int           *method
           )
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
    else if (strcmp(argv[i], "-implicit") == 0) {
      *method = 4;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = 3;
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

  /* Set default values */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
  int                n_zone    = 1;

  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /* Read args */

  _read_args(argc, argv, &n_vtx_seg, &n_part, &post, (int *) &method);

  /* Init */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Alloc for distributed mesh */
  int *dn_cell       = (int *) malloc(n_zone * sizeof(int));
  int *dn_face       = (int *) malloc(n_zone * sizeof(int));
  int *dn_vtx        = (int *) malloc(n_zone * sizeof(int));
  int *n_face_group  = (int *) malloc(n_zone * sizeof(int));
  int *dface_vtx_s   = (int *) malloc(n_zone * sizeof(int));
  int *dface_group_s = (int *) malloc(n_zone * sizeof(int));

  PDM_g_num_t  **dface_cell      = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t*));
  int          **dface_vtx_idx   = (int **)         malloc(n_zone * sizeof(int*));
  PDM_g_num_t  **dface_vtx       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t*));
  double       **dvtx_coord      = (double **)      malloc(n_zone * sizeof(double*));
  int          **dface_group_idx = (int **)         malloc(n_zone * sizeof(int*));
  PDM_g_num_t  **dface_group     = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t*));
  int          **dface_bnd_idx   = (int **)         malloc(n_zone * sizeof(int*));
  PDM_g_num_t  **dface_bnd       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t*));
  int          **dface_join_idx  = (int **)         malloc(n_zone * sizeof(int*));
  PDM_g_num_t  **dface_join      = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t*));
  int          **djoins_ids      = (int **)         malloc(n_zone * sizeof(int*));

  /* Initialize multipart */
  int *n_part_zones  = (int *) malloc(n_zone * sizeof(int));
  for (int i_zone = 0; i_zone < n_zone; i_zone++){
    n_part_zones[i_zone] = n_part;
  }
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone, n_part_zones, PDM_FALSE,
						   method, PDM_PART_SIZE_HOMOGENEOUS,
						   NULL, comm, PDM_OWNERSHIP_KEEP);
  PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_CUTHILL", NULL, "PDM_PART_RENUM_FACE_LEXICOGRAPHIC");
  if (n_zone > 1)
    PDM_multipart_set_reordering_options(mpart,  1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");

  /* Generate mesh */
  PDM_dcube_t **dcube = (PDM_dcube_t **) malloc(n_zone * sizeof(PDM_dcube_t *));
  PDM_dmesh_t **dmesh = (PDM_dmesh_t **) malloc(n_zone * sizeof(PDM_dmesh_t *));
  for (int i_zone = 0; i_zone < n_zone; i_zone++)
  {
    // Create a cube for this zone
    if (i_rank == 0) PDM_printf("Creating dcube for zone %d\n", i_zone);

    dcube[i_zone] = PDM_dcube_gen_init(comm, n_vtx_seg, length, i_zone, 0., 0., PDM_OWNERSHIP_KEEP);
    PDM_dcube_gen_dim_get(dcube[i_zone],
                         &n_face_group[i_zone],
                         &dn_cell[i_zone],
                         &dn_face[i_zone],
                         &dn_vtx[i_zone],
                         &dface_vtx_s[i_zone],
                         &dface_group_s[i_zone]);
    PDM_dcube_gen_data_get(dcube[i_zone],
                          &dface_cell[i_zone],
                          &dface_vtx_idx[i_zone],
                          &dface_vtx[i_zone],
                          &dvtx_coord[i_zone],
                          &dface_group_idx[i_zone],
                          &dface_group[i_zone]);
    /* Les faces groups du dcube sont : zmin, zmax, xmin, xmax, ymin, ymax
    Il faut les séparer en faces de bords et faces raccord, sachant que
    les zones sont alignées selon X */
    int n_bnd = 4;
    int n_jn  = 2;
    if (i_zone == 0){
      n_bnd++;
      n_jn-- ;
    }
    if (i_zone == n_zone-1){
      n_bnd++;
      n_jn-- ;
    }

    // Join numbering (left to right, increasing i_zone)
    printf("n_jn = %i \n", n_jn);
    if(n_jn > 0 ) {
      djoins_ids[i_zone] = (int *) malloc(n_jn * sizeof(int));
      if (i_zone == 0)
        djoins_ids[i_zone][0] = 0;
      else if (i_zone == n_zone-1)
        djoins_ids[i_zone][0] = 2*i_zone - 1;
      else {
        djoins_ids[i_zone][0] = 2*i_zone - 1;
        djoins_ids[i_zone][1] = 2*i_zone;
      }
    }

    dface_bnd_idx[i_zone]  = (int *) malloc((n_bnd+1) * sizeof(int));
    dface_join_idx[i_zone] = (int *) malloc((n_jn +1) * sizeof(int));
    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = 1;
    dface_bnd_idx[i_zone][0]  = 0;
    dface_join_idx[i_zone][0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++)
    {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      int group_size = dface_group_idx[i_zone][igroup+1] - dface_group_idx[i_zone][igroup];
      if (copy_to_bnd) //Its a boundary
        dface_bnd_idx[i_zone][i_bnd++] = group_size;
      else { //Its a join
        dface_join_idx[i_zone][i_jn++] = group_size;
      }
    }
    for (int i = 0; i < n_bnd; i++) {
      dface_bnd_idx[i_zone][i+1] = dface_bnd_idx[i_zone][i+1] + dface_bnd_idx[i_zone][i];
    }
    for (int i = 0; i < n_jn; i++) {
      dface_join_idx[i_zone][i+1] = dface_join_idx[i_zone][i+1] + dface_join_idx[i_zone][i];
    }

    // Second pass to copy
    dface_bnd [i_zone] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_zone][n_bnd] * sizeof(PDM_g_num_t));
    dface_join[i_zone] = (PDM_g_num_t *) malloc(dface_join_idx[i_zone][n_jn ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      if (copy_to_bnd){ //Its a boundary
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++)
          dface_bnd[i_zone][i_bnd++] = dface_group[i_zone][i];
      } else { //Its a join
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++)
          dface_join[i_zone][i_jn++] = dface_group[i_zone][i];
      }
    }

    // Store it in dmesh struct
    dmesh[i_zone] = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell[i_zone],
                                     dn_face[i_zone],
                                     0, // dn_edge
                                     dn_vtx[i_zone],
                                     comm);

    PDM_dmesh_vtx_coord_set(dmesh[i_zone],
                            dvtx_coord[i_zone],
                            PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh[i_zone],
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx[i_zone],
                               dface_vtx_idx[i_zone],
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh[i_zone],
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell[i_zone],
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh[i_zone],
                        PDM_BOUND_TYPE_FACE,
                        n_bnd,
                        dface_bnd[i_zone],
                        dface_bnd_idx[i_zone],
                        PDM_OWNERSHIP_USER);

    // PDM_dmesh_set(dmesh[i_zone],
    //               dvtx_coord[i_zone],
    //               dface_vtx_idx[i_zone],
    //               dface_vtx[i_zone],
    //               dface_cell[i_zone],
    //               dface_bnd_idx[i_zone],
    //               dface_bnd[i_zone],
    //               djoins_ids[i_zone],
    //               dface_join_idx[i_zone],
    //               dface_join[i_zone]);
    PDM_multipart_register_block(mpart, i_zone, dmesh[i_zone]);
  }

  /* Connection between zones */
  int n_total_joins = 2*(n_zone-1);
  int *join_to_opposite = (int *) malloc(n_total_joins*sizeof(int));
  for (int ijoin = 0; ijoin < n_total_joins; ijoin++){
    if (ijoin % 2 == 0)
      join_to_opposite[ijoin] = ijoin + 1;
    else
      join_to_opposite[ijoin] = ijoin - 1;
  }
  PDM_multipart_register_joins(mpart, n_total_joins, join_to_opposite);

  /* Run */
  PDM_multipart_compute(mpart);
  if (i_rank==0)
    PDM_printf("Partitioning done !\n");



  if (post == 1)
  {
    /* Prepare writer */
    int *geom_ids = (int *) malloc(n_zone * sizeof(int));
    PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                            PDM_WRITER_FMT_ASCII,
                                            PDM_WRITER_TOPO_CST,
                                            PDM_WRITER_OFF,
                                            "test_mpart_cube",
                                            "mpart",
                                            PDM_MPI_COMM_WORLD,
                                            PDM_IO_KIND_MPI_SIMPLE,
                                            1.,
                                            NULL);
    for (int i_zone = 0; i_zone < n_zone; i_zone++){
      geom_ids[i_zone] = PDM_writer_geom_create(id_cs,
                                                "mesh",
                                                n_part_zones[i_zone]); //total nb of part for this proc/zone
    }
    // Cell local id
    int id_var_cell_id = PDM_writer_var_create(id_cs, PDM_WRITER_OFF, PDM_WRITER_VAR_SCALAR, PDM_WRITER_VAR_ELEMENTS, "cell_id");
    // Global partition Id (ordred by proc / zone), staring at 1
    int id_var_gpart_id = PDM_writer_var_create(id_cs, PDM_WRITER_OFF, PDM_WRITER_VAR_SCALAR, PDM_WRITER_VAR_ELEMENTS, "gpart_id");
    // Local partition Id on the proc / zone, starting at 0
    int id_var_lpart_id = PDM_writer_var_create(id_cs, PDM_WRITER_OFF, PDM_WRITER_VAR_SCALAR, PDM_WRITER_VAR_ELEMENTS, "lpart_id");
    // Proc Id
    int id_var_proc_id = PDM_writer_var_create(id_cs, PDM_WRITER_OFF, PDM_WRITER_VAR_SCALAR, PDM_WRITER_VAR_ELEMENTS, "iproc");
    // Id of opposite proc
    int id_var_opp_proc_id = PDM_writer_var_create(id_cs, PDM_WRITER_OFF, PDM_WRITER_VAR_SCALAR, PDM_WRITER_VAR_VERTICES, "oppProc");
    // Id of opposite part (in the local numerotation of the opposite proc)
    int id_var_opp_part_id = PDM_writer_var_create(id_cs, PDM_WRITER_OFF, PDM_WRITER_VAR_SCALAR, PDM_WRITER_VAR_VERTICES, "oppPart");
    PDM_writer_step_beg(id_cs, 0.);

    /* Alloc for part meshes */
    int tn_part_proc = 0;
    for (int i_zone = 0; i_zone < n_zone; i_zone++)
      tn_part_proc += n_part_zones[i_zone];

    /* Get results */
    int **pface_vtx_n  = (int **) malloc(tn_part_proc * sizeof(int *));
    int **pcell_face_n = (int **) malloc(tn_part_proc * sizeof(int *));
    int  *pn_cell      = (int  *) malloc(tn_part_proc * sizeof(int  ));
    int  *pn_vtx       = (int  *) malloc(tn_part_proc * sizeof(int  ));
    int **comm_visu    = (int **) malloc(tn_part_proc * sizeof(int *));

    int ipartzone = 0;
    /* Write geometry */
    for (int i_zone = 0; i_zone < n_zone; i_zone++){
      for (int i_part = 0; i_part < n_part_zones[i_zone]; i_part++){
        if (i_rank==0) PDM_printf("Get zone %i part %i\n", i_zone, i_part);

        PDM_g_num_t* cell_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL,
                                        &cell_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        int *cell_face     = NULL;
        int *cell_face_idx = NULL;
        int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                         &cell_face,
                                                         &cell_face_idx,
                                                         PDM_OWNERSHIP_KEEP);


        int *face_cell     = NULL;
        int *face_cell_idx = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                            &face_cell,
                                            &face_cell_idx,
                                            PDM_OWNERSHIP_KEEP);
        assert(face_cell_idx == NULL);

        int *face_vtx     = NULL;
        int *face_vtx_idx = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &face_vtx,
                                            &face_vtx_idx,
                                            PDM_OWNERSHIP_KEEP);

        int *face_edge     = NULL;
        int *face_edge_idx = NULL;
        int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                         &face_edge,
                                                         &face_edge_idx,
                                                         PDM_OWNERSHIP_KEEP);

        PDM_g_num_t* face_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                         &edge_vtx,
                                                         &edge_vtx_idx,
                                                         PDM_OWNERSHIP_KEEP);
        assert(edge_vtx_idx == NULL);
        PDM_g_num_t* edge_ln_to_gn = NULL;
        int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      0,
                                                      i_part,
                                                      PDM_MESH_ENTITY_EDGE,
                                                      &edge_ln_to_gn,
                                                      PDM_OWNERSHIP_KEEP);
        assert(n_edge2 == n_edge);

        PDM_g_num_t* vtx_ln_to_gn = NULL;
        int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    0,
                                                    i_part,
                                                    PDM_MESH_ENTITY_VERTEX,
                                                    &vtx_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);
        int          pn_face_group        = 0;
        int         *face_group_idx      = NULL;
        int         *face_group          = NULL;
        PDM_g_num_t *face_group_ln_to_gn = NULL;
        PDM_multipart_bound_get(mpart,
                                0,
                                i_part,
                                PDM_BOUND_TYPE_FACE,
                                &pn_face_group,
                                &face_group_idx,
                                &face_group,
                                &face_group_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);

        int *vtx_part_bound_proc_idx = NULL;
        int *vtx_part_bound_part_idx = NULL;
        int *vtx_part_bound          = NULL;
        PDM_multipart_part_graph_comm_get(mpart,
                                          0,
                                          i_part,
                                          PDM_BOUND_TYPE_VTX,
                                          &vtx_part_bound_proc_idx,
                                          &vtx_part_bound_part_idx,
                                          &vtx_part_bound,
                                          PDM_OWNERSHIP_KEEP);

        int *face_part_bound_proc_idx = NULL;
        int *face_part_bound_part_idx = NULL;
        int *face_part_bound          = NULL;
        PDM_multipart_part_graph_comm_get(mpart,
                                          0,
                                          i_part,
                                          PDM_BOUND_TYPE_FACE,
                                          &face_part_bound_proc_idx,
                                          &face_part_bound_part_idx,
                                          &face_part_bound,
                                          PDM_OWNERSHIP_KEEP);
        double *vtx = NULL;
        PDM_multipart_part_vtx_coord_get(mpart,
                                         0,
                                         i_part,
                                         &vtx,
                                         PDM_OWNERSHIP_KEEP);

        // Store opposite proc & part for rendering
        comm_visu[ipartzone] = (int *) malloc(2*n_vtx * sizeof(int));
        for (int i=0; i < 2*n_vtx; i++)  {
          comm_visu[ipartzone][i] = -1;
        }

        // for (int ijoin = 0; ijoin < n_joins; ijoin++){
        //   for (int iface = face_join_idx[ijoin]; iface < face_join_idx[ijoin+1]; iface++){
        //     int face_lid = face_join[4*iface] - 1;
        //     for (int ivtx = face_vtx_idx[face_lid]; ivtx < face_vtx_idx[face_lid+1]; ivtx++){
        //       int vtx_lid = face_vtx[ivtx] - 1;
        //       comm_visu[ipartzone][2*vtx_lid]   = face_join[4*iface+1];
        //       comm_visu[ipartzone][2*vtx_lid+1] = face_join[4*iface+2];
        //     }
        //   }
        // }
        // Also do it for internal partitions
        for (int p = 0; p < n_rank; p++){
          for (int iface = face_part_bound_proc_idx[p]; iface < face_part_bound_proc_idx[p+1]; iface++){
            int face_lid = face_part_bound[4*iface]-1;
            for (int ivtx = face_vtx_idx[face_lid]; ivtx < face_vtx_idx[face_lid+1]; ivtx++){
              int vtx_lid = face_vtx[ivtx] - 1;
              comm_visu[ipartzone][2*vtx_lid]   = face_part_bound[4*iface+1];
              comm_visu[ipartzone][2*vtx_lid+1] = face_part_bound[4*iface+2]-1;
            }
          }
        }

        pface_vtx_n[ipartzone]  = (int *) malloc(sizeof(int) * n_face);
        pcell_face_n[ipartzone] = (int *) malloc(sizeof(int) * n_cell);
        for (int i = 0; i < n_cell; i++)
          pcell_face_n[ipartzone][i] = cell_face_idx[i+1] - cell_face_idx[i];
        for (int i = 0; i < n_face; i++)
          pface_vtx_n[ipartzone][i] = face_vtx_idx[i+1] - face_vtx_idx[i];

        PDM_writer_geom_coord_set(id_cs,
                                  geom_ids[i_zone],
                                  i_part,
                                  n_vtx,
                                  vtx,
                                  vtx_ln_to_gn,
                                  PDM_OWNERSHIP_USER);
        
        PDM_writer_geom_cell3d_cellface_add(id_cs,
                                            geom_ids[i_zone],
                                            i_part,
                                            n_cell,
                                            n_face,
                                            face_vtx_idx,
                                            pface_vtx_n[ipartzone],
                                            face_vtx,
                                            cell_face_idx,
                                            pcell_face_n[ipartzone],
                                            cell_face,
                                            cell_ln_to_gn);
        pn_cell[ipartzone] = n_cell;
        pn_vtx [ipartzone] = n_vtx;
        ipartzone++;
      }
    if (i_rank==0) PDM_printf("Write geometry for zone %i\n", i_zone);
    PDM_writer_geom_write(id_cs, geom_ids[i_zone]);
    }

    /* Write data */
    PDM_g_num_t *partzoneshift = PDM_compute_entity_distribution(comm, tn_part_proc);

    ipartzone = 0;
    for (int i_zone = 0; i_zone < n_zone; i_zone++){
      for (int i_part = 0; i_part < n_part_zones[i_zone]; i_part++){
        PDM_real_t *val_cell_id     = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_cell[ipartzone]);
        PDM_real_t *val_gpart_id    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_cell[ipartzone]);
        PDM_real_t *val_lpart_id    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_cell[ipartzone]);
        PDM_real_t *val_proc_id     = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_cell[ipartzone]);
        PDM_real_t *val_opp_proc_id = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_vtx [ipartzone]);
        PDM_real_t *val_opp_part_id = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_vtx [ipartzone]);
        for (int i=0; i < pn_cell[ipartzone]; i++) {
          val_cell_id[i] = i;
          val_gpart_id[i] = (PDM_real_t) (partzoneshift[i_rank] + ipartzone);
          val_lpart_id[i] = i_part;
          val_proc_id[i]  = i_rank;
        }
        for (int i=0; i < pn_vtx[ipartzone]; i++){
          val_opp_proc_id[i] = comm_visu[ipartzone][2*i  ];
          val_opp_part_id[i] = comm_visu[ipartzone][2*i+1];
        }
        PDM_writer_var_set(id_cs, id_var_cell_id    , geom_ids[i_zone], i_part, val_cell_id    );
        PDM_writer_var_set(id_cs, id_var_gpart_id   , geom_ids[i_zone], i_part, val_gpart_id   );
        PDM_writer_var_set(id_cs, id_var_lpart_id   , geom_ids[i_zone], i_part, val_lpart_id   );
        PDM_writer_var_set(id_cs, id_var_proc_id    , geom_ids[i_zone], i_part, val_proc_id    );
        PDM_writer_var_set(id_cs, id_var_opp_proc_id, geom_ids[i_zone], i_part, val_opp_proc_id);
        PDM_writer_var_set(id_cs, id_var_opp_part_id, geom_ids[i_zone], i_part, val_opp_part_id);

        free(val_cell_id    );
        free(val_gpart_id   );
        free(val_lpart_id   );
        free(val_proc_id    );
        free(val_opp_proc_id);
        free(val_opp_part_id);

        ipartzone++;
      }
    }
    free(partzoneshift);

    if (i_rank==0) PDM_printf("Write variables\n");
    PDM_writer_var_write(id_cs, id_var_cell_id    );
    PDM_writer_var_free (id_cs, id_var_cell_id    );
    PDM_writer_var_write(id_cs, id_var_gpart_id   );
    PDM_writer_var_free (id_cs, id_var_gpart_id   );
    PDM_writer_var_write(id_cs, id_var_lpart_id   );
    PDM_writer_var_free (id_cs, id_var_lpart_id   );
    PDM_writer_var_write(id_cs, id_var_proc_id    );
    PDM_writer_var_free (id_cs, id_var_proc_id    );
    PDM_writer_var_write(id_cs, id_var_opp_proc_id);
    PDM_writer_var_free (id_cs, id_var_opp_proc_id);
    PDM_writer_var_write(id_cs, id_var_opp_part_id);
    PDM_writer_var_free (id_cs, id_var_opp_part_id);

    PDM_writer_step_end(id_cs);

    for (int i_zone = 0; i_zone < n_zone; i_zone++){
      PDM_writer_geom_data_free(id_cs, geom_ids[i_zone]);
      PDM_writer_geom_free(id_cs, geom_ids[i_zone]);
    }
    free(geom_ids);
    PDM_writer_free(id_cs);
    for (int i_part = 0; i_part < tn_part_proc; i_part++){
      free(comm_visu[i_part]);
      free(pface_vtx_n[i_part]);
      free(pcell_face_n[i_part]);
    }
    free(comm_visu);
    free(pface_vtx_n);
    free(pcell_face_n);
    free(pn_cell);
    free(pn_vtx);
    if (i_rank==0) PDM_printf("Post treatment completed\n");
  }

  /* Free memory */
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(dface_bnd_idx[i_zone]);
    free(dface_bnd[i_zone]);
    free(dface_join_idx[i_zone]);
    free(dface_join[i_zone]);
    if(n_zone > 1) {
      free(djoins_ids[i_zone]);
    }
    PDM_dmesh_free(dmesh[i_zone]);
    PDM_dcube_gen_free(dcube[i_zone]);
  }
  free(dcube);
  free(dn_cell);
  free(dn_face);
  free(dn_vtx);
  free(n_face_group);
  free(dface_group_s);
  free(dface_vtx_s);
  free(dface_cell);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dvtx_coord);
  free(dface_group_idx);
  free(dface_group);
  free(dface_bnd_idx);
  free(dface_bnd);
  free(dface_join_idx);
  free(dface_join);
  free(djoins_ids);
  free(join_to_opposite);

  PDM_multipart_free(mpart);
  free(dmesh);
  free(n_part_zones);

  if (i_rank==0) PDM_printf("pdm_t_multipart run finalized\n");
  PDM_MPI_Finalize();
  return 0;
}
