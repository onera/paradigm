#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_extract_part.h"
#include "pdm_priv.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_array.h"
#include "pdm_generate_mesh.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_multipart.h"
#include "pdm_distrib.h"


static PDM_part_mesh_nodal_t *
_generate_mesh
(
 PDM_MPI_Comm         comm,
 PDM_Mesh_nodal_elt_t elt_type,
 int                  n_part
 )
{
  PDM_g_num_t      n_vtx_seg   = 5;
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  if (elt_type == PDM_MESH_NODAL_POLY_3D) {
    // Polyhedral mesh
    PDM_g_num_t ng_cell      = 0;
    PDM_g_num_t ng_face      = 0;
    PDM_g_num_t ng_vtx       = 0;
    int         dn_cell      = 0;
    int         dn_face      = 0;
    int         dn_edge      = 0;
    int         dn_vtx       = 0;
    int         n_face_group = 0;

    double      *dvtx_coord       = NULL;
    int         *dcell_face_idx   = NULL;
    PDM_g_num_t *dcell_face       = NULL;
    PDM_g_num_t *dface_cell       = NULL;
    int         *dface_vtx_idx    = NULL;
    PDM_g_num_t *dface_vtx        = NULL;
    int         *dface_group_idx  = NULL;
    PDM_g_num_t *dface_group      = NULL;

    PDM_poly_vol_gen(comm,
                     0.,
                     0.,
                     0.,
                     1.,
                     1.,
                     1.,
                     n_vtx_seg,
                     n_vtx_seg,
                     n_vtx_seg,
                     0,
                     0,
                     &ng_cell,
                     &ng_face,
                     &ng_vtx,
                     &n_face_group,
                     &dn_cell,
                     &dn_face,
                     &dn_vtx,
                     &dcell_face_idx,
                     &dcell_face,
                     &dface_cell,
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dvtx_coord,
                     &dface_group_idx,
                     &dface_group);

    /* Generate dmesh */
    PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                          dn_cell,
                                          dn_face,
                                          dn_edge,
                                          dn_vtx,
                                          comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);


    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_KEEP);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_KEEP);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               dcell_face,
                               dcell_face_idx,
                               PDM_OWNERSHIP_KEEP);



    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_KEEP);

    PDM_dmesh_compute_distributions(dmesh);

    int n_domain = 1;
    PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                  &n_part,
                                                  PDM_FALSE,
                                                  part_method,
                                                  PDM_PART_SIZE_HOMOGENEOUS,
                                                  NULL,
                                                  comm,
                                                  PDM_OWNERSHIP_KEEP);

    PDM_multipart_dmesh_set(mpart, 0, dmesh);

    PDM_multipart_compute(mpart);

    PDM_dmesh_free(dmesh);

    PDM_part_mesh_nodal_t *pmn = PDM_part_mesh_nodal_create(3, n_part, comm);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &cell_face_idx,
                                          &cell_face,
                                          PDM_OWNERSHIP_KEEP);

      int *face_vtx_idx = NULL;
      int *face_vtx     = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx_idx,
                                          &face_vtx,
                                          PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *cell_ln_to_gn = NULL;
      int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                                   0,
                                                   i_part,
                                                   PDM_MESH_ENTITY_CELL,
                                                   &cell_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *face_ln_to_gn = NULL;
      int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                   0,
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE,
                                                   &face_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  0,
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_USER);

      double *vtx_coord = NULL;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       0,
                                       i_part,
                                       &vtx_coord,
                                       PDM_OWNERSHIP_USER);

      PDM_part_mesh_nodal_coord_set(pmn,
                                    i_part,
                                    n_vtx,
                                    vtx_coord,
                                    vtx_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

      PDM_part_mesh_nodal_cell3d_cellface_add(pmn,
                                              i_part,
                                              n_cell,
                                              n_face,
                                              face_vtx_idx,
                                              face_vtx,
                                              face_ln_to_gn,
                                              cell_face_idx,
                                              cell_face,
                                              cell_ln_to_gn,
                                              PDM_OWNERSHIP_KEEP);
    }
    PDM_multipart_free(mpart);

    return pmn;
  }
  else {
    int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(elt_type);

    if (mesh_dimension == 2) {
      return PDM_generate_mesh_rectangle(comm,
                                         elt_type,
                                         1,
                                         NULL,
                                         0.,
                                         0.,
                                         0.,
                                         1.,
                                         1.,
                                         n_vtx_seg,
                                         n_vtx_seg,
                                         n_part,
                                         part_method);
    }
    else {
      return PDM_generate_mesh_parallelepiped(comm,
                                              elt_type,
                                              1,
                                              NULL,
                                              0.,
                                              0.,
                                              0.,
                                              1.,
                                              1.,
                                              1.,
                                              n_vtx_seg,
                                              n_vtx_seg,
                                              n_vtx_seg,
                                              n_part,
                                              part_method);
    }
  }
}

static void
_get_extracted
(
 PDM_MPI_Comm                   comm,
 PDM_part_mesh_nodal_elmts_t   *pmne,
 int                            n_part,
 int                            empty_extraction,
 int                          **out_n_extract,
 int                         ***out_extract_lnum
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int  *n_extract    = PDM_array_zeros_int(n_part);
  int **extract_lnum = NULL;
  PDM_malloc(extract_lnum, n_part, int *);

  if (empty_extraction) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      n_extract   [i_part] = 0;
      extract_lnum[i_part] = NULL;
    }
    *out_n_extract    = n_extract;
    *out_extract_lnum = extract_lnum;
    return;
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);

      PDM_g_num_t *ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                  sections_id[i_section],
                                                                  i_part,
                                                                  PDM_OWNERSHIP_KEEP);
      for (int i_elt = 0; i_elt < n_elt; i_elt++) {
        if (ln_to_gn[i_elt]%5 == 0) {
          n_extract[i_part]++;
        }
      }
    }


    PDM_malloc(extract_lnum[i_part], n_extract[i_part], int);

    n_extract[i_part] = 0;


    int i_cell = -1;
    for (int i_section = 0; i_section < n_section; i_section++) {

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);


      PDM_g_num_t *ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                  sections_id[i_section],
                                                                  i_part,
                                                                  PDM_OWNERSHIP_KEEP);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 sections_id[i_section],
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      for (int i_elt = 0; i_elt < n_elt; i_elt++) {
        if (ln_to_gn[i_elt]%5 == 0) {
          if (parent_num == NULL) {
            i_cell++;
          }
          else {
            i_cell = parent_num[i_elt];
          }
          extract_lnum[i_part][n_extract[i_part]++] = i_cell + 1;
        }
      }
    }

    // PDM_log_trace_array_int(extract_lnum[i_part], n_extract[i_part], "extract_lnum : ");
  }

  *out_n_extract    = n_extract;
  *out_extract_lnum = extract_lnum;
}


static void
_get_targets
(
 PDM_MPI_Comm                   comm,
 PDM_part_mesh_nodal_elmts_t   *pmne,
 int                            n_part_in,
 int                            n_part_out,
 int                            empty_extraction,
 int                          **n_target,
 PDM_g_num_t                 ***target_gnum
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_malloc(*n_target,    n_part_out, int          );
  PDM_malloc(*target_gnum, n_part_out, PDM_g_num_t *);

  if (empty_extraction) {
    for (int i_part = 0; i_part < n_part_out; i_part++) {
      (*n_target   )[i_part] = 0;
      (*target_gnum)[i_part] = NULL;
    }
    return;
  }

  PDM_g_num_t lmax_gnum = -1;

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_part = 0; i_part < n_part_in; i_part++) {
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);

      PDM_g_num_t *ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                  sections_id[i_section],
                                                                  i_part,
                                                                  PDM_OWNERSHIP_KEEP);
      for (int i_elt = 0; i_elt < n_elt; i_elt++) {
        lmax_gnum = PDM_MAX(lmax_gnum, ln_to_gn[i_elt]);
      }
    }
  }


  PDM_g_num_t max_gnum = -1;
  PDM_MPI_Allreduce(&lmax_gnum, &max_gnum, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);


  PDM_g_num_t gn_extract = PDM_MAX(1, max_gnum/5);

  PDM_g_num_t *distrib = PDM_compute_uniform_entity_distribution(comm, gn_extract);


  int n_target_tot = distrib[i_rank+1] - distrib[i_rank];

  int div = n_target_tot / n_part_out;
  int rem = n_target_tot % n_part_out;
  int idx = 0;
  for (int i_part = 0; i_part < n_part_out; i_part++) {
    (*n_target)[i_part] = div;
    if (i_part < rem) {
      (*n_target)[i_part]++;
    }

    PDM_malloc((*target_gnum)[i_part], (*n_target)[i_part], PDM_g_num_t);
    for (int i = 0; i < (*n_target)[i_part]; i++) {
      (*target_gnum)[i_part][i] = 5*(distrib[i_rank] + idx + 1);
      idx++;
    }
  }

  PDM_free(distrib);
}


static void
_check_groups
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                          n_part_in,
 PDM_part_mesh_nodal_elmts_t *extract_pmne,
 int                          n_part_out
)
{
  if (pmne == NULL) {
    return;
  }
  PDM_UNUSED(n_part_in);

  int n_group         = PDM_part_mesh_nodal_elmts_n_group_get(pmne);
  int extract_n_group = PDM_part_mesh_nodal_elmts_n_group_get(extract_pmne);
  CHECK(n_group == extract_n_group);

  if (n_group == 0) return;

  /* Check that all elements are assigned a group */
  for (int i_part = 0; i_part < n_part_out; i_part++) {

    int n_elt_tot = PDM_part_mesh_nodal_elmts_n_elmts_get(extract_pmne, i_part);

    int n_elt_in_all_groups = 0;
    for (int i_group = 0; i_group < n_group; i_group++) {
      int          n_group_elmt   = 0;
      int         *group_elmt     = NULL;
      PDM_g_num_t *group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(extract_pmne,
                                          i_part,
                                          i_group,
                                          &n_group_elmt,
                                          &group_elmt,
                                          &group_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);

      // log_trace("group %d : ", i_group);
      // PDM_log_trace_array_int(group_elmt, n_group_elmt, "group_elmt : ");

      n_elt_in_all_groups += n_group_elmt;
    }

    // log_trace("n_elt_in_all_groups = %d / %d\n", n_elt_in_all_groups, n_elt_tot);
    CHECK(n_elt_in_all_groups == n_elt_tot);
  }
}

// mpirun -n 2 pdm_t_unit_test -tc="[pdm_extract_part] - 2p - nodal"


MPI_TEST_CASE("[pdm_extract_part] - 2p - nodal", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /* Generate mesh */
  PDM_Mesh_nodal_elt_t    elt_type         = PDM_MESH_NODAL_PRISM6;
  int                     n_part_in        = 1;
  int                     n_part_out       = 1;
  PDM_split_dual_t        part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_extract_part_kind_t extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
  int                     empty_extraction = 0;


  /* Define subcases */
  SUBCASE("LOCAL, n_part_in 1, n_part_out 1, POLY_2D") {
    if (i_rank == 0) printf("LOCAL, n_part_in 1, n_part_out 1, POLY_2D\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
  }
  SUBCASE("LOCAL, n_part_in 1, n_part_out 1, POLY_2D, empty extraction") {
    if (i_rank == 0) printf("LOCAL, n_part_in 1, n_part_out 1, POLY_2D, empty extraction\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
  }
  SUBCASE("LOCAL, n_part_in 2, n_part_out 2, POLY_2D") {
    if (i_rank == 0) printf("LOCAL, n_part_in 2, n_part_out 2, POLY_2D\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
  }
  SUBCASE("LOCAL, n_part_in 2, n_part_out 2, POLY_2D, empty extraction") {
    if (i_rank == 0) printf("LOCAL, n_part_in 2, n_part_out 2, POLY_2D, empty extraction\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
  }
  SUBCASE("LOCAL, n_part_in 1, n_part_out 1, PRISM6") {
    if (i_rank == 0) printf("LOCAL, n_part_in 1, n_part_out 1, PRISM6\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
  }
  SUBCASE("LOCAL, n_part_in 1, n_part_out 1, PRISM6, empty extraction") {
    if (i_rank == 0) printf("LOCAL, n_part_in 1, n_part_out 1, PRISM6, empty extraction\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
  }
  SUBCASE("LOCAL, n_part_in 2, n_part_out 2, PRISM6") {
    if (i_rank == 0) printf("LOCAL, n_part_in 2, n_part_out 2, PRISM6\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
  }
  SUBCASE("LOCAL, n_part_in 2, n_part_out 2, PRISM6, empty extraction") {
    if (i_rank == 0) printf("LOCAL, n_part_in 2, n_part_out 2, PRISM6, empty extraction\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
  }
  SUBCASE("LOCAL, n_part_in 1, n_part_out 1, POLY_3D") {
    if (i_rank == 0) printf("LOCAL, n_part_in 1, n_part_out 1, POLY_3D\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
  }
  SUBCASE("LOCAL, n_part_in 1, n_part_out 1, POLY_3D, empty extraction") {
    if (i_rank == 0) printf("LOCAL, n_part_in 1, n_part_out 1, POLY_3D, empty extraction\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
  }
  SUBCASE("LOCAL, n_part_in 2, n_part_out 2, POLY_3D") {
    if (i_rank == 0) printf("LOCAL, n_part_in 2, n_part_out 2, POLY_3D\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
  }
  SUBCASE("LOCAL, n_part_in 2, n_part_out 2, POLY_3D, empty extraction") {
    if (i_rank == 0) printf("LOCAL, n_part_in 2, n_part_out 2, POLY_3D, empty extraction\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, HILBERT") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, HILBERT\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  }
#ifdef PDM_HAVE_PARMETIS
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PARMETIS") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PARMETIS\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PARMETIS;
  }
#endif
#ifdef PDM_HAVE_PTSCOTCH
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("REEQUILIBRATE, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_2D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_2D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, PRISM6, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_PRISM6;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 1, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 1;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 1, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 1;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 0;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
  SUBCASE("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH") {
    if (i_rank == 0) printf("FROM_TARGET, n_part_in 2, n_part_out 2, POLY_3D, empty extraction, PTSCOTCH\n");
    extract_kind     = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_in        = 2;
    n_part_out       = 2;
    elt_type         = PDM_MESH_NODAL_POLY_3D;
    empty_extraction = 1;
    part_method      = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  }
#endif





  PDM_geometry_kind_t geom_kind_parent = PDM_GEOMETRY_KIND_MAX;
  int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(elt_type);
  if (mesh_dimension == 3) {
    geom_kind_parent = PDM_GEOMETRY_KIND_VOLUMIC;
  }
  else if (mesh_dimension == 2) {
    geom_kind_parent = PDM_GEOMETRY_KIND_SURFACIC;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "1D not yet available\n");
  }
  PDM_part_mesh_nodal_t *pmn = _generate_mesh(pdm_comm, elt_type, n_part_in);

  PDM_part_mesh_nodal_elmts_t *pmne_parent = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                           geom_kind_parent);


  int  *n_extract    = NULL;
  int **extract_lnum = NULL;

  int          *n_target    = NULL;
  PDM_g_num_t **target_gnum = NULL;

  /* Extract part */
  PDM_bool_t compute_child_gnum = PDM_TRUE;
  PDM_extract_part_t *extrp = PDM_extract_part_create(mesh_dimension,
                                                      n_part_in,
                                                      n_part_out,
                                                      extract_kind,
                                                      part_method,
                                                      compute_child_gnum,
                                                      PDM_OWNERSHIP_KEEP,
                                                      pdm_comm);

  PDM_extract_part_part_nodal_set(extrp, pmn);



  if (extract_kind == PDM_EXTRACT_PART_KIND_FROM_TARGET) {
    /* PDM_EXTRACT_PART_KIND_FROM_TARGET */
    _get_targets(pdm_comm,
                 pmne_parent,
                 n_part_in,
                 n_part_out,
                 empty_extraction,
                 &n_target,
                 &target_gnum);

    for (int i_part = 0; i_part < n_part_out; i_part++) {
      PDM_extract_part_target_set(extrp,
                                  i_part,
                                  n_target   [i_part],
                                  target_gnum[i_part],
                                  NULL);
    }
  }
  else {
    /* LOCAL or REEQUILIBRATE */
    _get_extracted(pdm_comm,
                   pmne_parent,
                   n_part_in,
                   empty_extraction,
                   &n_extract,
                   &extract_lnum);

    for (int i_part = 0; i_part < n_part_in; i_part++) {
      PDM_extract_part_selected_lnum_set(extrp,
                                         i_part,
                                         n_extract   [i_part],
                                         extract_lnum[i_part]);
    }
  }

  // Compute extraction
  PDM_extract_part_compute(extrp);

  /* Check */
  PDM_part_mesh_nodal_t *extract_pmn = NULL;
  PDM_extract_part_part_mesh_nodal_get(extrp,
                                       &extract_pmn,
                                       PDM_OWNERSHIP_KEEP);

  for (int geom_kind = (int) geom_kind_parent; geom_kind < (int) PDM_GEOMETRY_KIND_CORNER; geom_kind++) {

    PDM_part_mesh_nodal_elmts_t *pmne         = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                              (PDM_geometry_kind_t) geom_kind);

    PDM_part_mesh_nodal_elmts_t *extract_pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(extract_pmn,
                                                                                              (PDM_geometry_kind_t) geom_kind);

    // Grouds
    _check_groups(pmne, n_part_in, extract_pmne, n_part_out);
  }



  /* Free memory */
  PDM_part_mesh_nodal_free(pmn);
  PDM_extract_part_free(extrp);

  if (extract_kind == PDM_EXTRACT_PART_KIND_FROM_TARGET) {
    for (int i_part = 0; i_part < n_part_out; i_part++) {
      PDM_free(target_gnum[i_part]);
    }
    PDM_free(n_target   );
    PDM_free(target_gnum);
  }
  else {
    for (int i_part = 0; i_part < n_part_in; i_part++) {
      PDM_free(extract_lnum[i_part]);
    }
    PDM_free(n_extract   );
    PDM_free(extract_lnum);
  }
}


/*
The subcases were generated using the following Python script :


n_subcase = 0

def write_subcase(dico):
  if dico["extract_kind"] == "LOCAL":
    dico["n_part_out"] = dico["n_part_in"]

  subcase = "%s, n_part_in %d, n_part_out %d, %s" % \
  (dico["extract_kind"], dico["n_part_in"], dico["n_part_out"], dico["elt_type"])

  if dico["extract_kind"] != "LOCAL":
    subcase += ", %s" % dico["part_method"]

  s = ""
  s += '  SUBCASE("%s") {\n' % subcase


  s += '    if (i_rank == 0) printf("%s\\n");\n' % subcase

  s += '    extract_kind = PDM_EXTRACT_PART_KIND_%s;\n' % dico["extract_kind"]
  s += '    n_part_in    = %d;\n' % dico["n_part_in"]
  s += '    n_part_out   = %d;\n' % dico["n_part_out"]
  s += '    elt_type     = PDM_MESH_NODAL_%s;\n' % dico["elt_type"]
  if dico["extract_kind"] != "LOCAL":
    s += '    part_method  = PDM_SPLIT_DUAL_WITH_%s;\n' % dico["part_method"]
  s += '  }\n'

  global n_subcase

  n_subcase += 1

  return s




if __name__ == "__main__":

  s = ""

  # LOCAL
  for elt_type in ["POLY_2D", "PRISM6", "POLY_3D"]:
    for n_part in [1, 2]:
      for empty_extraction in [0, 1]:
        s += write_subcase({
                           "extract_kind"    : "LOCAL",
                           "n_part_in"       : n_part,
                           "elt_type"        : elt_type,
                           "empty_extraction": empty_extraction,
                           })

  # REEQUILIBRATE & FROM_TARGET
  for part_method in ["HILBERT", "PARMETIS", "PTSCOTCH"]:

    if part_method != "HILBERT":
      s += "#ifdef PDM_HAVE_%s\n" % part_method

    for extract_kind in ["REEQUILIBRATE", "FROM_TARGET"]:

      for elt_type in ["POLY_2D", "PRISM6", "POLY_3D"]:
        for n_part_in in [1, 2]:
          for n_part_out in [1, 2]:
            for empty_extraction in [0, 1]:
              s += write_subcase({
                                 "extract_kind"    : extract_kind,
                                 "n_part_in"       : n_part_in,
                                 "n_part_out"      : n_part_out,
                                 "elt_type"        : elt_type,
                                 "part_method"     : part_method,
                                 "empty_extraction": empty_extraction,
                                 })

    if part_method != "HILBERT":
      s += "#endif\n"



  print(f"n_subcase = {n_subcase}")

  with open("extract_part_unit_test.c", "w") as out:
    out.write(s)
*/
