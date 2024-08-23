#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_extract_part.h"
#include "pdm_priv.h"

#include "pdm_array.h"
#include "pdm_generate_mesh.h"
#include "pdm_distrib.h"

static void
_get_extracted
(
 PDM_MPI_Comm                   comm,
 PDM_part_mesh_nodal_elmts_t   *pmne,
 int                            n_part,
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
 int                          **n_target,
 PDM_g_num_t                 ***target_gnum
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

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

  PDM_malloc(*n_target,    n_part_out, int          );
  PDM_malloc(*target_gnum, n_part_out, PDM_g_num_t *);

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
  int n_group         = PDM_part_mesh_nodal_elmts_n_group_get(pmne);
  int extract_n_group = PDM_part_mesh_nodal_elmts_n_group_get(extract_pmne);
  CHECK(n_group == extract_n_group);

  for (int i_part = 0; i_part < n_part_out; i_part++) {

    int n_elt_tot = PDM_part_mesh_nodal_elmts_n_elmts_get(extract_pmne, i_part);

    int n_elt_in_all_groups = 0;
    for (int i_group = 0; i_group < extract_n_group; i_group++) {
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

      n_elt_in_all_groups += n_group_elmt;
    }

    CHECK(n_elt_in_all_groups == n_elt_tot);
  }
}



MPI_TEST_CASE("[pdm_extract_part] - 2p - nodal", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /* Generate mesh */
  PDM_Mesh_nodal_elt_t elt_type    = PDM_MESH_NODAL_PRISM6;
  PDM_g_num_t          n_vtx_seg   = 5;
  int                  n_part_in   = 1;
  PDM_split_dual_t     part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(pdm_comm,
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
                                                                n_part_in,
                                                                part_method);

  PDM_part_mesh_nodal_elmts_t *pmne_vol = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                        PDM_GEOMETRY_KIND_VOLUMIC);


  int n_part_out = -1;
  PDM_extract_part_kind_t extract_kind;

  int  *n_extract    = NULL;
  int **extract_lnum = NULL;

  int          *n_target    = NULL;
  PDM_g_num_t **target_gnum = NULL;

  SUBCASE("Local") {
    if (i_rank == 0) printf("Local\n");
    extract_kind = PDM_EXTRACT_PART_KIND_LOCAL;
    n_part_out   = n_part_in;
  }

  SUBCASE("Reequilibrate - 1 part in, 1 part out") {
    if (i_rank == 0) printf("Reequilibrate - 1 part in, 1 part out\n");
    extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_out   = 1;
  }

  SUBCASE("Reequilibrate - 1 part in, 2 part out") {
    if (i_rank == 0) printf("Reequilibrate - 1 part in, 2 part out\n");
    extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    n_part_out   = 2;
  }

  SUBCASE("From target - 1 part in, 1 part out") {
    if (i_rank == 0) printf("From target - 1 part in, 1 part out\n");
    extract_kind = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_out   = 1;
  }

  SUBCASE("From target - 1 part in, 2 part out") {
    if (i_rank == 0) printf("From target - 1 part in, 2 part out\n");
    extract_kind = PDM_EXTRACT_PART_KIND_FROM_TARGET;
    n_part_out   = 2;
  }


  /* Extract part */
  PDM_bool_t compute_child_gnum = PDM_TRUE;
  PDM_extract_part_t *extrp = PDM_extract_part_create(3,
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
                 pmne_vol,
                 n_part_in,
                 n_part_out,
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
                   pmne_vol,
                   n_part_in,
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

  for (PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC; geom_kind < PDM_GEOMETRY_KIND_CORNER; geom_kind++) {

    PDM_part_mesh_nodal_elmts_t *pmne         = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                              geom_kind);

    PDM_part_mesh_nodal_elmts_t *extract_pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(extract_pmn,
                                                                                              geom_kind);

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
