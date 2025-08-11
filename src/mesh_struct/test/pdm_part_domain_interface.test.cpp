#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm.h"
#include "pdm_logging.h"
#include "pdm_doctest.h"
#include "pdm_domain_interface.h"
#include "pdm_domain_interface_priv.h"
#include "pdm_part_domain_interface.h"

#include <vector>

/*
 *  Use case
 *
 *            4
 *     5 +----+----+ 1
 *        \       /
 *         \     + 2
 *          \   /
 *           \ /
 *          3 +
 *
 *  5 <-> 1 (gnum=1)
 *  3 <-> 3 (gnum=2)
 *
 */
MPI_TEST_CASE("[PDM_part_domain_interface] PDM_part_domain_interface_to_domain_interface - butterfly", 2) {

  int debug_verbose = 1;

  int i_rank;
  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_interface = 1;
  int n_domain    = 1;
  int n_part      = 1;
  PDM_part_domain_interface_t *pdi = PDM_part_domain_interface_create(n_interface,
                                                                      n_domain,
                                                                     &n_part,
                                                                      PDM_DOMAIN_INTERFACE_MULT_YES,
                                                                      PDM_OWNERSHIP_USER,
                                                                      comm);

  // > Partition description
  std::vector<std::vector<int        >> g_n_vtx    = {      {4},     {3}};
  std::vector<std::vector<PDM_g_num_t>> g_vtx_gnum = {{3,1,5,4}, {2,1,3}};

  // > Partitioned domain interface description
  std::vector<            int         > itrf_g_pn       = {4, 3};
  std::vector<std::vector<PDM_g_num_t>> itrf_g_gnum     = {{ 2, 1,  1, 2}, { 2, 2,  1}};
  std::vector<std::vector<int        >> itrf_g_sens     = {{ 1, 1,  1, 1}, { 1, 1,  1}};
  std::vector<std::vector<int        >> itrf_g_sign     = {{-1, 1, -1, 1}, {-1, 1, -1}};
  std::vector<std::vector<int        >> itrf_g_trpt_idx = {{0,3,6,8,11},
                                                           {0,3,6,8  }};
  std::vector<std::vector<int        >> itrf_g_trpt     = {{0,0,0,   0,0,0, 1,0,2,
                                                            0,0,2,   0,0,1, 1,0,1,
                                                            0,0,1,   0,0,2,
                                                            0,0,0,   0,0,0, 1,0,2},
                                                           {1,0,2,   1,0,2, 0,0,0,
                                                            1,0,2,   1,0,2, 0,0,0,
                                                            1,0,1,   0,0,2}};
  std::vector<std::vector<int        >> itrf_g_dom      = {{0, 0,0,
                                                            0, 0,0,
                                                            0, 0,
                                                            0, 0,0},
                                                           {0, 0,0,
                                                            0, 0,0,
                                                            0, 0}};

  int          *n_vtx    = g_n_vtx   [i_rank].data();
  PDM_g_num_t *_vtx_gnum = g_vtx_gnum[i_rank].data();
  PDM_g_num_t **vtx_gnum = &_vtx_gnum;

  int          itrf_pn       = itrf_g_pn      [i_rank];
  PDM_g_num_t *itrf_gnum     = itrf_g_gnum    [i_rank].data();
  int         *itrf_sens     = itrf_g_sens    [i_rank].data();
  int         *itrf_sign     = itrf_g_sign    [i_rank].data();
  int         *itrf_trpt_idx = itrf_g_trpt_idx[i_rank].data();
  int         *itrf_trpt     = itrf_g_trpt    [i_rank].data();
  int         *itrf_dom      = itrf_g_dom     [i_rank].data();

  PDM_part_domain_interface_set(pdi,
                                PDM_BOUND_TYPE_VTX,
                                0,
                                0,
                                0,
                                itrf_pn,
                                itrf_gnum,
                                itrf_sign,
                                itrf_sens,
                                itrf_trpt,
                                itrf_trpt_idx,
                                itrf_dom);
  // PDM_part_domain_interface_set(pdi,
  //                               PDM_BOUND_TYPE_EDGE,
  //                               0,
  //                               0,
  //                               0,
  //                               0,
  //                               NULL,
  //                               NULL,
  //                               NULL,
  //                               NULL,
  //                               NULL,
  //                               NULL);
  // PDM_part_domain_interface_set(pdi,
  //                               PDM_BOUND_TYPE_FACE,
  //                               0,
  //                               0,
  //                               0,
  //                               0,
  //                               NULL,
  //                               NULL,
  //                               NULL,
  //                               NULL,
  //                               NULL,
  //                               NULL);


  PDM_domain_interface_t  *di = NULL;
  int                    **is_entity1_on_itrf = NULL;

  PDM_part_domain_interface_to_domain_interface(pdi,
                                                PDM_BOUND_TYPE_VTX,
                                               &n_part,
                                               &n_vtx,
                                               &vtx_gnum,
                                               &di,
                                               &is_entity1_on_itrf);

  if(debug_verbose == 1) {
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      PDM_log_trace_array_long(di->interface_ids_vtx[i_interface], 2 *di->interface_dn_vtx[i_interface], "interface_ids_vtx ::");
      PDM_log_trace_array_int (di->interface_dom_vtx[i_interface], 2 *di->interface_dn_vtx[i_interface], "interface_dom_vtx ::");
    }
  }

  std::vector<std::vector<int>> g_expctd_ids = {{5,1}, {3,3}};
  std::vector<std::vector<int>> g_expctd_dom = {{0,0}, {0,0}};
  int *expctd_ids = g_expctd_ids[i_rank].data();
  int *expctd_dom = g_expctd_dom[i_rank].data();
  MPI_CHECK_EQ_C_ARRAY(0, di->interface_ids_vtx[0], expctd_ids,  2);
  MPI_CHECK_EQ_C_ARRAY(0, di->interface_dom_vtx[0], expctd_dom,  2);

  free(is_entity1_on_itrf[0]);
  free(is_entity1_on_itrf);
  PDM_part_domain_interface_free(pdi);
  PDM_domain_interface_free(di);

}
