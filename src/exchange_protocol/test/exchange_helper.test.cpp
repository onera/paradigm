#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_doctest.h"
#include "pdm_exchange_helper.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include <functional>


MPI_TEST_CASE("[PDM_exchange_helper] - Create ", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  // int i_rank;

  PDM_exchange_helper_t* exch_helper = PDM_exchange_helper_create(pdm_comm, 10);

  // Dans les échanges c'est pas forcement bijectif : ex --> RMA
  // Dans l'utilisation avec sonics par exemple l'allocation est externe, et on aura pas le choix en RMA aussi (car synchro amont obligatoire)

  // PDM_exchange_helper_send_exch_add(exch_helper,
  //                                   cst_stride,
  //                                   s_data,
  //                                   tag,
  //                                   PDM_OWNERSHIP_KEEP);

  // PDM_exchange_helper_add_recv_exch(exch_helper,
  //                                   cst_stride,
  //                                   s_data, )

  // PDM_exchange_helper_add_collective_exch(exch_helper,
  //                                         cst_stride,
  //                                         s_data, )


  PDM_exchange_helper_free(exch_helper);
}
