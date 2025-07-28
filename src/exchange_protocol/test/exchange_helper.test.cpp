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


MPI_TEST_CASE("[PDM_exchange_helper] - Create ", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_exchange_helper_t* exch_helper = PDM_exchange_helper_create(pdm_comm, 10);

  std::vector<std::vector<int>> send_buffer = {{1 ,  2,  3,  4, 5},
                                               {10, 20, 30, 40   }};

  std::vector<std::vector<int>> send_idx = {{0, 3, 5},
                                            {0, 1, 4}};

  std::vector<int> send_n  (n_rank  );
  std::vector<int> recv_n  (n_rank  );
  std::vector<int> recv_idx(n_rank+1);

  for(int i = 0; i < n_rank; ++i ) {
    send_n[i] = send_idx[i_rank][i+1] - send_idx[i_rank][i];
  }

  PDM_MPI_Alltoall(send_n.data(), 1, PDM_MPI_INT,
                   recv_n.data(), 1, PDM_MPI_INT, pdm_comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i ) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  std::vector<int> recv_buffer(recv_idx[n_rank]);


  // Scenario 1 :
  //   - Simple exch of data (two-way)
  // PDM_exchange_helper_exch(exch_helper,
  //                          PDM_MPI_COMM_KIND_COLLECTIVE,
  //                          1,
  //                          sizeof(int),
  //                          send_idx   [i_rank].data(),
  //                          send_n             .data(),
  //                          send_buffer[i_rank].data(),
  //                          recv_idx           .data(),
  //                          recv_n             .data(),
  //                          recv_buffer        .data());

  int request_id = PDM_exchange_helper_exch_init(exch_helper,
                                                 PDM_MPI_COMM_KIND_P2P,
                                                 1,
                                                 sizeof(int),
                                                 send_idx   [i_rank].data(),
                                                 send_n             .data(),
                                                 send_buffer[i_rank].data(),
                                                 recv_idx           .data(),
                                                 recv_n             .data(),
                                                 recv_buffer        .data());


  PDM_exchange_helper_exch_start(exch_helper,
                                 request_id);

  PDM_exchange_helper_exch_wait(exch_helper,
                                request_id);

  // PDM_log_trace_array_int(recv_buffer.data(), recv_idx[n_rank], "recv_buffer ::");

  // Changement des buffers !
  for(int i = 0; i < static_cast<int>(send_buffer[i_rank].size()); ++i) {
    send_buffer[i_rank][i] += 10;
  }

  PDM_exchange_helper_exch_start(exch_helper,
                                 request_id);

  PDM_exchange_helper_exch_wait(exch_helper,
                                request_id);

  // PDM_log_trace_array_int(recv_buffer.data(), recv_idx[n_rank], "recv_buffer ::");


  PDM_exchange_helper_exch_free(exch_helper, request_id);

  // Scenario 2 :
  //   - Simple exch of data (two-way)


  // Dans les échanges c'est pas forcement bijectif : ex --> RMA
  // Dans l'utilisation avec sonics par exemple l'allocation est externe, et on aura pas le choix en RMA aussi (car synchro amont obligatoire)

  // PDM_exchange_helper_exch_add(exch_helper,
  //                              cst_stride,
  //                              s_data,
  //                              tag,
  //                               PDM_OWNERSHIP_KEEP);

  // PDM_exchange_helper_add_recv_exch(exch_helper,
  //                                   cst_stride,
  //                                   s_data, )

  // PDM_exchange_helper_add_collective_exch(exch_helper,
  //                                         cst_stride,
  //                                         s_data, )


  PDM_exchange_helper_free(exch_helper);
}
