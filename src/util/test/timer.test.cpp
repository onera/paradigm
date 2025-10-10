#include <stddef.h>
#include <unistd.h>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"



MPI_TEST_CASE("[pdm_timer] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_timer_t* timer = PDM_timer_create(pdm_comm);

  PDM_timer_start(timer, "compute"  , 1);

  sleep(2.);

  PDM_timer_start(timer, "sub_step1", 1);

  if(i_rank == 0) {
    sleep(1.);
  }

  PDM_timer_end  (timer, "sub_step1", 1);

  PDM_timer_start(timer, "sub_step2", 1);

  if(i_rank == 1) {
    sleep(3.);
  }
  PDM_timer_end  (timer, "sub_step2", 1);

  PDM_timer_end  (timer, "compute"  , 1);

  // char filename[999];
  // sprintf(filename, "profiling_%i.json", i_rank);
  // PDM_timer_dump_json(timer, filename);

  // sprintf(filename, "debug.log", i_rank);
  // PDM_timer_gather_dump(timer, filename);

  if(i_rank == 0) {
    PDM_timer_print(timer, 0);
    PDM_timer_print(timer, 1);
  }
  // PDM_timer_log(timer, 0);
  // PDM_timer_log(timer, 1);


  // get_time_from_path(timer, this_part, path)
  // orig = get_current_root  (timer, ....)


  PDM_timer_free(timer);

}
