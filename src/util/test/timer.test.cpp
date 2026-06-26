#include <stddef.h>
#include <unistd.h>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_timer.h"



MPI_TEST_CASE("[pdm_timer] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  /*
   * Example with mulitple level (depth > 1)
   */
  PDM_timer_t* timer = PDM_timer_create(pdm_comm);

  PDM_timer_start(timer, "compute"  , 1);

  sleep(2.);

  PDM_timer_start(timer, "sub_step1", 1);

  // New Level
  PDM_timer_start(timer, "sub_sub_step1", 1);

  if(i_rank == 0) {
    sleep(1.);
  }

  PDM_timer_end  (timer, "sub_sub_step1", 1);
  PDM_timer_end  (timer, "sub_step1", 1);

  PDM_timer_start(timer, "sub_step2", 1);

  if(i_rank == 1) {
    sleep(3.);
  }
  PDM_timer_end  (timer, "sub_step2", 1);

  PDM_timer_end  (timer, "compute"  , 1);

  /*
   * Multiple dump possibilities
   */
  // 1 file per proc - json (to be post-treat after)
  char filename[999];
  sprintf(filename, "profiling_%i.json", i_rank);
  PDM_timer_dump_json(timer, filename);

  // 1 file - All gather
  const char* filename_gather = "profiling_gather.log";
  PDM_timer_gather_dump(timer, filename_gather);

  // Same but in stdout
  PDM_timer_gather_dump(timer, NULL);

  if(i_rank == 0) {
    // Print timer - Hierarchical view
    PDM_timer_print(timer, 0);
    // Print timer - Flat view
    PDM_timer_print(timer, 1);
  }
  // PDM_timer_log(timer, 0);
  // PDM_timer_log(timer, 1);

  PDM_timer_gather_dump_json(timer, "profiling_gather.json");

  remove(filename);

  if(i_rank == 0) {
    remove("profiling_gather.log");
    remove("profiling_gather.json");
  }

  PDM_timer_free(timer);
}


MPI_TEST_CASE("[pdm_timer] - get", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_timer_t *timer = PDM_timer_create(pdm_comm);

  PDM_timer_start(timer, "root", 1);
  {
    PDM_timer_start(timer, "left", 1);
    {
      PDM_timer_start(timer, "child", 1);
      sleep(1.);
      PDM_timer_end  (timer, "child", 1);
    }
    {
      PDM_timer_start(timer, "child", 1);
      sleep(3.);
      PDM_timer_end  (timer, "child", 1);
    }
    PDM_timer_end(timer, "left", 1);
  }
  {
    PDM_timer_start(timer, "right", 1);
    {
      PDM_timer_start(timer, "child", 1);
      sleep(2.);
      PDM_timer_end  (timer, "child", 1);
    }
    PDM_timer_end(timer, "right", 1);
  }
  PDM_timer_end(timer, "root", 1);


  long   n_call;
  double duration;

  n_call = PDM_timer_get(timer, "/root/left/child", &duration);
  CHECK(n_call == 2);
  CHECK(PDM_ABS(duration - 4.) < 1e-3);

  n_call = PDM_timer_get(timer, "/root/right/child", &duration);
  CHECK(n_call == 1);
  CHECK(PDM_ABS(duration - 2.) < 1e-3);

  PDM_timer_free(timer);
}