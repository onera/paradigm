#include "doctest/extensions/doctest_mpi.h"
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_io_utils.h"


MPI_TEST_CASE("[PDM_io_utils] - PDM_io_utils_file_extension", 1) {

  const char *extension;

  extension = PDM_io_utils_file_extension("toto.tata");
  CHECK(strcmp(extension, "tata") == 0);

  extension = PDM_io_utils_file_extension("titi.tutu.tata");
  CHECK(strcmp(extension, "tata") == 0);
}


MPI_TEST_CASE("[PDM_io_utils] - PDM_io_utils_file_name_from_path", 1) {

  const char *name;

  name = PDM_io_utils_file_name_from_path("foo.bar");
  CHECK(strcmp(name, "foo.bar") == 0);

  name = PDM_io_utils_file_name_from_path("/titi/foo.bar");
  CHECK(strcmp(name, "foo.bar") == 0);

  name = PDM_io_utils_file_name_from_path("tata/titi/foo.bar");
  CHECK(strcmp(name, "foo.bar") == 0);

}


MPI_TEST_CASE("[PDM_io_utils] - PDM_io_utils_diff_files", 1) {

  int expected;

  const char *filename1 = PDM_MESH_DIR"box.mesh";
  const char *filename2 = PDM_MESH_DIR"sphere.stl";
  const char *filename3 = "404/not.found";

  int diff;
  SUBCASE("Indentical") {
    expected = 0;
    diff = PDM_io_utils_diff_files(filename1, filename1);
  }
  SUBCASE("Different") {
    expected = 1;
    diff = PDM_io_utils_diff_files(filename1, filename2);
  }
  SUBCASE("Invalid file 1") {
    expected = -1;
    diff = PDM_io_utils_diff_files(filename3, filename2);
  }
  SUBCASE("Invalid file 2") {
    expected = -2;
    diff = PDM_io_utils_diff_files(filename1, filename3);
  }

  CHECK(diff == expected);
}