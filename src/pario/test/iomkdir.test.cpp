#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "doctest/extensions/doctest_mpi.h"
#include "pdm_io.h"


#define PDM_DIRNAME_TESTCASE    "pdm_io_mkdir_test_case"
#define PDM_MAX_STR_LEN_TESTCASE 1024

// =============================================================================
MPI_TEST_CASE("[PDM_io] - PDM_io_mkdir", 1) {

  std::string path = PDM_DIRNAME_TESTCASE;

  // ---------------------------------------------------------------------------
  SUBCASE("Relative_path - <dirname>") {
    int         err            = -1;
    bool        is_dir_created = false;
    std::string path_rel       = path;
    struct stat sb;

    err            = PDM_io_mkdir(path_rel.c_str());
    is_dir_created = stat(path_rel.c_str(), &sb) == 0;

    CHECK(err            == 0);
    CHECK(is_dir_created == true);

    if (is_dir_created == true) {
      err = rmdir(path_rel.c_str());
    }

    CHECK(err != -1);
  }

  // ---------------------------------------------------------------------------
  SUBCASE("Relative_path - ./<dirname>") {
    int         err            = -1;
    bool        is_dir_created = false;
    std::string path_rel       = "./" + path;
    struct stat sb;

    err            = PDM_io_mkdir(path_rel.c_str());
    is_dir_created = stat(path_rel.c_str(), &sb) == 0;

    CHECK(err            == 0);
    CHECK(is_dir_created == true);

    if (is_dir_created == true) {
      err = rmdir(path_rel.c_str());
    }

    CHECK(err != -1);
  }

  // ---------------------------------------------------------------------------
  SUBCASE("absolute_path - /<pwd>/<dirname>") {
    char  path_pwd[PDM_MAX_STR_LEN_TESTCASE];
    char *err_ptr = NULL;

    err_ptr = getcwd(path_pwd, PDM_MAX_STR_LEN_TESTCASE);

    CHECK(err_ptr != NULL);

    int         err            = -1;
    bool        is_dir_created = false;
    std::string path_abs       = (std::string)path_pwd + "/" + path;
    struct stat sb;

    err            = PDM_io_mkdir(path_abs.c_str());
    is_dir_created = stat(path_abs.c_str(), &sb) == 0;

    CHECK(err            == 0);
    CHECK(is_dir_created == true);

    if (is_dir_created == true) {
      err = rmdir(path_abs.c_str());
    }

    CHECK(err != -1);
  }

  // ---------------------------------------------------------------------------
  SUBCASE("Relative_chained_path  - <dirname>/<dirname>_lv2") {
    int         err              = -1;
    bool        is_dir_created   = false;
    bool        is_child_created = false;
    std::string path_rel         = path;
    std::string path_child       = path + "/" + path + "_lv2";
    struct stat sb;

    err              = PDM_io_mkdir(path_child.c_str());
    is_dir_created   = stat(path_rel.c_str(), &sb) == 0;
    is_child_created = stat(path_child.c_str(), &sb) == 0;

    CHECK(err              == 0);
    CHECK(is_dir_created   == true);
    CHECK(is_child_created == true);

    if (is_child_created == true) {
      err = rmdir(path_child.c_str());
    }
    CHECK(err != -1);

    if (is_dir_created == true) {
      err = rmdir(path_rel.c_str());
    }
    CHECK(err != -1);
  }

  // ---------------------------------------------------------------------------
  SUBCASE("Relative_chained_path_with_revert  - <dirname>/../<dirname>_lv2") {
    int         err              = -1;
    bool        is_dir_created   = false;
    bool        is_child_created = false;
    std::string path_rel         = path;
    std::string path_child       = path + "/../" + path + "_lv2";
    struct stat sb;

    err              = PDM_io_mkdir(path_child.c_str());
    is_dir_created   = stat(path_rel.c_str(), &sb) == 0;
    is_child_created = stat(path_child.c_str(), &sb) == 0;

    CHECK(err              == 0);
    CHECK(is_dir_created   == true);
    CHECK(is_child_created == true);

    if (is_child_created == true) {
      err = rmdir(path_child.c_str());
    }
    CHECK(err != -1);

    if (is_dir_created == true) {
      err = rmdir(path_rel.c_str());
    }
    CHECK(err != -1);
  }

}

#undef PDM_DIRNAME_TESTCASE
#undef PDM_MAX_STR_LEN_TESTCASE
