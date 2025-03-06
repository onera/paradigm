#include <ext/alloc_traits.h>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_doctest.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"


MPI_TEST_CASE("[PDM_sphere_surf_gen] - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;
  PDM_g_num_t nv        =     n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;

  PDM_sphere_surf_gen(pdm_comm,
                      nu,
                      nv,
                      x_center,
                      y_center,
                      z_center,
                      radius,
                      &dvtx_coord,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &distrib_vtx,
                      &distrib_face);

  free(dvtx_coord   );
  free(dface_vtx_idx);
  free(dface_vtx    );
  free(distrib_vtx  );
  free(distrib_face );

}


MPI_TEST_CASE("[PDM_sphere_surf_gen] - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;
  PDM_g_num_t nv        =     n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;

  PDM_sphere_surf_gen(pdm_comm,
                      nu,
                      nv,
                      x_center,
                      y_center,
                      z_center,
                      radius,
                      &dvtx_coord,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &distrib_vtx,
                      &distrib_face);

  free(dvtx_coord   );
  free(dface_vtx_idx);
  free(dface_vtx    );
  free(distrib_vtx  );
  free(distrib_face );

}

MPI_TEST_CASE("[PDM_sphere_surf_gen_nodal] - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;
  PDM_g_num_t nv        =     n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  PDM_dmesh_nodal_t *dmn = NULL;

  PDM_sphere_surf_gen_nodal(pdm_comm,
                            nu,
                            nv,
                            x_center,
                            y_center,
                            z_center,
                            radius,
                            &dmn);
  PDM_DMesh_nodal_free(dmn);

}

MPI_TEST_CASE("[PDM_sphere_surf_gen_nodal] - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;
  PDM_g_num_t nv        =     n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  PDM_dmesh_nodal_t *dmn = NULL;

  PDM_sphere_surf_gen_nodal(pdm_comm,
                            nu,
                            nv,
                            x_center,
                            y_center,
                            z_center,
                            radius,
                            &dmn);
  PDM_DMesh_nodal_free(dmn);

}



MPI_TEST_CASE("[PDM_sphere_surf_icosphere_gen] - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;

  PDM_sphere_surf_icosphere_gen(pdm_comm,
                                nu,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dvtx_coord,
                                &dface_vtx_idx,
                                &dface_vtx,
                                &distrib_vtx,
                                &distrib_face);

  free(dvtx_coord   );
  free(dface_vtx_idx);
  free(dface_vtx    );
  free(distrib_vtx  );
  free(distrib_face );

}


MPI_TEST_CASE("[PDM_sphere_surf_icosphere_gen] - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;

  PDM_sphere_surf_icosphere_gen(pdm_comm,
                                nu,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dvtx_coord,
                                &dface_vtx_idx,
                                &dface_vtx,
                                &distrib_vtx,
                                &distrib_face);

  free(dvtx_coord   );
  free(dface_vtx_idx);
  free(dface_vtx    );
  free(distrib_vtx  );
  free(distrib_face );

}

MPI_TEST_CASE("[PDM_sphere_surf_icosphere_gen_nodal] - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  PDM_dmesh_nodal_t *dmn = NULL;

  PDM_sphere_surf_icosphere_gen_nodal(pdm_comm,
                                      nu,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn);
  PDM_DMesh_nodal_free(dmn);

}

MPI_TEST_CASE("[PDM_sphere_surf_icosphere_gen_nodal] - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_vtx_seg = 10;
  PDM_g_num_t nu        = 2 * n_vtx_seg;

  double length   = 1.;
  double x_center = 0.;
  double y_center = 0.;
  double z_center = 0.;
  double radius   = 0.8*length;


  PDM_dmesh_nodal_t *dmn = NULL;

  PDM_sphere_surf_icosphere_gen_nodal(pdm_comm,
                                      nu,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn);
  PDM_DMesh_nodal_free(dmn);

}
