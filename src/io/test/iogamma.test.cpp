#include <iostream>
#include <stdio.h>

#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm_mem_tool.h"
#include "pdm_predicate.h"
#include "pdm_generate_mesh.h"
#include "pdm_part_mesh_nodal.h"

// =======================================================================================
#define FILL_TRIA3_VTX_COORD(vtx_coord) \
vtx_coord[0] =  0.0; vtx_coord[1] = 0.0; vtx_coord[2] = 0.0; \
vtx_coord[3] =  1.5; vtx_coord[4] = 1.5; vtx_coord[5] = 0.0; \
vtx_coord[6] = -1.5; vtx_coord[7] = 1.5; vtx_coord[8] = 0.0;

#define FILL_QUAD4_VTX_COORD(vtx_coord) \
vtx_coord[0] = -1.5; vtx_coord[1] = 0.0; vtx_coord[2] = 0.0; \
vtx_coord[3] =  1.5; vtx_coord[4] = 0.0; vtx_coord[5] = 0.0; \
vtx_coord[6] =  1.5; vtx_coord[7] = 1.5; vtx_coord[8] = 0.0; \
vtx_coord[6] = -1.5; vtx_coord[7] = 1.5; vtx_coord[8] = 0.0;


// =======================================================================================
MPI_TEST_CASE("[PDM_io_gamma_reader] - element orientation", 1) {

	double *vtx_coord = NULL;

	// -----------------------------------------------------------------------------------
	SUBCASE("tria3 - no correction") {
		int n_vtx = 3;
		PDM_malloc(vtx_coord, 3*n_vtx, double);

		FILL_TRIA3_VTX_COORD(vtx_coord)

		double *p1 = vtx_coord;
		double *p2 = vtx_coord + 3;
		double *p3 = vtx_coord + 6;

		double surf = PDM_predicate_orient2d(p1, p2, p3);

		CHECK(surf > 0);
	}

	// -----------------------------------------------------------------------------------
	SUBCASE("tria3 - correction") {
		int n_vtx = 3;
		PDM_malloc(vtx_coord, 3*n_vtx, double);

		FILL_TRIA3_VTX_COORD(vtx_coord)

		double *p1   = vtx_coord;
		double *p2   = vtx_coord + 3;
		double *p3   = vtx_coord + 6;
		double  surf = PDM_predicate_orient2d(p1, p3, p2);

		CHECK(surf < 0);
	}

	// -----------------------------------------------------------------------------------
	SUBCASE("quad4 - no correction") {
		int n_vtx = 4;
		PDM_malloc(vtx_coord, 3*n_vtx, double);

		FILL_QUAD4_VTX_COORD(vtx_coord);

		double *p1   = vtx_coord;
		double *p2   = vtx_coord + 3;
		double *p3   = vtx_coord + 6;
		double *p4   = vtx_coord + 9;
		double  surf = PDM_predicate_orient2d_quad(p1, p2, p3, p4);

		CHECK(surf > 0);
	}

	// -----------------------------------------------------------------------------------
	SUBCASE("quad4 - correction") {
		int n_vtx = 4;
		PDM_malloc(vtx_coord, 3*n_vtx, double);

		FILL_QUAD4_VTX_COORD(vtx_coord);

		double *p1   = vtx_coord;
		double *p2   = vtx_coord + 3;
		double *p3   = vtx_coord + 6;
		double *p4   = vtx_coord + 9;
		double  surf = PDM_predicate_orient2d_quad(p1, p4, p3, p2);

		CHECK(surf < 0);
	}


	PDM_free(vtx_coord);
}


MPI_TEST_CASE("[PDM_io_gamma_reader] - Read and write", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

	const char *filename_in_2d = PDM_MESH_DIR"mixed_elements_2d.mesh";
	const char *filename_in_3d = PDM_MESH_DIR"mixed_elements_3d.mesh";

	char filename_out[999];

	PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

	const char *filename_in;
	int n_part;
	int dim;
	SUBCASE("2D") {
		dim         = 2;
		filename_in = filename_in_2d;
		SUBCASE("n_part = 1") {
			n_part = 1;
		}
		SUBCASE("n_part = 2") {
			n_part = 2;
		}
	}
	SUBCASE("3D") {
		dim         = 3;
		filename_in = filename_in_3d;
		SUBCASE("n_part = 1") {
			n_part = 1;
		}
		SUBCASE("n_part = 2") {
			n_part = 2;
		}
	}

	// Load and split mesh
	PDM_part_mesh_nodal_t *mesh = PDM_generate_mesh_nodal_from_file(pdm_comm,
		                                                              n_part,
																																	part_method,
																																	filename_in);

  // Re-write mesh
	sprintf(filename_out, "unit_test_gamma_io_dim_%d_n_part_%d.mesh", dim, n_part);
	PDM_part_mesh_nodal_dump_gamma(mesh, filename_out);

	// Check
	// TODO

	remove(filename_out);

	PDM_part_mesh_nodal_free(mesh);
}

#undef FILL_TRIA3_VTX_COORD
#undef FILL_QUAD4_VTX_COORD