#include <iostream>

#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm_mem_tool.h"
#include "pdm_predicate.h"

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

#undef FILL_TRIA3_VTX_COORD
#undef FILL_QUAD4_VTX_COORD