/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_check.h"
#include "pdm.h"
#include "pdm_mem_tool.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_check.h"


/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \brief Find unconnected vertex in a mesh connectivity
 *
 * If a vertex is not cited in the face->vtx connectivity, the function
 * removes it from the mesh to ensure contiguity
 *
 * \param [in, out] n_vtx                 Number of vertices
 * \param [in, out] l_face_vtx            Size of face->vtx connectivity
 * \param [in, out] face_vtx              Face->vtx connectivity
 * \param [in     ] nb_vtx_last_problem   last unconnected vertex
 *
 * \return the first unconnected vertex behind \ref nb_vtx_last_problem
 */


static PDM_g_num_t
_unconnected_vertex_find
(
PDM_g_num_t* n_vtx,
PDM_g_num_t* l_face_vtx,
PDM_g_num_t* face_vtx,
PDM_g_num_t  nb_vtx_last_problem
)
{

  // Déclarations
  PDM_g_num_t *check_som = NULL;
  PDM_g_num_t index_courant ;
  PDM_g_num_t nb_som_problem = 0 ;

  // Allocation
  check_som = (PDM_g_num_t*) calloc((*n_vtx), sizeof(PDM_g_num_t));

  // check_som à 1 pour tous les sommets cités dans la connectivité face->som
  for (PDM_g_num_t i=0; i < *l_face_vtx; i++) {
    index_courant = face_vtx[i] - 1;
    check_som[index_courant] = 1 ;
  }

  // fix trailing zeros in first pass
  PDM_g_num_t _nb_vtx_last_problem = nb_vtx_last_problem;
  if (nb_vtx_last_problem == *n_vtx) {
    while (check_som[_nb_vtx_last_problem-1] == 0 && _nb_vtx_last_problem > 0) {
      _nb_vtx_last_problem--;
    }
  }

  // si check_som est 0 quelquepart, c'est qu'il y a un trou dans la connectivité
  for (PDM_g_num_t i=_nb_vtx_last_problem; i > 0; i--) {
    if (check_som[i-1] == 0) {
      nb_som_problem = i ;
      break ;
    }
  }

  PDM_free(check_som) ;

  return nb_som_problem ;
}


/**
 * \brief Remove unconnected vertex in a mesh connectivity
 *
 * If a vertex is not cited in the face->vtx connectivity, the function
 * removes it from the mesh to ensure contiguity
 *
 * \param [in, out] n_vtx                 Number of vertices
 * \param [in     ] nb_vtx_problem        unconnected vertex
 * \param [in     ] coords                Coordinates
 * \param [in, out] l_face_vtx            Size of face->vtx connectivity
 * \param [in, out] face_vtx              Face->vtx connectivity
 *
 */

static void
_unconnected_vertex_remove
(
PDM_g_num_t* n_vtx,
PDM_g_num_t  nb_vtx_problem,
double*      coords,
PDM_g_num_t* l_face_vtx,
PDM_g_num_t* face_vtx
)
{

  PDM_g_num_t index_courant ;

  // Coordonnees
  for (PDM_g_num_t i=nb_vtx_problem; i < *n_vtx; i++) {
    coords[3*(i-1)  ] = coords[3*i  ] ; //x(n-1) <- x(n)
    coords[3*(i-1)+1] = coords[3*i+1] ; //y(n-1) <- y(n)
    coords[3*(i-1)+2] = coords[3*i+2] ; //z(n-1) <- z(n)
  }

  // Connectivité faces->sommets
  for (PDM_g_num_t i=0; i < *l_face_vtx; i++) {
    index_courant = face_vtx[i] ;
    if (index_courant > nb_vtx_problem) {
      face_vtx[i] = index_courant - 1 ;
    }
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_mesh_check_unconnected_vertex
(
PDM_g_num_t* nb_vtx,
PDM_g_num_t* l_face_vtx,
PDM_g_num_t* face_vtx,
double*      coords,
int*         nb_holes
)
{

  PDM_g_num_t nb_som_problem = 0 ; // (absolute) number of problematic vertex
  PDM_g_num_t nb_som_last_problem = *nb_vtx ; // problematic vertex from the previous iteration
  *nb_holes = 0 ;               // total number of holes in the mesh

  do {

    nb_som_problem = _unconnected_vertex_find // Find the number of a problematic vertex
                     (nb_vtx,
                      l_face_vtx,
                      face_vtx,
                      nb_som_last_problem) ;

    if (nb_som_problem != 0) {
      _unconnected_vertex_remove                // if found, fix it
        (nb_vtx,
         nb_som_problem,
         coords,
         l_face_vtx,
         face_vtx) ;

      (*nb_holes)++ ;
      nb_som_last_problem = nb_som_problem ;
    }

  } while (nb_som_problem != 0) ;

}


/*----------------------------------------------------------------------------*/

