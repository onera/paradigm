#ifndef __PDM_H__
#define __PDM_H__

#include <limits.h>
#include "pdm_config.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if defined (__uxpv__)
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
#endif

#ifdef PDM_LONG_G_NUM
#define PDM_FMT_G_NUM "%ld"
#else
#define PDM_FMT_G_NUM "%d"
#endif

#define PDM_FMT_L_NUM "%d"

#define PDM_MAX_CHAR_LENGTH 100

#define PDM_UNUSED(x) (void)(x)

/**
 * \brief Interface to for hexa manipulation
 * \param [in]       Hex config flags
 * \param [in]       Flag to test
 * \param [out]      An integer, 1 if the flag is set, 0 if unset
 */

#define PDM_HASFLAG( config, flag )   ((config & flag) == flag)
/**
 * \brief Interface to for hexa manipulation
 * \param Arg:       Hex config flags
 * \param Arg:       Flag to test
 * \param [out]      The flags as integer
 */
#define PDM_SETFLAG( config, flag )   ( config |= flag)

/*
 * \brief Interface to for hexa manipulation
 * \param [in]       Hex config flags
 * \param [in]       flag to test
 * \param [out]      The flags as integer
 */
#define PDM_UNSETFLAG( config, flag ) ( config &= ~flag)

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

#if defined(_WIN64)
# define __int64 long long
#endif


/**
 * \enum PDM_g_num_t
 * \brief Long int in pdm
 *
 */

#ifdef PDM_LONG_G_NUM
typedef  long PDM_g_num_t;
#define  PDM__MPI_G_NUM MPI_LONG
#define  PDM__PDM_MPI_G_NUM PDM_MPI_LONG
#define  PDM_G_NUM_MAX LONG_MAX
#else
typedef int PDM_g_num_t;
#define  PDM__MPI_G_NUM MPI_INT
#define  PDM__PDM_MPI_G_NUM PDM_MPI_INT
#define  PDM_G_NUM_MAX INT_MAX
#endif

typedef  double PDM_real_t;
#define  PDM__MPI_REAL MPI_DOUBLE
#define  PDM__PDM_MPI_REAL PDM_MPI_DOUBLE

/**
 * \enum PDM_l_num_t
 * \brief Long int in pdm
 *
 */

typedef int PDM_l_num_t;
#define  PDM__MPI_L_NUM MPI_INT
#define  PDM__PDM_MPI_L_NUM PDM_MPI_INT


/**
 * \enum PDM_data_t
 * \brief Type of data
 *
 */

typedef enum {

  PDM_INT    = 0,  /*!< Integer */
  PDM_DOUBLE = 1,   /*!< Double */
  PDM_WRONG_DATA = -1   /*!< Double */

} PDM_data_t;

/**
 * \enum PDM_stride_t
 * \brief The notion of "stride" represents the number of field components. There are 3 modes:
 *    - PDM_STRIDE_CST_INTERLACED   : The number of components is constant for each element. The field is stored according to this pattern (c_1,1 ... c_s,1 ... c_1,n ... c_s,n), whhere 's' is the stride and 'n' the number of field elements
 *    - PDM_STRIDE_CST_INTERLEAVED  : The number of components is constant for each element. The field is stored according to this pattern (c_1,1 ... c_1,n ... c_s,1 ... c_s,n), whhere 's' is the stride and 'n' the number of field elements
 *    - PDM_STRIDE_VAR_INTERLACED   : The number of components is variable for each element. The field is stored according to this pattern (c_1,1 ... c_s1,1 ... c_1,n ... c_sn,n), whhere 's_i' is the 'i' element stride and 'n' the number of field elements
 *
 */

typedef enum {

  PDM_STRIDE_CST_INTERLACED = 0, /*!< Constant stride interlaced */
  PDM_STRIDE_CST_INTERLEAVED= 1, /*!< Constant stride interleaved */
  PDM_STRIDE_VAR_INTERLACED = 2  /*!< Variable stride interlaced */

} PDM_stride_t;


/**
 * \enum PDM_mpi_comm_t
 * \brief Framework used for MPI communications
 *
 */

typedef enum {

  PDM_MPI_COMM_KIND_P2P        = 0, /*!< Peer-to-peer (MPI_Isend/MPI_Irecv) */
  PDM_MPI_COMM_KIND_COLLECTIVE = 1, /*!< Collective communications (MPI_Ialltoall, ...) */
  PDM_MPI_COMM_KIND_WIN_RMA    = 2, /*!< RMA windows (MPI_Put, MPI_GET, ...) */
  /** \cond INTERNAL_USE */
  PDM_MPI_COMM_KIND_INVALID    = 3  /*!< Initial state, not for users */
  /** \endcond */

} PDM_mpi_comm_kind_t;


typedef enum {
  PDM_EXCHANGE_DIRECTION_SEND = 0, /*!< Trigger a MPI_Send like communication */
  PDM_EXCHANGE_DIRECTION_RECV = 1, /*!< Trigger a MPI_Recv like communication */
} PDM_exchange_direction_t;


/**
 * \enum PDM_bool_t
 * \brief Bool type
 *
 */

typedef enum {

  PDM_FALSE = 0, /*!< False */
  PDM_TRUE  = 1  /*!< True  */

} PDM_bool_t;


/**
 * \enum PDM_mesh_entities_t
 * \brief Mesh entities
 *
 */

typedef enum {

  PDM_MESH_ENTITY_CELL = 0,  /*!< Cell entity  */
  PDM_MESH_ENTITY_FACE = 1,  /*!< Face entity  */
  PDM_MESH_ENTITY_EDGE = 2,  /*!< Edge entity  */
  PDM_MESH_ENTITY_VTX  = 3,  /*!< Vertex entity  */
  PDM_MESH_ENTITY_MAX  = 4   /*!<  */

} PDM_mesh_entities_t;


/**
 * \enum PDM_mesh_nature_t
 * \brief Mesh natures
 *
 */

typedef enum {

  PDM_MESH_NATURE_NODAL_SHARED   = 0,  /*!< Shared PDM_mesh_nodal  */
  PDM_MESH_NATURE_MESH_SETTED    = 1   /*!< PDm_surface_mesh  */

} PDM_mesh_nature_t;


/**
 * \enum PDM_ownership_t
 * \brief Bool type
 *
 */

typedef enum {

  PDM_OWNERSHIP_KEEP                 = 0, /*!< paradigm will free results */
  PDM_OWNERSHIP_USER                 = 1, /*!< Ownership is gives to user  */
  PDM_OWNERSHIP_UNGET_RESULT_IS_FREE = 2, /*!< Free all memory that not be getted by user */
  PDM_OWNERSHIP_BAD_VALUE            = 3  /*!< Wrong  */

} PDM_ownership_t;

/**
 * \enum PDM_connectivity_type_t
 * \brief Mesh connectivity
 *
 */

typedef enum {
  PDM_CONNECTIVITY_TYPE_CELL_ELMT   = 0,    /*!< cell->element connectivity                                     */
  PDM_CONNECTIVITY_TYPE_CELL_CELL   = 1,    /*!< cell->cell connectivity                                        */
  PDM_CONNECTIVITY_TYPE_CELL_FACE   = 2,    /*!< cell->face connectivity                                        */
  PDM_CONNECTIVITY_TYPE_CELL_EDGE   = 3,    /*!< cell->edge connectivity                                        */
  PDM_CONNECTIVITY_TYPE_CELL_VTX    = 4,    /*!< cell->vertex connectivity                                      */
  PDM_CONNECTIVITY_TYPE_FACE_ELMT   = 5,    /*!< face->element connectivity                                     */
  PDM_CONNECTIVITY_TYPE_FACE_CELL   = 6,    /*!< face->cell connectivity                                        */
  PDM_CONNECTIVITY_TYPE_FACE_FACE   = 7,    /*!< face->face connectivity                                        */
  PDM_CONNECTIVITY_TYPE_FACE_EDGE   = 8,    /*!< face->edge connectivity                                        */
  PDM_CONNECTIVITY_TYPE_FACE_VTX    = 9,    /*!< face->vertex connectivity                                      */
  PDM_CONNECTIVITY_TYPE_EDGE_ELMT   = 10,   /*!< edge->element connectivity                                     */
  PDM_CONNECTIVITY_TYPE_EDGE_CELL   = 11,   /*!< edge->cell connectivity                                        */
  PDM_CONNECTIVITY_TYPE_EDGE_FACE   = 12,   /*!< edge->face connectivity                                        */
  PDM_CONNECTIVITY_TYPE_EDGE_EDGE   = 13,   /*!< edge->edge connectivity                                        */
  PDM_CONNECTIVITY_TYPE_EDGE_VTX    = 14,   /*!< edge->vertex connectivity                                      */
  PDM_CONNECTIVITY_TYPE_VTX_ELMT    = 15,   /*!< vertex->element connectivity                                   */
  PDM_CONNECTIVITY_TYPE_VTX_CELL    = 16,   /*!< vertex->cell connectivity                                      */
  PDM_CONNECTIVITY_TYPE_VTX_FACE    = 17,   /*!< vertex->face connectivity                                      */
  PDM_CONNECTIVITY_TYPE_VTX_EDGE    = 18,   /*!< vertex->edge connectivity                                      */
  PDM_CONNECTIVITY_TYPE_VTX_VTX     = 19,   /*!< vertex->vertex connectivity                                    */
  PDM_CONNECTIVITY_TYPE_ELMT_CELL   = 20,   /*!< element->cell connectivity                                     */
  PDM_CONNECTIVITY_TYPE_ELMT_FACE   = 21,   /*!< element->face connectivity                                     */
  PDM_CONNECTIVITY_TYPE_ELMT_EDGE   = 22,   /*!< element->edge connectivity                                     */
  PDM_CONNECTIVITY_TYPE_ELMT_VTX    = 23,   /*!< element->vertex connectivity                                   */
  PDM_CONNECTIVITY_TYPE_MAX         = 24    /*!< enumerator with the maximal integer refering to a connectivity */
} PDM_connectivity_type_t;

typedef enum {
  PDM_BOUND_TYPE_ELMT   = 0,
  PDM_BOUND_TYPE_CELL   = 1,
  PDM_BOUND_TYPE_FACE   = 2,
  PDM_BOUND_TYPE_EDGE   = 3,
  PDM_BOUND_TYPE_VTX    = 4,
  PDM_BOUND_TYPE_MAX    = 5,
} PDM_bound_type_t;

typedef enum {
  PDM_GEOMETRY_KIND_VOLUMIC  = 0,
  PDM_GEOMETRY_KIND_SURFACIC = 1,
  PDM_GEOMETRY_KIND_RIDGE    = 2,
  PDM_GEOMETRY_KIND_CORNER   = 3,
  PDM_GEOMETRY_KIND_MAX      = 4,
} PDM_geometry_kind_t;

typedef enum {

  PDM_MESH_LOCATION_OCTREE,
  PDM_MESH_LOCATION_DBBTREE,
  PDM_MESH_LOCATION_LOCATE_ALL_TGT,
  PDM_MESH_LOCATION_DOCTREE

} PDM_mesh_location_method_t;

typedef enum {

  PDM_DOCTREE_LOCAL_TREE_OCTREE,
  PDM_DOCTREE_LOCAL_TREE_LINEAR_OCTREE,
  PDM_DOCTREE_LOCAL_TREE_KDTREE,

} PDM_doctree_local_tree_t;

typedef enum {
  PDM_TREE_SOLICITATION_BOXES_POINTS,
  PDM_TREE_SOLICITATION_BOXES_BOXES
} PDM_tree_solicitation_t;


typedef enum {
  PDM_VTX_KIND_NONE               = 0x0000000, /*=* empty flag, all values set to false */
  PDM_VTX_KIND_NULL               = 0xFFFFFFF, /*=* unsignificant flags value           */
  PDM_VTX_KIND_ON_VOLUME          = 0x0000001, /*=*  0                                  */
  PDM_VTX_KIND_ON_SURFACE         = 0x0000002, /*=*  1                                  */
  PDM_VTX_KIND_ON_RIDGE           = 0x0000004, /*=*  2                                  */
  PDM_VTX_KIND_ON_CORNER          = 0x0000008, /*=*  3                                  */
} PDM_vtx_kind;

typedef enum {
  PDM_REDUCE_OP_MIN,
  PDM_REDUCE_OP_MAX,
  PDM_REDUCE_OP_SUM//,
  // PDM_REDUCE_OP_MEAN
} PDM_reduce_op_t;


typedef enum {

  PDM_DOMAIN_INTERFACE_MULT_NO  = 0,  /*!< Each interface involves only 2 domains */
  PDM_DOMAIN_INTERFACE_MULT_YES = 1,  /*!< Each interface involves several domains */

} PDM_domain_interface_mult_t;


/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_SPLIT_DUAL_WITH_PARMETIS = 1, /*!< Use the <a href="https://github.com/KarypisLab/ParMETIS">ParMETIS</a> graph partitioning library */
  PDM_SPLIT_DUAL_WITH_PTSCOTCH = 2, /*!< Use the <a href="https://gitlab.inria.fr/scotch/scotch">PT-Scotch</a> graph partitioning library */
  PDM_SPLIT_DUAL_WITH_HILBERT  = 3, /*!< Use in-house method based on the <a href="https://en.wikipedia.org/wiki/Hilbert_curve">Hilbert space-filling</a> curve */
  PDM_SPLIT_DUAL_WITH_IMPLICIT = 4  /*!< Split into contiguous chunks of global IDs */
} PDM_split_dual_t;

/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_FIELD_KIND_SCALAR     = 1,
  PDM_FIELD_KIND_COORDS     = 2,
  PDM_FIELD_KIND_VECTOR     = 3,
  PDM_FIELD_KIND_TENSOR_SYM = 4
} PDM_field_kind_t;


/**
 * \enum PDM_iso_surface_kind_t
 * \brief Type of iso-surface
 *
 */
typedef enum {
  PDM_ISO_SURFACE_KIND_FIELD    = 0, /*!< (\c Isosurface.FIELD in Python).     Isosurface follows given partitioned (or block-distributed) field */
  PDM_ISO_SURFACE_KIND_PLANE    = 1, /*!< (\c Isosurface.PLANE in Python).     Isosurface follows the plane of equation (4 coefficients) : \n \f$\scriptsize\phi(x,y,z) = \texttt{coeff[0]} \cdot x + \texttt{coeff[1]} \cdot y + \texttt{coeff[2]} \cdot z - \texttt{coeff[3]}\f$                                                                                                                                                                           */
  PDM_ISO_SURFACE_KIND_SPHERE   = 2, /*!< (\c Isosurface.SPHERE in Python).    Isosurface follows the sphere of equation (4 coefficients) : \n \f$\scriptsize\phi(x,y,z) = (x - \texttt{coeff[0]})^2 + (y - \texttt{coeff[1]})^2 + (z - \texttt{coeff[2]})^2 - \texttt{coeff[3]}^2\f$                                                                                                                                                                        */
  PDM_ISO_SURFACE_KIND_ELLIPSE  = 3, /*!< (\c Isosurface.ELLIPSE in Python).   Isosurface follows the ellipse of equation (7 coefficients) : \n \f$\scriptsize\phi(x,y,z) = \left(\frac{x - \texttt{coeff[0]}}{\texttt{coeff[3]}}\right)^2 + \left(\frac{y - \texttt{coeff[1]}}{\texttt{coeff[4]}}\right)^2 + \left(\frac{z - \texttt{coeff[2]}}{\texttt{coeff[5]}}\right)^2 - \texttt{coeff[6]}^2\f$                                                        */
  PDM_ISO_SURFACE_KIND_QUADRIC  = 4, /*!< (\c Isosurface.QUADRIC in Python).   Isosurface follows the quadric of equation (10 coefficients) : \n \f$\scriptsize\phi(x,y,z) = \texttt{coeff[6]} \left(\frac{x - \texttt{coeff[0]}}{\texttt{coeff[3]}}\right)^2 + \texttt{coeff[7]} \left(\frac{y - \texttt{coeff[1]}}{\texttt{coeff[4]}}\right)^2 + \texttt{coeff[8]} \left(\frac{z - \texttt{coeff[2]}}{\texttt{coeff[5]}}\right)^2 - \texttt{coeff[9]}^2\f$ */
  PDM_ISO_SURFACE_KIND_HEART    = 5, /*!< (\c Isosurface.HEART in Python).     Isosurface follows a heart-shaped function                          */
  PDM_ISO_SURFACE_KIND_FUNCTION = 6, /*!< (\c Isosurface.FUNCTION in Python).  Isosurface follows a user-defined function                          */
  PDM_ISO_SURFACE_KIND_MAX      = 7  /*!< Enumerator to list the number of iso-surface kinds available */
} PDM_iso_surface_kind_t;

typedef enum {
  PDM_CELL_TO_VTX_INTERP_KIND_IDW  = 0, /*!< Inverse Distance Weighting    */
  PDM_CELL_TO_VTX_INTERP_KIND_RBF  = 1, /*!< Radial Basis Function         */
  PDM_CELL_TO_VTX_INTERP_KIND_LSQ  = 2, /*!< Least Square                  */
  PDM_CELL_TO_VTX_INTERP_KIND_USER = 3, /*!< User, we must define callback */
} PDM_cell_to_vtx_interp_kind_t;


typedef enum {

  PDM_EXTRACT_PART_KIND_LOCAL         = 0, /*!< The entities of interest are selected locally and are extracted in place */
  PDM_EXTRACT_PART_KIND_REEQUILIBRATE = 1, /*!< The entities of interest are selected locally and are redistributed according to the chosen partitioning method */
  PDM_EXTRACT_PART_KIND_FROM_TARGET   = 2, /*!< The requested target entities are fetched from their location in the input mesh */
} PDM_extract_part_kind_t;

typedef enum {
  PDM_MESH_INTERSECTION_KIND_PREPROCESS    = 0, /*! Only candidate selection and load balancing of elementary intersection tasks */
  PDM_MESH_INTERSECTION_KIND_WEIGHT        = 1, /*! Resulting is only relative weight of A inside B                              */
  PDM_MESH_INTERSECTION_KIND_UNMERGED_POLY = 2, /*! Resulting is the union of disjoint intersection polyhedra                    */
  PDM_MESH_INTERSECTION_KIND_MESH          = 3, /*! Resulting is new mesh                                                        */
} PDM_mesh_intersection_kind_t;

typedef enum {

  PDM_CLUSTERING_KIND_OCTREE                = 0, /*!< Clustering by parallel octree */
  PDM_CLUSTERING_KIND_OCTREE_WITH_SMOOTHING = 1, /*!< Clustering by parallel octree with smoothing */

} PDM_clustering_kind_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Helper to get entity type according to a connectivity
 *
 */
PDM_mesh_entities_t
PDM_connectivity_type_to_entity_type
(
 PDM_connectivity_type_t   connectivity_type
);

/**
 * \brief Helper to get entity type according to a bound
 *
 */
PDM_mesh_entities_t
PDM_bound_type_to_entity_type
(
 PDM_bound_type_t   bound_type
);

/**
 * \brief Helper to get bound type according to an entity
 *
 */
PDM_bound_type_t
PDM_entity_type_to_bound_type
(
 PDM_mesh_entities_t  entity_type
);

/**
 * \brief Helper to get entity type according to a geometry kind
 *
 */
PDM_mesh_entities_t
PDM_geometry_kind_to_entity_type
(
 PDM_geometry_kind_t   geom_kind
);


/**
 * \brief Helper to get dimension type according to a geometry kind
 *
 */
int
PDM_geometry_kind_to_dimension
(
 PDM_geometry_kind_t   geom_kind
);

/**
 * \brief Helper to get geometry kind according to an entity type
 *
 */
PDM_geometry_kind_t
PDM_entity_type_to_geometry_kind
(
 PDM_mesh_entities_t   entity_type
);

/**
 * \brief Helper to get entity type according to a dimension
 *
 */
PDM_mesh_entities_t
PDM_dimension_to_entity_type
(
  const int dim
);

/**
 * \brief Helper to get geometry kind according to a dimension
 *
 */
PDM_geometry_kind_t
PDM_dimension_to_geometry_kind
(
  const int dim
);

/**
 * \brief Convert a pair of entity types to its corresponding connectivity type
 *
 * \param [in]  connectivity_type   Connectivity type
 * \param [out] entity_type1        First entity type
 * \param [out] entity_type2        Second entity type
 */
PDM_connectivity_type_t
PDM_entity_pair_to_connectivity_type
(
  PDM_mesh_entities_t entity_type1,
  PDM_mesh_entities_t entity_type2
);

/**
 * \brief Convert a connectivity type to its corresponding pair of entity types
 *
 * \param [in]  connectivity_type   Connectivity type
 * \param [out] entity_type1        First entity type
 * \param [out] entity_type2        Second entity type
 *
 * \return 1 if is a valid connectivity type, 0 else
 */
int
PDM_connectivity_type_to_entity_pair
(
  PDM_connectivity_type_t  connectivity_type,
  PDM_mesh_entities_t     *entity_type1,
  PDM_mesh_entities_t     *entity_type2
);

/**
 * \brief Finalize PDM
 *
 * This function frees all allocated global variables
 *
 */

void
PDM_Finalize
(
void
);


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  /* __PDM_H__ */
