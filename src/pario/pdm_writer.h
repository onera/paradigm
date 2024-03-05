/*
 * \file
 */

#ifndef __PDM_WRITER_H__
#define __PDM_WRITER_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"

/*=============================================================================
 * Definitions des macro
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Statut
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_OFF,
  PDM_WRITER_ON

} PDM_writer_status_t;

/*----------------------------------------------------------------------------
 * Type de topologie
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_TOPO_CST ,
  PDM_WRITER_TOPO_DEFORMABLE,
  PDM_WRITER_TOPO_VARIABLE

} PDM_writer_topology_t;

/*----------------------------------------------------------------------------
 * Type d'elements géometriques (It's the same than the type defined into PDM_Mesh_nodal)
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_POINT = PDM_MESH_NODAL_POINT,
  PDM_WRITER_BAR2 = PDM_MESH_NODAL_BAR2,
  PDM_WRITER_TRIA3 = PDM_MESH_NODAL_TRIA3,
  PDM_WRITER_QUAD4 = PDM_MESH_NODAL_QUAD4,
  PDM_WRITER_POLY_2D = PDM_MESH_NODAL_POLY_2D,
  PDM_WRITER_TETRA4 = PDM_MESH_NODAL_TETRA4,
  PDM_WRITER_PYRAMID5 = PDM_MESH_NODAL_PYRAMID5,
  PDM_WRITER_PRISM6 = PDM_MESH_NODAL_PRISM6,
  PDM_WRITER_HEXA8 = PDM_MESH_NODAL_HEXA8,
  PDM_WRITER_POLY_3D  = PDM_MESH_NODAL_POLY_3D

} PDM_writer_elt_geom_t;


typedef struct _PDM_writer_t PDM_writer_t;

typedef struct _PDM_writer_geom_t PDM_writer_geom_t;

typedef struct _PDM_writer_var_t PDM_writer_var_t;

typedef struct _PDM_writer_cst_global_var_t PDM_writer_cst_global_var_t;


/**
 * \struct PDM_writer_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_t_t instance
 *
 */

typedef void (*PDM_writer_fct_t) (PDM_writer_t *pw);


/**
 * \struct PDM_writer_geom_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_t_geom_t instance
 *
 */

typedef void (*PDM_writer_geom_fct_t) (PDM_writer_geom_t *geom);


/**
 * \struct PDM_writer_var_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_t_var_t instance
 *
 */

typedef void (*PDM_writer_var_fct_t) (PDM_writer_var_t *var);


/*----------------------------------------------------------------------------
 * Format du fichie de sortie
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_FMT_BIN,
  PDM_WRITER_FMT_ASCII

} PDM_writer_fmt_fic_t;

/*----------------------------------------------------------------------------
 * Types des variables
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_VAR_CST           = 0,
  PDM_WRITER_VAR_SCALAR        = 1,
  PDM_WRITER_VAR_VECTOR        = 3,
  PDM_WRITER_VAR_TENSOR_SYM    = 6,
  PDM_WRITER_VAR_TENSOR_ASYM   = 9

} PDM_writer_var_dim_t;

/*----------------------------------------------------------------------------
 * Localisation des variables
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_VAR_VERTICES,
  PDM_WRITER_VAR_ELEMENTS,
  PDM_WRITER_VAR_PARTICLES

} PDM_writer_var_loc_t;


/*=============================================================================
 * Variables globales
 *============================================================================*/

/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/

/**
 *
 * \brief Create a structure for parallel writing of geometry and associated variables
 *
 * \param [in] fmt                  Output format
 * \param [in] fmt_fic              Binary or ASCII
 * \param [in] topologie            Indicates whether the mesh is mobile
 * \param [in] st_reprise           Finalizes previous outputs before restart
 * \param [in] rep_sortie           Output repository
 * \param [in] nom_sortie           Output filename
 * \param [in] pdm_mpi_com          MPI communicator
 * \param [in] acces                Access type
 * \param [in] prop_noeuds_actifs   Amount of active nodes:
 *                                    *  -1 : all active
 *                                    *   1 : one process per node
 *                                    * 0 < val < 1 : one process per active node
 * \param [in] options              Complementary options for the format structured as
 *                                  ("name_1 = val_1 : ... : name_n = val_n")
 *
 * \return   Pointer to \ref PDM_writer_t object
 *
 */

PDM_writer_t *
PDM_writer_create
(
const char                   *fmt,
const PDM_writer_fmt_fic_t    fmt_fic,
const PDM_writer_topology_t  topologie,
const PDM_writer_status_t     st_reprise,
const char                   *rep_sortie,
const char                   *nom_sortie,
const PDM_MPI_Comm            pdm_mpi_comm,
const PDM_io_kind_t          acces,
const double                  prop_noeuds_actifs,
const char                   *options
);


/**
 * \brief Free a writer structure
 *
 * \param [in] cs    Pointer to \ref PDM_writer_t object
 *
 */

void
PDM_writer_free
(
 PDM_writer_t *cs
);


/**
 * \brief Begin a time step
 *
 * \param [in] cs             Pointer to \ref PDM_writer_t object
 * \param [in] physical_time  Time
 */

void
PDM_writer_step_beg
(
 PDM_writer_t  *cs,
 const double   physical_time
);


/**
 * \brief Is there a open time step?
 *
 * \param [in] cs             Pointer to \ref PDM_writer_t object
 *
 * \return   Flag which indicates if a time step is open
 */

int
PDM_writer_is_open_step
(
 PDM_writer_t  *cs
);


/**
 * \brief End a time step
 *
 * \param [in] cs             Pointer to \ref PDM_writer_t object
 *
 */

void
PDM_writer_step_end
(
 PDM_writer_t  *cs
);


/**
 * \brief Create a new geometry in the writer structure
 *
 * \param [in]  cs                Pointer to \ref PDM_writer_t object
 * \param [in]  nom_geom          Name of the geometry
 * \param [in]  n_part            Number of partitions
 *
 * \return   Identifier of the geometry within the writer instance
 *
 */

int
PDM_writer_geom_create
(
 PDM_writer_t               *cs,
 const char                 *nom_geom,
 const int                   n_part
);


/**
 * \brief Create a geometry from a nodal mesh structure
 *
 * \param [in]  cs                Pointer to \ref PDM_writer_t object
 * \param [in]  nom_geom          Name of the geometry
 * \param [in]  mesh              Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 * \return   Identifier of the geometry within the writer instance
 *
 */

int
PDM_writer_geom_create_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const char                *nom_geom,
 PDM_part_mesh_nodal_t     *mesh
);


/**
 * \brief Set a geometry from a nodal mesh structure
 *
 * \param [in]  cs                Pointer to \ref PDM_writer_t object
 * \param [in]  id_geom           Geometry identifier
 * \param [in]  mesh              Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */

void
PDM_writer_geom_set_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const int                  id_geom,
 PDM_part_mesh_nodal_t     *mesh
);


/**
 * \brief Define the coordinates of the current partition
 *
 * \param [in] cs        Pointer to \ref PDM_writer_t object
 * \param [in] id_geom   Geometry identifier
 * \param [in] id_part   Partition identifier
 * \param [in] n_som     Number of vertices
 * \param [in] coords    Coordinates (size = 3 * \p n_som)
 * \param [in] numabs    Vertex global numbering (size = \p n_som)
 * \param [in] ownership Ownership
 *
 */

void
PDM_writer_geom_coord_set
(
 PDM_writer_t      *cs,
 const int          id_geom,
 const int          id_part,
 const int          n_som,
 const PDM_real_t  *coords,
 const PDM_g_num_t *numabs,
 const PDM_ownership_t ownership
);


/**
 * \brief Definition of the coordinates of the vertices in the current partition from a parent set
 *
 *
 * \param [in] cs               Pointer to \ref PDM_writer_t object
 * \param [in] id_geom          Geometry identifier
 * \param [in] id_part          Partition identifier
 * \param [in] n_som            Number of vertices
 * \param [in] n_som_parent     Number of parent vertices
 * \param [in] numabs           Vertex global numbering (size = \p n_som)
 * \param [in] num_parent       Vertex parent local numbering (size = \p n_som)
 * \param [in] coords_parent    Coordinates of parent vertices (size = 3 * \p n_som_parent)
 * \param [in] numabs_parent    Vertex parent global numbering (size = \p n_som_parent)
 * \param [in] ownership        Ownership
 *
 */

void
PDM_writer_geom_coord_from_parent_set
(
 PDM_writer_t      *cs,
 const int          id_geom,
 const int          id_part,
 const int          n_som,
 const int          n_som_parent,
 const PDM_g_num_t *numabs,
 const int         *num_parent,
 const PDM_real_t  *coords_parent,
 const PDM_g_num_t *numabs_parent,
 const PDM_ownership_t ownership
);


/**
 * \brief Add a block of elements of a given type
 *
 * \param [in] cs             Pointer to \ref PDM_writer_t object
 * \param [in] id_geom        Geometry identifier
 * \param [in] t_elt          Element type
 * \param [in] owner          Block ownership
 *
 * \return   Block identifier
 *
 */

int
PDM_writer_geom_bloc_add
(
 PDM_writer_t                *cs,
 const int                    id_geom,
 const PDM_writer_elt_geom_t  t_elt,
 const PDM_ownership_t        owner
);


/**
 * \brief Set in the given geometry a block of elements of a given type
 *
 * \param [in] cs                  Pointer to \ref PDM_writer_t object
 * \param [in] id_geom             Geometry identifier
 * \param [in] id_bloc             Block identifier
 * \param [in] id_part             Partition identifier
 * \param [in] n_elt               Number of elements
 * \param [in] connec              Element->Vertex connectivity
 * \param [in] numabs              Element global numbering
 *
 */

void
PDM_writer_geom_bloc_std_set
(
 PDM_writer_t  *cs,
 const int      id_geom,
 const int      id_bloc,
 const int      id_part,
 const int      n_elt,
 PDM_l_num_t   *connec,
 PDM_g_num_t   *numabs
);


/**
 * \brief Add a block of polygons to the current partition
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 * \param [in] id_bloc         Block identifier
 * \param [in] id_part         Partition identifier
 * \param [in] n_elt           Number of elements
 * \param [in] connec_idx      Index of the Element->Vertex connectivity (size = \p n_elt+1)
 * \param [in] connec          Element->Vertex connectivity (size = \p connec_idx[\p n_elt])
 * \param [in] numabs          Element global numbering (size = \p n_elt)
 *
 */

void
PDM_writer_geom_bloc_poly2d_set
(
PDM_writer_t        *cs,
const int            id_geom,
const int            id_bloc,
const int            id_part,
const PDM_l_num_t    n_elt,
      PDM_l_num_t   *connec_idx,
      PDM_l_num_t   *connec,
      PDM_g_num_t   *numabs
);


/**
 * \brief Add a block of polyhedra to the current partition
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 * \param [in] id_bloc         Block identifier
 * \param [in] id_part         Partition identifier
 * \param [in] n_elt           Number of elements
 * \param [in] n_face          Number of faces
 * \param [in] facsom_idx      Index of the Face->Vertex connectivity (size = \p n_face + 1)
 * \param [in] facsom          Face->Vertex connectivity (size = \p facsom_idx[\p n_face])
 * \param [in] cellfac_idx     Index of the Cell->Face connectivity (size = \p n_elt+1)
 * \param [in] cellfac         Cell->Face connectivity (size = \p cellfac_idx[\p n_elt])
 * \param [in] numabs          Cell global numbering (size = \p n_elt)
 *
 */

void
PDM_writer_geom_bloc_poly3d_set
(
PDM_writer_t        *cs,
const int            id_geom,
const int            id_bloc,
const int            id_part,
const PDM_l_num_t    n_elt,
const PDM_l_num_t    n_face,
      PDM_l_num_t   *facsom_idx,
      PDM_l_num_t   *facsom,
      PDM_l_num_t   *cellfac_idx,
      PDM_l_num_t   *cellfac,
      PDM_g_num_t   *numabs
);


/**
 *
 * \brief Add 3D cells described in terms of faces
 *
 * This function determines element types and creates
 * blocks grouping elements of the same type.
 *
 * \param [in]  cs              Pointer to \ref PDM_writer_t object
 * \param [in]  id_geom         Geometry identifier
 * \param [in]  id_part         Partition identifier
 * \param [in]  n_cell          Number of 3D cells
 * \param [in]  n_face          Number of faces
 * \param [in]  face_som_idx    Index of the Face->Vertex connectivity (size = \p n_face + 1)
 * \param [in]  face_som_nb     Number of vertices per face (optional)
 * \param [in]  face_som        Face->Vertex connectivity (size = \p face_som_idx[\p n_face])
 * \param [in]  cell_face_idx   Index of the Cell->Face connectivity (size = \p n_cell + 1)
 * \param [in]  cell_face_nb    Number of faces per cell (optional)
 * \param [in]  cell_face       Cell->Face connectivity (size = \p cell_face_idx[\p n_cell])
 * \param [in]  numabs          Cell global numbering (size = \p n_cell)
 *
 */

void
PDM_writer_geom_cell3d_cellface_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_cell,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_l_num_t  *cell_face_idx,
 PDM_l_num_t  *cell_face_nb,
 PDM_l_num_t  *cell_face,
 PDM_g_num_t  *numabs
 );


/**
 *
 * \brief Add 2D cells described in terms of faces
 *
 * This function determines element types and creates
 * blocks grouping elements of the same type.
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 * \param [in] id_part         Partition identifier
 * \param [in] n_cell          Number of 2D cells
 * \param [in] n_face          Number of faces
 * \param [in] face_som_idx    Index of the Face->Vertex connectivity (unused)
 * \param [in] face_som_nb     Number of vertices per face (unused)
 * \param [in] face_som        Face->Vertex connectivity (size = 2 * \p n_face)
 * \param [in] cell_face_idx   Index of the Cell->Face connectivity (size = \p n_cell + 1)
 * \param [in] cell_face_nb    Number of faces per cell (optional)
 * \param [in] cell_face       Cell->Face connectivity (size = \p cell_face_idx[\p n_cell])
 * \param [in] numabs          Cell global numbering (size = \p n_cell)
 *
 */

void
PDM_writer_geom_cell2d_cellface_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_cell,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_l_num_t  *cell_face_idx,
 PDM_l_num_t  *cell_face_nb,
 PDM_l_num_t  *cell_face,
 PDM_g_num_t  *numabs
);


/**
 *
 * \brief Add faces described in nodal fashion
 *
 * This function determines element types and creates
 * blocks grouping elements of the same type.
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 * \param [in] id_part         Partition identifier
 * \param [in] n_face          Number of faces
 * \param [in] face_som_idx    Index of the Face->Vertex connectivity (size = \p n_face + 1)
 * \param [in] face_som_nb     Number of vertices per face (optional)
 * \param [in] face_som        Face->Vertex connectivity (size = \p face_som_idx[\p n_face])
 * \param [in] numabs          Face global numbering (size = \p n_face)
 *
 */

void
PDM_writer_geom_faces_facesom_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_g_num_t  *numabs
);


/**
 * \brief Write current geometry
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 *
 */

void
PDM_writer_geom_write
(
 PDM_writer_t *cs,
 const int     id_geom
 );


/**
 * \brief Free data describing the current geometry
 *        Indirections on absolute numbering are retained
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 *
 */

void
PDM_writer_geom_data_free
(
 PDM_writer_t *cs,
 const int     id_geom
 );


/**
 * \brief Free data describing the current geometry
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 *
 */

void
PDM_writer_geom_free
(
 PDM_writer_t *cs,
 const int     id_geom
);


/**
 * \brief Create a global constant variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] nom_var         Variable name
 * \param [in] val_var         Variable value
 *
 * \return  Variable identifier
 *
 */

int
PDM_writer_cst_global_var_create
(
 PDM_writer_t               *cs,
 const char                 *nom_var,
 const double                val_var
);


/**
 * \brief Set a global constant variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_var          Variable identifier
 * \param [in] val_var         Variable value
 *
 * \return Variable identifier
 *
 */

void
PDM_writer_cst_global_var_set
(
 PDM_writer_t               *cs,
 const int                   id_var,
 const double                val_var
);



/**
 * \brief Create a variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] st_dep_temps    Indicates whether the variable is time dependent
 * \param [in] dim             Variable's dimension
 * \param [in] loc             Variable's location
 * \param [in] nom_var         Name of the variable
 *
 * \return  Variable identifier
 *
 */

int
PDM_writer_var_create
(
 PDM_writer_t               *cs,
 const PDM_writer_status_t   st_dep_tps,
 const PDM_writer_var_dim_t  dim,
 const PDM_writer_var_loc_t  loc,
 const char                 *nom_var
);


/**
 * \brief Variable name mapping
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] public_name     Public variable name
 * \param [in] pivate_name     Private variable name
 *
 */

void
PDM_writer_name_map_add
(
 PDM_writer_t *cs,
 const char   *public_name,
 const char   *private_name
);


/**
 * \brief Write variable values
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_var          Variable identifier
 *
 */

void
PDM_writer_var_write
(
 PDM_writer_t *cs,
 const int     id_var
);


/**
 * \brief Update variable values
 *
 * \warning the values defined for the elements must be defined in the order in which the blocks are defined!
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_var          Variable identifier
 * \param [in] id_geom         Geometry identifier
 * \param [in] id_part         Partition identifier
 * \param [in] val             Variable values
 *
 */

void
PDM_writer_var_set
(
 PDM_writer_t     *cs,
 const int         id_var,
 const int         id_geom,
 const int         id_part,
 const PDM_real_t *val
);


/**
 * \brief Free variable data arrays
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_var          Variable identifier
 *
 */

void
PDM_writer_var_data_free
(
 PDM_writer_t *cs,
 const int     id_var
);


/**
 * \brief Free variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_var          Variable identifier
 *
 */

void
PDM_writer_var_free
(
 PDM_writer_t *cs,
 const int     id_var
);


/**
 * \brief Add a writer format
 *
 * Define a new format writer
 *
 * \param [in] name            Name
 * \param [in] create_fct      Customize \ref PDM_writer_t_create function for the new format  (or NULL)
 * \param [in] free_fct        Customize \ref PDM_writer_t_free function for the new format (or NULL)
 * \param [in] beg_step_fct    Customize \ref PDM_writer_t_step_beg function for the new format (or NULL)
 * \param [in] end_step_fct    Customize \ref PDM_writer_t_step_end function for the new format (or NULL)
 * \param [in] geom_create_fct Customize \ref PDM_writer_t_geom_create function for the new format (or NULL)
 * \param [in] geom_write_fct  Customize \ref PDM_writer_t_geom_write function for the new format
 * \param [in] geom_free_fct   Customize \ref PDM_writer_t_geom_free function for the new format (or NULL)
 * \param [in] var_create_fct  Customize \ref PDM_writer_t_var_create function for the new format (or NULL)
 * \param [in] var_write_fct   Customize \ref PDM_writer_t_var_write function for the new format
 * \param [in] var_free_fct    Customize \ref PDM_writer_t_var_free function for the new format (or NULL)
 *
 */

void
PDM_writer_fmt_add
(
 const char                  *name,            /*!< Name                                                     */
 const PDM_writer_fct_t       create_fct,      /*!< Customize \ref PDM_writer_t_create function for the format */
 const PDM_writer_fct_t       free_fct,        /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_fct_t       beg_step_fct,    /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_fct_t       end_step_fct,    /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_geom_fct_t  geom_create_fct, /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_geom_fct_t  geom_write_fct,  /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_geom_fct_t  geom_free_fct,   /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_var_fct_t   var_create_fct,  /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_var_fct_t   var_write_fct,   /*!< Customize \ref PDM_writer_t_free function for the format   */
 const PDM_writer_var_fct_t   var_free_fct     /*!< Customize \ref PDM_writer_t_free function for the format   */
);


/**
 * \brief Free format
 *
 */

void
PDM_writer_fmt_free
(
 void
);


/**
 * \brief Reset data describing the current geometry
 *
 * \param [in] cs              Pointer to \ref PDM_writer_t object
 * \param [in] id_geom         Geometry identifier
 *
 */

void
PDM_writer_geom_data_reset
(
 PDM_writer_t *cs,
 const int     id_geom
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_H__ */
