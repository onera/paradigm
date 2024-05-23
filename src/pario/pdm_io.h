/*
 * \file
 */

#ifndef __PDM_IO_H__
#define __PDM_IO_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

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
 * Types des données
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_T_INT,
  PDM_IO_T_LONG,
  PDM_IO_T_DOUBLE,
  PDM_IO_T_FLOAT,
  PDM_IO_T_CHAR

} PDM_io_type_t;

/*----------------------------------------------------------------------------
 * Types de suffixe
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_SUFF_AUTO,
  PDM_IO_SUFF_MAN

} PDM_io_suff_t;

/*----------------------------------------------------------------------------
 * Types d'entrees/sorties paralleles
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_KIND_MPIIO_EO,
  PDM_IO_KIND_MPIIO_IP,
  PDM_IO_KIND_MPI_SIMPLE,
  PDM_IO_KIND_SEQ

} PDM_io_kind_t;

/*----------------------------------------------------------------------------
 * Mode d'acces lecture, ecriture, lecture/ecriture
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_MOD_READ,
  PDM_IO_MOD_WRITE,
  PDM_IO_MOD_APPEND

} PDM_io_mod_t;


/*----------------------------------------------------------------------------
 * Format du fichier
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_FMT_TXT,
  PDM_IO_FMT_BIN

} PDM_io_fmt_t;

/*----------------------------------------------------------------------------
 * Presence ou non de l'entete
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_BACKUP_ON,
  PDM_IO_BACKUP_OFF

} PDM_io_backup_t;

/*----------------------------------------------------------------------------
 * Type de rangement des donnees (par bloc ou entrelace)
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_BIGENDIAN,
  PDM_IO_LITTLEENDIAN,
  PDM_IO_NATIVE

} PDM_io_endian_t;

/*----------------------------------------------------------------------------
 * Mode de positionnement dans le fichier
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_SEEK_SET,   /* Position a partir du debut du fichier */
  PDM_IO_PAR_SEEK_CUR,   /* Position a partir de la position courante */
  PDM_IO_PAR_SEEK_END    /* Position a partir de la fin du fichier */

} PDM_io_seek_t;

/*----------------------------------------------------------------------------
 * Type decrivant un fichier de type parallele io
 *----------------------------------------------------------------------------*/

typedef struct _PDM_io_file_t PDM_io_file_t;

/*=============================================================================
 * Variables globales
 *============================================================================*/


/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/

/**
 * \brief Return the file name (or NULL if no file)
 *
 * \param [in]  fichier   Pointer to \ref PDM_io_file_t object
 *
 * \return   Name of the file
 *
 */

const char* PDM_io_file_name_get
(
 PDM_io_file_t *fichier
 );


/**
 * \brief Open a file for parallel access
 *
 * \param [in]  nom                 File name
 * \param [in]  fmt                 Text of Binary format
 * \param [in]  suff_t              Type of suffix (manual or automatic)
 * \param [in]  suff_u              Suffix (if manual)
 * \param [in]  s_backup            Activates the backup of a pre-existing file in write mode
 * \param [in]  acces               Type (parallel with MPI-IO, parallel without MPI-IO, sequential)
 * \param [in]  mode                Access mode (read, write, read/write)
 * \param [in]  endian              Endian type (big, little or native)
 * \param [in]  comm                Communicator associated to the file
 * \param [in]  prop_noeuds_actifs  Proportion of active nodes
 * \param [out] unite               Unit of the file
 * \param [out] ierr                Indicates whether the file is of type PDM_io or not (for read-only opening only)
 *
 */

void PDM_io_open
(
 const char             *nom,
 const PDM_io_fmt_t      fmt,
 const PDM_io_suff_t     suff_t,
 const char             *suff_u,
 const PDM_io_backup_t   s_backup,
 const PDM_io_kind_t     acces,
 const PDM_io_mod_t      mode,
 const PDM_io_endian_t   endian,
 PDM_MPI_Comm            comm,
 double                  prop_noeuds_actifs,
 PDM_io_file_t         **unite,
 PDM_l_num_t            *ierr
);


/**
 * \brief Set the file position indicator
 *
 * \param [in] fichier         Pointer to \ref PDM_io_file_t object
 * \param [in] offset          Address
 * \param [in] seek            Origin type
 *
 */

void PDM_io_seek
(
 PDM_io_file_t    *fichier,
 const PDM_g_num_t    offset,
 const PDM_io_seek_t  seek
);


/**
 * \brief Return the current file position
 *
 * \param [in] fichier         Pointer to \ref PDM_io_file_t object
 *
 * \return   Current position in file
 *
 */

PDM_g_num_t
PDM_io_tell
(
 PDM_io_file_t   *fichier
);


/**
 * \brief Global read: the master process alone accesses the
 *        file and redistributes the information to all the communicator's processes
 *
 * \param [in]  fichier         Pointer to \ref PDM_io_file_t object
 * \param [in]  taille_donnee   Size of a unit piece of data
 * \param [in]  n_donnees       Amount of data to be read
 * \param [out] donnees         Read data
 *
 */

void PDM_io_global_read
(
 PDM_io_file_t  *fichier,
 const PDM_l_num_t  taille_donnee,
 const PDM_g_num_t  n_donnees,
 void              *donnees
 );


/**
 * \brief Global write: The master process has sole access to the file
 *
 * \param [in]  fichier         Pointer to \ref PDM_io_file_t object
 * \param [in]  taille_donnee   Size of a unit piece of data
 * \param [in]  n_donnees       Amount of data to write
 * \param [in]  donnees         Data to write
 *
 */

void PDM_io_global_write
(
 PDM_io_file_t  *fichier,
 const PDM_l_num_t  taille_donnee,
 const PDM_g_num_t  n_donnees,
 const void        *donnees
);


/**
 * \brief Parallel reading of data blocks followed by
 *        redistribution of the data according to indirection
 *
 * \param [in]  fichier          Pointer to \ref PDM_io_file_t object
 * \param [in]  t_n_composantes  Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
 * \param [in]  n_composantes    Number of components for each piece of data
 * \param [in]  taille_donnee    Unit size of a piece of data
 * \param [in]  n_donnees        Number of data items to be read
 * \param [in]  indirection      Indirection of data redistribution
 * \param [out] donnees          Read data
 *
 */

void PDM_io_par_interlaced_read
(
 PDM_io_file_t             *fichier,
 const PDM_stride_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t            *indirection,
 void                         *donnees
 );


/**
 * \brief Parallel reading of data blocks
 *        The blocks must be arranged in ascending order
 *        according to the numbering of the processes
 *
 * \param [in]  fichier          Pointer to \ref PDM_io_file_t object
 * \param [in]  t_n_composantes  Component size type (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
 * \param [in]  n_composantes    Number of components for each data item
 * \param [in]  taille_donnee    Unit size of a piece of data
 * \param [in]  n_donnees        Number of data items to be read
 * \param [in]  debut_bloc       Relative address of start of block
 * \param [out] donnees          Read data
 *
 */

void PDM_io_par_block_read
(
 PDM_io_file_t             *fichier,
 const PDM_stride_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t             debut_bloc,
 void                         *donnees
);


/**
 * \brief Data sorted according to indirection, then parallel write of data blocks
 *
 * \param [in] fichier           Pointer to \ref PDM_io_file_t object
 * \param [in] t_n_composantes   Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
 * \param [in] n_composantes     Number of components for each data item
 * \param [in] taille_donnee     Unit size of the data
 * \param [in] n_donnees         Number of data items to be written
 * \param [in] indirection       Data redistribution direction
 * \param [in] donnees           Data to be written
 *
 */

void PDM_io_par_interlaced_write
(
 PDM_io_file_t             *fichier,
 const PDM_stride_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t            *indirection,
 const void                   *donnees
);


/**
 * \brief Parallel writing of data blocks
 *        Blocks must be arranged in ascending order according
 *        to numbering of the processes
 *
 * \param [in] fichier           Pointer to \ref PDM_io_file_t object
 * \param [in] t_n_composantes   Type of component sizes (PDM_STRIDE_CST_)INTERLACED or PDM_STRIDE_VAR_INTERLACED)
 * \param [in] n_composantes     Number of components for each data item
 * \param [in] taille_donnee     Unit size of the data
 * \param [in] n_donnees         Number of data to read
 * \param [in] debut_bloc        Relative address of start of block
 * \param [in] donnees           Data to be written
 *
 */

void PDM_io_par_block_write
(
 PDM_io_file_t             *fichier,
 const PDM_stride_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t             debut_bloc,
 const void                   *donnees
);


/**
 * \brief Closing the file without destroying the PDM_io structure
 *        associated with unit
 *
 * \param [in] fichier           Pointer to \ref PDM_io_file_t object
 *
 */

void PDM_io_close
(
 PDM_io_file_t   *fichier
);


/**
 * \brief Free of the PDM_io structure associated with the unit
 *
 * \param [in] fichier           Pointer to \ref PDM_io_file_t object
 *
 */

void PDM_io_free
(
 PDM_io_file_t   *fichier
);


/**
 * \brief Returns the cumulative files access time
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 * \param [out] t_cpu             CPU time
 * \param [out] t_elapsed         Elapsed time
 *
 */

void PDM_io_get_timer_fichier
(
 PDM_io_file_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Returns the cumulative time for data swap
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 * \param [out] t_cpu             CPU time
 * \param [out] t_elapsed         Elapsed time
 *
 */

void PDM_io_timer_swap_endian_get
(
 PDM_io_file_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Returns the cumulative time for data distribution
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 * \param [out] t_cpu             CPU time
 * \param [out] t_elapsed         Elapsed time
 *
 */

void PDM_io_timer_distrib_get
(
 PDM_io_file_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Returns the total cumulative time
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 * \param [out] t_cpu             CPU time
 * \param [out] t_elapsed         Elapsed time
 *
 */

void PDM_io_timer_total_get
(
 PDM_io_file_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Shows file information
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 *
 */

void PDM_io_dump
(
 PDM_io_file_t   *fichier
);


/**
 * \brief Returns the file communicator
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 * \param [out] pdm_mpi_comm      MPI communicator
 *
 */

void PDM_io_comm_get
(
 PDM_io_file_t *fichier,
 PDM_MPI_Comm     *pdm_mpi_comm
);


/**
 * \brief Activate endian swap
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 *
 */

void PDM_io_swap_endian_on
(
 PDM_io_file_t   *fichier
);


/**
 * \brief Deactivate endian swap
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 *
 */

void PDM_io_swap_endian_off
(
 PDM_io_file_t   *fichier
);


/**
 * \brief Swap endian pour conversion little endian <-> big endian
 *
 * \param [in]  taille_donnee   Size of a unit piece of data
 * \param [in]  n_donnee        Amount of data
 * \param [in]  donnees         Data
 * \param [out] resultats       Result
 *
 */

void PDM_io_swap_endian
(
 const size_t   taille_donnee,
 const size_t   n_donnees,
 const void    *donnees,
 void          *resultats
 );


/**
 * \brief Defines the format of the individual data for text output
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_file_t object
 * \param [in]  n_char_fmt        Number of characters in the format
 * \param [in]  data_type         Type of data
 * \param [in]  fmt               Format
 *
 */

void PDM_io_fmt_data_set
(
 PDM_io_file_t    *fichier,
 const PDM_l_num_t    n_char_fmt,
 const PDM_io_type_t  data_type,
 const char          *fmt
);


/**
 * \brief Create a directory
 *
 * \param [in] path   Path to new directory
 *
 * \return 0 if successful, -1 else
 *
 */

int PDM_io_mkdir
(
 const char* path
);


/**
 * \brief Calculating the total size of a data field
 *
 * \param [in]  fichier          Pointer to \ref PDM_io_file_t object
 * \param [in]  t_n_composantes  Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
 * \param [in]  n_composantes    Number of components for each data
 * \param [in]  n_donnees        Number of data
 * \param [in]  indirection      Data redistribution direction
 *
 * \return   Total size of a data field
 *
 */

PDM_g_num_t
PDM_io_n_data_get
(
 PDM_io_file_t             *fichier,
 const PDM_stride_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t            *indirection
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_IO_H__ */
