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

  PDM_IO_ACCES_MPIIO_EO,
  PDM_IO_ACCES_MPIIO_IP,
  PDM_IO_ACCES_MPI_SIMPLE,
  PDM_IO_ACCES_SEQ

} PDM_io_acces_t;

/*----------------------------------------------------------------------------
 * Mode d'acces lecture, ecriture, lecture/ecriture
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_MODE_LECTURE,
  PDM_IO_MODE_ECRITURE,
  PDM_IO_MODE_AJOUT

} PDM_io_mode_t;

/*----------------------------------------------------------------------------
 * Nombre de composantes par donnee constant ou variable
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_N_COMPOSANTE_CONSTANT,
  PDM_IO_N_COMPOSANTE_VARIABLE

} PDM_io_n_composantes_t;

/*----------------------------------------------------------------------------
 * Type de rangement des donnees (par bloc ou entrelace)
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_IO_RANGEMENT_BLOC,
  PDM_IO_RANGEMENT_ENTRELACE

} PDM_io_rangement_t;

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

typedef struct _PDM_io_fichier_t PDM_io_fichier_t;

/*=============================================================================
 * Variables globales
 *============================================================================*/


/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/

/**
 * \brief Return the file name (or NULL if no file)
 *
 * \param [in]  fichier   Pointer to \ref PDM_io_fichier_t object
 *
 * \return   Name of the file
 *
 */

const char* PDM_io_get_nom_fichier
(
 PDM_io_fichier_t *fichier
 );


/**
 * \brief Ouverture d'un fichier pour acces parallele
 *
 * \param [in]  nom             Nom du fichier
 * \param [in]  fmt             Fichier text ou binaire
 * \param [in]  suff_t          Type de suffixe (manuel ou automatique)
 * \param [in]  suff_u          Suffixe (si suffixe manuel)
 * \param [in]  s_backup        Active le backup d'un fichier preexistant en mode ecriture
 * \param [in]  accesio         Type (parallele avec mpiio, parallele sans mpiio, sequentiel)
 * \param [in]  mode            Mode d'acces (lecture, ecriture, lecture/ecriture)
 * \param [in]  pdm_mpi_comm    Communicateur lie au fichier
 * \param [out] unite           Unite du fichier
 * \param [out] ierr            Indique si le fichier est de type PDM_io ou non (uniquement pour une ouverture en lecture)
 *
 */
// void PROCF (pdm_io_open_cf, PDM_IO_OPEN_CF)
// (const char            *nom,
//  const PDM_l_num_t  *l_nom,
//  const int             *fmt,
//  const int             *suff_t,
//  const char            *suff_u,
//  const PDM_l_num_t  *l_suff_u,
//  const int             *s_backup,
//  const int             *acces,
//  const int             *mode,
//  const int             *endian,
//  PDM_MPI_Fint              *comm,
//  double                *prop_noeuds_actifs,
//  PDM_l_num_t        *unite,
//  PDM_l_num_t        *ierr
//  ARGF_SUPP_CHAINE
//  );

void PDM_io_open
(
 const char             *nom,
 const PDM_io_fmt_t      fmt,
 const PDM_io_suff_t     suff_t,
 const char             *suff_u,
 const PDM_io_backup_t   s_backup,
 const PDM_io_acces_t    acces,
 const PDM_io_mode_t     mode,
 const PDM_io_endian_t   endian,
 PDM_MPI_Comm            comm,
 double                  prop_noeuds_actifs,
 PDM_io_fichier_t      **unite,
 PDM_l_num_t            *ierr
);


/**
 * \brief Set the file position indicator
 *
 * \param [in] fichier         Pointer to \ref PDM_io_fichier_t object
 * \param [in] offset          Adress
 * \param [in] seek            Origin type
 *
 */
// void PROCF (pdm_io_seek, PDM_IO_SEEK)
// (
// const PDM_l_num_t   *unite,
// const PDM_g_num_t   *offset,
// const PDM_io_seek_t *seek
// );

void PDM_io_seek
(
 PDM_io_fichier_t    *fichier,
 const PDM_g_num_t    offset,
 const PDM_io_seek_t  seek
);


/**
 * \brief Return the current file position
 *
 * \param [in] fichier         Pointer to \ref PDM_io_fichier_t object
 * \param [in] offset          Adress
 *
 */
// void PROCF (pdm_io_tell, PDM_IO_TELL)
// (
// const PDM_l_num_t    *unite,
//       PDM_g_num_t    *offset
// );

PDM_g_num_t
PDM_io_tell
(
 PDM_io_fichier_t   *fichier
);


/**
 * \brief Lecture globale : Le processus maitre accede seul au fichier et redistribue
 * l'information a l'ensemble des processus du communicateur
 *
 * \param [in]  fichier         Pointer to \ref PDM_io_fichier_t object
 * \param [in]  taille_donnee   Taille unitaire de la donnee
 * \param [in]  n_donnees       Nombre de donnees a lire
 * \param [out] donnees         Donnees lues
 *
 */
// void PROCF (pdm_io_lecture_globale, PDM_IO_LECTURE_GLOBALE)
// (const PDM_l_num_t *unite,
//  const PDM_l_num_t *taille_donnee,
//  const PDM_g_num_t *n_donnees,
//  void                 *donnees
// );

void PDM_io_lecture_globale
(
 PDM_io_fichier_t  *fichier,
 const PDM_l_num_t  taille_donnee,
 const PDM_g_num_t  n_donnees,
 void              *donnees
 );


/**
 * \brief Ecriture globale : Le processus maitre accede seul au fichier
 *
 * \param [in]  fichier         Pointer to \ref PDM_io_fichier_t object
 * \param [in]  taille_donnee   Taille unitaire de la donnee
 * \param [in]  n_donnees       Nombre de donnees a ecrire
 * \param [in]  donnees         Donnees ecrites
 *
 */
// void PROCF (pdm_io_ecriture_globale, PDM_IO_ECRITURE_GLOBALE)
// (const PDM_l_num_t *unite,
//  const PDM_l_num_t *taille_donnee,
//  const PDM_g_num_t *n_donnees,
//  const void           *donnees
// );

void PDM_io_ecriture_globale
(
 PDM_io_fichier_t  *fichier,
 const PDM_l_num_t  taille_donnee,
 const PDM_g_num_t  n_donnees,
 const void        *donnees
);


/**
 * \brief Lecture parallele de blocs de donnees suivie d'une redistribution des
 * des donnees suivant l'indirection
 *
 * \param [in]  fichier          Pointer to \ref PDM_io_fichier_t object
 * \param [in]  t_n_composantes  Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
 * \param [in]  n_composantes    Nombre de composantes pour chaque donnee
 * \param [in]  taille_donnee    Taille unitaire de la donnee
 * \param [in]  n_donnees        Nombre de donnees a lire
 * \param [in]  indirection      Indirection de redistribition des donnees
 * \param [out] donnees          Donnees lues
 *
 */
// void PROCF (pdm_io_lec_par_entrelacee, PDM_IO_LEC_PAR_ENTRELACEE)
// (const PDM_l_num_t  *unite,
//  const int             *t_n_composantes,
//  const PDM_l_num_t  *n_composantes,
//  const PDM_l_num_t  *taille_donnee,
//  const PDM_l_num_t  *n_donnees,
//  const PDM_g_num_t *indirection,
//  void                  *donnees
//  );

void PDM_io_lec_par_entrelacee
(PDM_io_fichier_t   *fichier,
 const PDM_io_n_composantes_t t_n_composantes,
 const PDM_l_num_t          *n_composantes,
 const PDM_l_num_t           taille_donnee,
 const PDM_l_num_t           n_donnees,
 const PDM_g_num_t         *indirection,
 void                          *donnees
 );


/**
 * \brief Lecture parallele de blocs de donnees
 * Les blocs doivent etre rangÃ©s par ordre croissant suivant la numÃ©rotation
 * des processus
 *
 * \param [in]  fichier          Pointer to \ref PDM_io_fichier_t object
 * \param [in]  t_n_composantes  Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
 * \param [in]  n_composantes    Nombre de composantes pour chaque donnee
 * \param [in]  taille_donnee    Taille unitaire de la donnee
 * \param [in]  debut_bloc       Adresse relative du debut de bloc
 * \param [in]  n_donnees        Nombre de donnees a lire
 * \param [out] donnees          Donnees lues
 *
 */
// void PROCF (pdm_io_lec_par_bloc, PDM_IO_LEC_PAR_BLOC)
// (const PDM_l_num_t  *unite,
//  const int             *t_n_composantes,
//  const PDM_l_num_t  *n_composantes,
//  const PDM_l_num_t  *taille_donnee,
//  const PDM_l_num_t  *n_donnees,
//  const PDM_g_num_t *debut_bloc,
//  void                  *donnees
// );

void PDM_io_lec_par_bloc
(
 PDM_io_fichier_t             *fichier,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t             debut_bloc,
 void                         *donnees
);


/**
 * \brief Tri des donnees suivant l'indirection puis ecriture parallele des blocs de
 * donnees
 *
 * \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [in] t_n_composantes   Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
 * \param [in] n_composantes     Nombre de composantes pour chaque donnee
 * \param [in] taille_donnee     Taille unitaire de la donnee
 * \param [in] n_donnees         Nombre de donnees a ecrire
 * \param [in] indirection       Indirection de redistribition des donnees
 * \param [in] donnees           Donnees a ecrire
 *
 */
// void PROCF (pdm_io_ecr_par_entrelacee, PDM_IO_ECR_PAR_ENTRELACEE)
// (const PDM_l_num_t  *unite,
//  const int             *t_n_composantes,
//  const PDM_l_num_t  *n_composantes,
//  const PDM_l_num_t  *taille_donnee,
//  const PDM_l_num_t  *n_donnees,
//  const PDM_g_num_t *indirection,
//  const void            *donnees
// );

void PDM_io_ecr_par_entrelacee
(
 PDM_io_fichier_t             *fichier,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t            *indirection,
 const void                   *donnees
);


/**
 * \brief Ecriture parallele de blocs de donnees
 * Les blocs doivent etre rangés par ordre croissant suivant la numérotation
 * des processus
 *
 * \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [in] t_n_composantes   Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
 * \param [in] n_composantes     Nombre de composantes pour chaque donnee
 * \param [in] taille_donnee     Taille unitaire de la donnee
 * \param [in] debut_bloc        Adresse relative du debut de bloc
 * \param [in] n_donnees         Nombre de donnees a lire
 * \param [in] donnees           Donnees a ecrire
 *
 */
// void PROCF (pdm_io_ecr_par_bloc, PDM_IO_ECR_PAR_BLOC)
// (const PDM_l_num_t  *unite,
//  const int             *t_n_composantes,
//  const PDM_l_num_t  *n_composantes,
//  const PDM_l_num_t  *taille_donnee,
//  const PDM_l_num_t  *n_donnees,
//  const PDM_g_num_t *debut_bloc,
//  const void            *donnees
// );

void PDM_io_ecr_par_bloc
(
 PDM_io_fichier_t             *fichier,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             taille_donnee,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t             debut_bloc,
 const void                   *donnees
);


/**
 * \brief Fermeture du fichier sans destruction de la structure PDM_io associee a
 * l'unite
 *
 * \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
 *
 */
// void PROCF (pdm_io_close, PDM_IO_CLOSE)
// (const PDM_l_num_t *unite
// );

void PDM_io_close
(
 PDM_io_fichier_t   *fichier
);


/**
 * \brief Destruction de la structure PDM_io associee a l'unite
 *
 * \param [in] fichier           Pointer to \ref PDM_io_fichier_t object
 *
 */
// void PROCF (pdm_io_detruit, PDM_IO_DETRUIT)
// (const PDM_l_num_t *unite);

void PDM_io_detruit
(
 PDM_io_fichier_t   *fichier
);


/**
 * \brief Retourne le temps cumule d'acces aux fichiers
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [out] t_cpu             Temps CPU
 * \param [out] t_elapsed         Temps elapsed
 *
 */
// void PROCF (pdm_io_get_timer_fichier, PDM_IO_GET_TIMER_FICHIER)
// (const PDM_l_num_t *unite,
//  double *t_cpu,
//  double *t_elapsed
// );

void PDM_io_get_timer_fichier
(
 PDM_io_fichier_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Retourne le temps cumule pour le swap des donnees
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [out] t_cpu             Temps CPU
 * \param [out] t_elapsed         Temps elapsed
 *
 */
// void PROCF (pdm_io_get_timer_swap_endian, PDM_IO_GET_TIMER_SWAP_ENDIAN)
// (const PDM_l_num_t *unite,
//  double               *t_cpu,
//  double               *t_elapsed
// );

void PDM_io_get_timer_swap_endian
(
 PDM_io_fichier_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Retourne le temps cumule pour la distribution des donnees
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [out] t_cpu             Temps CPU
 * \param [out] t_elapsed         Temps elapsed
 *
 */
// void PROCF (pdm_io_get_timer_distrib, PDM_IO_GET_TIMER_DISTRIB)
// (const PDM_l_num_t *unite,
//  double               *t_cpu,
//  double               *t_elapsed
// );

void PDM_io_get_timer_distrib
(
 PDM_io_fichier_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Retourne le temps cumule total
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [out] t_cpu             Temps CPU
 * \param [out] t_elapsed         Temps elapsed
 *
 */
// void PROCF (pdm_io_get_timer_total, PDM_IO_GET_TIMER_TOTAL)
// (const PDM_l_num_t *unite,
//  double               *t_cpu,
//  double               *t_elapsed
// );

void PDM_io_get_timer_total
(
 PDM_io_fichier_t *fichier,
 double           *t_cpu,
 double           *t_elapsed
);


/**
 * \brief Affiche les informations sur le fichier
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 *
 */
// void PROCF (pdm_io_dump, PDM_IO_DUMP)
// (const PDM_l_num_t *unite
// );

void PDM_io_dump
(
 PDM_io_fichier_t   *fichier
);


/**
 * \brief Retourne le communicateur du fichier
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [out] pdm_mpi_comm      Communicateur MPI
 *
 */
// void PROCF (pdm_io_get_comm, PDM_IO_GET_COMM)
// (PDM_l_num_t *unite,
//  PDM_MPI_Fint       *pdm_mpi_comm
// );

void PDM_io_get_comm
(
 PDM_io_fichier_t *fichier,
 PDM_MPI_Comm     *pdm_mpi_comm
);


/**
 * \brief Active le swap endian
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 *
 */
// void PROCF (pdm_io_swap_endian_on, PDM_IO_SWAP_ENDIAN_ON)
// (
// PDM_l_num_t *unite
// );

void PDM_io_swap_endian_on
(
 PDM_io_fichier_t   *fichier
);


/**
 * \brief Désactive le swap endian
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 *
 */
// void PROCF (pdm_io_swap_endian_off, PDM_IO_SWAP_ENDIAN_OFF)
// (
// PDM_l_num_t *unite
// );

void PDM_io_swap_endian_off
(
 PDM_io_fichier_t   *fichier
);


/**
 * \brief Swap endian pour conversion little endian <-> big endian
 *
 * \param [in]  taille_donnee   Taille unitaire de la donnee
 * \param [in]  n_donnee        Nombre de donnees
 * \param [in]  donnees         Donnees
 * \param [out] resultats       Resultat
 *
 */
//  void PROCF (pdm_io_swap_endian, PDM_IO_SWAP_ENDIAN)
//  (
//   const int         *taille_donnee,
//   const PDM_g_num_t *n_donnees,
//   const void        *donnees,
//   void              *resultats
// );

void PDM_io_swap_endian
(
 const size_t   taille_donnee,
 const size_t   n_donnees,
 const void    *donnees,
 void          *resultats
 );


/**
 * \brief Définit le format de la donnée indviduelle pour la sortie text
 *
 * \param [in]  fichier           Pointer to \ref PDM_io_fichier_t object
 * \param [in]  n_char_fmt        Nombre de caractères du format
 * \param [in]  data_type         Type de donnees
 * \param [in]  fmt               Format
 *
 */
// void PROCF (pdm_io_fmt_donnee_set_cf, PDM_IO_FMT_DONNEE_SET_CF)
// (
//  const PDM_l_num_t *unite,
//  const PDM_l_num_t *n_char_fmt,
//  const PDM_l_num_t *data_type,
//  const char           *fmt,
//  const PDM_l_num_t *l_fmt
//  ARGF_SUPP_CHAINE
// );

void PDM_io_fmt_donnee_set
(
 PDM_io_fichier_t    *fichier,
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
 * \brief Calcul de la taille totale d'un champ de donnees
 *
 * \param [in]  fichier          Pointer to \ref PDM_io_fichier_t object
 * \param [in]  t_n_composantes  Type de tailles composantes (PDM_IO_N_COMPOSANTE_CONSTANT ou PDM_IO_N_COMPOSANTE_VARIABLE)
 * \param [in]  n_composantes    Nombre de composantes pour chaque donnee
 * \param [in]  n_donnees        Nombre de donnees
 * \param [in]  indirection      Indirection de redistribition des donnees
 * \param [out] donnees          Donnees
 *
 * \return   Taille totale d'un champ de donnees
 *
 */
// void PROCF (pdm_io_n_donnees_get, PDM_IO_N_DONNEES_GET)
// (const PDM_l_num_t  *unite,
//  const int             *t_n_composantes,
//  const PDM_l_num_t  *n_composantes,
//  const PDM_l_num_t  *n_donnees,
//  const PDM_g_num_t *indirection,
//        PDM_g_num_t *t_n_donnees
//  );

PDM_g_num_t
PDM_io_n_donnees_get
(
 PDM_io_fichier_t             *fichier,
 const PDM_io_n_composantes_t  t_n_composantes,
 const PDM_l_num_t            *n_composantes,
 const PDM_l_num_t             n_donnees,
 const PDM_g_num_t            *indirection
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_IO_H__ */
