/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_file_seq.h"
#include "pdm.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Type decrivant un fichier de type MSG (MPI_IO)
 *----------------------------------------------------------------------------*/

struct _PDM_file_seq_t {

  FILE               *fichier;     /* Pointeur sur le fichier C */
  char               *nom;     /* Nom du fichier */
  PDM_file_seq_mode_t  mode;    /* Mode */

};

/*============================================================================
 * Variables globales locales
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  PDM_file_par_sek -> pdm_mpi_file_seek_t
 *----------------------------------------------------------------------------*/

static int _std_file_seek[3] = {SEEK_SET,
                                      SEEK_CUR,
                                      SEEK_END};

/*============================================================================
 * Definitions des fonctions locales
 *============================================================================*/

/*============================================================================
 * Definitions des fonctions publiques
 *============================================================================*/

PDM_file_seq_t*
PDM_file_seq_open
(
  const char                *nom,
  const PDM_file_seq_mode_t  mode
)
{

  PDM_file_seq_t *fichier;
  PDM_malloc(fichier, 1, PDM_file_seq_t);

  PDM_malloc(fichier->nom, strlen(nom) + 1, char);
  strcpy(fichier->nom, nom);
  fichier->mode = mode;

  switch (mode) {
  case FICHIER_SEQ_MODE_LECTURE:
    fichier->fichier = fopen(nom, "r");
    break;
  case FICHIER_SEQ_MODE_ECRITURE:
    fichier->fichier = fopen(nom, "w");
    break;
  case FICHIER_SEQ_MODE_AJOUT:
    fichier->fichier = fopen(nom, "a");
    break;
  default:
    PDM_error("Unknow file mode");
  }

  if (fichier->fichier == NULL) {
    PDM_error("Error when open file %s", nom);
  }

  return fichier;

}


PDM_g_num_t
PDM_file_seq_write
(
       PDM_file_seq_t *fichier,
 const size_t          taille_donnee,
 const PDM_g_num_t     n_donnees,
       void           *donnees
)
{

  if (fichier->mode == FICHIER_SEQ_MODE_LECTURE) {
    PDM_error("Write forbidden for file '%s' opened in read mode", fichier->nom);
  }

  size_t _n_donnees = (size_t) n_donnees ;
  size_t _n_donnees_ecrites = fwrite(donnees, taille_donnee, _n_donnees, fichier->fichier);
  PDM_g_num_t n_donnees_ecrites = (PDM_g_num_t) _n_donnees_ecrites;

  return n_donnees_ecrites;

}


PDM_g_num_t
PDM_file_seq_read
(
        PDM_file_seq_t *fichier,
  const size_t          taille_donnee,
  const PDM_g_num_t     n_donnees,
        void           *donnees
)
{

  if (fichier->mode == FICHIER_SEQ_MODE_ECRITURE ||
      fichier->mode == FICHIER_SEQ_MODE_AJOUT) {
    PDM_error("Read forbidden for file '%s' opened in write/append mode", fichier->nom);
  }

  size_t _n_donnees = (size_t) n_donnees ;
  size_t _n_donnees_lues = fread(donnees, taille_donnee,
                                   _n_donnees, fichier->fichier);
  PDM_g_num_t n_donnees_lues = (PDM_g_num_t) _n_donnees_lues;

  return n_donnees_lues;

}


void
PDM_file_seq_seek
(
 PDM_file_seq_t     *fichier,
 long                offset,
 PDM_file_seq_seek_t whence
)
{
  fseek(fichier->fichier, offset, _std_file_seek[whence]);
}


long
PDM_file_seq_tell
(
 PDM_file_seq_t *fichier
)
{
  return ftell(fichier->fichier);
}


void
PDM_file_seq_close
(
 PDM_file_seq_t *fichier
)
{
  PDM_free(fichier->nom);
  fclose(fichier->fichier);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
