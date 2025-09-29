!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2022  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_io

  use iso_c_binding
  use pdm

  implicit none

  !
  ! Parametres
  ! ----------

  !
  ! Types de sufixe

  integer (kind = pdm_l_num_s), parameter :: pdm_io_suff_auto        = 0 ! Suffixe detemine automatiquement
  integer (kind = pdm_l_num_s), parameter :: pdm_io_suff_man         = 1 ! Suffixe fourni par l'utilisateur

  !
  ! Endianness

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_BIGENDIAN        = 0 ! Contenu en bigendian
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_LITTLEENDIAN     = 1 ! Contenu en little endian
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_NATIVE           = 3 ! Contenu natif machine

  !
  ! Types de données

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_INT    = 0 ! Type de donnee int
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_LONG   = 1 ! Type de donnee long
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_DOUBLE = 2 ! Type de donnee double
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_FLOAT  = 3 ! Type de donnee float
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_T_CHAR   = 4 ! Type de donnee char

  !
  ! Types d'entrees/sorties paralleles

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_KIND_MPIIO_EO   = 0 ! Acces parallele avec MPIIO (explicit offset)
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_KIND_MPIIO_IP   = 1 ! Acces parallele avec MPIIO (individual pointer)
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_KIND_MPI_SIMPLE = 2 ! Acces parallele sans MPI
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_KIND_SEQ        = 3 ! Acces 1 fichier par processus

  !
  ! Mode d'acces lecture, ecriture, lecture/ecriture

  integer (kind = pdm_l_num_s), parameter :: PDM_IO_MOD_READ  = 0    ! Acces en lecture
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_MOD_WRITE = 1    ! Acces en ecriture
  integer (kind = pdm_l_num_s), parameter :: PDM_IO_MOD_APPEND    = 2    ! Acces en lecture/ecriture

  !
  ! Indique si le fichier contient une entete IOCEDRE

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_ENTETE_ON  = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_ENTETE_OFF = 1

  !
  ! Indique le format du fichier

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_FMT_TXT = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_FMT_BIN = 1

  !
  ! Active ou non le backup d'un fichier preexistant

  integer(kind = pdm_l_num_s), parameter :: PDM_IO_BACKUP_ON  = 0
  integer(kind = pdm_l_num_s), parameter :: PDM_IO_BACKUP_OFF = 1

  interface






  !>
  !! \brief Initialise une phase d'écriture parallèle de tableaux de données associées
  !! aux numéros de variables PDM
  !! Chaque tableau a ses propres caractéristiques :
  !!         - taille de données
  !!         - nombre de donnée
  !!         - indirection (numérotation absolue)
  !!
  !! \param [in] unite              Unite du fichier
  !! \param [in] t_rangement        Type de rangement
  !! \param [in] num_var_cedre_max  Numéro max de variable PDM
  !! \param [in] n_partition_local  Nombre de partitions locales
  !!
  !!

  subroutine PDM_io_array_write_beg (unite,             &
                                   t_rangement,       &
                                   num_var_cedre_max, &
                                   n_partition_local) &
  bind (c, name='PDM_io_array_write_beg')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: unite
    integer(c_int), value :: t_rangement
    integer(c_int), value :: num_var_cedre_max
    integer(c_int), value :: n_partition_local

  end subroutine PDM_io_array_write_beg


  !>
  !! \brief Definition d'une variable en ecriture
  !!
  !! \param [in] num_var_cedre          Numéro de variable PDM
  !! \param [in] num_indirection_cedre  Numéro d'indirection PDM
  !! \param [in] t_n_composantes        Type de tailles composantes (PDM_
  !! \param [in] n_composantes          Nombre de composantes pour chaque
  !! \param [in] taille_donnee          Taille unitaire de la donnnee
  !!

  subroutine PDM_io_array_write_var_def (num_var_cedre,         &
                                     num_indirection_cedre, &
                                     t_n_composantes,       &
                                     n_composantes,         &
                                     taille_donnee)         &
  bind (c, name='PDM_io_array_write_var_def')
    use iso_c_binding
    implicit none

    integer(c_int), value :: num_var_cedre
    integer(c_int), value :: num_indirection_cedre
    integer(c_int), value :: t_n_composantes
    integer(c_int), value :: n_composantes
    integer(c_int), value :: taille_donnee

  end subroutine PDM_io_array_write_var_def


  !>
  !! \brief Finalise une phase d'écriture parallèle de tableaux de données associées
  !! aux numéros de variables PDM. Cette fonction déclenche réellement
  !! les écritures
  !!

  subroutine PDM_io_array_write_end () &
  bind (c, name='PDM_io_array_write_end')
    use iso_c_binding
    implicit none

  end subroutine PDM_io_array_write_end


  !>
  !! \brief Initialise une phase de lecture parallèle de tableaux de données associées
  !! aux numéros de variables PDM
  !! Chaque tableau a ses propres caractéristiques :
  !!         - taille de données
  !!         - nombre de donnée
  !!         - indirection (numérotation absolue)
  !!
  !! \param [in] unite              Unite du fichier
  !! \param [in] t_rangement        Type de rangement
  !! \param [in] num_var_cedre_max  Numéro max de variable PDM
  !! \param [in] n_partition_local  Nombre de partitions locales
  !!

  subroutine PDM_io_array_read_beg (unite,             &
                                   t_rangement,       &
                                   num_var_cedre_max, &
                                   n_partition_local) &
  bind (c, name='PDM_io_array_read_beg')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: unite
    integer(c_int), value :: t_rangement
    integer(c_int), value :: num_var_cedre_max
    integer(c_int), value :: n_partition_local

  end subroutine PDM_io_array_read_beg


  !>
  !! \brief Definition d'une variable en ecriture
  !!
  !! \param [in] num_var_cedre          Numéro de variable PDM
  !! \param [in] num_indirection_cedre  Numéro d'indirection PDM
  !! \param [in] t_n_composantes        Type de tailles composantes (PDM_STRIDE_CST_INTERLACED ou PDM_STRIDE_VAR_INTERLACED)
  !! \param [in] n_composantes          Nombre de composantes pour chaque donnee
  !! \param [in] taille_donnee          Taille unitaire de la donnnee
  !!

  subroutine PDM_io_array_read_var_def (num_var_cedre,         &
                                     num_indirection_cedre, &
                                     t_n_composantes,       &
                                     n_composantes,         &
                                     taille_donnee)         &
  bind (c, name='PDM_io_array_read_var_def')
    use iso_c_binding
    implicit none

    integer(c_int), value :: num_var_cedre
    integer(c_int), value :: num_indirection_cedre
    integer(c_int), value :: t_n_composantes
    integer(c_int), value :: n_composantes
    integer(c_int), value :: taille_donnee

  end subroutine PDM_io_array_read_var_def


  !>
  !! \brief Finalise une phase de lecture parallèle de tableaux de données associées
  !! aux numéros de variables PDM. Cette fonction déclenche réellement
  !! les écritures
  !!

  subroutine PDM_io_array_read_end () &
  bind (c, name='PDM_io_array_read_end')
    use iso_c_binding
    implicit none

  end subroutine PDM_io_array_read_end

  end interface

contains


subroutine PDM_io_open(nom,                &
                       fmt,                &
                       suff_t,             &
                       suff_u,             &
                       s_backup,           &
                       acces,              &
                       mode,               &
                       endian,             &
                       comm,               &
                       prop_noeuds_actifs, &
                       unite,              &
                       ierr)
  ! Open a file for parallel access
  implicit none

  character (len=*)             :: nom                ! File name
  integer,          intent(in)  :: fmt                ! ASCII or binary file
  integer,          intent(in)  :: suff_t             ! Suffix type (manual/automatic)
  character (len=*)             :: suff_u             ! Suffix (if manual)
  integer,          intent(in)  :: s_backup           ! Enable backup of preexisting file in writing mode
  integer,          intent(in)  :: acces              ! Access type (parallel with/without MPI-IO, serial)
  integer,          intent(in)  :: mode               ! Access mode (read, write, read & write)
  integer,          intent(in)  :: endian             ! Endian type (little or big)
  integer,          intent(in)  :: comm               ! MPI communicator
  real(8),          intent(in)  :: prop_noeuds_actifs ! Proportion of active nodes
  type(c_ptr),      intent(out) :: unite              ! PDM_io_file_t instance
  integer,          intent(out) :: ierr               ! Error code (indicates whether the file is of type PDM_io_file_t or not (for read-only opening only))

  integer(c_int)                :: c_fmt
  integer(c_int)                :: c_suff_t
  integer(c_int)                :: c_s_backup
  integer(c_int)                :: c_acces
  integer(c_int)                :: c_mode
  integer(c_int)                :: c_endian
  type(c_ptr)                   :: c_comm
  real(c_double)                :: c_prop_noeuds_actifs
  integer(c_int)                :: c_ierr

  interface
    subroutine PDM_io_open_c (nom,                &
                              fmt,                &
                              suff_t,             &
                              suff_u,             &
                              s_backup,           &
                              acces,              &
                              mode,               &
                              endian,             &
                              comm,               &
                              prop_noeuds_actifs, &
                              unite,              &
                              ierr)               &
    bind (c, name='PDM_io_open')
      use iso_c_binding
      implicit none

      character(c_char)     :: nom(*)
      integer(c_int), value :: fmt
      integer(c_int), value :: suff_t
      character(c_char)     :: suff_u(*)
      integer(c_int), value :: s_backup
      integer(c_int), value :: acces
      integer(c_int), value :: mode
      integer(c_int), value :: endian
      type(c_ptr), value    :: comm
      real(c_double), value :: prop_noeuds_actifs
      type(c_ptr)           :: unite
      integer(c_int)        :: ierr

    end subroutine PDM_io_open_c
  end interface

  c_fmt      = fmt
  c_suff_t   = suff_t
  c_s_backup = s_backup
  c_acces    = acces
  c_mode     = mode
  c_endian   = endian

  c_comm = PDM_MPI_Comm_f2c(comm)

  c_prop_noeuds_actifs = prop_noeuds_actifs

  call PDM_io_open_c(trim(nom)//C_NULL_CHAR,    &
                     c_fmt,                     &
                     c_suff_t,                  &
                     trim(suff_u)//C_NULL_CHAR, &
                     c_s_backup,                &
                     c_acces,                   &
                     c_mode,                    &
                     c_endian,                  &
                     c_comm,                    &
                     c_prop_noeuds_actifs,      &
                     unite,                     &
                     c_ierr)

  ierr = c_ierr

end subroutine PDM_io_open


subroutine PDM_io_seek(fichier, &
                        offset,  &
                        seek)
  ! Set the file position indicator
  implicit none

  type(c_ptr),          intent(in) :: fichier ! PDM_io_file_t instance
  integer(pdm_g_num_s), intent(in) :: offset  ! Address
  integer,              intent(in) :: seek    ! Origin type

  interface
    subroutine PDM_io_seek_c(fichier, &
                              offset,  &
                              seek)    &
    bind (c, name='PDM_io_seek')
      use iso_c_binding
      implicit none
      type(c_ptr),     value :: fichier
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: offset
#else
      integer(c_int),  value :: offset
#endif
      integer(c_int),  value :: seek
    end subroutine PDM_io_seek_c
  end interface

  call PDM_io_seek_c(fichier, &
                      offset,  &
                      seek)

end subroutine PDM_io_seek


subroutine PDM_io_tell(fichier, &
                       offset)
  ! Return the current file position
  implicit none

  type(c_ptr),          intent(in)  :: fichier ! PDM_io_file_t instance
  integer(pdm_g_num_s), intent(out) :: offset  ! Current position in file

  interface
    function PDM_io_tell_c (fichier) &
    result (offset)                  &
    bind (c, name='PDM_io_tell')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
#ifdef PDM_LONG_G_NUM
      integer(c_long)    :: offset
#else
      integer(c_int)     :: offset
#endif
    end function PDM_io_tell_c
  end interface

  offset = PDM_io_tell_c(fichier)

end subroutine PDM_io_tell



subroutine PDM_io_global_read(fichier,       &
                              taille_donnee, &
                              n_donnees,     &
                              donnees)
  ! Global read: the master process alone accesses the file and redistributes the information to all the communicator's processes
  implicit none

  type(c_ptr),          intent(in)  :: fichier       ! PDM_io_file_t instance
  integer,              intent(in)  :: taille_donnee ! Size of a unit piece of data
  integer(pdm_g_num_s), intent(in)  :: n_donnees     ! Amount of data to be read
  type(c_ptr)                       :: donnees       ! Read data

  interface
    subroutine PDM_io_global_read_c(fichier,       &
                                    taille_donnee, &
                                    n_donnees,     &
                                    donnees)       &
    bind (c, name='PDM_io_global_read')
      use iso_c_binding
      implicit none

      type(c_ptr), value          :: fichier
      integer(c_int),  value :: taille_donnee
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: n_donnees
#else
      integer(c_int),  value :: n_donnees
#endif
      type(c_ptr), value        :: donnees

    end subroutine PDM_io_global_read_c
  end interface

  call PDM_io_global_read_c(fichier,       &
                            taille_donnee, &
                            n_donnees,     &
                            donnees)

end subroutine PDM_io_global_read



subroutine PDM_io_par_interlaced_read(fichier,         &
                                      t_n_composantes, &
                                      n_composantes,   &
                                      taille_donnee,   &
                                      n_donnees,       &
                                      indirection,     &
                                      donnees)
  ! Parallel reading of data blocks followed by redistribution of the data according to indirection
  implicit none

  type(c_ptr),          intent(in) :: fichier          ! PDM_io_file_t instance
  integer(pdm_l_num_s), intent(in) :: t_n_composantes  ! Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
  integer(pdm_l_num_s), pointer    :: n_composantes(:) ! Number of components for each piece of data
  integer(pdm_l_num_s), intent(in) :: taille_donnee    ! Unit size of a piece of data
  integer(pdm_l_num_s), intent(in) :: n_donnees        ! Number of data items to be read
  integer(pdm_g_num_s), pointer    :: indirection(:)   ! Indirection of data redistribution
  type(c_ptr)                      :: donnees          ! Read data

  integer(c_int)                   :: c_t_n_composantes
  integer(c_int)                   :: c_taille_donnee
  integer(c_int)                   :: c_n_donnees
  type(c_ptr)                      :: c_n_composantes
  type(c_ptr)                      :: c_indirection

  interface
    subroutine PDM_io_par_interlaced_read_c(fichier,         &
                                            t_n_composantes, &
                                            n_composantes,   &
                                            taille_donnee,   &
                                            n_donnees,       &
                                            indirection,     &
                                            donnees)         &
    bind (c, name='PDM_io_par_interlaced_read')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: fichier
      integer(c_int), value :: t_n_composantes
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: taille_donnee
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
      type(c_ptr),    value :: donnees

    end subroutine PDM_io_par_interlaced_read_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  c_indirection = C_NULL_PTR
  if (associated(indirection)) then
    c_indirection   = c_loc(indirection  )
  endif


  call PDM_io_par_interlaced_read_c(fichier,           &
                                    c_t_n_composantes, &
                                    c_n_composantes,   &
                                    c_taille_donnee,   &
                                    c_n_donnees,       &
                                    c_indirection,     &
                                    donnees)

end subroutine PDM_io_par_interlaced_read



subroutine PDM_io_par_block_read(fichier,         &
                                 t_n_composantes, &
                                 n_composantes,   &
                                 taille_donnee,   &
                                 n_donnees,       &
                                 debut_bloc,      &
                                 donnees)
  ! Parallel reading of data blocks. The blocks must be arranged in ascending order according to the numbering of the processes
  implicit none

  type(c_ptr),          intent(in) :: fichier          ! PDM_io_file_t instance
  integer(pdm_l_num_s), intent(in) :: t_n_composantes  ! Component size type (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
  integer(pdm_l_num_s), pointer    :: n_composantes(:) ! Number of components for each data item
  integer(pdm_l_num_s), intent(in) :: taille_donnee    ! Unit size of a piece of data
  integer(pdm_l_num_s), intent(in) :: n_donnees        ! Number of data items to be read
  integer(pdm_g_num_s), intent(in) :: debut_bloc       ! Relative address of start of block
  type(c_ptr)                      :: donnees          ! Read data

  integer(c_int)                   :: c_t_n_composantes
  integer(c_int)                   :: c_taille_donnee
  integer(c_int)                   :: c_n_donnees
  type(c_ptr)                      :: c_n_composantes
#ifdef PDM_LONG_G_NUM
  integer(c_long)                  :: c_debut_bloc
#else
  integer(c_int)                   :: c_debut_bloc
#endif

  interface
    subroutine PDM_io_par_block_read_c(fichier,         &
                                       t_n_composantes, &
                                       n_composantes,   &
                                       taille_donnee,   &
                                       n_donnees,       &
                                       debut_bloc,      &
                                       donnees)         &
    bind (c, name='PDM_io_par_block_read')
      use iso_c_binding
      implicit none

      type(c_ptr),     value :: fichier
      integer(c_int),  value :: t_n_composantes
      type(c_ptr),     value :: n_composantes
      integer(c_int),  value :: taille_donnee
      integer(c_int),  value :: n_donnees
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: debut_bloc
#else
      integer(c_int),  value :: debut_bloc
#endif
      type(c_ptr),     value :: donnees

    end subroutine PDM_io_par_block_read_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees
  c_debut_bloc      = debut_bloc

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  call PDM_io_par_block_read_c(fichier,           &
                               c_t_n_composantes, &
                               c_n_composantes,   &
                               c_taille_donnee,   &
                               c_n_donnees,       &
                               c_debut_bloc,      &
                               donnees)

end subroutine PDM_io_par_block_read


subroutine PDM_io_global_write(fichier,       &
                               taille_donnee, &
                               n_donnees,     &
                               donnees)

  ! Global write: The master process has sole access to the file
  implicit none

  type(c_ptr),          intent(in) :: fichier       ! PDM_io_file_t instance
  integer,              intent(in) :: taille_donnee ! Size of a unit piece of data
  integer(pdm_g_num_s), intent(in) :: n_donnees     ! Amount of data to write
  type(c_ptr),          intent(in) :: donnees       ! Data to write

  interface
    subroutine PDM_io_global_write_c(fichier,       &
                                     taille_donnee, &
                                     n_donnees,     &
                                     donnees)       &
    bind (c, name='PDM_io_global_write')
      use iso_c_binding
      implicit none
      type(c_ptr),     value :: fichier
      integer(c_int),  value :: taille_donnee
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: n_donnees
#else
      integer(c_int),  value :: n_donnees
#endif
      type(c_ptr),     value :: donnees
    end subroutine PDM_io_global_write_c
  end interface

  call PDM_io_global_write_c(fichier,       &
                             taille_donnee, &
                             n_donnees,     &
                             donnees)

end subroutine PDM_io_global_write



subroutine PDM_io_par_interlaced_write(fichier,         &
                                       t_n_composantes, &
                                       n_composantes,   &
                                       taille_donnee,   &
                                       n_donnees,       &
                                       indirection,     &
                                       donnees)
  ! Data sorted according to indirection, then parallel write of data blocks
  implicit none

  type(c_ptr),          intent(in) :: fichier          ! PDM_io_file_t instance
  integer(pdm_l_num_s), intent(in) :: t_n_composantes  ! Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
  integer(pdm_l_num_s), pointer    :: n_composantes(:) ! Number of components for each data item
  integer(pdm_l_num_s), intent(in) :: taille_donnee    ! Unit size of the data
  integer(pdm_l_num_s), intent(in) :: n_donnees        ! Number of data items to be written
  integer(pdm_g_num_s), pointer    :: indirection(:)   ! Data redistribution direction
  type(c_ptr),          intent(in) :: donnees          ! Data to be written

  integer(c_int)                   :: c_t_n_composantes
  integer(c_int)                   :: c_taille_donnee
  integer(c_int)                   :: c_n_donnees
  type(c_ptr)                      :: c_n_composantes
  type(c_ptr)                      :: c_indirection

  interface
    subroutine PDM_io_par_interlaced_write_c (fichier,         &
                                            t_n_composantes, &
                                            n_composantes,   &
                                            taille_donnee,   &
                                            n_donnees,       &
                                            indirection,     &
                                            donnees)         &
    bind (c, name='PDM_io_par_interlaced_write')
      use iso_c_binding
      implicit none
      type(c_ptr),    value :: fichier
      integer(c_int), value :: t_n_composantes
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: taille_donnee
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
      type(c_ptr),    value :: donnees
    end subroutine PDM_io_par_interlaced_write_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  c_indirection = C_NULL_PTR
  if (associated(indirection)) then
    c_indirection   = c_loc(indirection  )
  endif


  call PDM_io_par_interlaced_write_c(fichier,           &
                                     c_t_n_composantes, &
                                     c_n_composantes,   &
                                     c_taille_donnee,   &
                                     c_n_donnees,       &
                                     c_indirection,     &
                                     donnees)

end subroutine PDM_io_par_interlaced_write



subroutine PDM_io_par_block_write(fichier,         &
                                  t_n_composantes, &
                                  n_composantes,   &
                                  taille_donnee,   &
                                  n_donnees,       &
                                  debut_bloc,      &
                                  donnees)
  ! Parallel writing of data blocks. Blocks must be arranged in ascending order according to numbering of the processes
  implicit none

  type(c_ptr),          intent(in) :: fichier          ! PDM_io_file_t instance
  integer(pdm_l_num_s), intent(in) :: t_n_composantes  ! Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
  integer(pdm_l_num_s), pointer    :: n_composantes(:) ! Number of components for each data item
  integer(pdm_l_num_s), intent(in) :: taille_donnee    ! Unit size of the data
  integer(pdm_l_num_s), intent(in) :: n_donnees        ! Number of data to read
  integer(pdm_g_num_s), intent(in) :: debut_bloc       ! Relative address of start of block
  type(c_ptr),          intent(in) :: donnees          ! Data to be written

  integer(c_int)                   :: c_t_n_composantes
  integer(c_int)                   :: c_taille_donnee
  integer(c_int)                   :: c_n_donnees
  type(c_ptr)                      :: c_n_composantes
#ifdef PDM_LONG_G_NUM
  integer(c_long)                  :: c_debut_bloc
#else
  integer(c_int)                   :: c_debut_bloc
#endif

  interface
    subroutine PDM_io_par_block_write_c(fichier,         &
                                        t_n_composantes, &
                                        n_composantes,   &
                                        taille_donnee,   &
                                        n_donnees,       &
                                        debut_bloc,      &
                                        donnees)         &
    bind (c, name='PDM_io_par_block_write')
      use iso_c_binding
      implicit none
      type(c_ptr),     value :: fichier
      integer(c_int),  value :: t_n_composantes
      type(c_ptr),     value :: n_composantes
      integer(c_int),  value :: taille_donnee
      integer(c_int),  value :: n_donnees
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: debut_bloc
#else
      integer(c_int),  value :: debut_bloc
#endif
      type(c_ptr),     value :: donnees
    end subroutine PDM_io_par_block_write_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_taille_donnee   = taille_donnee
  c_n_donnees       = n_donnees
  c_debut_bloc      = debut_bloc

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  call PDM_io_par_block_write_c(fichier,           &
                                c_t_n_composantes, &
                                c_n_composantes,   &
                                c_taille_donnee,   &
                                c_n_donnees,       &
                                c_debut_bloc,      &
                                donnees)

end subroutine PDM_io_par_block_write



subroutine PDM_io_get_timer_fichier(fichier,   &
                                    t_cpu,     &
                                    t_elapsed)
  ! Returns the cumulative files access time
  implicit none

  type(c_ptr), intent(in)  :: fichier   ! PDM_io_file_t instance
  real(8),     intent(out) :: t_cpu     ! CPU time
  real(8),     intent(out) :: t_elapsed ! Elapsed time

  real(c_double)           :: c_t_cpu
  real(c_double)           :: c_t_elapsed

  interface
    subroutine PDM_io_get_timer_fichier_c(fichier,   &
                                          t_cpu,     &
                                          t_elapsed) &
    bind (c, name='PDM_io_get_timer_fichier')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed
    end subroutine PDM_io_get_timer_fichier_c
  end interface

  call PDM_io_get_timer_fichier_c(fichier,     &
                                  c_t_cpu,     &
                                  c_t_elapsed)

  t_cpu     = c_t_cpu
  t_elapsed = c_t_elapsed

end subroutine PDM_io_get_timer_fichier



subroutine PDM_io_timer_swap_endian_get(fichier,   &
                                        t_cpu,     &
                                        t_elapsed)
  ! Returns the cumulative time for data swap
  implicit none

  type(c_ptr), intent(in)  :: fichier   ! PDM_io_file_t instance
  real(8),     intent(out) :: t_cpu     ! CPU time
  real(8),     intent(out) :: t_elapsed ! Elapsed time

  real(c_double)           :: c_t_cpu
  real(c_double)           :: c_t_elapsed

  interface
    subroutine PDM_io_timer_swap_endian_get_c(fichier,   &
                                              t_cpu,     &
                                              t_elapsed) &
    bind (c, name='PDM_io_timer_swap_endian_get')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed
    end subroutine PDM_io_timer_swap_endian_get_c
  end interface

  call PDM_io_timer_swap_endian_get_c(fichier,     &
                                      c_t_cpu,     &
                                      c_t_elapsed)

  t_cpu     = c_t_cpu
  t_elapsed = c_t_elapsed

end subroutine PDM_io_timer_swap_endian_get



subroutine PDM_io_timer_distrib_get(fichier,   &
                                    t_cpu,     &
                                    t_elapsed)
  ! Returns the cumulative time for data distribution
  implicit none

  type(c_ptr), intent(in)  :: fichier   ! PDM_io_file_t instance
  real(8),     intent(out) :: t_cpu     ! CPU time
  real(8),     intent(out) :: t_elapsed ! Elapsed time

  real(c_double)           :: c_t_cpu
  real(c_double)           :: c_t_elapsed

  interface
    subroutine PDM_io_timer_distrib_get_c(fichier,   &
                                          t_cpu,     &
                                          t_elapsed) &
    bind (c, name='PDM_io_timer_distrib_get')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed
    end subroutine PDM_io_timer_distrib_get_c
  end interface

  call PDM_io_timer_distrib_get_c(fichier,     &
                                  c_t_cpu,     &
                                  c_t_elapsed)

  t_cpu     = c_t_cpu
  t_elapsed = c_t_elapsed

end subroutine PDM_io_timer_distrib_get



subroutine PDM_io_timer_total_get(fichier,   &
                                  t_cpu,     &
                                  t_elapsed)
  ! Returns the total cumulative time
  implicit none

  type(c_ptr), intent(in)  :: fichier   ! PDM_io_file_t instance
  real(8),     intent(out) :: t_cpu     ! CPU time
  real(8),     intent(out) :: t_elapsed ! Elapsed time

  real(c_double)       :: c_t_cpu
  real(c_double)       :: c_t_elapsed

  interface
    subroutine PDM_io_timer_total_get_c(fichier,   &
                                        t_cpu,     &
                                        t_elapsed) &
    bind (c, name='PDM_io_timer_total_get')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
      real(c_double)     :: t_cpu
      real(c_double)     :: t_elapsed
    end subroutine PDM_io_timer_total_get_c
  end interface

  call PDM_io_timer_total_get_c(fichier,     &
                                c_t_cpu,     &
                                c_t_elapsed)

  t_cpu     = c_t_cpu
  t_elapsed = c_t_elapsed

end subroutine PDM_io_timer_total_get



subroutine PDM_io_swap_endian(taille_donnee, &
                              n_donnees,     &
                              donnees,       &
                              resultats)
  ! Swap endian pour conversion little endian <-> big endian
  implicit none

  integer, intent(in) :: taille_donnee ! Size of a unit piece of data
  integer, intent(in) :: n_donnees     ! Amount of data
  type(c_ptr)         :: donnees       ! Data
  type(c_ptr)         :: resultats     ! Result

  integer(c_size_t)   :: c_taille_donnee
  integer(c_size_t)   :: c_n_donnees

  interface
    subroutine PDM_io_swap_endian_c(taille_donnee, &
                                    n_donnees,     &
                                    donnees,       &
                                    resultats)     &
    bind (c, name='PDM_io_swap_endian')
      use iso_c_binding
      implicit none

      integer(c_size_t), value :: taille_donnee
      integer(c_size_t), value :: n_donnees
      type(c_ptr),       value :: donnees
      type(c_ptr)              :: resultats

    end subroutine PDM_io_swap_endian_c
  end interface

  c_taille_donnee = taille_donnee
  c_n_donnees     = n_donnees

  call PDM_io_swap_endian_c(c_taille_donnee, &
                            c_n_donnees,     &
                            donnees,         &
                            resultats)

end subroutine PDM_io_swap_endian



subroutine PDM_io_swap_endian_on(fichier)
  ! Activate endian swap
  implicit none

  type(c_ptr), intent(in) :: fichier ! PDM_io_file_t instance

  interface
    subroutine PDM_io_swap_endian_on_c(fichier) &
    bind (c, name='PDM_io_swap_endian_on')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
    end subroutine PDM_io_swap_endian_on_c
  end interface

  call PDM_io_swap_endian_on_c(fichier)

end subroutine PDM_io_swap_endian_on



subroutine PDM_io_swap_endian_off(fichier)
  ! Deactivate endian swap
  implicit none

  type(c_ptr), intent(in) :: fichier ! PDM_io_file_t instance

  interface
    subroutine PDM_io_swap_endian_off_c(fichier) &
    bind (c, name='PDM_io_swap_endian_off')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
    end subroutine PDM_io_swap_endian_off_c
  end interface

  call PDM_io_swap_endian_off_c(fichier)

end subroutine PDM_io_swap_endian_off



subroutine PDM_io_fmt_data_set(fichier,    &
                               n_char_fmt, &
                               data_type,  &
                               fmt)
  ! Defines the format of the individual data for text output
  implicit none

  type(c_ptr), intent(in) :: fichier    ! PDM_io_file_t instance
  integer,     intent(in) :: n_char_fmt ! Number of characters in the format
  integer,     intent(in) :: data_type  ! Type of data
  character(len=*)        :: fmt        ! Format

  integer(c_int)          :: c_n_char_fmt
  integer(c_int)          :: c_data_type

  interface
    subroutine PDM_io_fmt_data_set_c(fichier,    &
                                     n_char_fmt, &
                                     data_type,  &
                                     fmt)        &
    bind (c, name='PDM_io_fmt_data_set')
      use iso_c_binding
      implicit none
      type(c_ptr),    value :: fichier
      integer(c_int), value :: n_char_fmt
      integer(c_int), value :: data_type
      character(c_char)     :: fmt(*)
    end subroutine PDM_io_fmt_data_set_c
  end interface

  c_n_char_fmt = n_char_fmt
  c_data_type  = data_type

  call PDM_io_fmt_data_set_c(fichier,      &
                             c_n_char_fmt, &
                             c_data_type,  &
                             trim(fmt)//C_NULL_CHAR)

end subroutine PDM_io_fmt_data_set



subroutine PDM_io_mkdir(path, &
                        code)
  ! Create a directory
  implicit none

  character(len=*), intent(in)  :: path ! Path to new directory
  integer,          intent(out) :: code ! 0 if successful, -1 else

  interface
    function PDM_io_mkdir_c(path) &
    result (code)                 &
    bind (c, name='PDM_io_mkdir')
      use iso_c_binding
      implicit none
      character(c_char) :: path(*)
      integer(c_int)    :: code
    end function PDM_io_mkdir_c
  end interface

  code = PDM_io_mkdir_c(trim(path)//C_NULL_CHAR)

end subroutine PDM_io_mkdir



subroutine PDM_io_n_data_get(fichier,         &
                             t_n_composantes, &
                             n_composantes,   &
                             n_donnees,       &
                             indirection,     &
                             taille)
  ! Calculate the total size of a data field
  implicit none

  type(c_ptr),          intent(in)  :: fichier          ! PDM_io_file_t instance
  integer,              intent(in)  :: t_n_composantes  ! Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
  integer(pdm_l_num_s), pointer     :: n_composantes(:) ! Number of components for each data
  integer(pdm_l_num_s), intent(in)  :: n_donnees        ! Number of data
  integer(pdm_g_num_s), pointer     :: indirection(:)   ! Data redistribution direction
  integer(pdm_g_num_s), intent(out) :: taille           ! Total size of a data field

  integer(c_int)                    :: c_t_n_composantes
  type(c_ptr)                       :: c_n_composantes
  integer(c_int)                    :: c_n_donnees
  type(c_ptr)                       :: c_indirection
#ifdef PDM_LONG_G_NUM
  integer(c_long)                   :: c_taille
#else
  integer(c_int)                    :: c_taille
#endif

  interface
    function PDM_io_n_data_get_c(fichier,         &
                                 t_n_composantes, &
                                 n_composantes,   &
                                 n_donnees,       &
                                 indirection)     &
    result (taille)                               &
    bind (c, name='PDM_io_n_data_get')
      use iso_c_binding
      implicit none
      type(c_ptr),    value :: fichier
      integer(c_int), value :: t_n_composantes
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
#ifdef PDM_LONG_G_NUM
      integer(c_long)       :: taille
#else
      integer(c_int)        :: taille
#endif
    end function PDM_io_n_data_get_c
  end interface

  c_t_n_composantes = t_n_composantes
  c_n_donnees       = n_donnees

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  c_indirection = C_NULL_PTR
  if (associated(indirection)) then
    c_indirection   = c_loc(indirection)
  endif


  c_taille = PDM_io_n_data_get_c(fichier,           &
                                 c_t_n_composantes, &
                                 c_n_composantes,   &
                                 c_n_donnees,       &
                                 c_indirection)

  taille = c_taille

end subroutine PDM_io_n_data_get


!>
!! \brief Ajoute une partie des donnees dans un tableau associés à une variable
!! PDM
!!
!! \param [in] num_var_cedre          Numéro de variable PDM
!! \param [in] i_part                 indice de partition
!! \param [in] n_composantes          Nombre de composantes pour chaque donnee
!! \param [in] n_donnees              Nombre de donnees a lire
!! \param [in] indirection            Indirection de redistribition des donnees
!! \param [in] donnees                Donnees a écrire
!!
!!

subroutine PDM_io_array_write_data_append (num_var_cedre, &
                                         i_part,        &
                                         n_composantes, &
                                         n_donnees,     &
                                         indirection,   &
                                         donnees)
  use iso_c_binding
  implicit none

  integer(pdm_l_num_s), intent(in) :: num_var_cedre
  integer(pdm_l_num_s), intent(in) :: i_part
  integer(pdm_l_num_s), pointer    :: n_composantes(:)
  integer(pdm_l_num_s), intent(in) :: n_donnees
  integer(pdm_g_num_s), pointer    :: indirection(:)
  type(c_ptr), value               :: donnees

  type(c_ptr)                      :: c_n_composantes
  type(c_ptr)                      :: c_indirection

  interface
    subroutine PDM_io_array_write_data_append_c (num_var_cedre, &
                                               i_part,        &
                                               n_composantes, &
                                               n_donnees,     &
                                               indirection,   &
                                               donnees)       &
    bind (c, name='PDM_io_array_write_data_append')
      use iso_c_binding
      implicit none

      integer(c_int), value :: num_var_cedre
      integer(c_int), value :: i_part
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
      type(c_ptr),    value :: donnees

    end subroutine PDM_io_array_write_data_append_c
  end interface

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  c_indirection = C_NULL_PTR
  if (associated(indirection)) then
    c_indirection   = c_loc(indirection  )
  endif


  call PDM_io_array_write_data_append_c (num_var_cedre,   &
                                       i_part,          &
                                       c_n_composantes, &
                                       n_donnees,       &
                                       c_indirection,   &
                                       donnees)

end subroutine PDM_io_array_write_data_append


!>
!! \brief Ajoute une partie des donnees dans un tableau associés à une variable PDM
!!
!! \param [in] num_var_cedre          Numéro de variable PDM
!! \param [in] i_part                 indice de partition
!! \param [in] n_composantes          Nombre de composantes pour chaque donnee
!! \param [in] n_donnees              Nombre de donnees a lire
!! \param [in] indirection            Indirection de redistribition des donnees
!! \param [in] donnees                Donnees a écrire
!!
!!

subroutine PDM_io_array_read_data_append (num_var_cedre, &
                                         i_part,        &
                                         n_composantes, &
                                         n_donnees,     &
                                         indirection,   &
                                         donnees)
  use iso_c_binding
  implicit none

  integer(pdm_l_num_s), intent(in) :: num_var_cedre
  integer(pdm_l_num_s), intent(in) :: i_part
  integer(pdm_l_num_s), pointer    :: n_composantes(:)
  integer(pdm_l_num_s), intent(in) :: n_donnees
  integer(pdm_g_num_s), pointer    :: indirection(:)
  type(c_ptr), value               :: donnees

  type(c_ptr)                      :: c_n_composantes
  type(c_ptr)                      :: c_indirection

  interface
    subroutine PDM_io_array_read_data_append_c (num_var_cedre, &
                                               i_part,        &
                                               n_composantes, &
                                               n_donnees,     &
                                               indirection,   &
                                               donnees)       &
    bind (c, name='PDM_io_array_read_data_append')
      use iso_c_binding
      implicit none

      integer(c_int), value :: num_var_cedre
      integer(c_int), value :: i_part
      type(c_ptr),    value :: n_composantes
      integer(c_int), value :: n_donnees
      type(c_ptr),    value :: indirection
      type(c_ptr),    value :: donnees

    end subroutine PDM_io_array_read_data_append_c
  end interface

  c_n_composantes = C_NULL_PTR
  if (associated(n_composantes)) then
    c_n_composantes = c_loc(n_composantes)
  endif

  c_indirection = C_NULL_PTR
  if (associated(indirection)) then
    c_indirection   = c_loc(indirection  )
  endif


  call PDM_io_array_read_data_append_c (num_var_cedre,   &
                                       i_part,          &
                                       c_n_composantes, &
                                       n_donnees,       &
                                       c_indirection,   &
                                       donnees)

end subroutine PDM_io_array_read_data_append



subroutine PDM_io_close(fichier)
  ! Close a file without destroying the PDM_io structure associated with unit
  implicit none

  type(c_ptr), intent(in) :: fichier ! PDM_io_file_t instance

  interface
    subroutine PDM_io_close_c(fichier) &
      bind (c, name='PDM_io_close')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
    end subroutine PDM_io_close_c
  end interface

  call PDM_io_close_c(fichier)
end subroutine PDM_io_close



subroutine PDM_io_free(fichier)
  ! Free of the PDM_io structure associated with the unit
  implicit none

  type(c_ptr), intent(in) :: fichier ! PDM_io_file_t instance

  interface
    subroutine PDM_io_free_c(fichier) &
      bind (c, name='PDM_io_free')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
    end subroutine PDM_io_free_c
  end interface

  call PDM_io_free_c(fichier)

end subroutine PDM_io_free



subroutine PDM_io_dump(fichier)
  ! Shows file information
  implicit none

  type(c_ptr), intent(in) :: fichier ! PDM_io_file_t instance

  interface
    subroutine PDM_io_dump_c(fichier) &
      bind (c, name='PDM_io_dump')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
    end subroutine PDM_io_dump_c
  end interface

  call PDM_io_dump_c(fichier)

end subroutine PDM_io_dump



subroutine PDM_io_comm_get(fichier, &
                           f_comm)
  ! Returns the file communicator
  implicit none

  type(c_ptr), intent(in)  :: fichier ! PDM_io_file_t instance
  integer,     intent(out) :: f_comm  ! MPI communicator

  interface
    subroutine PDM_io_comm_get_c(fichier, &
                                 f_comm)  &
    bind (c, name='PDM_io_comm_get')
      use iso_c_binding
      implicit none
      type(c_ptr), value :: fichier
      integer(c_int)     :: f_comm
    end subroutine PDM_io_comm_get_c
  end interface

  call PDM_io_comm_get_c(fichier, &
                         f_comm)

end subroutine PDM_io_comm_get


end module pdm_io
