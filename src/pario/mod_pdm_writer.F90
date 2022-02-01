module pdm_writer

  use pdm

  implicit none

  !
  ! Statut

  integer, parameter :: PDM_WRITER_OFF  = 0
  integer, parameter :: PDM_WRITER_ON   = 1

  !
  ! Type de topologie

  integer, parameter :: PDM_WRITER_TOPO_CONSTANTE  = 0
  integer, parameter :: PDM_WRITER_TOPO_DEFORMABLE = 1
  integer, parameter :: PDM_WRITER_TOPO_VARIABLE   = 2

  !
  ! Type d'elements géometriques

  integer, parameter :: PDM_WRITER_POINT    = 0
  integer, parameter :: PDM_WRITER_BAR2     = 1
  integer, parameter :: PDM_WRITER_TRIA3    = 2
  integer, parameter :: PDM_WRITER_QUAD4    = 3
  integer, parameter :: PDM_WRITER_POLY_2D  = 4
  integer, parameter :: PDM_WRITER_TETRA4   = 5
  integer, parameter :: PDM_WRITER_PYRAMID5 = 6
  integer, parameter :: PDM_WRITER_PRISM6   = 7
  integer, parameter :: PDM_WRITER_HEXA8    = 8
  integer, parameter :: PDM_WRITER_POLY_3D  = 9

  !
  ! Format de sortie

  integer, parameter ::  PDM_WRITER_FMT_ENSIGHT = 0

  !
  ! Format du fichier

  integer, parameter ::  PDM_WRITER_FMT_BIN   = 0
  integer, parameter ::  PDM_WRITER_FMT_ASCII = 1

  !
  ! Dimension géométrique de la sortie

  integer, parameter ::  PDM_WRITER_DIM_2 = 0
  integer, parameter ::  PDM_WRITER_DIM_3 = 1

  !
  ! Dim des variables

  integer, parameter ::  PDM_WRITER_VAR_CSTE         = 0
  integer, parameter ::  PDM_WRITER_VAR_SCALAIRE     = 1
  integer, parameter ::  PDM_WRITER_VAR_VECTEUR      = 3
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_SYM  = 6
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_ASYM = 9

  !
  ! Localisation des variables

  integer, parameter ::  PDM_WRITER_VAR_SOMMETS      = 0
  integer, parameter ::  PDM_WRITER_VAR_ELEMENTS     = 1
  integer, parameter ::  PDM_WRITER_VAR_PARTICULES   = 2


  contains

  !>
  !!
  !! \brief Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet
  !!
  !! \param [out] cs                   Pointer to \ref PDM_writer object
  !! \param [in]  fmt                  Format de sortie
  !! \param [in]  fmt_fic              Binary or ASCII
  !! \param [in]  topologie            Indique le maillage est mobile ou non
  !! \param [in]  st_reprise           Complete les sorties des calculs precedents en reprise
  !! \param [in]  rep_sortie           Repertoire de sortie
  !! \param [in]  nom_sortie           Nom de la sortie
  !! \param [in]  pdm_mpi_com          Communicateur MSG
  !! \param [in]  acces                Type d'acces
  !! \param [in]  prop_noeuds_actifs   Proportion des noeuds actifs dans les acces au fichier
  !!                                     *  -1 : tous les processus actifs
  !!                                     *   1 : un processus par noeud
  !!                                     * 0 < val < 1 : un processus par noeud actif
  !! \param [in]  options              Options complementaires propres au format sous
  !!                                 la forme ("nom_1 = val_1 : ... : nom_n = val_n")
  !!
  !!

  subroutine PDM_writer_create (cs,                 &
                                fmt,                &
                                fmt_fic,            &
                                topologie,          &
                                st_reprise,         &
                                rep_sortie,         &
                                nom_sortie,         &
                                f_comm,             &
                                acces,              &
                                prop_noeuds_actifs, &
                                options)
    use iso_c_binding
    implicit none

    type(c_ptr)        :: cs
    character(len = *) :: fmt
    integer            :: fmt_fic
    integer            :: topologie
    integer            :: st_reprise
    character(len = *) :: rep_sortie
    character(len = *) :: nom_sortie
    integer            :: f_comm
    integer            :: acces
    double precision   :: prop_noeuds_actifs
    character(len = *) :: options

    integer(c_int)     :: c_fmt_fic
    integer(c_int)     :: c_topologie
    integer(c_int)     :: c_st_reprise
    integer(c_int)     :: c_comm
    integer(c_int)     :: c_acces
    real(c_double)     :: c_prop_noeuds_actifs

    interface
      function PDM_writer_create_c (fmt,                &
                                    fmt_fic,            &
                                    topologie,          &
                                    st_reprise,         &
                                    rep_sortie,         &
                                    nom_sortie,         &
                                    pdm_mpi_comm,       &
                                    acces,              &
                                    prop_noeuds_actifs, &
                                    options)            &
      result (cs)                                       &
      bind (c, name='PDM_writer_create')
        use iso_c_binding
        implicit none

        type(c_ptr)           :: cs
        character(c_char)     :: fmt
        integer(c_int), value :: fmt_fic
        integer(c_int), value :: topologie
        integer(c_int), value :: st_reprise
        character(c_char)     :: rep_sortie
        character(c_char)     :: nom_sortie
        integer(c_int), value :: pdm_mpi_comm
        integer(c_int), value :: acces
        real(c_double), value :: prop_noeuds_actifs
        character(c_char)     :: options

      end function PDM_writer_create_c
    end interface

    c_comm               = PDM_MPI_Comm_f2c(f_comm)
    c_fmt_fic            = fmt_fic
    c_topologie          = topologie
    c_st_reprise         = st_reprise
    c_acces              = acces
    c_prop_noeuds_actifs = prop_noeuds_actifs

    cs = PDM_writer_create_c (fmt//C_NULL_CHAR,        &
                              c_fmt_fic,               &
                              c_topologie,             &
                              c_st_reprise,            &
                              rep_sortie//C_NULL_CHAR, &
                              nom_sortie//C_NULL_CHAR, &
                              c_comm,                  &
                              c_acces,                 &
                              c_prop_noeuds_actifs,    &
                              options//C_NULL_CHAR)

  end subroutine PDM_writer_create



  !>
  !! \brief Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
  !!
  !! \param [in] cs    Pointer to \ref PDM_writer object
  !!
  !!

  subroutine PDM_writer_free (cs) &
  bind (c, name='PDM_writer_free')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cs

  end subroutine PDM_writer_free



  !>
  !! \brief Debut d'increment
  !!
  !! \param [in] cs             Pointer to \ref PDM_writer object
  !! \param [in] physical_time  Temps
  !!

  subroutine PDM_writer_step_beg (cs,            &
                                  physical_time)
    use iso_c_binding
    implicit none

    type(c_ptr), value :: cs
    double precision   :: physical_time

    real(c_double)     :: c_physical_time

    interface
      subroutine PDM_writer_step_beg_c (cs,            &
                                        physical_time) &
      bind (c, name="PDM_writer_step_beg")

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        real(c_double), value :: physical_time

      end subroutine PDM_writer_step_beg_c
    end interface

    c_physical_time = physical_time

    call PDM_writer_step_beg_c (cs,            &
                                c_physical_time)

  end subroutine PDM_writer_step_beg



  !>
  !! \brief Fin d'increment
  !!
  !! \param [in] cs             Pointer to \ref PDM_writer object
  !!
  !!

  subroutine PDM_writer_step_end (cs) &
  bind (c, name="PDM_writer_step_end")

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cs

  end subroutine PDM_writer_step_end



  !>
  !! \brief Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
  !!
  !! \param [in]  cs                Pointer to \ref PDM_writer object
  !! \param [in]  nom_geom          Nom de l'objet geometrique
  !! \param [in]  st_decoup_poly2d  Active le decoupage des polygones
  !! \param [in]  st_decoup_poly3d  Active le decoupage des polyedres
  !!
  !! \return   Identificateur de l'objet geom dans cs
  !!
  !!




end module pdm_writer
