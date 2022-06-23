program pdm_io_fortran
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Test du bon fonctionnement des procedures :
  ! PDM_io_open
  ! PDM_io_global_write
  ! PDM_io_par_block_write
  ! PDM_io_par_interlaced_write
  ! PDM_io_global_read
  ! PDM_io_par_block_read
  ! PDM_io_par_interlaced_read
  ! PDM_io_close
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_c_binding
  use mpi
  use pdm
  use pdm_io
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  character(80)                        :: name
  type(c_ptr)                          :: id
  integer(4), pointer                  :: ptr_int(:)
  real(8), pointer                     :: ptr_r8 (:)
  complex(8), pointer                  :: ptr_c8 (:)
  character(len=:), pointer            :: buffer=>null()
  type(c_ptr)                          :: cptr
  integer(4)                           :: rank,size,iErr
  integer(kind=pdm_g_num_s)            :: s_data
  integer(kind=pdm_g_num_s)            :: n_data
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call MPI_INIT(iErr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, iErr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(character(len=80) :: buffer)
  allocate(ptr_int(1:1))
  allocate(ptr_r8 (1:1))
  allocate(ptr_c8 (1:1))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Ecriture
  if( rank==0 )print '("Tests d''ecriture")'
  
  name="toto.dat"
  call PDM_io_open(                    &
  &    nom=trim(name)                 ,& !> nom       <-- Nom du fichier
  &    fmt=PDM_IO_FMT_BIN             ,& !> fmt       <-- Fichier text ou binaire
  &    suff_t=PDM_IO_SUFF_MAN         ,& !> suff_t    <-- Type de suffixe (manuel ou automatique)
  &    suff_u=""                      ,& !> suff_u    <-- Suffixe (si suffixe manuel)
  &    s_backup=PDM_IO_BACKUP_OFF     ,& !> s_backup  <-- Active le backup d'un fichier preexistant en mode ecriture
  &    acces=PDM_IO_KIND_MPIIO_EO     ,& !> accesio   <-- Type (parallele avec mpiio, parallele sans mpiio, sequentiel)
  &    mode=PDM_IO_MOD_WRITE          ,& !> mode      <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
 !&    endian=PDM_IO_LITTLEENDIAN     ,& !> PDM_IO_LITTLEENDIAN
  &    endian=PDM_IO_BIGENDIAN        ,& !> PDM_IO_BIGENDIAN           
  &    comm=MPI_COMM_WORLD            ,& !> msg_comm  <-- Communicateur lie au fichier
  &    prop_noeuds_actifs=-1d0        ,& !> - 1: tous les rangs participent  +1: seuls les rangs maitres de chaque noeud participent
  &    unite=id                       ,& !> unite     --> Unite du fichier
  &    ierr=iErr                       ) !> ierr      --> Indique si le fichier est de type PDM_io ou non Utiliser uniquement pour une ouverture en lecture
  
  !>>> Ecriture Globale
  if( rank==0 )print '(3x,"Ecriture Globale")'
  
  ptr_int(1)=10+rank
  n_data=1
  s_data=4
  call PDM_io_global_write(id,s_data,n_data,c_loc(ptr_int))
  if( rank==0 )print '(6x,"ptr_int(1)=",i0)',ptr_int(1)
  
  ptr_r8(1)=1d0
  n_data=1
  s_data=8
  call PDM_io_global_write(id,s_data,n_data,c_loc(ptr_r8))
  if( rank==0 )print '(6x,"ptr_r8 (1)=",e22.15)',ptr_r8(1)
  
  ptr_c8(1)=(1d0,1d0)
  n_data= 1
  s_data=16
  call PDM_io_global_write(id,s_data,n_data,c_loc(ptr_c8))
  if( rank==0 )print '(6x,"ptr_c8 (1)=",e22.15,1x,e22.15)',ptr_c8(1)
  
  write(buffer,'("Ecriture Global depuis le rank: ",i3)')rank
  !buffer(80:80)=C_NULL_CHAR
  n_data=80
  s_data= 1
  call PDM_io_global_write(id,s_data,n_data,c_loc(buffer))
  if( rank==0 )print '(6x,"buffer=""",a,"""")',buffer
  !<<< Ecriture Globale
  
  
  !>>> Ecriture Blocs
  block
    integer                              :: iRank
    integer                              :: iLine,nLines
    integer, pointer                     :: nLinesRank(:)
    character(80)                        :: ligne
    character(80), pointer               :: lignes(:)
    integer(4), pointer                  :: iTab(:)
    integer(kind=pdm_g_num_s)            :: shift
    
    if( rank==0 )print '(3x,"Ecriture Blocs")'
    
    nLines=4+rank !> blocs de taille variable
    allocate(lignes(1:nLines))
    do iLine=1,nLines
      write(ligne,'("Ecriture Bloc Rank: ",i3,2x,"ligne: ",i3)')rank,iLine
      !ligne(80:80)=C_NEW_LINE
      lignes(iLine)=ligne
    enddo
    
    allocate(nLinesRank(0:size-1))
    call mpi_allgather(                   &
    &    nLines     , 1,mpi_integer      ,&
    &    nLinesRank , 1,mpi_integer      ,&
    &    MPI_COMM_WORLD                  ,&
    &    iErr                             )
    shift=sum([(nLinesRank(iRank),iRank=0,rank-1)])+1    !> debut_bloc
    
    allocate(iTab(1:1)) ; iTab(1)=1                    !> <=
    call PDM_io_par_block_write(                    &
    &    fichier=id                                ,&  !> unite
    &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,&  !> t_n_composantes
    &    n_composantes=iTab                        ,&  !> n_composantes
    &    taille_donnee=80                          ,&  !> taille_donnee
    &    n_donnees=nLines                          ,&  !> n_donnees
    &    debut_bloc=shift                          ,&  !> debut_bloc
    &    donnees=c_loc(lignes)                      )  !> donnees
    deallocate(iTab)
    
    deallocate(nLinesRank)
    deallocate(lignes)
  end block
  !<<< Ecriture Blocs
  
  !>>> Ecriture Entrelacee
  block
    integer                              :: iLine,nLines
    character(80)                        :: ligne
    character(80), pointer               :: lignes(:)
    integer(kind=pdm_g_num_s), pointer   :: indirection(:)
    integer, pointer                     :: iTab(:)
    
    if( rank==0 )print '(3x,"Ecriture Entrelacee")'
    
    nLines=4
    allocate(lignes(1:nLines))
    do iLine=1,nLines
      write(ligne,'("Ecriture Entrelacee Rank: ",i3,2x,"ligne: ",i3)')rank,iLine
      !ligne(80:80)=C_NEW_LINE
      lignes(iLine)=ligne
    enddo
    
    !> indirection
    allocate(indirection(1:nLines))
    do iLine=1,nLines
      indirection(iLine)=rank+(iLine-1)*nLines+1 !> Rangement Rank/iLine
    enddo
    
    allocate( iTab(1:1) ) ; iTab(1)=1                 !> <=
    call PDM_io_par_interlaced_write(               &
    &    fichier        =id                        ,& !> unite
    &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,& !> t_n_composantes
    &    n_composantes  =iTab                      ,& !> *n_composantes
    &    taille_donnee  =80                        ,& !> taille_donnee
    &    n_donnees      =nLines                    ,& !> n_donnees
    &    indirection    =indirection               ,& !> *indirection
    &    donnees        =c_loc(lignes)              ) !> *donnees
    deallocate(iTab)
    
    deallocate(indirection)
    
    deallocate(lignes)
  end block
  !<<< Ecriture Entrelacee
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Lecture
  if( rank==0 )print '(/"Tests de lecture")'
  
  call PDM_io_open(                    &
  &    nom=trim(name)                 ,& !> nom       <-- Nom du fichier
  &    fmt=PDM_IO_FMT_BIN             ,& !> fmt       <-- Fichier text ou binaire
  &    suff_t=PDM_IO_SUFF_MAN         ,& !> suff_t    <-- Type de suffixe (manuel ou automatique)
  &    suff_u=""                      ,& !> suff_u    <-- Suffixe (si suffixe manuel)
  &    s_backup=PDM_IO_BACKUP_OFF     ,& !> s_backup  <-- Active le backup d'un fichier preexistant en mode ecriture
  &    acces=PDM_IO_KIND_MPIIO_EO     ,& !> accesio   <-- Type (parallele avec mpiio, parallele sans mpiio, sequentiel)
  &    mode=PDM_IO_MOD_READ           ,& !> mode      <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
 !&    endian=PDM_IO_LITTLEENDIAN     ,& !> PDM_IO_LITTLEENDIAN
  &    endian=PDM_IO_BIGENDIAN        ,& !> PDM_IO_BIGENDIAN           
  &    comm=MPI_COMM_WORLD            ,& !> msg_comm  <-- Communicateur lie au fichier
  &    prop_noeuds_actifs=-1d0        ,& !> - 1: tous les rangs participent  +1: seuls les rangs maitres de chaque noeud participent
  &    unite=id                       ,& !> unite     --> Unite du fichier
  &    ierr=iErr                       ) !> ierr      --> Indique si le fichier est de type PDM_io ou non Utiliser uniquement pour une ouverture en lecture
  
  !>>> Lecture Globale
  if( rank==0 )print '(3x,"Lecture Globale")'
  
  n_data=1
  s_data=4
  call PDM_io_global_read(id,s_data,n_data,c_loc(ptr_int))
  if( rank==0 )print '(6x,"ptr_int(1)=",i0)',ptr_int(1)
  
  n_data=1
  s_data=8
  call PDM_io_global_read(id,s_data,n_data,c_loc(ptr_r8))
  if( rank==0 )print '(6x,"ptr_r8 (1)=",e22.15)',ptr_r8(1)
  
  n_data= 1
  s_data=16
  call PDM_io_global_read(id,s_data,n_data,c_loc(ptr_c8))
  if( rank==0 )print '(6x,"ptr_c8 (1)=",e22.15,1x,e22.15)',ptr_c8(1)
  
  n_data=80
  s_data= 1
  call PDM_io_global_read(id, s_data, n_data, c_loc(buffer))
  if( rank==0 )print '(6x,"buffer: """,a,"""")',buffer  
  !<<< Lecture Globale
  
  !>>> Lecture Blocs
  block
    integer                              :: iRank
    integer                              :: iLine,nLines
    integer, pointer                     :: nLinesRank(:)
    character(80)                        :: ligne
    character(80), pointer               :: lignes(:)
    integer(kind=pdm_g_num_s)            :: shift
    integer, pointer                     :: iTab(:)
    
    if( rank==0 )print '(/3x,"Lecture Blocs")'
    
    nLines=4+rank !> blocs de taille variable
    
    allocate(nLinesRank(0:size-1))
    call mpi_allgather(                   &
    &    nLines     , 1,mpi_integer      ,&
    &    nLinesRank , 1,mpi_integer      ,&
    &    MPI_COMM_WORLD                  ,&
    &    iErr                             )
    shift=sum([(nLinesRank(iRank),iRank=0,rank-1)])+1 !> debut_bloc
    
    allocate(lignes(1:nLines))
    
    allocate(iTab(1:1)) ; iTab(1)=1                    !> <=
    call PDM_io_par_block_read(                     &
    &    fichier=id                                ,&  !> unite
    &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,&  !> t_n_composantes
    &    n_composantes=iTab                        ,&  !> n_composantes
    &    taille_donnee=80                          ,&  !> taille_donnee
    &    n_donnees=nLines                          ,&  !> n_donnees
    &    debut_bloc=shift                          ,&  !> debut_bloc
    &    donnees=c_loc(lignes)                      )  !> donnees
    deallocate(iTab)
        
    do iRank=0,size-1
      if( rank==iRank )then
        print '("rank= ",i0)',rank
        do iLine=1,nLines
          print '(3x,"Lu sur rank",i0,": """,a,"""")',rank,lignes(iLine)
        enddo
      endif
      call mpi_barrier(MPI_COMM_WORLD,iErr)
    enddo
    
    deallocate(lignes)
    lignes=>null()
    
  end block
  !<<< Lecture par bloc

  !>>> Lecture Entrelacee
  block
    integer                              :: iRank
    integer                              :: iLine,nLines
    character(80)                        :: ligne
    character(80), pointer               :: lignes(:)
    integer, pointer                     :: indirection(:)
    integer, pointer                     :: iTab(:)
    type(c_ptr)                          :: cptr
    
    if( rank==0 )print '(/3x,"Lecture Entrelacee")'
    
    nLines=4
    
    !> indirection
    allocate(indirection(1:nLines))
    do iLine=1,nLines
      indirection(iLine)=rank+(iLine-1)*nLines+1 !> Rangement Rank/iLine
    enddo
    
    allocate(lignes(1:nLines))
    
    allocate( iTab(1:1) ) ; iTab(1)=1                 !> <=
    call PDM_io_par_interlaced_read(                &
    &    fichier        =id                        ,& !> unite
    &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,& !> t_n_composantes
    &    n_composantes  =iTab                      ,& !> *n_composantes
    &    taille_donnee  =80                        ,& !> taille_donnee
    &    n_donnees      =nLines                    ,& !> n_donnees
    &    indirection    =indirection               ,& !> *indirection
    &    donnees        =c_loc(lignes)              ) !> *donnees
    deallocate(iTab)
    
    deallocate(indirection)
        
    do iRank=0,size-1
      if( rank==iRank )then
        print '("rank= ",i0)',rank
        do iLine=1,nLines
          print '(3x,"Lu sur rank",i0,": """,a,"""")',rank,lignes(iLine)
        enddo
      endif
      call mpi_barrier(MPI_COMM_WORLD,iErr)
    enddo
    
    deallocate(lignes)
    lignes=>null()
  end block
  !<<< Lecture Entrelacee
  
  call PDM_io_close(id)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Lecture Fortran sur proc0

  block
    integer       :: i,iUnit
    character(1)  :: char
    integer       :: cpt,iChar
    character(80) :: ligne
    integer       :: Reason
    if( rank==0 )then
      print '(/"Lecture Fortran sur rank 0")'
      
      open(newunit=iUnit                        ,&
      &    file=trim(name)                      ,&
      &    access='stream'                      ,&
      &    form='unformatted'                   ,&
      &    convert="BIG_ENDIAN"                 ,&
      &    action='read'                         )
      
      read(iUnit)ptr_int(1) ; print '(3x,"ptr_int(1)=",i0          )',ptr_int(1)
      read(iUnit)ptr_r8 (1) ; print '(3x,"ptr_r8 (1)=",e22.15      )',ptr_r8 (1)
      read(iUnit)ptr_c8 (1) ; print '(3x,"ptr_c8 (1)=",2(e22.15,1x))',ptr_c8 (1)
      read(iUnit)buffer     ; print '(3x,"buffer: """,a,""""       )',buffer
      
      cpt=0
      iChar=0
      lecture: do
        read(iUnit,iostat=Reason)char
        if( Reason>0 )then
         stop "probleme de lecture"
        elseif( Reason<0 )then
         exit lecture
        endif
        cpt=cpt+1
        
        iChar=iChar+1
        !if( char==C_NEW_LINE )then
        if( iChar==80 )then
          ligne(iChar:iChar)=char
          print '(3x,"lecture fortran  iChar=",i0," ligne: """,a,"""")',iChar,ligne
          ligne=""
          iChar=0
        else
          ligne(iChar:iChar)=char
        endif
      enddo lecture
      print '("cpt=",i0)',cpt
      
      close(iUnit)
    endif
  end block
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(ptr_int)
  deallocate(ptr_r8 )
  deallocate(ptr_c8 )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call MPI_FINALIZE(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


end program pdm_io_fortran
