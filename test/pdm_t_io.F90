program pdm_io_fortran
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
  integer, pointer                     :: ptr_int(:)
  real(8), pointer                     :: ptr_r8 (:)
  complex(8), pointer                  :: ptr_c8 (:)
  character(len=:), pointer            :: buffer=>null()
  type(c_ptr)                          :: cptr
  integer                              :: rank,size,iErr
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

  name="toto.dat"
  
  call PDM_io_open(                    &
  &    nom=trim(name)                 ,& !> nom       <-- Nom du fichier
  &    fmt=PDM_IO_FMT_BIN             ,& !> fmt       <-- Fichier text ou binaire
  &    suff_t=PDM_IO_SUFF_MAN         ,& !> suff_t    <-- Type de suffixe (manuel ou automatique)
  &    suff_u=""                      ,& !> suff_u    <-- Suffixe (si suffixe manuel) 
  &    s_backup=PDM_IO_BACKUP_OFF     ,& !> s_backup  <-- Active le backup d'un fichier preexistant en mode ecriture 
  &    acces=PDM_IO_KIND_MPIIO_EO     ,& !> accesio   <-- Type (parallele avec mpiio, parallele sans mpiio, sequentiel) 
  &    mode=PDM_IO_MOD_WRITE          ,& !> mode      <-- Mode d'acces (lecture, ecriture, lecture/ecriture) 
  &    endian=PDM_IO_LITTLEENDIAN     ,& !> PDM_io_littleendian            
  &    comm=MPI_COMM_WORLD            ,& !> msg_comm  <-- Communicateur lie au fichier 
  &    prop_noeuds_actifs=-1d0        ,& !> - 1: tous les rangs participent  +1: seuls les rangs maitres de chaque noeud participent
  &    unite=id                       ,& !> unite     --> Unite du fichier 
  &    ierr=iErr                       ) !> ierr      --> Indique si le fichier est de type PDM_io ou non Utiliser uniquement pour une ouverture en lecture
  
  !> Ecritrue Globale
  ptr_int(1)=1
  print '("0 ptr_int(1)=",i0)',ptr_int(1)
  call PDM_io_global_write(id,4,1,c_loc(ptr_int))
  print '("1 ptr_int(1)=",i0)',ptr_int(1)

  ptr_r8(1)=1d0
  call PDM_io_global_write(id,1,4,c_loc(ptr_int))

  ptr_c8(1)=(1d0,1d0)
  call PDM_io_global_write(id,1,4,c_loc(ptr_int))
  
  buffer="Bonjour !" ; buffer(80:80)=C_NULL_CHAR
  call PDM_io_global_write(id,1,80,c_loc(buffer))

  call PDM_io_close(id)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Lecture
  call PDM_io_open(                    &
  &    nom=trim(name)                 ,& !> nom       <-- Nom du fichier
  &    fmt=PDM_IO_FMT_BIN             ,& !> fmt       <-- Fichier text ou binaire
  &    suff_t=PDM_IO_SUFF_MAN         ,& !> suff_t    <-- Type de suffixe (manuel ou automatique)
  &    suff_u=""                      ,& !> suff_u    <-- Suffixe (si suffixe manuel) 
  &    s_backup=PDM_IO_BACKUP_OFF     ,& !> s_backup  <-- Active le backup d'un fichier preexistant en mode ecriture 
  &    acces=PDM_IO_KIND_MPIIO_EO     ,& !> accesio   <-- Type (parallele avec mpiio, parallele sans mpiio, sequentiel) 
  &    mode=PDM_IO_MOD_READ           ,& !> mode      <-- Mode d'acces (lecture, ecriture, lecture/ecriture) 
  &    endian=PDM_IO_LITTLEENDIAN     ,& !> PDM_io_littleendian            
  &    comm=MPI_COMM_WORLD            ,& !> msg_comm  <-- Communicateur lie au fichier 
  &    prop_noeuds_actifs=-1d0        ,& !> - 1: tous les rangs participent  +1: seuls les rangs maitres de chaque noeud participent
  &    unite=id                       ,& !> unite     --> Unite du fichier 
  &    ierr=iErr                       ) !> ierr      --> Indique si le fichier est de type PDM_io ou non Utiliser uniquement pour une ouverture en lecture
  
  call PDM_io_global_read(id,1,4,cptr)
  call c_f_pointer(cptr=cptr, fptr=ptr_int, shape=[1])
  if( rank==0 )print '("ptr_int(1)=",i0)',ptr_int(1)
  
  
  call PDM_io_global_read(id, 1, 8, cptr)
  call c_f_pointer(cptr=cptr, fptr=ptr_r8 , shape=[1])
  if( rank==0 )print '("ptr_r8 (1)=",e22.15)',ptr_r8(1)

  call PDM_io_global_read(id, 1,16, cptr)
  call c_f_pointer(cptr=cptr, fptr=ptr_c8 , shape=[1])
  if( rank==0 )print '("ptr_r8 (1)=",e22.15,1x,e22.15)',ptr_c8(1)
  
  call PDM_io_global_read(id, 1, 80, cptr)
  call c_f_pointer(cptr=cptr, fptr=buffer, shape=[80])
  if( rank==0 )print '("buffer: ",a)',buffer


  !if( rank==0 )then
  !  
  !  open(newunit=iUnit                        ,&
  !  &    file=trim(name)                      ,&
  !  &    access='stream'                      ,&
  !  &    form='unformatted'                   ,&
  ! !&    convert="BIG_ENDIAN"                 ,&
  !  &    action='read'                         )
  !  
  !  read(iUnit) ptr_int(1) ; print '("ptr(1)=",i0)',ptr(1)
  !  read(iUnit) ptr_r8 (1) ; print '("ptr(1)=",i0)',ptr(1)
  !  read(iUnit) ptr_c8 (1) ; print '("ptr(1)=",i0)',ptr(1)
  !  close (iUnit)
  !endif

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(buffer)
  deallocate(ptr_int)
  deallocate(ptr_r8 )
  deallocate(ptr_c8 )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call MPI_FINALIZE(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  
end program pdm_io_fortran
