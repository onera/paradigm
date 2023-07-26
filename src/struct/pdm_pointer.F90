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

module pdm_array

  use pdm
  use pdm_fortran
  use iso_c_binding

  implicit none

  type PDM_array_t

    integer               :: type      = -1
    integer               :: length    = -1
    type(c_ptr)           :: cptr   
    integer               :: ownership = PDM_OWNERSHIP_KEEP

  end type PDM_array_t


   interface PDM_array_part_get
     module procedure PDM_array_part_get_int_
     module procedure PDM_array_part_get_int_2_
     module procedure PDM_array_part_get_int_3_
! #ifdef PDM_LONG_G_NUM
!     module procedure PDM_array_part_get_g_num
! #endif
!     module procedure PDM_array_part_get_double
!     module procedure PDM_array_part_get_double_2
!     module procedure PDM_array_part_get_double_3
!     module procedure PDM_array_part_get_complex8
!     module procedure PDM_array_part_get_complex4
!     module procedure PDM_array_part_get_real4
   end interface

  interface PDM_array_create
    module procedure PDM_array_create_
!     module procedure PDM_array_create_int_2
!     module procedure PDM_array_create_int_3
! #ifdef PDM_LONG_G_NUM
!     module procedure PDM_array_part_create_g_num
! #endif
!     module procedure PDM_array_create_get_double
!     module procedure PDM_array_create_get_double_2
!     module procedure PDM_array_create_get_double_3
!     module procedure PDM_array_create_get_complex8
!     module procedure PDM_array_create_get_complex4
!     module procedure PDM_array_create_get_real4
  end interface


  private :: &
!              PDM_array_create_type, &
              PDM_array_part_get_int_, &
              PDM_array_part_get_int_2_, &
              PDM_array_part_get_int_3_, &
!              PDM_array_part_get_double, &
!              PDM_array_part_get_double_2, &
!              PDM_array_part_get_double_3, &
!              PDM_array_part_get_complex8, &
!              PDM_array_part_get_complex4, &
!              PDM_array_part_get_real4, &
! #ifdef PDM_LONG_G_NUM
!              PDM_array_part_get_g_num, &
! #endif
             PDM_array_create_ 
!              PDM_array_part_create_int_2, &
!              PDM_array_part_create_int_3, &
! #ifdef PDM_LONG_G_NUM
!              PDM_array_part_create_g_num, &
! #endif
!              PDM_array_create_get_double, &
!              PDM_array_create_get_double_2, &
!              PDM_array_create_get_double_3, &
!              PDM_array_create_get_complex8, &
!              PDM_array_create_get_complex4, &
!              PDM_array_create_get_real4
  contains


  !>
  !! \brief Initialize a \ref PDM_array_t object
  !!
  !! \param [out]  pa      \ref PDM_array_t object
  !! \param [in]   type    Data type of pointers
  !! \param [in]   length  Length of array
  !!

  subroutine PDM_array_create_ (pa,        &
                               type,      &
                               length,    &
                               cptr,      &
                               ownership)
    use iso_c_binding
    implicit none

    type(PDM_array_t), pointer :: pa
    integer, intent(in)        :: type
    integer, intent(in)        :: length
    type(c_ptr), intent(in)    :: cptr
    integer, intent(in)        :: ownership

    if (associated(pa)) then
      print*, "Error PDM_array_create : pa is already associated ! "
      call exit
    endif

    allocate(pa)

    pa%type      = type
    pa%length    = length
    pa%cptr      = cptr
    pa%ownership = ownership

  end subroutine PDM_array_create_


  !>
  !! \brief Free a \ref PDM_array_t object
  !!
  !! \param [in, out]  pa      \ref PDM_array_t object
  !!

  subroutine PDM_array_free_ (pa)
    use iso_c_binding
    implicit none

    type(PDM_array_t), pointer  :: pa

    if (.not. associated(pa)) then
      print*, "Error PDM_array_free : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%ownership .eq. PDM_OWNERSHIP_KEEP) then
      call pdm_fortran_free_c(pa%cptr)
    endif

    deallocate(pa)
    pa => null()

  end subroutine PDM_array_free_


  !>
  !! \brief Get an array
  !!
  !! Maps a Fortran pointer onto a C pointer
  !!
  !! \param [in]       pa         Array of \ref PDM_array_t
  !! \param [in]       stride1    Optional : dimension 1 of multi-dimension array
  !! \param [in]       stride2    Optional : dimension 2 of multi-dimension array
  !! \param [in]       stride3    Optional : dimension 3 of multi-dimension array
  !! \param [in, out]  pointer_f  Pointer to an integer array
  !!

  subroutine PDM_array_part_get_int_ (pa,        &
                                      pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_array_t),     pointer :: pa
    integer(pdm_l_num_s),  pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_array_part_get_int : wrong type"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_array_part_get_int_


  subroutine PDM_array_part_get_int_2_ (pa,        &
                                        stride1,   & 
                                        pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer(pdm_l_num_s),  pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_array_part_get_int : wrong type"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_array_part_get_int_2_


  subroutine PDM_array_part_get_int_3_ (pa,        &
                                        stride1,   & 
                                        stride2,   & 
                                        pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    integer(pdm_l_num_s),  pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_array_part_get_int : wrong type"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_array_part_get_int_3_

!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in, out]  pointer_f  Pointer to a g_num array
!   !!

! #ifdef PDM_LONG_G_NUM
!   subroutine PDM_array_part_get_g_num (pa,        &
!                                                i_part,    &
!                                                pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     integer(pdm_g_num_s),      pointer :: pointer_f(:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (pa%type .ne. PDM_TYPE_G_NUM) then
!       print *, "PDM_array_part_get_g_num : wrong type"
!       stop
!     end if

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_g_num : wrong i_part"
!       stop
!     end if


!     call c_f_pointer(pa%cptr(i_part+1),     &
!                      pointer_f,             &
!                      [pa%length(i_part+1)])

!   end subroutine PDM_array_part_get_g_num
! #endif


!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in, out]  pointer_f  Pointer to a double array
!   !!

!   subroutine PDM_array_part_get_double (pa,        &
!                                                 i_part,    &
!                                                 pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     double precision,          pointer :: pointer_f(:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (pa%type .ne. PDM_TYPE_DOUBLE) then
!       print *, "PDM_array_part_get_double : wrong type"
!       stop
!     end if

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_double : wrong i_part"
!       stop
!     end if


!     call c_f_pointer(pa%cptr(i_part+1),     &
!                      pointer_f,             &
!                      [pa%length(i_part+1)])

!   end subroutine PDM_array_part_get_double


!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in]       t_stride   Type of stride
!   !! \param [in]       stride     Stride
!   !! \param [in, out]  pointer_f  Pointer to a double array
!   !!

!   subroutine PDM_array_part_get_double_2 (pa,        &
!                                                   i_part,    &
!                                                   t_stride,  &
!                                                   stride,    &
!                                                   pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     integer, intent(in)                :: t_stride
!     integer, intent(in)                :: stride
!     double precision,          pointer :: pointer_f(:,:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (pa%type .ne. PDM_TYPE_DOUBLE) then
!       print *, "PDM_array_part_get_double_2 : wrong type"
!       stop
!     end if

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_double_2 : wrong i_part"
!       stop
!     end if

!     select case (t_stride)
!     case (PDM_STRIDE_CST_INTERLACED)
!       call c_f_pointer(pa%cptr(i_part+1),     &
!                        pointer_f,             &
!                        [stride,pa%length(i_part+1)/stride])
!     case (PDM_STRIDE_CST_INTERLEAVED)
!       call c_f_pointer(pa%cptr(i_part+1),     &
!                        pointer_f,             &
!                        [pa%length(i_part+1)/stride,stride])
!     case default
!       print *, "PDM_array_part_get_double_2 : wrong stride type"
!       stop
!     end select

!   end subroutine PDM_array_part_get_double_2


!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in]       t_stride   Type of stride
!   !! \param [in]       stride 1   Stride 1
!   !! \param [in]       stride 2   Stride 2
!   !! \param [in, out]  pointer_f  Pointer to a double array
!   !!

!   subroutine PDM_array_part_get_double_3 (pa,        &
!                                                   i_part,    &
!                                                   t_stride,  &
!                                                   stride1,   &
!                                                   stride2,   &
!                                                   pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     integer, intent(in)                :: t_stride
!     integer, intent(in)                :: stride1
!     integer, intent(in)                :: stride2
!     double precision,          pointer :: pointer_f(:,:,:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (pa%type .ne. PDM_TYPE_DOUBLE) then
!       print *, "PDM_array_part_get_double_3 : wrong type"
!       stop
!     end if

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_double_3 : wrong i_part"
!       stop
!     end if

!     select case (t_stride)
!     case (PDM_STRIDE_CST_INTERLACED)
!       call c_f_pointer(pa%cptr(i_part+1),     &
!                        pointer_f,             &
!                        [stride1,stride2,pa%length(i_part+1)/(stride1*stride2)])
!     case (PDM_STRIDE_CST_INTERLEAVED)
!       call c_f_pointer(pa%cptr(i_part+1),     &
!                        pointer_f,             &
!                        [pa%length(i_part+1)/(stride1*stride2),stride1,stride2])
!     case default
!       print *, "PDM_array_part_get_double_3 : wrong stride type"
!       stop
!     end select

!   end subroutine PDM_array_part_get_double_3


!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in, out]  pointer_f  Pointer to a real4 array
!   !!

!   subroutine PDM_array_part_get_real4 (pa,        &
!                                                i_part,    &
!                                                pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     real (kind=4),             pointer :: pointer_f(:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (associated(pa)) then
!       if (pa%type .ne. PDM_TYPE_REAL4) then
!         print *, "PDM_array_part_set_double : wrong type"
!         stop
!       end if

!       if (i_part .ge. size(pa%cptr)) then
!         print *, "PDM_array_part_set_double : wrong i_part"
!         stop
!       end if


!       call c_f_pointer(pa%cptr(i_part+1),     &
!                        pointer_f,             &
!                        [pa%length(i_part+1)])
!     endif

!   end subroutine PDM_array_part_get_real4

!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in, out]  pointer_f  Pointer to a complex4 array
!   !!

!   subroutine PDM_array_part_get_complex4 (pa,        &
!                                                i_part,    &
!                                                pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     complex (kind=4),             pointer :: pointer_f(:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (pa%type .ne. PDM_TYPE_COMPLEX4) then
!       print *, "PDM_array_part_get_complex4 : wrong type"
!       stop
!     end if

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_complex4 : wrong i_part"
!       stop
!     end if


!     call c_f_pointer(pa%cptr(i_part+1),     &
!                      pointer_f,             &
!                      [pa%length(i_part+1)])

!   end subroutine PDM_array_part_get_complex4


!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in, out]  pointer_f  Pointer to a complex8 array
!   !!

!   subroutine PDM_array_part_get_complex8 (pa,        &
!                                                i_part,    &
!                                                pointer_f)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     complex (kind=8),             pointer :: pointer_f(:)

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (pa%type .ne. PDM_TYPE_COMPLEX8) then
!       print *, "PDM_array_part_get_complex8 : wrong type"
!       stop
!     end if

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_complex8 : wrong i_part"
!       stop
!     end if


!     call c_f_pointer(pa%cptr(i_part+1),     &
!                      pointer_f,             &
!                      [pa%length(i_part+1)])

!   end subroutine PDM_array_part_get_complex8



!   !>
!   !! \brief Get a partition
!   !!
!   !! Maps a Fortran pointer onto a C pointer
!   !!
!   !! \param [in]       pa         Array of \ref PDM_array_t
!   !! \param [in]       i_part     Id of partition
!   !! \param [in, out]  pointer_f  Pointer to a complex8 array
!   !!

!   subroutine PDM_array_part_get_cptr (pa,        &
!                                                i_part,    &
!                                                pointer_c,     &
!                                                length)
!     use iso_c_binding
!     implicit none

!     type(PDM_array_t), pointer  :: pa
!     integer, intent(in)                :: i_part
!     type(c_ptr)                        :: pointer_c
!     integer                            :: length

!     if (.not. associated(pa)) then
!       print*, "Error PDM_array_part_get : 'pa' pointer is not associated "
!       call exit
!     endif

!     if (i_part .ge. size(pa%cptr)) then
!       print *, "PDM_array_part_get_complex8 : wrong i_part"
!       stop
!     end if

!     pointer_c = pa%cptr(i_part+1)
!     length    = pa%length(i_part+1)

!   end subroutine PDM_array_part_get_cptr



end module pdm_array
