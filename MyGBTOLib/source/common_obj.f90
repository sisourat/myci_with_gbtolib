! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
module common_obj_gbl
use precisn_gbl
use utils_gbl, only: xermsg
use const_gbl, only: level2, level3, nuc_nam_len

   !> \class <darray_1d>
   !> 1D array of reals
   type darray_1d
      !> The array values.
      real(kind=cfp), allocatable :: a(:)
      !> Dimension of the array.
      integer :: d1 = 0
   end type darray_1d

   !> \class <array_2d>
   !> 2D array of reals
   type darray_2d
      !> The array values.
      real(kind=cfp), allocatable :: a(:,:)
      !> First dimension of the array.
      integer :: d1 = 0
      !> Second dimension of the array.
      integer :: d2 = 0
   contains
      !> Allocates the space for the array a. The values d1, d2 must be set before. Returns a non-zero number if allocation was unsucessfull.
      procedure :: init => init_darray_2d
      !> Deallocates space. Returns a non-zero number if the array has not been allocated before.
      procedure :: final => final_darray_2d
   end type darray_2d

   !> \class <iarray_1d>
   !> 1D array of integers
   type iarray_1d
      integer, allocatable :: a(:)
      !> Dimension of the array.
      integer :: d1 = 0
   end type iarray_1d

   !> \class <nucleus_type>
   type nucleus_type
      !> Cartesian coordinates of the center of the nucleus.
      real(kind=cfp) :: center(1:3)
      !> Charge on the nucleus. In theory this is integer, but it may be useful to allow non-integer values for some experiments.
      !> If this "nucleus" represents the center for the basis of the continuum, then Z must be set to 0.
      real(kind=cfp) :: charge
      !> ID of the nucleus, i.e. its number.
      integer :: nuc
      !> Nucleus name, i.e. atomic symbol.
      character(len=nuc_nam_len) :: name
   contains
      !> This function checks that if nuc = 0 then the corresponding distance from the origin is 0, i.e. nuc = 0 is reserved for CMS functions. 
      procedure :: check => check_nucleus_type
      !> Print the nuclear data.
      procedure :: print => print_nucleus_type
      !> Returns .true. if the nucleus has all attributes of the scattering centre: nuc = 0; Z = 0.0_cfp; nam = 'sc'; center = 0.0_cfp
      procedure :: is_continuum => is_scattering_center
   end type nucleus_type

   interface resize_array
      module procedure resize_2d_array_cfp
      module procedure resize_3d_array_cfp
      module procedure resize_2d_array_int
      module procedure resize_3d_array_int
   end interface resize_array

   !private routines
   private print_nucleus_type, check_nucleus_type, is_scattering_center

contains

   subroutine print_orbital_table(energies,num_sym,n_orbs,n_sym,convert_to_eV)
      use sort_gbl, only: cfp_sort_float_int_1d
      use phys_const_gbl, only: to_eV
      implicit none
      integer, intent(in) :: n_orbs, n_sym, num_sym(2,n_orbs)
      real(kind=cfp) :: energies(n_orbs)
      logical, intent(in) :: convert_to_eV

      integer, allocatable :: permutation(:), nob(:)
      integer :: err, i, num, sym

         allocate(permutation(n_orbs),nob(n_sym),stat=err)
         if (err .ne. 0) call xermsg ('common_obj', 'print_orbital_table', 'Memory allocation failed.',err,1)

         nob = 0

         do i=1,n_orbs
            permutation(i) = i
         enddo !i

         !sort the orbitals according to their energies (if required)
         write(level2,'(/,10X,"Sorting orbital energies...")')
         call cfp_sort_float_int_1d(n_orbs,energies,permutation)
         write(level2,'("...done")')
   
         !convert the energies from H to eV
         if (convert_to_eV) energies = energies*to_eV
   
         write(level2,'(/,5X,"List of Energy-sorted molecular orbitals (number.symmetry) follows:",/)')
         write(level2,'(5X,"|No |NUM.SYMM|ENERGY[eV]|",8(" ",i1," |"))') (i,i=1,n_sym)
   
         do i=1,n_orbs
            num = num_sym(1,permutation(i))
            sym = num_sym(2,permutation(i))
            nob(sym) = nob(sym) + 1
            if (i > 1) then
               if (energies(i)*energies(i-1) < 0.0_cfp) then
                  write(level2,'(5X,25("-"))')
               endif
            endif
            write(level2,'(5X,"|",i3,"|",i6,".",i1,"|",F10.3,"|",8(i3,"|"))') i,num,sym, energies(i), nob(1:n_sym)
         enddo

   end subroutine print_orbital_table

   function init_darray_2d(this)
      implicit none
      class(darray_2d) :: this
      integer :: init_darray_2d

      integer :: err
      
         init_darray_2d = 0

         if (this%d1 .le. 0) then 
            init_darray_2d = 1
            return
         endif

         if (this%d2 .le. 0) then
            init_darray_2d = 2
            return
         endif

         if (allocated(this%a)) deallocate(this%a)

         allocate(this%a(this%d1,this%d2),stat=err)
         if (err .ne. 0) init_darray_2d = 3
      
   end function init_darray_2d

   function final_darray_2d(this)
      implicit none
      class(darray_2d) :: this
      integer :: final_darray_2d

         final_darray_2d = 0

         if (allocated(this%a)) then
            deallocate(this%a)
            this%d1 = 0
            this%d2 = 0
         else
            final_darray_2d = 1
            return
         endif

   end function final_darray_2d

   !> \todo Sort out the nuc == 0 question.
   function check_nucleus_type(this)
      implicit none
      class(nucleus_type) :: this
      integer :: check_nucleus_type

         check_nucleus_type = 0

         !nuc = 0 is reserved only for nuclei centered on the CMS. todo This is not true and it should be probably removed.
!         if (dot_product(this%center,this%center) == 0.0_cfp) then
!            if (this%nuc .ne. 0) then
!               check_nucleus_type = 1
!               call this%print
!               call xermsg('common_obj','check_nucleus_type','nuc = 0 must be used for nuclei centered on the CMS.',1,1)
!            endif
!         endif

         if (this%nuc .eq. 0 .and. dot_product(this%center,this%center) .ne. 0.0_cfp) then
            check_nucleus_type = 2
            call this%print
            call xermsg ('common_obj', 'check_nucleus_type', &
                         'nuc = 0 (i.e. nucleus at the CMS), but the distance from the origin is not 0.', 2, 1)
         endif

   end function check_nucleus_type

   subroutine print_nucleus_type(this)
      implicit none
      class(nucleus_type) :: this

         write(level2,'("Nucleus name: ",a)') adjustl(this%name)
         write(level2,'("Nucleus number: ",i5)') this%nuc
         write(level2,'("Nucleus Z: ",f20.15)') this%charge
         write(level2,'("Nucleus center: ",3f25.15)') this%center(:)

   end subroutine print_nucleus_type

   function is_scattering_center(this)
      implicit none
      class(nucleus_type) :: this
      logical :: is_scattering_center

      real(kind=cfp) :: r

         is_scattering_center = .false.

         r = dot_product(this%center,this%center)
         if (this%name .eq. 'sc' .and. this%nuc .eq. 0 .and. r .eq. 0.0_cfp) is_scattering_center = .true.

   end function

   function resize_2d_array_cfp(array,d1,d2)
      implicit none
      integer, intent(in) :: d1, d2
      real(kind=cfp), allocatable :: array(:,:)
      integer :: resize_2d_array_cfp

      logical :: do_allocation
      integer :: err

         resize_2d_array_cfp = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2),stat=err)
            resize_2d_array_cfp = err
         endif

   end function resize_2d_array_cfp

   function resize_3d_array_cfp(array,d1,d2,d3)
      implicit none
      integer, intent(in) :: d1, d2, d3
      real(kind=cfp), allocatable :: array(:,:,:)
      integer :: resize_3d_array_cfp

      logical :: do_allocation
      integer :: err

         resize_3d_array_cfp = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2 .or. size(array,3) < d3) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2,d3),stat=err)
            resize_3d_array_cfp = err
         endif

   end function resize_3d_array_cfp

   function resize_2d_array_int(array,d1,d2)
      implicit none
      integer, intent(in) :: d1, d2
      integer, allocatable :: array(:,:)
      integer :: resize_2d_array_int

      logical :: do_allocation
      integer :: err

         resize_2d_array_int = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2),stat=err)
            resize_2d_array_int = err
         endif

   end function resize_2d_array_int

   function resize_3d_array_int(array,d1,d2,d3)
      implicit none
      integer, intent(in) :: d1, d2, d3
      integer, allocatable :: array(:,:,:)
      integer :: resize_3d_array_int

      logical :: do_allocation
      integer :: err

         resize_3d_array_int = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2 .or. size(array,3) < d3) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2,d3),stat=err)
            resize_3d_array_int = err
         endif

   end function resize_3d_array_int

   subroutine resize_copy_2d_array(array,d1,d2)
       implicit none
       real(kind=cfp), allocatable :: array(:,:)
       integer, intent(in) :: d1,d2

       real(kind=cfp), allocatable :: tmp(:,:)
       integer :: err, d1_al, d2_al, d1_old, d2_old

          d1_al = d1
          d2_al = d2

          d1_old = size(array,1)
          d2_old = size(array,2)

          if (.not.(allocated(array))) then
             allocate(array(d1_al,d2_al),stat=err)
             if (err .ne. 0) call xermsg('common_obj','resize_copy_2d_array','Memory allocation 1 failed.',err,1)
             array = 0.0_cfp
          else
             d1_al = max(d1_al,d1_old)
             d2_al = max(d2_al,d2_old)

             call move_alloc(array,tmp)

             allocate(array(d1_al,d2_al),stat=err)
             if (err .ne. 0) call xermsg('common_obj','resize_copy_2d_array','Memory allocation 2 failed.',err,1)
             array = 0.0_cfp

             array(1:d1_old,1:d2_old) = tmp(1:d1_old,1:d2_old)
             deallocate(tmp)
          endif

   end subroutine resize_copy_2d_array

end module common_obj_gbl
