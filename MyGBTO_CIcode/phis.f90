module phis
!!reads only 1e matrix elements from scatci_integrals_vp
implicit none

double precision, dimension(:,:), allocatable :: vpMO, vpMOasymp
integer :: nmos

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine freephis
if(allocated(vpMO)) deallocate(vpMO)
end subroutine freephis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phis_vp
use iso_fortran_env 
implicit none

integer :: i, j, stat
double precision :: buf1,buf2,buf3

 
 open(unit=11,file='moints_vp_asymp')
 read(11,*)nmos
 write(*,*)nmos
 allocate(vpMOasymp(nmos,nmos))
 vpMOasymp(:,:) = 0d0
 do
   read(11,*,iostat=stat) i,j,buf1,buf2,buf3
   if (stat == iostat_end) exit
   vpMOasymp(i,j) = buf3
   vpMOasymp(j,i) = buf3
 end do
 close(11)


 open(unit=11,file='moints_vp')
 read(11,*)nmos
 write(*,*)nmos
 allocate(vpMO(nmos,nmos))
 vpMO(:,:) = 0d0
 do
   read(11,*,iostat=stat) i,j,buf1,buf2,buf3
   if (stat == iostat_end) exit
   vpMO(i,j) = buf3
   vpMO(j,i) = buf3
!   if(abs(buf3)<1d-09) then
!     write(*,*)"WARNING: sign(vpMO(:,:),vpMOasymp(:,:)) - vpMOasymp(:,:) might NOT be OK"
!   endif
 end do
 close(11)

! do i = 1, 2 !nmos
!   do j = i,2 ! nmos
!     write(*,*)j,i,vpMO(j,i),vpMOasymp(j,i)
!   enddo
! enddo
!write(*,*)"1e integrals read"
vpMO(:,:) = sign(vpMO(:,:),vpMOasymp(:,:)) !- vpMOasymp(:,:)
write(*,*)"1e integrals read"

end subroutine phis_vp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module phis
