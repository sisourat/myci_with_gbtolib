program coll_coupling
! From CIcode
  use SlaterDeterminant
  use diagonalizer
  use phis

! From Prema
  USE GlobalMod, ONLY : get_unit,stdlog
  USE ukrmol_interface_gbl, ONLY : START_MPI, FINALIZE_MPI, GET_INTEGRAL, GET_KINETIC_ENERGY_INTEGRAL, LPQRS
  USE PhisMod, ONLY : GBTO_INIT, NBAS, NOCC

  IMPLICIT NONE

integer :: ndim
double precision, dimension(:,:), allocatable :: istavec, fstavec, matcsfs
double precision, dimension(:), allocatable :: istaeig, fstaeig

type csf
integer :: ndets
double precision, dimension(:), allocatable :: coeffs
type(Sdeterminant), dimension(:), allocatable :: dets 
end type csf

integer :: nicsfs, nialpha, nibeta
integer :: nfcsfs, nfalpha, nfbeta
type(csf), dimension(:), allocatable :: icsfs, fcsfs

character(2) :: ca, cb

character(50) :: fista, ffsta
double precision :: prttol, fanomat, mat, mat1e, mat2e
integer :: nista, nfsta, nmainfsta

integer :: n1, n2, n3, iout

integer :: i1, j1

! Other
  INTEGER(2) i,j,k,l
  INTEGER ivec,iounit,iou,ninit,nvec,nchan,ind, nint2e
  REAL(8) time1,time2,iover,gtotal,E0,E1,E2
  LOGICAL isthere
  CHARACTER(4) gprefix
  CHARACTER(11) flnamep

  write(6,'(/,"Starting MPI -- dummy in this version")')
  call START_MPI

!  call GBTO_INIT(.true.)
!  write(6,'(/,"GBTOlib initialized.")')

  call phis_vp

  read(*,*)nista, nfsta, fista, ffsta

  open(unit=22,file=fista)
! reads CSF
   read(22,*)nialpha,nibeta
   read(22,*)nicsfs
   allocate(icsfs(nicsfs),istavec(nicsfs,nista),istaeig(nista))

  do i = 1, nicsfs
    read(22,*)icsfs(i)%ndets
!    write(*,*)i,icsfs(i)%ndets
    allocate(icsfs(i)%dets(icsfs(i)%ndets),icsfs(i)%coeffs(icsfs(i)%ndets))
    do j = 1, icsfs(i)%ndets
      icsfs(i)%dets(j)%nalpha = nialpha
      icsfs(i)%dets(j)%nbeta = nibeta
      allocate(icsfs(i)%dets(j)%alpha(nialpha))
      allocate(icsfs(i)%dets(j)%beta(nibeta))
      read(22,*)icsfs(i)%coeffs(j)
      read(22,*)ca,(icsfs(i)%dets(j)%alpha(k),k=1,nialpha)
      read(22,*)cb,(icsfs(i)%dets(j)%beta(k),k=1,nibeta)
    enddo
  enddo

   do i = 1, nista
     read(22,*)ca,istaeig(i)
     do j = 1, nicsfs
       read(22,*)istavec(j,i)
     enddo
!       write(*,*)sum(istavec(:,i)**2)
   enddo
  close(22)

  open(unit=23,file=ffsta)
! reads CSF
  read(23,*)nfalpha,nfbeta
  read(23,*)nfcsfs!, nmainfsta
  allocate(fcsfs(nfcsfs),fstavec(nfcsfs,nfsta),fstaeig(nfsta))
  do i = 1, nfcsfs
    read(23,*)fcsfs(i)%ndets
    allocate(fcsfs(i)%dets(fcsfs(i)%ndets),fcsfs(i)%coeffs(fcsfs(i)%ndets))
!    write(*,*)i,fcsfs(i)%ndets!fcsfs(i)%coeffs(j),ca,(fcsfs(i)%dets(j)%alpha(k),k=1,nfalpha),cb,(fcsfs(i)%dets(j)%beta(k),k=1,nfbeta)
    do j = 1, fcsfs(i)%ndets
      fcsfs(i)%dets(j)%nalpha = nfalpha
      fcsfs(i)%dets(j)%nbeta = nfbeta
      allocate(fcsfs(i)%dets(j)%alpha(nfalpha))
      allocate(fcsfs(i)%dets(j)%beta(nfbeta))
      read(23,*)fcsfs(i)%coeffs(j),ca,(fcsfs(i)%dets(j)%alpha(k),k=1,nfalpha),cb,(fcsfs(i)%dets(j)%beta(k),k=1,nfbeta)
      !write(*,*)fcsfs(i)%coeffs(j),ca,(fcsfs(i)%dets(j)%alpha(k),k=1,nfalpha),cb,(fcsfs(i)%dets(j)%beta(k),k=1,nfbeta)
    enddo
  enddo

   do i = 1, nfsta
     read(23,*)ca,fstaeig(i)
   write(*,*)ca,fstaeig(i)
     do j = 1, nfcsfs
       read(23,*)fstavec(j,i)
     enddo
   enddo

  close(23)

  allocate(matcsfs(nfcsfs,nicsfs))
  matcsfs(:,:) = 0d0

   do i = 1, nicsfs
    do j = 1, nfcsfs
     do i1 = 1, icsfs(i)%ndets
       do j1 = 1, fcsfs(j)%ndets
          call matvp2dets(icsfs(i)%dets(i1),fcsfs(j)%dets(j1),mat1e)
          mat = mat1e ! no need for coll coupling + mat2e
          matcsfs(j,i) = matcsfs(j,i) + icsfs(i)%coeffs(i1)*fcsfs(j)%coeffs(j1)*mat
       enddo
     enddo

    enddo
   enddo


!    write(10,*)nista, nfsta
    do i = 1, nista
      do j = 1, nfsta
       mat = 0d0
       do i1 = 1, nicsfs
        do j1 = 1, nfcsfs
            mat = mat + fstavec(j1,j)*istavec(i1,i)*matcsfs(j1,i1)
        enddo
       enddo
      write(10,*)istaeig(i),fstaeig(j),mat
     enddo

  enddo

 deallocate(icsfs,istavec,istaeig)
 deallocate(fcsfs,fstavec,fstaeig)
 deallocate(matcsfs)

 call FINALIZE_MPI
 call freephis

end program coll_coupling
