PROGRAM myCI

! From CIcode
  use SlaterDeterminant
  use diagonalizer

! From Prema
  USE GlobalMod, ONLY : get_unit,stdlog
  USE ukrmol_interface_gbl, ONLY : START_MPI, FINALIZE_MPI, GET_INTEGRAL, GET_KINETIC_ENERGY_INTEGRAL, LPQRS
  USE PhisMod, ONLY : GBTO_INIT, NBAS, NOCC

  IMPLICIT NONE

! From CIcode
  integer :: ndim
  double precision, dimension(:,:), allocatable :: cimat, cimat1e, cimat2e
  double precision, dimension(:), allocatable :: eigci, tempa

  type csf
   integer :: ndets
   double precision, dimension(:), allocatable :: coeffs
   type(Sdeterminant), dimension(:), allocatable :: dets 
  end type csf

  integer :: ncsfs, nalpha, nbeta
  type(csf), dimension(:), allocatable :: csfs
  double precision :: mat, mat1e, mat2e, ovl

  type(diag) :: mydiag
  character(8) :: cityp
  character(1) :: ca, cb

  double precision :: prttol
  integer :: nsta, nmaincsfs

  character(50) :: fileigvec, filinput

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

  call GBTO_INIT(.true.)
  write(6,'(/,"GBTOlib initialized.")')

! From CIcode
  mydiag%typ = 'exact'
  CALL getarg(1, filinput)
  write(*,*)filinput

  open(unit=21,file=filinput)
   read(21,*)nsta, prttol
!   write(*,*)nsta, prttol
   read(21,*)fileigvec
!   write(*,*)fileigvec

  open(unit=22,file=fileigvec)

! reads CSF
   read(21,*)nalpha,nbeta
   read(21,*)ncsfs !,nmaincsfs
   write(22,*)nalpha,nbeta
   write(22,*)ncsfs !,nmaincsfs
   allocate(csfs(ncsfs))
   do i = 1, ncsfs
     read(21,*)csfs(i)%ndets
     write(22,*)csfs(i)%ndets
!     write(*,*)i
     allocate(csfs(i)%dets(csfs(i)%ndets),csfs(i)%coeffs(csfs(i)%ndets))
     do j = 1, csfs(i)%ndets
       csfs(i)%dets(j)%nalpha = nalpha
       csfs(i)%dets(j)%nbeta = nbeta
       allocate(csfs(i)%dets(j)%alpha(nalpha))
       allocate(csfs(i)%dets(j)%beta(nbeta))
       read(21,*)csfs(i)%coeffs(j),ca,(csfs(i)%dets(j)%alpha(k),k=1,nalpha),cb,(csfs(i)%dets(j)%beta(k),k=1,nbeta)
       write(22,*)csfs(i)%coeffs(j)
       write(22,'(a,500(i4,1X))')ca,(csfs(i)%dets(j)%alpha(k),k=1,nalpha)
       write(22,'(a,500(i4,1X))')cb,(csfs(i)%dets(j)%beta(k),k=1,nbeta)
     enddo
   enddo

!! creates CI matrix

 ndim = ncsfs
 allocate(cimat(ndim,ndim),cimat1e(ndim,ndim),cimat2e(ndim,ndim),eigci(ndim),tempa(ndim))

 cimat(:,:) = 0d0
 cimat1e(:,:) = 0d0
 cimat2e(:,:) = 0d0

 do i = 1, ndim
   do j = 1, ndim 
     
    do i1 = 1, csfs(i)%ndets
      do j1 = 1, csfs(j)%ndets
         call mat2dets(csfs(i)%dets(i1),csfs(j)%dets(j1),mat1e,mat2e,ovl)
         cimat(j,i) = cimat(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*(mat1e+mat2e)
         cimat1e(j,i) = cimat1e(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*mat1e
         cimat2e(j,i) = cimat2e(j,i) + csfs(i)%coeffs(i1)*csfs(j)%coeffs(j1)*mat2e
!         write(*,'(4(i3,1X),3(f20.15,1X))')j,i,j1,i1,csfs(i)%coeffs(i1),csfs(j)%coeffs(j1),mat1e,mat2e
!         cimat(j,i) = i*j*1d0
      enddo
    enddo
!         write(*,'(2(i3,1X),2(f20.15,1X))')j,i,cimat1e(j,i),cimat2e(j,i)
 
   enddo
enddo

!! diagonalizes the CI matrix

 write(*,*)"DIAG"
 call eig(mydiag,ndim,cimat,eigci)
 write(*,*)"NSTATES=",ndim

  do i = 1, nsta
    write(22,*)
    write(22,*)'E=',eigci(i)
    write(*,*)i,eigci(i)
    do j = 1, ndim
!     if(abs(cimat(j,i))<prttol) cycle
      write(22,*)cimat(j,i),"     ",j
    enddo
     write(22,*)
!     tempa(:) = matmul(cimat1e(:,:),cimat(:,i))
!     write(*,*)" ONE ELECTRON ENERGY = ", dot_product(cimat(:,i),tempa)
!     tempa(:) = matmul(cimat2e(:,:),cimat(:,i))
!     write(*,*)" TWO ELECTRON ENERGY = ", dot_product(cimat(:,i),tempa)
  enddo
 close(21)
 close(22)

 deallocate(csfs,cimat,cimat1e,cimat2e,eigci,tempa)
 call FINALIZE_MPI

END PROGRAM myCI
