program Prop
use propdyn
use general, only : lenmax
implicit none

  character(len=lenmax) :: input, finput
  character(len=lenmax) :: fout

  double precision :: spac

  double precision :: t1, t2, tmec, tdyn, ttot, bproj, vproj
  integer :: i, j, k, ista

!  call getarg(1,finput)
!  call getarg(2,input)

!  read(input,*)ista

ista = 1

  open(unit=10,file='matcoll')
    read(10,*)bproj, vproj
    write(*,*)bproj, vproj
    read(10,*)ntotsta, ntime
    write(*,*)ntotsta
    allocate(mcoup(ntime,ntotsta,ntotsta),tgrid(ntime),esta(ntotsta),mat(ntotsta,ntotsta))

!    do i = 1, ntotsta
!     read(10)esta(i)
!    enddo

    do i = 1, ntime
      read(10,*)tgrid(i)
      tgrid(i) = tgrid(i) / vproj
      do j = 1, ntotsta
        do k = 1, ntotsta
          read(10,*)esta(j),esta(k),spac
          mcoup(i,j,k) = spac
        enddo
      enddo
!      write(*,*)tgrid(i),real(mcoup(i,1,1))
    enddo
  close(10)

  allocate(rmat2intrp(ntime,ntotsta,ntotsta),cmat2intrp(ntime,ntotsta,ntotsta))
  allocate(matintrp(ntotsta,ntotsta))
  allocate(psi(ntotsta))
  
  psi(:) = 0d0
  psi(ista) = 1d0
   
  call cpu_time(t1)
  call dyn
  call cpu_time(t2)
  tdyn = t2-t1
  write(*,*)'DYN takes',tdyn
  write(*,'(5000(f12.6,1X))')(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2)
  write(100,'(5000(f12.6,1X))')bproj,(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2),vproj

  deallocate(mcoup,matintrp,psi,rmat2intrp,cmat2intrp)

end program Prop
