!! Nico 11.04.13
!! module SlaterDeterminant
!! defines a determinant, tools for dealing with dets
!! includes compare2dets, countunpairel

module SlaterDeterminant
implicit none

type Sdeterminant
 integer :: nunpairel, nalpha, nbeta
 integer, dimension(:), allocatable :: alpha, beta
end type Sdeterminant

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compare2dets(deta,detb,cperm,ndiffalpha,orb1alpha,orb2alpha,ndiffbeta,orb1beta,orb2beta)
! returns the number of differents orbitals 
! and the latters (returns -1 if more than 2 differents orbitals)
implicit none

type(Sdeterminant), intent(in) :: deta, detb
type(Sdeterminant) :: det1, det2
integer, intent(out) :: ndiffalpha, ndiffbeta
integer, dimension(2), intent(out) :: orb1alpha, orb2alpha, orb1beta, orb2beta
double precision, intent(out) :: cperm
integer :: tmp, permut
character(1) :: same

integer :: i, j

det1 = deta
det2 = detb

ndiffalpha = 0
ndiffbeta = 0

orb1alpha(:) = 0
orb2alpha(:) = 0
orb1beta(:) = 0
orb2beta(:) = 0

permut=0
!write(*,*)"compare alpha",det1%alpha,det2%alpha
!write(*,*)"compare beta",det1%beta,det2%beta

! check alpha electrons

do i = 1, det1%nalpha
  same = 'f'

  do j = 1, det2%nalpha
     if(det1%alpha(i)==det2%alpha(j)) then
       same='t'
       if(i/=j) then
        tmp = det2%alpha(i)
        det2%alpha(i) = det2%alpha(j)
        det2%alpha(j) = tmp
        permut = permut + 1
        cycle
       endif
     endif
  enddo

 if(same=='f') then
   ndiffalpha = ndiffalpha + 1
   if(ndiffalpha>2) then
     cperm = 0d0
     ndiffalpha = -10
     return
   endif
 endif
enddo

 j = 0
do i = 1, det1%nalpha
  if(det1%alpha(i)/=det2%alpha(i)) then
   j = j + 1
   orb1alpha(j) = det1%alpha(i)
   orb2alpha(j) = det2%alpha(i)
  endif
enddo

! check beta electrons

do i = 1, det1%nbeta
  same = 'f'

  do j = 1, det2%nbeta
     if(det1%beta(i)==det2%beta(j)) then
       same='t'
       if(i/=j) then
        tmp = det2%beta(i)
        det2%beta(i) = det2%beta(j)
        det2%beta(j) = tmp
        permut = permut + 1
        cycle
       endif
     endif
  enddo

 if(same=='f') then
   ndiffbeta = ndiffbeta + 1
   if(ndiffalpha+ndiffbeta>2) then
     cperm = 0d0
     ndiffbeta = -10
     return
   endif
 endif

enddo

 j = 0
do i = 1, det1%nbeta
  if(det1%beta(i)/=det2%beta(i)) then
   j = j + 1
   orb1beta(j) = det1%beta(i)
   orb2beta(j) = det2%beta(i)
  endif
enddo

!write(*,*)ndiffalpha,ndiffbeta

     cperm = (-1d0)**permut

end subroutine compare2dets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mat2dets(det1,det2,mat1e,mat2e,ovl)
  USE GlobalMod, ONLY : get_unit,stdlog
  USE ukrmol_interface_gbl, ONLY : START_MPI, FINALIZE_MPI, GET_INTEGRAL, GET_KINETIC_ENERGY_INTEGRAL, LPQRS
  USE PhisMod, ONLY : GBTO_INIT, NBAS, NOCC
!! returns the Hamiltonian matrix element between det1 and det2
!! using the Slater Condon rules
!! spin-integration is taken care here
implicit none

type(Sdeterminant), intent(in) :: det1, det2
double precision :: mat1e, mat2e, ovl, h1e, g2e, cperm
integer :: ndiffalpha, ndiffbeta, ndiff
integer, dimension(2) :: orb1alpha, orb2alpha, orb1beta, orb2beta
integer :: ie, je, ke, le, inde, intindex

integer :: i, j

mat1e = 0d0
mat2e = 0d0
ovl = 0d0
h1e = 0d0
g2e = 0d0

call compare2dets(det1,det2,cperm,ndiffalpha,orb1alpha,orb2alpha,ndiffbeta,orb1beta,orb2beta)
ndiff = ndiffalpha + ndiffbeta

if(ndiff==0) then
  ovl = 1d0
  do i = 1, det1%nalpha
      ie =  det1%alpha(i)
      h1e = h1e + GET_INTEGRAL(ie,ie,0,0,0) !h1eMO(ie,ie)
    do j = i+1, det2%nalpha
      je =  det2%alpha(j)
!      inde = intindex(ie,ie,je,je)
      g2e = g2e + GET_INTEGRAL(ie,ie,je,je,0) !int2e(inde)
!      inde = intindex(ie,je,ie,je)
      g2e = g2e - GET_INTEGRAL(ie,je,ie,je,0) ! int2e(inde)
    enddo
  enddo

  do i = 1, det1%nbeta
      ie =  det1%beta(i)
!      h1e = h1e + h1eMO(ie,ie)
      h1e = h1e + GET_INTEGRAL(ie,ie,0,0,0) !h1eMO(ie,ie)
    do j = i+1, det2%nbeta
      je =  det2%beta(j)
!      inde = intindex(ie,ie,je,je)
!      g2e = g2e + int2e(inde)
      g2e = g2e + GET_INTEGRAL(ie,ie,je,je,0) !int2e(inde)
!      inde = intindex(ie,je,ie,je)
!      g2e = g2e - int2e(inde)
      g2e = g2e - GET_INTEGRAL(ie,je,ie,je,0) ! int2e(inde)
    enddo
  enddo

  do i = 1, det1%nalpha
      ie =  det1%alpha(i)
    do j = 1, det2%nbeta
      je =  det2%beta(j)
!      inde = intindex(ie,ie,je,je)
!      g2e = g2e + int2e(inde)
      g2e = g2e + GET_INTEGRAL(ie,ie,je,je,0) !int2e(inde)
    enddo
  enddo

 mat1e = cperm*h1e
 mat2e = cperm*g2e

elseif(ndiff==1) then

   if(ndiffalpha==1 .and. ndiffbeta==0) then
    ie=orb1alpha(1)
    je=orb2alpha(1)
    h1e = cperm*GET_INTEGRAL(ie,je,0,0,0)
!    h1e = cperm*h1eMO(orb1alpha(1),orb2alpha(1))
   elseif(ndiffalpha==0 .and. ndiffbeta==1) then
    ie=orb1beta(1)
    je=orb2beta(1)
    h1e = cperm*GET_INTEGRAL(ie,je,0,0,0)
!    h1e = cperm*h1eMO(orb1beta(1),orb2beta(1))
   else
    write(*,*)"error in applying Slater Condon rules, I stop"
    write(*,*)ndiffalpha,ndiffbeta
    stop
   endif

  do i = 1, det1%nalpha
    do j = 1, ndiffalpha
       if(det1%alpha(i)/=orb1alpha(j)) then
         ie = orb1alpha(j)
         je = orb2alpha(j)
         ke = det1%alpha(i)
         le = det1%alpha(i)
!        inde = intindex(orb1alpha(j),orb2alpha(j),det1%alpha(i),det1%alpha(i))
!        g2e = g2e + cperm*int2e(inde)
         g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)
!        inde = intindex(orb1alpha(j),det1%alpha(i),det1%alpha(i),orb2alpha(j))
!        g2e = g2e - cperm*int2e(inde)
         ie = orb1alpha(j)
         je = det1%alpha(i)
         ke = det1%alpha(i)
         le = orb2alpha(j)
         g2e = g2e - cperm*GET_INTEGRAL(ie,je,ke,le,0)
       endif
    enddo
    do j = 1, ndiffbeta
!        inde = intindex(orb1beta(j),orb2beta(j),det1%alpha(i),det1%alpha(i))
!        g2e = g2e + cperm*int2e(inde)
         ie = orb1beta(j)
         je = orb2beta(j)
         ke = det1%alpha(i)
         le = det1%alpha(i)
         g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)
    enddo
  enddo

  do i = 1, det1%nbeta
    do j = 1, ndiffbeta
       if(det1%beta(i)/=orb1beta(j)) then
!        inde = intindex(orb1beta(j),orb2beta(j),det1%beta(i),det1%beta(i))
!        g2e = g2e + cperm*int2e(inde)
         ie = orb1beta(j)
         je = orb2beta(j)
         ke = det1%beta(i)
         le = det1%beta(i)
         g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)

!        inde = intindex(orb1beta(j),det1%beta(i),det1%beta(i),orb2beta(j))
!        g2e = g2e - cperm*int2e(inde)
         ie = orb1beta(j)
         je = det1%beta(i)
         ke = det1%beta(i)
         le = orb2beta(j)
         g2e = g2e - cperm*GET_INTEGRAL(ie,je,ke,le,0)
       endif
    enddo
    do j = 1, ndiffalpha
!        inde = intindex(orb1alpha(j),orb2alpha(j),det1%beta(i),det1%beta(i))
!        g2e = g2e + cperm*int2e(inde)
         ie = orb1alpha(j)
         je = orb2alpha(j)
         ke = det1%beta(i)
         le = det1%beta(i)
         g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)
    enddo
  enddo

    mat1e = h1e
    mat2e = g2e

elseif(ndiff==2) then 

 if (ndiffalpha == 2) then

!    inde = intindex(orb1alpha(1),orb2alpha(1),orb1alpha(2),orb2alpha(2))
!    g2e = g2e + cperm*int2e(inde)
     ie = orb1alpha(1)
     je = orb2alpha(1)
     ke = orb1alpha(2)
     le = orb2alpha(2)
     g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)

!    inde = intindex(orb1alpha(1),orb2alpha(2),orb1alpha(2),orb2alpha(1))
!    g2e = g2e - cperm*int2e(inde)
     ie = orb1alpha(1)
     je = orb2alpha(2)
     ke = orb1alpha(2)
     le = orb2alpha(1)
     g2e = g2e - cperm*GET_INTEGRAL(ie,je,ke,le,0)

 elseif (ndiffbeta == 2) then
!    inde = intindex(orb1beta(1),orb2beta(1),orb1beta(2),orb2beta(2))
!    g2e = g2e + cperm*int2e(inde)
     ie = orb1beta(1)
     je = orb2beta(1)
     ke = orb1beta(2)
     le = orb2beta(2)
     g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)

!    inde = intindex(orb1beta(1),orb2beta(2),orb1beta(2),orb2beta(1))
!    g2e = g2e - cperm*int2e(inde)
     ie = orb1beta(1)
     je = orb2beta(2)
     ke = orb1beta(2)
     le = orb2beta(1)
     g2e = g2e - cperm*GET_INTEGRAL(ie,je,ke,le,0)
 else
!    inde = intindex(orb1alpha(1),orb2alpha(1),orb1beta(1),orb2beta(1))
!    g2e = g2e + cperm*int2e(inde)
     ie = orb1alpha(1)
     je = orb2alpha(1)
     ke = orb1beta(1)
     le = orb2beta(1)
     g2e = g2e + cperm*GET_INTEGRAL(ie,je,ke,le,0)
 endif

    mat1e = h1e 
    mat2e = g2e

else
 mat1e = 0d0
 mat2e = 0d0
endif

end subroutine mat2dets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matvp2dets(det1,det2,mat1e)
  use phis
  USE GlobalMod, ONLY : get_unit,stdlog
  USE ukrmol_interface_gbl, ONLY : START_MPI, FINALIZE_MPI, GET_INTEGRAL, GET_KINETIC_ENERGY_INTEGRAL, LPQRS
  USE PhisMod, ONLY : GBTO_INIT, NBAS, NOCC
!! returns the Vp matrix element between det1 and det2
!! using the Slater Condon rules
!! spin-integration is taken care here
implicit none

type(Sdeterminant), intent(in) :: det1, det2
double precision :: mat1e, h1e, cperm
integer :: ndiffalpha, ndiffbeta, ndiff
integer, dimension(2) :: orb1alpha, orb2alpha, orb1beta, orb2beta
integer :: ie, je, ke, le, inde, intindex

integer :: i, j

mat1e = 0d0
h1e = 0d0

call compare2dets(det1,det2,cperm,ndiffalpha,orb1alpha,orb2alpha,ndiffbeta,orb1beta,orb2beta)
ndiff = ndiffalpha + ndiffbeta

if(ndiff==0) then
  do i = 1, det1%nalpha
      ie =  det1%alpha(i)
      h1e = h1e + vpMO(ie,ie)
  enddo

  do i = 1, det1%nbeta
      ie =  det1%beta(i)
      h1e = h1e + vpMO(ie,ie)
  enddo
 mat1e = cperm*h1e

elseif(ndiff==1) then

   if(ndiffalpha==1 .and. ndiffbeta==0) then
    h1e = cperm*vpMO(orb1alpha(1),orb2alpha(1))
   elseif(ndiffalpha==0 .and. ndiffbeta==1) then
    h1e = cperm*vpMO(orb1beta(1),orb2beta(1))
   else
    write(*,*)"error in applying Slater Condon rules, I stop"
    write(*,*)ndiffalpha,ndiffbeta
    stop
   endif

    mat1e = h1e

else
 mat1e = 0d0
endif

end subroutine matvp2dets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SlaterDeterminant
