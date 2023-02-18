PROGRAM main

  USE GlobalMod, ONLY : get_unit,stdlog
  USE ukrmol_interface_gbl, ONLY : START_MPI, FINALIZE_MPI, GET_INTEGRAL, GET_KINETIC_ENERGY_INTEGRAL
  USE PhisMod, ONLY : GBTO_INIT

  IMPLICIT NONE

! Other
  INTEGER i,j,ivec,iounit,iou,ninit,nvec,nchan,ind
  REAL(8) time1,time2,iover,gtotal,E0,E1,E2
  LOGICAL isthere
  CHARACTER(4) gprefix
  CHARACTER(11) flnamep


  write(6,'(/,"Starting MPI -- dummy in this version")')
  call START_MPI

  call GBTO_INIT(.true.)
  write(6,'(/,"GBTOlib initialized.")')

  write(*,*)GET_INTEGRAL(1,1,0,0,0),GET_KINETIC_ENERGY_INTEGRAL(1,1)
!  write(*,*)GET_INTEGRAL(1,1,1,1,0)

  call FINALIZE_MPI

END PROGRAM main
