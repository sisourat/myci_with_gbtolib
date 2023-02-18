MODULE PhisMod

  USE precisn_gbl, ONLY : wp
  USE ukrmol_interface_gbl, ONLY : READ_PHIS_INTS, GET_ENERSYMOCC
  USE GlobalMod, ONLY : get_unit

  !# DEFINE ARRAY OF POINT GROUP NAMES AND FOR EACH GROUP NAMES OF SYMMETRIES
  !(8x8) ARRAY OF STRINGS)

  IMPLICIT NONE

  ! SCF
  INTEGER(4) pgroup,nirep
  INTEGER(4) nbas,nocc
  INTEGER(4),DIMENSION(8,8) :: mt

  !> MO basis data 
  integer(4),allocatable,public :: sym(:)
  real(wp),allocatable,public :: epsilon(:),occ(:)
  logical,allocatable,public :: iscont(:),iscont_ext(:)
  logical,public :: no_cont = .false., &
                    zeroFock = .false., &
                    zeroFock_occ = .false., &
                    override_scatci_class = .true.
 
  INTEGER irep

  INTEGER nholes,nparticles,ncorep,ncontp
  INTEGER(2),ALLOCATABLE :: holes(:),particles(:)
  INTEGER(2),ALLOCATABLE :: corep(:),contp(:)

  REAL(8) :: max_peps = 1.d9

  REAL(8) :: GSE_correction = 0.d0

  CHARACTER(21), PARAMETER :: fl_contclass = "moints.classification"

CONTAINS

  SUBROUTINE GBTO_INIT(printflag)
    LOGICAL,INTENT(in) :: printflag
   !*******************************
    integer i,j,ind,mi,mj,sm
    real(8) max_offd
    logical isthere
    !GBTO parameters and auxiliary variables
    INTEGER :: nfti = -1
    INTEGER :: nfta = 6
    INTEGER :: nint1e = 0, &
               nint2e = 0

100  FORMAT (A5,T9,A3,T17,A3,T28,A,T44,A,1x,A)
110  FORMAT (I4,T8,I3,T16,F4.1,T24,F18.10,2x,L3,2x,L3)
120  FORMAT (3x,I1,1x,I1,1x,I1,1x,I1,1x,I1,1x,I1,1x,I1,1x,I1,1x,I1)

   !Executables

      call READ_PHIS_INTS(nfti,nfta,nint1e,nint2e,nirep)

      call GET_ENERSYMOCC(pgroup,nirep,nbas,nocc, &
        epsilon,sym,occ,iscont,mt(1:8,1:8))

      allocate(iscont_ext(nbas))

      iscont_ext(1:nbas) = iscont(1:nbas)
      inquire(file=fl_contclass,exist=isthere)
      if (isthere) then
        call get_unit(ind)
        open(ind,file=fl_contclass,status='old')
        do i = 1,nbas
          read(ind,*) j,sm,mi,mj
          if (sm.ne.sym(j)) &
            write(6,'(/,"WARNING, possible moints.canonical/", &
              A21," mismatch")') fl_contclass
          if (mj.eq.0) then
            iscont_ext(i) = .false.
          else if (mj.eq.1) then
            iscont_ext(i) = .true.
          end if
        end do
        if (override_scatci_class) iscont = iscont_ext
      end if

      if (no_cont) then
        iscont = .false.
        iscont_ext = .false.
      end if

      if (printflag) then
        write(6,*)
        write(6,'(A,I1)') 'POINT GROUP: ', pgroup
        write(6,'(A,I1)') 'NO. OF IRREPS: ', nirep
        write(6,'(A,I4)') 'NO. OF BASIS FUNCTIONS: ', nbas
        write(6,'(A,I4)') 'NO. OF OCCUPIED MOs: ', nocc
      end if

      if (printflag) then
        write(6,*)
        write(6,*) "Symmetry table:"
        write(6,*) "..............."
        write(6,120) 0,(i,i=1,nirep)
        do i = 1,nirep
          write(6,120) i,(mt(i,j), j=1,nirep)
        end do

        write(6,*)
        write(6,*) '                 RHF data:'
        write(6,'(50("-"))')
        write(6,100) 'NO.', 'SYM', 'OCC', 'ORBITAL ENERGY', 'CORE', 'CONT'
        write(6,'(50("-"))')
        do i = 1,nbas
          write(6,110) i,sym(i),occ(i),epsilon(i),.not.(iscont(i)),iscont_ext(i)
        end do
      end if

  END SUBROUTINE GBTO_INIT

END MODULE PhisMod
