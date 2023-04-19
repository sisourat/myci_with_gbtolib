MODULE GlobalMod

  IMPLICIT NONE

  INTEGER,PARAMETER :: stderr = 0, &
                       stdout = 6
  INTEGER           :: stdlog

  INTEGER :: buffer_size = 1048576

  REAL(8),PARAMETER :: infty = huge(1.d10)

  ! Holes on the initially ionized atom
  INTEGER,ALLOCATABLE :: hole_in(:)
  LOGICAL,ALLOCATABLE :: confC1IN(:)
  LOGICAL use_hole_in

  ! Initial vector
  INTEGER,DIMENSION(2) :: i2h
  INTEGER ih
  CHARACTER*10 infile

  ! Pre-diag main subspace for state selection
  LOGICAL :: diagC1 = .false.
  REAL(8),ALLOCATABLE :: vectorC1(:)

  ! Re-run options
  LOGICAL dump_ham
  INTEGER saveexp
  REAL(8),ALLOCATABLE :: mat_init(:,:)

  ! Lanczos
  LOGICAL lancfin,lancfull,FDMfin,FDMinit
  INTEGER :: max_full = 70000

  CHARACTER(9) :: phifile="Hivec.dat"
  CHARACTER(12) :: phifile_short="QHQspect.dat"
  CHARACTER(10) :: Pphifile="PHivec.dat"

  CHARACTER(11) :: fn_sym = "syminfo.dat"
  CHARACTER(9):: fn_diag = 'final.dia'
  CHARACTER(10) :: fn_offd = 'final.offc'
  CHARACTER(12) :: fn_C12 = "confC12.dump"
  CHARACTER(7) :: fn_stiel="LANCTRM"
!  CHARACTER(8) :: fn_pstiel="LANCTRMp"
  CHARACTER(7) :: fn_tmat = "TMATEVC"
  CHARACTER(13) :: fn_pkrvbase = "pkrvbase.dat"
  CHARACTER(12) :: fn_SIefin = "efin_inp.dat"
  CHARACTER(13) :: fn_SIinput = "stieltjes.inp"
  CHARACTER(6) :: fn_rband0 = "RBAND0"
  CHARACTER(5) :: fn_rband1 = "RBAND"
  CHARACTER(5) :: fn_rband2 = "pkrvb"
  CHARACTER(8) :: fn_lancdump = "LANCDUMP"
  CHARACTER(11) :: fn_chanvec = "chanvec.dat"
  CHARACTER(12) :: fn_chanconf = "chanconf.dat"

  CHARACTER(11) :: fn_chanE = "chan_en.dat"

  !Matrix-vector
  LOGICAL :: hamfly = .false.

  INTEGER :: maxMV = 500
    ! maximum number of vectors to be multiplied at the same time
    ! (for OPENMP version)

CONTAINS

  SUBROUTINE get_unit(iounit)
    INTEGER,INTENT(out) :: iounit
    !*****************************************
    !* Returns number of an available I/O unit
    !*****************************************
    logical unitok,unitop
    integer i
    !Executables
    iounit = -1
    i = 10
    do while (iounit.eq.-1)
      inquire(unit=i,exist=unitok,opened=unitop)
      if (unitok.and.(.not.unitop)) then
        iounit = i
      else
        i = i+1
      end if
    end do
  END SUBROUTINE get_unit

  INTEGER FUNCTION min2(i,j)
    INTEGER,INTENT(in) :: i,j
    if (i.le.j) then
      min2=i
    else
      min2=j
    end if
  END FUNCTION min2 

  INTEGER FUNCTION max2(i,j)
    INTEGER,INTENT(in) :: i,j
    if (i.ge.j) then
      max2=i
    else
      max2=j
    end if
  END FUNCTION max2

  LOGICAL FUNCTION InArray(m,arr)
    INTEGER(2),INTENT(in) :: m
    INTEGER(2),INTENT(in) :: arr(:)
    !*******************************
    !* Checks whether m is in arr(:)
    !*******************************
    integer i
    !Executables
    InArray = .false.
    do i = 1,size(arr)
      if (m == arr(i)) then
        InArray = .true.
        return
      end if
    end do
  END FUNCTION InArray

END MODULE GlobalMod
