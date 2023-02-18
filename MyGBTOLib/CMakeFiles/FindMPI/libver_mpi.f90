      program mpi_ver
      use mpi_f08
      implicit none
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpilibver_str
      integer(kind=MPI_INTEGER_KIND) :: ierror, reslen
      call MPI_GET_LIBRARY_VERSION(mpilibver_str, reslen, ierror)
      print *, mpilibver_str
      end program mpi_ver
