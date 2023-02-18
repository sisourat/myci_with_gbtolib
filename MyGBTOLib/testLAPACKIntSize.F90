program mkl_interface
    use iso_fortran_env, only: INT, real32
    real(real32) :: A(3,3) = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), (/ 3, 3 /))
    real(real32) :: M(3,3)
    integer(INT) :: P(3) = (/ 3, 1, 2 /), Q(2) = (/ 3, 3 /)
    integer(INT) :: one = 1, two = 2, three = 3
    M(:,:) = A(:,:)
    call slaswp(three, M, three, one, two, Q, one)
    write (*, '(9I0)') nint(M(:,:) - A(P,:))
end program mkl_interface
