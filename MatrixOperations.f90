Module MatrixOperations
use GlobalVariables


contains



function inv(A) result(Ainv)
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.

    ! Source: Fortran Wiki
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
        stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
        stop 'Matrix inversion failed!'
    end if
end function inv

subroutine MatrixMultiply(A,B,C)
    implicit none

    real(dp), dimension(:,:), intent(in)            :: A, B
    real(dp), dimension(:,:), intent(out)           :: C
    integer                                         :: i,j, m, k, n
    real(dp)                                        :: alpha, beta

    m = size(A(:,1))
    n = size(B(1,:))
    k = size(A(1,:))
    alpha = 1.0_dp
    beta = 0.0_dp
   
    
    call DGEMM('N','N',m,n,k,alpha,A,M,B,k,beta,C,M)
    
end subroutine MatrixMultiply









end Module