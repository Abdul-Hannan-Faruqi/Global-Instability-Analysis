! Author: Abdul Hannan Faruqi
! Last modified: Apr 01, 2020 20:22

module data
    implicit none
    SAVE
    complex*16 :: alpha = 1.0, beta = 0.0
    integer:: n=3, m=2
end module

Program Subspace Iteration

USE data
implicit none
integer::i, j, K, ok, row, col, iter, ierr
character(len=50):: e_str, filename
complex*16, dimension(:,:), allocatable:: A, B, Th, evec, R, Ax, Bx, C, Phi
complex*16, dimension(:), allocatable:: eval, WORK1, work2, eigs, tau
double precision, dimension(:), allocatable:: w, w2
integer, dimension(:), allocatable::IPIV
complex*16 DUMMY(1,1)
real :: v, tol= 0.00000001, e = 1.0
real:: zlange, Th_norm

!filename = 'M.dat'
!open(unit=100, file=filename, status = 'old', action = 'read', iostat = ierr, iomsg = e_str)
!print *, e_str
!read(100, *) x


allocate(A(n,n))
allocate(B(n,n))
allocate(R(n,m))
allocate(Th(n,m))
allocate(Ax(m,m))
allocate(Bx(m,m))
allocate(C(m,m))
allocate(w(n))
allocate(w2(2*m))
allocate(work1(m*m))
allocate(work2(2*m))
allocate(eval(m))
allocate(eigs(m))
allocate(tau(m))
allocate(IPIV(m))
allocate(evec(m,m))
allocate(Phi(n,m))

eigs=0.0

do i=1,n
    do j=1, m
        if (i.eq.j) then
            R(i,j) = 1
        else
            R(i,j)= 0
        end if
    end do
end do

do i=1,n
    do j=1,i
        A(i,j) = AINT(100*rand(0))
        A(j,i) = A(i,j)
    end do
end do

print *, "A = "
do i=1,n
    print *, A(i,:)
end do

do i=1,n
    do j = 1,n
        if ((i.eq.j).and.(mod(i,3).ne.1)) then
            B(i,j) = (1)
        else
            B(i,j) = (0)
        endif
    end do
end do

iter = 0
do
    iter = iter+1
!   Multiplication
    call zgemm('N', 'N', n, m, n, alpha, A, n, R, n, beta, Th, n)
    Th_norm = zlange("F", n, m, Th, n, w)
    Th = Th/Th_norm
!   QR decomposition
    call zgeqrf(n, m, Th, n, tau, work1, m*m, ok)
    call zungqr(n, m, m, Th, n, tau, work1, m*m, ok)
!   Project A and B onto the subspace spanned by Th
    call Project(A, Th, Ax)
    call Project(B, Th, Bx)
!   Calculating inv(Bx)
    call zgetrf(m, m, Bx, m, IPIV, ok)
    call zgetri(m, Bx, m, IPIV, work1, m*m, ok)
!   C = inv(Bx)*Ax
    call zgemm('N', 'N', m, m, m, alpha, Bx, m, Ax, m, beta, C, m)
!   Eigenvalues and right eigenvectors of C
    call ZGEEV('N', 'V', m, C, m, eval, DUMMY, 1, evec, m, WORK2, 2*m, W2, ok)
!   Calculate new R
    call zgemm('N', 'N', n, m, m, alpha, Th, n, evec, m, beta, Phi, n)
    call zgemm('N', 'N', n, m, n, alpha, B, n, Phi, n, beta, R, n)
    e = abs(maxval(real(eval)-real(eigs)))
    eigs = eval
    if ((e.lt.tol)) then
        exit
    end if
end do
print *, "error = ", e
print *, "Iterations = ", iter
print *, "Eigenvalues = ", eval
print *, "Eigenvectors = "
do i=1,n
    print *, R(i,:)
end do
!
!close(100)
end program

subroutine Project(Mat, Th, Ax)
    USE data
    implicit none
    complex*16 Mat(n,n), Th(n,m), Ax(m,m), Int_arr(m,n)

    call zgemm('T', 'N', m, n, n, alpha, Th, n, Mat, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, n, n, alpha, Int_arr, m, Th, n, beta, Ax, m)
end subroutine
