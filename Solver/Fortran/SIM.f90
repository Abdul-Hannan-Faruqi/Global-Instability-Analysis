! Author: Abdul Hannan Faruqi
! Last modified: Apr 04, 2020 0827

module data_mod
    implicit none
    SAVE
    complex*16 :: alpha = (1.0,0.0), beta = (0.0,0.0)
    integer:: n, m=20
    complex*16::shift_par
end module

Program Eigenvalue
USE data_mod
implicit none
integer::i, j, K, ok, row, col, iter, ierr, x
character(len=50):: e_str, filename
complex*16, dimension(:,:), allocatable:: A, B, Th, evec, R, Ax, Bx, C, Phi
complex*16, dimension(:), allocatable:: eval, WORK1, work2, eigs, tau
double precision, dimension(:), allocatable:: w, w2
integer, dimension(:), allocatable::IPIV
integer, dimension(:,:), allocatable::eye
complex*16 DUMMY(1,1)
real(8) :: val, tol= 0.00000001, e = 1.0
double precision:: Th_norm, zlange
real::sqrt, start_time, stop_time

filename = 'M.dat'
open(unit=100, file=filename, status = 'old', action = 'read', iostat = ierr, iomsg = e_str)
if (ierr.eq.0) then
    print *, "File opened successfully"
else
    print *, "Error!"
endif
read(100, *) n
print *, n

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
allocate(eye(n,n))

eigs=0.0
eye = 0.0
do i=1,n
    eye(i,i)=1
end do

do
    read(100,*, iostat = ierr) row, col, val
    if (ierr.ne.0)    exit
    A(row+1,col+1) = val
enddo
close(100)
print *, "File closed"

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
    do j = 1,n
        if ((i.eq.j).and.(mod(i,3).ne.1)) then
            B(i,j) = (1)
        else
            B(i,j) = (0)
        endif
    end do
end do

iter = 0
call cpu_time(start_time)
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
    if (iter.eq.50) then
        shift_par = minval(real(eigs))
        A = A - shift_par*eye
    end if
    if ((e.lt.tol)) then
        exit
    end if
end do
call cpu_time(stop_time)
   
print *, "Execution time= ", stop_time - start_time, "s"

print *, "Iterations = ", iter
print *, "error = ", e
print *, "Shift = ", shift_par
eigs = eigs + shift_par
print *, "Eigenvalues = ", eigs

end program

subroutine Project(Mat, Th, Mat_x)
    USE data_mod
    implicit none
    integer::i
    complex*16 Mat(n,n), Th(n,m), Mat_x(m,m), Int_arr(m,n)

    call zgemm('T', 'N', m, n, n, alpha, Th, n, Mat, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Mat_x, m)
end subroutine