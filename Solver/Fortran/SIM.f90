! Author: Abdul Hannan Faruqi
! Last modified: Apr 04, 2020 0827

module data_mod
    implicit none
    SAVE
    complex*16 :: alpha = (1.0,0.0), beta = (0.0,0.0)
    integer:: n, m=15
    complex*16::shift_par=0.0
end module

Program Eigenvalue
USE data_mod
implicit none
integer::i, j, K, ok, row, col, iter, ierr, ierr2, sh_switch= 0, d1, d2, d3
character(len=50):: e_str, filename, e_str2
complex*16, dimension(:,:), allocatable:: A, B, Th, evec, Ax, Bx, C, Phi, Int_arr
complex*16, dimension(:), allocatable:: eval, WORK1, work2, eigs, tau
double precision, dimension(:), allocatable:: w, w2
integer, dimension(:), allocatable::IPIV, IPIV2
integer, dimension(:,:), allocatable::eye
complex*16 DUMMY(1,1)
real(8) :: val, tol= 0.000001, e = 1.0
double precision:: Th_norm, zlange
real::sqrt, start_time, stop_time
character(len=6):: num

filename = 'Mn.dat'
open(unit=10000, file=filename, status = 'old', action = 'read', iostat = ierr, iomsg = e_str)
if (ierr.eq.0) then
    print *, "File opened successfully"
else
    print *, "Error!"
endif
read(10000, *) n
print *, n

allocate(A(n,n))
allocate(B(n,n))
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
allocate(IPIV2(n))
allocate(evec(m,m))
allocate(Phi(n,m))
allocate(Int_arr(m,n))
allocate(eye(n,n))

eigs=0.0
eye = 0.0
do i=1,n
    eye(i,i)=1
end do

do
    read(10000,*, iostat = ierr) row, col, val
    if (ierr.ne.0)    exit
    A(row+1,col+1) = val
enddo
close(10000)
print *, "File closed"

open(unit=20000, file= "Eigs.dat", status = 'new', action = 'write', iostat = ierr2, iomsg = e_str2)

do i=1,n
    do j=1, m
        if (i.eq.j) then
            Th(i,j) = 1
        else
            Th(i,j)= 0
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
!   Pre-multiplication
    call zgesv(n, m, A, n, ipiv2, Th, n, ok)

!   Normalize Th
    Th_norm = zlange("F", n, m, Th, n, w)
    Th = Th/Th_norm

!   QR decomposition
!    if(mod(iter,5).eq.0) then
	call zgeqrf(n, m, Th, n, tau, work1, m*m, ok)
    	call zungqr(n, m, m, Th, n, tau, work1, m*m, ok)
!    endif

!   Project A and B onto the subspace spanned by Th
!    call Project(A, Th, Ax)   | Replace subroutine to prevent 
!    call Project(B, Th, Bx)   | Stack Overflow
    call zgemm('T', 'N', m, n, n, alpha, Th, n, A, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Ax, m)
    call zgemm('T', 'N', m, n, n, alpha, Th, n, B, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Bx, m)

!   Calculating inv(Bx)
    call zgetrf(m, m, Bx, m, IPIV, ok)
    call zgetri(m, Bx, m, IPIV, work1, m*m, ok)

!   C = inv(Bx)*Ax
    call zgemm('N', 'N', m, m, m, alpha, Bx, m, Ax, m, beta, C, m)

!   Eigenvalues and right eigenvectors of C
    call ZGEEV('N', 'V', m, C, m, eval, DUMMY, 1, evec, m, WORK2, 2*m, W2, ok)

!   Calculate new Th
    call zgemm('N', 'N', n, m, m, alpha, Th, n, evec, m, beta, Phi, n)
    call zgemm('N', 'N', n, m, n, alpha, B, n, Phi, n, beta, Th, n)

    e = abs(maxval(real(eval)-real(eigs)))
    eigs = eval
    if ((e.lt.0.0001).and.(sh_switch.eq.0)) then
        shift_par = minval(real(eigs))
        A = A - shift_par*eye
        sh_switch = 1
    end if

    if (mod(iter,10).eq.0)  then
        filename="Iterations\Iter000.dat"
        d3=iter/100
        d2=(iter-100*d3)/10
        d1=iter-d2*10-d3*100
        filename(16:16)=char(d3+48)
        filename(17:17)=char(d2+48)
        filename(18:18)=char(d1+48)
        open(unit=iter, file=filename, status = 'new', action = 'write')
        write(iter, *) "Error =", e
        write(iter, *) "Eigenvalues = " 
        do i=1,m
            write(iter, *) 1/(eigs(i)+shift_par)
	enddo
        close(iter)
    endif

    if ((e.lt.tol)) then
        exit
    end if
    print *, "Iterations = ", iter
    print *, "error = ", e
end do
call cpu_time(stop_time)
   
print *, "Execution time= ", stop_time - start_time, "s"

print *, "Iterations = ", iter
print *, "error = ", e
print *, "Shift = ", shift_par
eigs = eigs + shift_par
print *, "Eigenvalues = ", 1/eigs

do i=1,m
    write(20000, *) 1/eigs(i)
enddo
close(20000)
end program

subroutine Project(Mat, Th, Mat_x)
    USE data_mod
    implicit none
    integer::i
    complex*16 Mat(n,n), Th(n,m), Mat_x(m,m), Int_arr(m,n)

    call zgemm('T', 'N', m, n, n, alpha, Th, n, Mat, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Mat_x, m)
end subroutine