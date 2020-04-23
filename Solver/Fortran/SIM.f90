! Author: Abdul Hannan Faruqi
! Last modified: Apr 04, 2020 0827

module data_mod
    implicit none
    SAVE
    complex*16 :: alpha = (1.0,0.0), beta = (0.0,0.0)
    integer:: n, m=10
    complex*16::shift_par=0.0
    type :: Eigenpair
        complex*16 :: eval
        complex*16, dimension(:), allocatable:: evec
    end type Eigenpair 
end module

Program Eigenvalue
USE data_mod
implicit none

integer::i, j, K, ok, row, col, iter, ierr, ierr2, sh_switch= 0, d1, d2, d3, d4
character(len=200):: e_str, filename, e_str2
complex*16, dimension(:,:), allocatable:: A, A_copy, B, Th, evec, Ax, Bx, C, Phi, Int_arr, R, AR, BR
complex*16, dimension(:), allocatable:: eval, WORK1, lwork, eigs, tau, eval_num, eval_den
double precision, dimension(:), allocatable:: w, rwork
integer, dimension(:), allocatable::IPIV, IPIV2
integer, dimension(:,:), allocatable::eye
complex*16 DUMMY(1,1)
real(8) :: val, tol= 1e-6, e = 1.0, err
double precision:: Th_norm, zlange
real::sqrt, start_time, stop_time, rand


filename = 'FJ_matrix.csv'
open(unit=10000, file=filename, status = 'old', action = 'read', iostat = ierr, iomsg = e_str)
if (ierr.eq.0) then
    print *, "File opened successfully"
else
    print *, "Error!"
endif
read(10000, *) n
print *, n

allocate(A(n,n))
allocate(A_copy(n,n))
allocate(B(n,n))
allocate(Th(n,m))
allocate(AR(n,m))
allocate(BR(n,m))
allocate(R(n,m))
allocate(Ax(m,m))
allocate(Bx(m,m))
allocate(C(m,m))
allocate(w(n))
allocate(rwork(8*m))
allocate(work1(m*m))
allocate(lwork(2*m))
allocate(eval(m))
allocate(eval_num(m))
allocate(eval_den(m))
allocate(evec(m,m))
allocate(eigs(m))
allocate(tau(m))
allocate(IPIV(m))
allocate(IPIV2(n))
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
open(unit=50000, file= "Error.dat", status = 'new', action = 'write')
 

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
!   Pre-multiplication
    A_copy = A
    Th = R
    call zgesv(n, m, A_copy, n, ipiv2, Th, n, ok)

!   Normalize Th
    Th_norm = zlange("F", n, m, Th, n, w)
    Th = Th/Th_norm 

!   QR decomposition
!    if(mod(iter,5).eq.0) then
	call zgeqrf(n, m, Th, n, tau, work1, m*m, ok)
   	call zungqr(n, m, m, Th, n, tau, work1, m*m, ok)
!       call zunmqr('R', 'N', n, n, m, Th, n, tau, eye, n, w, n, ok)
!    endif

!   Project A and B onto the subspace spanned by Th
!    call Project(A, Th, Ax)   | Replace subroutine to prevent 
!    call Project(B, Th, Bx)   | Stack Overflow
    call zgemm('T', 'N', m, n, n, alpha, Th, n, A, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Ax, m)
    call zgemm('T', 'N', m, n, n, alpha, Th, n, B, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Bx, m)

!   Calculating inv(Bx)
!    call zgetrf(m, m, Bx, m, IPIV, ok)
!    call zgetri(m, Bx, m, IPIV, work1, m*m, ok)
    
!   C = inv(Bx)*Ax
!    call zgemm('N', 'N', m, m, m, alpha, Bx, m, Ax, m, beta, C, m)

!   Eigenvalues and right eigenvectors of C
!    call ZGEEV('N', 'V', m, C, m, eval, DUMMY, 1, evec, m, WORK2, 2*m, W2, ok)

    call zggev('N', 'V', m, Ax, m, Bx, m, eval_num, eval_den, DUMMY, 1, evec, m, lwork, 2*m, rwork, ok) 

!   Calculate new Th
    call zgemm('N', 'N', n, m, m, alpha, Th, n, evec, m, beta, Phi, n)
    call zgemm('N', 'N', n, m, n, alpha, B, n, Phi, n, beta, R, n)

    do i=1,m
        if (real(eval_den(i)).ne.0)    then
            eval(i) = eval_num(i)/eval_den(i)
        else 
            eval(i)= 0
        endif
    enddo
    if (iter.gt.9)    call sort(eval, evec)

    e = abs(maxval(abs(eval-eigs)))

    eigs = eval

!    if ((e.lt.1e-3).and.(sh_switch.eq.0)) then
!        shift_par = minval(real(eigs))
!        A = A - shift_par*eye
!        sh_switch = 1
!    end if

    if (mod(iter,10).eq.0)  then
        filename="Iterations_40x40\Iter0000.dat"
        d4=iter/1000
        d3=(iter-1000*d4)/100
        d2=(iter-1000*d4-100*d3)/10
        d1=iter-d2*10-d3*100-d4*1000
        filename(22:22)=char(d4+48)
        filename(23:23)=char(d3+48)
        filename(24:24)=char(d2+48)
        filename(25:25)=char(d1+48)
        open(unit=iter, file=filename, status = 'new', action = 'write')
        write(iter, *) "Shift =", shift_par
        write(iter, *) "Eigenvalues = " 
        do i=1,m
            write(iter, *) eigs(i)+shift_par
	enddo
        close(iter)
    endif
       
    write(50000, *) e

    if ((e.lt.tol).or.(iter.eq.100)) then
        exit
    end if
    print *, "Iteration = ", iter
    print *, "error = ", e
end do
close(50000)
call cpu_time(stop_time)
  
print *, "Execution time= ", stop_time - start_time, "s"
print *, "Iterations = ", iter
print *, "error = ", e
print *, "Shift = ", shift_par
eigs = eigs + shift_par
print *, "Eigenvalues = "

do i=1,m
    if (real(eval_den(i)).ne.0)    then
        write(20000, *) i, eval_num(i)/eval_den(i)
        write(*, *) i, eval_num(i)/eval_den(i)
    else   
        write(20000, *) i, eval_num(i), eval_den(i)
        write(*, *) i, eval_num(i), eval_den(i)
    endif
enddo
close(20000)
call zgemm('N', 'N', n, m, n, alpha, A, n, R, n, beta, AR, n)

end program

subroutine sort(eval, evec)
    USE data_mod
    implicit none
    integer::i, j, ptr
    complex*16::mine
    complex*16 eval(m)
    complex*16 evec(m,m)
    type (Eigenpair), dimension(:), allocatable :: Eigp
    allocate(Eigp(m))
    
    do i=1,m
        allocate(Eigp(i)%evec(m))
        ptr=i
        mine= eval(ptr)
        do j=i+1,m
            if (real(eval(j)).lt.real(mine))    then
                mine = eval(j)       
                ptr = j
            endif
        enddo
        eval(ptr) = eval(i)
        Eigp(i)%eval = mine
        Eigp(i)%evec(:) = evec(:,ptr)
    enddo
    eval = Eigp(:)%eval
  
    do i=1,m
        evec(:,i) = Eigp(i)%evec
    enddo

end subroutine sort
    

subroutine Project(Mat, Th, Mat_x)
    USE data_mod
    implicit none
    integer::i
    complex*16 Mat(n,n), Th(n,m), Mat_x(m,m), Int_arr(m,n)

    call zgemm('T', 'N', m, n, n, alpha, Th, n, Mat, n, beta, Int_arr, m)
    call zgemm('N', 'N', m, m, n, alpha, Int_arr, m, Th, n, beta, Mat_x, m)
end subroutine