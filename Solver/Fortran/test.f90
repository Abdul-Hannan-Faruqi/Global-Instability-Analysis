Program Eigenvalue

implicit none

integer:: m
integer::i, j, K, ok, row, col, ierr, ierr2, x
character(len=50):: e_str, filename, e_str2
complex*16, dimension(:,:), allocatable:: evec, B, A
complex*16, dimension(:), allocatable:: eval_num, eval_den, lwork
double precision, dimension(:), allocatable:: rwork
complex*16 DUMMY(1,1)
real(8) :: val
integer, dimension(:), allocatable::IPIV
real::sqrt, start_time, stop_time

filename = 'FJ_matrix.csv'
open(unit=100, file=filename, status = 'old', action = 'read', iostat = ierr, iomsg = e_str)
if (ierr.eq.0) then
    print *, "File opened successfully"
else
    print *, "Error!"
endif
read(100, *) m
print *, m

open(unit=200, file= "Eigs_test.dat", status = 'new', action = 'write', iostat = ierr2, iomsg = e_str2)
print *, "ierr2 = ", ierr2

allocate(B(m,m))
allocate(A(m,m))
allocate(rwork(8*m))
allocate(lwork(2*m))
allocate(eval_num(m))
allocate(eval_den(m))
allocate(evec(m,m))

do
    read(100,*, iostat = ierr, iomsg=e_str) row, col, val
    if (ierr.ne.0)    exit
    A(row+1,col+1) = val
enddo
close(100)
print *, "File closed"

do i=1,m
    do j = 1,m
        if ((i.eq.j).and.(mod(i,3).ne.1)) then
            B(i,j) = (1)
        else
            B(i,j) = (0)
        endif
    end do
end do

print *, "B formed"

call cpu_time(start_time)
call zggev('N', 'V', m, A, m, B, m, eval_num, eval_den, DUMMY, 1, evec, m, lwork, 2*m, rwork, ok) 
call cpu_time(stop_time)
   
print *, "Execution time= ", stop_time - start_time, "s"

do i=1,m
    if (real(eval_den(i)).ne.0)    then
        write(200, *) i, eval_num(i)/eval_den(i)
    else
        write(200,*) i, eval_num(i), eval_den(i)
    endif
enddo

end program