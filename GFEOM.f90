
program GFEOM

use matrixmod

!definition of the precsion
implicit none
integer, parameter :: i1=selected_int_kind(1)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)


!Physical system varibales
real(kind=r8)::mu,U,beta 
!Kind system physical variables ... 
real(kind=r8)::Ku(2,2),Kc(2,2),R(2,2),av(2)
real(kind=r8)::M_mat(2,2),v(2),Q(2,2)
real(kind=r8)::E(2),e_r(2),e_i(2),R_temp(2,2) 
real(kind=r8)::FR(2,2),L_sys(2,2),L1(2,2),L2(2,2),L3(2,2)     
!Dummy variable 
integer(kind=i4)::i,j


!------------------------------------------
!
! PHYSICAL CONSTANT
!
!-----------------------------------------
U=10
mu=1
beta=500 


!K_uncorr
Ku(1,1)=0

Ku(1,2)=0

Ku(2,1)=0

Ku(2,2)=U   
!K_correlate

Kc(1,1)=0

Kc(1,2)=U

Kc(2,1)=0

Kc(2,2)=0

!M_matrix 

M_mat(1,1)=0

M_mat(1,2)=0

M_mat(2,1)=1

M_mat(2,2)=0 

!v know vector
v(1)=1
v(2)=0
!trial to check 
E(1)=0
E(2)=U
!Let us diagonalize the matrix Ku and definition of R  
call eigengenleft(Ku,e_r,e_i,R_temp,2)   
R=transpose(R_temp) 
write(*,*) "The Spectrum of the Theory" 
do i=1,2
	write(*,*) i,e_r(i),e_i(i) 
enddo
E=e_r
write(*,*) "the lefteigenvector store by column" 
call mprint(R_temp,2)
write(*,*) "the matrix Q"  
!Let us generate the ( Q=f(ei)-f(e2) )  / (e1-e2) and we should do it in a subroutine. 
call generate_Q(Q,Kc,E,beta,mu,2)
call generate_FR(FR,R,E,beta,mu,2) 
call mprint(Q,2)
write(*,*) "the matrix FR"
call mprint(FR,2) 
L1=R 
L2=matmul(FR,M_mat) 
L3=matmul(Q,matmul(R,M_mat)) 
L_sys=L1-L2-L3
write(*,*) "check the matmul if it works with vector" 
call mprint(M_mat,2) 
write(*,*) "chekc the vector matrix multiplication"
call matvect(M_mat,v,e_i,2) 
do i=1,2
	write(*,*) i,e_i(i) 
enddo










end program 
