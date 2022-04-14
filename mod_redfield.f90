Module mod_redfield
!! JPCA_117_6196
implicit none
integer, parameter :: nquant=2
real*8, parameter :: gama=0.2d0
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11
complex*16,parameter :: iota = (0.d0,1.d0)

complex*16 g_p(nquant,nquant,nquant,nquant),g_p_dot(nquant,nquant,nquant,nquant),g_p_dot2(nquant,nquant,nquant,nquant)
complex*16 R_diss1(nquant,nquant),R_coh(nquant,nquant),R_diss_t(nquant,nquant)
real*8 E_excit(nquant),pi,wave_to_j,U_exc(nquant,nquant),R_diss(nquant,nquant)
real*8 tot_tim,dt 
integer n_t,n_p
!! Common Variables
real*8 omg_c,s01,s02,V_coup,V_exothermicity,V_barrier,gamma,V_reorg,lambda
real*8 temperature,mass,x_cr,beta
complex*16 :: rho(nquant,nquant),rho1(nquant,nquant)
contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  real*8 c_0,c_e

  pi=dacos(-1.d0)

  wave_to_j=2*pi*clight*hbar

  open(10,file="redfield.inp")
  read(10,*) mass
  read(10,*) omg_c
  read(10,*) V_coup
  read(10,*) lambda
  read(10,*) V_exothermicity
  read(10,*) temperature
  read(10,*) tot_tim
  read(10,*) dt
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif

  !---------------------------------------------------------- 
E_excit(1)=-707.106781186548
E_excit(2)=707.106781186548

E_excit=E_excit*wave_to_j
!write(*,*) E_excit
!stop
  mass=mass*au2kg
  omg_c=omg_c*(2*pi*clight)
!  write(*,*) 1/omg_c,hbar*omg_c
!  stop
  V_exothermicity=V_exothermicity*wave_to_j
  V_coup=V_coup*wave_to_j
!  V_reorg=V_reorg*wave_to_j

   U_exc(1,1)=0.382683432365090
   U_exc(1,2)=-0.92387953251128
   U_exc(2,1)=-0.92387953251128
   U_exc(2,2)=-0.38268343236509


n_t=nint(tot_tim/dt)

rho=0d0
rho(1,1)=1d0

rho=matmul(U_exc,matmul(rho,transpose(U_exc)))
!write(*,*) rho
!stop
end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i
  real*8 log_gama
  real*8 time(n_t)
  real*8 k_FGR,k_TST,k_M,k_red
  complex*16,dimension(nquant,nquant):: k1

  open(100,file="red_rho_ex.out")
  open(101,file="red_rho.out")

!  beta=1.d0/(kb*temperature)
  beta=1.d0/(hbar*omg_c)

R_diss=0d0  
do i=1,n_t

time(i)=(i-1)*dt

 write(100,'(21f15.7)') time(i)*omg_c,real(rho(:,:))
rho1=(matmul(transpose(U_exc),matmul(rho,(U_exc))))
 write(101,'(21f15.7)') time(i)*omg_c,real(rho1(:,:))

 call compute_R_matrix(dt,i)

R_diss=R_diss+R_diss1*dt

call calculate_k(rho,k1)

  rho=rho+dt*(k1)
 
enddo

  close(100)
  close(101)

end subroutine main
!---------------------------------------------------------- 

subroutine compute_R_matrix(dt,n_p)
  implicit none
  real*8,intent(in) :: dt
  integer,intent(in) :: n_p
  integer,parameter :: n_omg=20000
  integer i,j,k,l,i1,i2,n,m,n_t1
  real*8 t,tmin,tmax,dt1
  real*8 w,wmin,wmax,dw,omg(n_omg),wt,c_p(nquant,nquant,nquant,nquant),l_exit(nquant,nquant,nquant,nquant)
  complex*16 su1,su2,fac(n_omg),su1_dot,su1_dot2,A(nquant),F(nquant),kai(nquant,nquant)
  real*8 t1,t2

  call cpu_time(t1)

V_reorg=gama*hbar/(pi**2*omg_c**2)*(2d0*omg_c**3)
  

  wmin=0.d0
  wmax=20*omg_c

  dw=(wmax-wmin)/dfloat(n_omg-1)
  
do i=1,n_omg
    w=wmin+i*dw
!    write(11,*)w/(2*pi*clight),drude_density(w)/wave_to_j!,drude_density(w)/wave_to_j
  enddo

  tmin=0.d0

i1=n_p
    t=tmin+(i1-1)*dt
    
    su1=0.d0
    su1_dot=0.d0
    su1_dot2=0.d0

    do i2=1,n_omg
      w=wmin+(i2)*dw
      wt=w*t

su1=su1+drude_density(w)/(w*w)*((1-dcos(wt))/dtanh(0.5d0*beta*hbar*w)+iota*(dsin(wt)-wt))
su1_dot=su1_dot+drude_density(w)/(w*w)*((dsin(wt)*w)/dtanh(0.5d0*beta*hbar*w)+iota*(w*dcos(wt)-w))
su1_dot2=su1_dot2+drude_density(w)/(w*w)*((dcos(wt)*w**2)/dtanh(0.5d0*beta*hbar*w)-iota*(w**2*dsin(wt)))


    enddo


    su1=su1/(hbar) * dw
    su1_dot=su1_dot/(hbar) * dw
    su1_dot2=su1_dot2/(hbar) * dw

l_exit=0d0
c_p=0d0
do i=1,nquant
do j=1,nquant
do k=1,nquant
do l=1,nquant

do n=1,nquant
c_p(i,j,k,l)=c_p(i,j,k,l)+(U_exc(n,i)*U_exc(n,j)*U_exc(n,k)*U_exc(n,l))

l_exit(i,j,k,l)=l_exit(i,j,k,l)+(U_exc(n,i)*U_exc(n,j)*U_exc(n,k)*U_exc(n,l))*V_reorg

enddo

g_p(i,j,k,l)=c_p(i,j,k,l)*su1
g_p_dot(i,j,k,l)=c_p(i,j,k,l)*su1_dot
g_p_dot2(i,j,k,l)=c_p(i,j,k,l)*su1_dot2


enddo
enddo
enddo
enddo


do i=1,nquant
A(i)=exp(-iota*E_excit(i)*t/hbar-g_p(i,i,i,i))
F(i)=exp(-iota*(E_excit(i)-2d0*l_exit(i,i,i,i))*t/hbar-conjg(g_p(i,i,i,i)))
enddo

do i=1,nquant
do j=1,nquant

kai(i,j)=exp(2d0*(g_p(i,i,j,j)+iota*l_exit(i,i,j,j)*t/hbar))*(g_p_dot2(j,i,i,j)-(g_p_dot(j,i,i,i)-g_p_dot(j,i,j,j)&
&-2d0*iota*l_exit(j,i,j,j)/hbar)*(g_p_dot(i,j,i,i)-g_p_dot(i,j,j,j)-2d0*iota*l_exit(i,j,j,j)/hbar))

R_diss1(i,j)=2d0*real(kai(i,j)*A(i)*conjg(F(j)))


enddo
enddo


!R_diss=real(2d0*R_diss1)



do i=1,nquant
do j=1,nquant
R_coh(i,j)=(real(g_p_dot(i,i,i,i)+g_p_dot(j,j,j,j)-2d0*g_p_dot(i,i,j,j))+iota*aimag(g_p_dot(i,i,i,i)-g_p_dot(j,j,j,j)))
enddo
enddo



  call cpu_time(t2)

end subroutine compute_R_matrix

!-----------------------------------------------------------------  
subroutine calculate_k(rho,ki)
implicit none
  complex*16,intent(in)::rho(nquant,nquant)
  complex*16,intent(out),dimension(nquant,nquant):: ki
  integer :: i,j,k,l

ki=0d0
do i=1,nquant

   do k=1,nquant
      ki(i,i)=ki(i,i)+((R_diss(i,k)*rho(k,k))-(R_diss(k,i)*rho(i,i)))
   enddo


do j=1,nquant

if (i.ne.j) then
 
  do l=1,nquant
      ki(i,j)=ki(i,j)-0.5d0*(R_diss(l,i)+R_diss(l,j))*rho(i,j)
   enddo

ki(i,j)=ki(i,j)-R_coh(i,j)*rho(i,j)
ki(i,j)=ki(i,j)-iota*((E_excit(i)-E_excit(j))/hbar)*rho(i,j)


!write(*,*) i,j,ki(i,j)
endif
enddo
enddo
!stop

end subroutine calculate_k

!---------------------------------------------------------------
subroutine rk4(rho,dt)
  implicit none
  complex*16,intent(inout)::rho(nquant,nquant)
  real*8,intent(in) :: dt
  complex*16,dimension(nquant,nquant):: k1,k2,k3,k4
  integer :: i,j,k

!k1
!call calculate_k(rho,k1)
!write(*,*) k1(1,1),k1(2,2)
!write(*,*) k1(1,2),k1(2,1)


!  call compute_R_matrix(dt/2,2*n_p-1)
!rho=rho+0.5d0*dt*k1

!k2
!call calculate_k(rho,k2)
!rho=rho+0.5d0*dt*k2

!k3
!call calculate_k(rho,k3)


! call compute_R_matrix(dt,n_p)
!rho=rho+dt*k3

!k4
call calculate_k(rho,k1)


  rho=rho+dt*(k1)

end subroutine rk4
!----------------------------------------------------------------

function density(w)
  implicit none
  real*8 density,w

!  density=0.5*V_reorg*omg*omg * gamma*w/((w*w-omg*omg)**2+(gamma*w)**2)

 ! density=0.5*V_reorg* gamma*w/(w*w+gamma**2)

end function density
!-----------------------------------------------------------------  

function drude_density(w)
  implicit none
  real*8 drude_density,w
  real*8 gama_d,lambda_D

  lambda_D=V_reorg!/4.d0
!  gama_D=omg**2/gamma
  !gama_D=106.15*(2*pi*clight)
! gama_D=gama

!  drude_density=0.5*lambda_D* gama_D*w/(w*w+gama_D**2)
  drude_density=gama*hbar*(w**3)*exp(-w/omg_c)/(pi*omg_c**2)

end function drude_density
!-----------------------------------------------------------------  

End Module mod_redfield
