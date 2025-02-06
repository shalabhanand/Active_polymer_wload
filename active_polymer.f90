! This code provides a 3D simulation of an active polar filament which has a load at the front.
! n= no of monomers in a chain, f_act=active force strength,rho=factor which controls the bending stiffness of the load, k_spring=spring constant, k_bend=bending stiffness constant, k_clamp=bending stiffness related to the load, eps=epsilon of the LJ potential, sig=diameter of a monomer, zeta=friction coefficient, xboxl,yboxl,zboxl are the dimension of the simulation box in x, y and z direction respectively, delt=simulation time-step, time_sim=total simulation time
module polymer_modeling
implicit none
save
integer, parameter :: n = 201
real(8), parameter :: f_act = 20.0, rho=0.0
real(8), parameter :: k_spring=1000, k_bend=100, k_clamp=rho*k_bend
real(8), parameter :: PI=3.1415926,eps=1.0, sig=1.0d0,zeta=1.0
real(8), parameter :: xboxl=n+100,yboxl=n+100,zboxl=n+100, temp=0.1
real(8), parameter :: factor=1.0, delt=0.00001,a=0.5d0*(delt**2), b1=0.5d0*delt
real(8), parameter :: rc = (2.0d0**(1.0d0/6.0d0))*sig, rc2=rc**2
integer(kind=8), parameter :: time_sim=nint(500/delt)
integer(kind=8), parameter :: sav_data=nint(5/delt),time_cut=nint(0/delt)
real(8) :: sigma(n), fric(n),b(n),rc, rc2
integer :: iseed
real(8) :: fx_act(n),fy_act(n), fz_act(n)
real(8) :: v_lj,v_spring,sumu, v_bend,vtot,w,ke,etot
real(8) :: u,p,g,u1
real(8) :: x(n),y(n),z(n)
real(8) :: fx(n),fy(n),fz(n),x0(n),y0(n),z0(n)
real(8) :: frx(n),fry(n),frz(n)
real(8) :: f1x(n),f1y(n),f1z(n),f_bx(n),f_by(n),f_bz(n)
end module

program main
use polymer_modeling
implicit none
integer(kind=8)::i,j , k, jseed
integer*4 :: timeArray( 3 )
character(10) :: time, t1, t2, date

do i=1, n
if (mod(i,n) .eq. 0) then
sigma(i) = factor*sig
fric(i)=(factor**2)*zeta     !!! We consider friction coefficient to be proportional to the surface area, thus it is factored by the factor squared (prop. to 4pi*r^2)
b(i) = fric(i) * temp
else
fric(i) = zeta
sigma(i) = sig
b(i) = fric(i) * temp
endif
end do

!!!! Seed value for random number generator!!!!!!!!
ISEED = - 10895

open(32,file='fil_conf.dat')

call init()
call force_lj()
call force_spring()
call force_bend()
call force_bd()
call force_sites()
do i = 1, time_sim
call update_position()
call force_lj()
call force_spring()
call force_bend()
call force_bd()
call force_sites()
if (i .gt. time_cut) then
if ( mod(i,sav_data) .eq. 0) then
do k = 1, n
write(32,*) real(x(k)), real(y(k)), real(z(k))
end do
endif
endif
end do
stop
end


subroutine force_spring()
use polymer_modeling
implicit none
integer :: i , j, k
real(8) :: xi,yi,zi,r2i,r6i,vij,wij,r2,w1,r0
real(8) ::rx,ry,rz,f1xi,f1yi,f1zi,r12,f1xij,f1yij,f1zij
        do i= 1,n
        f1x(i) = 0.0
        f1y(i) = 0.0
        f1z(i) = 0.0
        end do
w1 = 0.0
v_spring = 0.0
        do i = 1, n-1
!        if (mod(i,n) .ne. 0) then
        xi = x(i)
        yi = y(i)
        zi = z(i)
        f1xi = f1x(i)
        f1yi = f1y(i)
        f1zi = f1z(i)
                j = i+1
                rx=xi-x(j)
                rx=rx-xboxl*nint(rx/xboxl)

                ry=yi-y(j)
                ry=ry-yboxl*nint(ry/yboxl)

                rz=zi-z(j)
                rz=rz-zboxl*nint(rz/zboxl)

                r2=rx**2+ry**2+rz**2
                r0=0.5*(sigma(i)+sigma(j))
                u = k_spring*(sqrt(r2)-r0)**2
                v_spring = v_spring + u
                w1 = k_spring*(1-r0/sqrt(r2))
                f1xij = -rx*w1
                f1yij = -ry*w1
                f1zij = -rz*w1
                f1x(j) = f1x(j) - f1xij
                f1y(j) = f1y(j) - f1yij
                f1z(j) = f1z(j) - f1zij
!print*, 0.5*v_spring , i, j
        f1x(i) = f1xi + f1xij   
        f1y(i) = f1yi + f1yij
        f1z(i) = f1zi + f1zij
!        endif
        end do
v_spring = 0.5*v_spring
return
end

subroutine force_bend()
use polymer_modeling
implicit none
integer :: i
real(8) :: s, temp_bend
v_bend = 0.0
        do i = 1, n
        f_bx(i) = 0.0
        f_by(i) = 0.0
        f_bz(i) = 0.0
        end do
do i = 1, n-2
        if (mod(i+2,n) .eq. 0) then
        temp_bend=k_clamp
        else
        temp_bend=k_bend
        endif
        s = (x(i+2)-2*x(i+1)+x(i))**2+(y(i+2)-2*y(i+1)+y(i))**2+(z(i+2)-2*z(i+1)+z(i))**2
        v_bend = v_bend + temp_bend*s
        f_bx(i) = f_bx(i) - temp_bend*(x(i+2)-2*x(i+1)+x(i))
        f_by(i) = f_by(i) - temp_bend*(y(i+2)-2*y(i+1)+y(i))
        f_bz(i) = f_bz(i) - temp_bend*(z(i+2)-2*z(i+1)+z(i))
        f_bx(i+1) = f_bx(i+1) +2* temp_bend*(x(i+2)-2*x(i+1)+x(i))
        f_by(i+1) = f_by(i+1) +2* temp_bend*(y(i+2)-2*y(i+1)+y(i))
        f_bz(i+1) = f_bz(i+1) +2* temp_bend*(z(i+2)-2*z(i+1)+z(i))
        f_bx(i+2) = f_bx(i+2) - temp_bend*(x(i+2)-2*x(i+1)+x(i))
        f_by(i+2) = f_by(i+2) - temp_bend*(y(i+2)-2*y(i+1)+y(i))
        f_bz(i+2) = f_bz(i+2) - temp_bend*(z(i+2)-2*z(i+1)+z(i))
end do
v_bend = 0.5 * v_bend
return
end

subroutine force_bd()
use polymer_modeling
implicit none
integer :: i , j
real(8) :: ran2
        do i = 1 , n
        p = 2 * sqrt( 6 * b(i)/delt)
        frx(i) = p * ( ran2(iseed) - 0.5 )
        fry(i) = p * ( ran2(iseed) - 0.5 )
        frz(i) = p * ( ran2(iseed) - 0.5 )
        end do
return
end

subroutine init()
use polymer_modeling
implicit none
integer :: i , j
real(8) :: ran2 , t_prime
integer :: ix , iy , iz
real(8) :: mot_sumx,mot_sumy,mot_sumz,sumvx,sumvy,sumvz
real(8) :: rx,ry,rz,r,x_cap,y_cap,z_cap
i = 0
do ix = 50 ,xboxl-1,3
        do iy = 50 ,n+49
                do iz = 20 ,n +19
                  i = i+1
                  if ( i .le. n) then
                  x(i) = ix
                  y(i) = iy
                  z(i) = iz
                  endif
                end do
        end do
end do
do i=1, n
x0(i) = x(i) -xboxl*floor(x(i)/xboxl)
y0(i) = y(i) -yboxl*floor(y(i)/yboxl)
z0(i) = z(i) -zboxl*floor(z(i)/zboxl)
end do
return
end

subroutine force_sites()
use polymer_modeling
implicit none
integer :: i , j, act_mon
real(8) :: rm, fxxx, fyyy, fzzz, rx, ry, rz
do i = 1, n
fx_act(i) = 0
fy_act(i) = 0
fz_act(i) = 0
end do
do i = 1, n
if ((mod(i,n) .ne. 0) .and. (mod(i,n) .ne. (n-1)) ) then
rx = (x(i+1)-x(i))
ry = (y(i+1)-y(i))
rz = (z(i+1)-z(i))
rx=rx-xboxl*nint(rx/xboxl)
ry=ry-yboxl*nint(ry/yboxl)
rz=rz-zboxl*nint(rz/zboxl)
rm=sqrt(rx**2+ry**2+rz**2)
fxxx = f_act*rx/rm
fyyy = f_act*ry/rm
fzzz = f_act*rz/rm
fx_act(i) = fx_act(i) + fxxx/2
fy_act(i) = fy_act(i) + fyyy/2
fz_act(i) = fz_act(i) + fzzz/2
fx_act(i+1) = fx_act(i+1) + fxxx/2
fy_act(i+1) = fy_act(i+1) + fyyy/2
fz_act(i+1) = fz_act(i+1) + fzzz/2
endif
end do
return
end

subroutine update_position()
use polymer_modeling
implicit none
integer :: i, j,k
real(8) :: x_cap, y_cap, z_cap, rr
do i=1,n
x(i)=x(i)+delt*(fx(i)+f1x(i)+f_bx(i)+fx_act(i)+frx(i))/fric(i)
y(i)=y(i)+delt*(fy(i)+f1y(i)+f_by(i)+fy_act(i)+fry(i))/fric(i)
z(i)=z(i)+delt*(fz(i)+f1z(i)+f_bz(i)+fz_act(i)+frz(i))/fric(i)
x0(i) = x(i) -xboxl*floor(x(i)/xboxl)
y0(i) = y(i) -yboxl*floor(y(i)/yboxl)
z0(i) = z(i) -zboxl*floor(z(i)/zboxl)
end do
return
end


subroutine force_lj()
use polymer_modeling
implicit none
integer :: i ,j , k
real(8) :: xi,yi,zi,r2i,r6i,vij,wij,fij,fxij,fyij,fzij,r2,ss
real(8) :: rx,ry,rz,f1xi,f1yi,f1zi,r12,fxi,fyi,fzi,ss2
do i = 1,n
fx(i) = 0.0
fy(i) = 0.0
fz(i) = 0.0
end do
v_lj = 0.0
w = 0.0

do i = 1 , n-1
xi = x(i)
yi = y(i)
zi = z(i)
fxi = fx(i)
fyi = fy(i)
fzi = fz(i)
do j = i+1 , n

   rx=xi-x(j)
   rx=rx-xboxl*nint(rx/xboxl)

   ry=yi-y(j)
   ry=ry-yboxl*nint(ry/yboxl)

   rz=zi-z(j)
   rz=rz-zboxl*nint(rz/zboxl)

   r2=rx**2+ry**2+rz**2
   
ss=0.5*(sigma(i)+sigma(j)); ss2=ss**2
rc=(2.0d0**(1.0d0/6.0d0))*ss
rc2=rc**2
if ( r2 .lt. rc2 ) then
r2i = ss2 / r2
r6i = r2i**3
vij = r6i * ( r6i - 1.0 ) +0.25d0
wij = r6i * ( r6i - 0.5 )
v_lj = v_lj + vij
w = w + wij
fij = wij / r2
fxij = rx * fij
fyij = ry * fij
fzij = rz * fij
fxi = fxi + fxij
fyi = fyi + fyij
fzi = fzi + fzij
fx(j) = fx(j) - fxij
fy(j) = fy(j) - fyij
fz(j) = fz(j) - fzij
endif
end do
fx(i) =  fxi
fy(i) =  fyi
fz(i) =  fzi
end do
v_lj = 4.0 *eps* v_lj
do i = 1, n
fx(i) = 48 *eps* fx(i) 
fy(i) = 48 *eps* fy(i) 
fz(i) = 48 *eps* fz(i) 
end do

return
end

FUNCTION ran2(idum)
      IMPLICIT NONE
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8  ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1& 
        ,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791&
        ,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then 

        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue  
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END  FUNCTION
