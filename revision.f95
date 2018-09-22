! PROGRAM TO CALCULATE THE ENERGY DISTRIBUTION OF 1-3 DECAY TYPE PARTICLES
! /** revision.f90 SOUMYANANDA GOSWAMI ; Uses gfortran compiler, fortran 95 type layout on linux mint 17.2 */
! /**REFER BARGER & PHILLIPS COLLLIDER PHYSICS,REVISED ED, CHAPTER 11-MC SIMULATIONS --ALGO COPYRIGHTED*/

! Main program starts
program main

!all variables declared explicitly
implicit none ! no prev defined datatype

 ! three arrays for different purposes, index for three particles and error array for standard deviation
! bin is used to store the number of events in that energy interval, row index is for the energy range binning, column index is for the particle label s,e,gnu (gnu is neutrino)
! error calculates the n shot stdev
! nb is used to sort events into appropriate bin numbers

		real*8, dimension(:,:) :: bin(25,3) , error(25,3) ! Range 0-25 GeV, 3 particles
		integer, dimension(:) ::	nb(3) ! used to sort events into appropriate bin numbers

		real*8	:: mc,ms,me,mg,dx,dy,dz,gamm,w1, se1, sx1, sy1, sz1, ee1, ex1, ey1, ez1, ge1, gx1, gy1, gz1,r

		integer :: m,n,i,j,s,t,x,lr ! s,t are dummy variables n is for number of iterations in monte carlo
	

! OPEN FILE FOR WRITING DATA OF BINNED EVENTS
	OPEN(20,ACCESS='SEQUENTIAL',STATUS='REPLACE',FILE='bin.txt',IOSTAT=m)
	write(20,*)"  j 		           s               e             nu" ! FILE HEADER

! calculate energy distributions from s - k, e, nu
! gnu decay in flight

	mc = 1.87d0 ! mass of charm quark
	ms = 0.50d0 ! mass of strange quark
	me= 0.00d0
	mg=0.00d0
	dz = 20.0d0! Initial momentum of incident charm in GeV/c
	n = 10000 ! iterations in monte carlo

	print 1, mc, ms, dz, n
	1 format (' mc, ms, dz, n =  ', 2f7.2,f8.2,i6)

! initialise the arrays
do  i = 1,25
	do  j = 1,3
	error(i,j) = 0.0d0
	 bin(i,j) = 0.00d0
	end do
end do
gamm=0.0d0

! implement event generation
do i = 1,n
		call cdec (mc, ms, me, mg, w1, se1, sx1, sy1, sz1, ee1, ex1, ey1, ez1, ge1, gx1, gy1, gz1)
		call boost (mc, 0.0d0, 0.0d0, dz, se1, sx1, sy1, sz1)
		call boost (mc, 0.0d0, 0.0d0, dz, ee1, ex1, ey1, ez1)
		call boost (mc, 0.0d0, 0.0d0, dz, ge1, gx1, gy1, gz1)

! possible acceptance cuts can be inserted here
! sort each event into appropriate bin numbers nb(i) ,

		nb(1) = int( se1 + 1.0d0)
		nb(2) = int( ee1+ 1.0d0)
		nb(3) = int( ge1 + 1.0d0)
	
	do j=1,3
		bin(nb(j),j) = bin(nb(j),j) + (w1/real(n))
		error(nb(j),j) = error(nb(j),j) + (w1**2)/(real(n))
	end do

gamm = gamm + (w1/real(n)) ! implementing monte carlo for decay width

end do
write(*,*)gamm ! write decay width
do  i = 1,25
		do  j = 1,3
			error(i,j) = error(i,j) - bin(i,j)**2 ! n shot standard deviation
			error(i,j) = sqrt(error(i,j)/real(n))/gamm
			bin(i,j) = (bin(i,j)/gamm) !*(mc/2.0d0)
		end do
		! x=int(2*j/1.87d0)
	write(20, *) i, (bin(i,j), error(i,j), j = 1,3) ! write events to a file 
	!7 format (I4, 6e12.4)
end do

! check file integrity
t=0;
do 
	read(20,*,iostat=m)
	if(m<0) exit
	t=t+1
end do

	rewind(20)
	close (20)
!close file


!open gnuplot command file
OPEN(30,ACCESS='SEQUENTIAL',STATUS='REPLACE',FILE='gp.txt',IOSTAT=s,position='rewind')
write(30,*)"set terminal postscript eps enhanced"
write(30,*)"set output 'binned.eps' "
write(30,*)"set title 'ENERGY DISTRIBUTION PLOT' "
!write(30,*)"set xlabel '2*ENERGY/m_{c} ' "
write(30,*)"set xlabel 'ENERGY in GeV' "
write(30,*)"set xtics 0,1,25 "
write(30,*)"set xrange [0:25] "
write(30,*)"show mxtics"
! write(30,*)"set ylabel '1/ {/Symbol G} d {/Symbol G} /dx_{i}'  "
write(30,*)"set ylabel 'weight distribution w (Barger and Phillips)'  "
write(30,*)"plot for [col=2:4] 'bin.txt' using 1:col smooth bezier title columnheader "
!write(30,*)"plot for [col=2:4] 'bin.txt' using 1:col w p title columnheader "

	rewind(30)
	close(30)
!close file


!subroutine and function declaration
contains

!three square formula
real*8 function alam(a,b,c)
real*8,intent(in) ::a,b,c
alam = a**2 + b**2 + c**2 - 2*a*b - 2*b*c - 2*c*a
return
end function alam


! boost to lab frame
subroutine boost (mc,dx,dy,dz,ee,ex,ey,ez)

implicit none

! boost e-mom from d-rest to lab; pd=dx.dy.dz
! rotate z-axis to d-direction. boost. rotate back
! define boost and angles; avoid theta = 0 ambiguities;
! ct = cos(theta). cps = cos(phi)sin(theta). etc.

real*8 :: g,bg,ct,cps,sps,rx,ry,rz
real*8,intent(in):: mc,dx,dy,dz
real*8,intent(inout):: ee, ex, ey, ez

	bg = sqrt(dx**2 + dy**2 + dz**2)/mc ! beta*gamma
	g = sqrt(1 + bg**2) ! gamma
	ct = dz/(mc*bg)
	cps = dy/(mc*bg)
	sps = dx/(mc*bg)
	rz = ex*sps + ey*cps + ez*ct
	ee = ee*g + rz*bg
	rz = ee*bg/g + rz/g
	rx = ex*(1-sps**2)-ey*sps*cps-ez*sps*ct+rz*sps
	ry =-ex*sps*cps+ey*(1-cps**2)-ez*cps*ct+rz*cps
	ez =-ex*sps*ct-ey*cps*ct+ez*(1-ct**2)+rz*ct
	ey = ry
	ex = rx

return
end subroutine boost

! begin charm decay kinematics
!subroutine cdec (mc, ms, me, mg,w, se, sx, sy,sz, ee,ex,ey,ez,ge,gx,gy,gz,r)
subroutine cdec (mc, ms, me, mg,w, se, sx, sy,sz, ee,ex,ey,ez,ge,gx,gy,gz)

implicit none

real*8 :: pmax2,p1sq,p1,e1,m23,pesq,pe,costh,sinth,phi,ce,cz,pcsq,msq,ct,st,cp,sp,b,g,bg,Gf,pie,mw,gw
real*8,intent(in)	:: 	mc,ms,me,mg
real*8,intent(inout)	:: 	w, se, sx, sy, sz, ee, ex, ey, ez, ge, gx, gy, gz
! declared above are masses, four momenta and constants 

Gf = 1.166e-5
pie = 3.14159265358979323d0
mw = 81.0d0
gw = 2.5d0
dz = 20.0d0
w=0.0d0

! generates c to s + e + gnu (denoted g) stands for neutrino in c rest frame.
! weight w

! start in c rest frame
pmax2 = alam(mc**2, ms**2, (me+mg)**2)/(4.0d0*mc**2)

	 !CALL RANDOM_SEED()

	p1sq = pmax2*rand()
p1 = sqrt(p1sq)
e1 = sqrt(p1sq + ms**2)
m23 = sqrt(mc**2 + ms**2 - 2*mc*e1) ! Jackson frame 

! now work in eg rest frame; assume boosted along z-axis.

pesq = alam(m23**2, me**2, mg**2)/(4.0d0*m23**2)
pe = sqrt(pesq)
ee = (m23**2 + me**2 - mg**2)/(2.0d0*m23)
	 !call RANDOM_SEED()
	costh = 1 - 2*rand()
	sinth = sqrt(1 - costh**2)
	!call RANDOM_SEED()
	phi = 2*pie*rand()

ex = pe*sinth*sin(phi)
ey = pe*sinth*cos(phi)
ez = pe*costh

ge = m23 - ee
gx = - ex
gy = - ey
gz = - ez

pcsq = alam(m23**2, mc**2, ms**2)/(4.0d0*m23**2)

	ce = (mc**2 - ms**2 + m23**2)/(2*m23)
	cz = sqrt(pcsq)

se = ce - m23
sz = cz

! matrix element squared and averaged
! Gf in units gev**(-2)
msq = 64*(Gf**2)*(ce*ee - cz*ez)*(se*ge - sz*gz)

! for very heavy mc, include w-propagator, too
msq = (msq*(mw**4))/((m23**2 - mw**2)**2 + (mw*gw)**2)
w = (msq*p1*pmax2*pe)/(64.0d0*(pie**3)*mc*e1*m23)
! go to 100
! boost back to c rest frame.
! b. g. bg = beta. gamma. beta*gamma
g = ce/mc
bg = - cz/mc
b = bg/g

ee = ee*g + ez*bg
ge = ge*g + gz*bg
ez = ee*b + ez/g
gz = ge*b + gz/g
! 100 continue

! randomize orientation in c-restframe
	ct = 1 - 2*rand()
	st = sqrt(1 - ct**2)

	phi = 2*pie*rand()

	cp = cos (phi)
	sp = sin (phi)

ez = ez*ct - ey*st
gz = gz*ct - gy*st
ey = (ez*st + ey)/ct
gy = (gz*st + gy)/ct

ex = ex*cp - ey*sp
gx = gx*cp - gy*sp

ey = (ex*sp + ey)/cp
gy = (gx*sp + gy)/cp

se = mc - ee - ge
sx = -ex - gx
sy = -ey - gy
sz = -ez - gz

return
end subroutine cdec

! end decay dynamics subroutine 

end program main

! end main

! compilation : gfortran -o revision revision.f90
! run: ./revision
! generate plot: gnuplot 'gp.txt'
