! PROGRAM TO CALCULATE THE ENERGY DISTRIBUTION OF 1-3 DECAY TYPE PARTICLES
! /** present.f90 SOUMYANANDA GOSWAMI ; Uses gfortran compiler, fortran 95 type layout on Ubuntu 16.04 */
! /**REFER BARGER & PHILLIPS COLLLIDER PHYSICS,REVISED ED, CHAPTER 11-MC SIMULATIONS --ALGO COPYRIGHTED*/

! Main program starts


PROGRAM charmdecay
implicit none
	real*8, dimension(:,:) :: bin(25,3) , err(25,3)
	integer,dimension(:)   :: nb(3) ! Range 0-25 GeV, 3 particles
	real*8	:: mc,ms,me,mg,dx,dy,dz,gamma,sd,var,w, se, sx, sy, sz
	real*8  :: ee, ex, ey, ez, ge, gx, gy, gz
	integer :: m,n,i,j,t,s !
	real*8 ,external:: alam

	! OPEN FILE FOR WRITING DATA OF BINNED EVENTS
		OPEN(20,ACCESS='SEQUENTIAL',STATUS='REPLACE',FILE='bin.txt',IOSTAT=m)
		write(20,*)"  j 		           s               e             nu" ! FILE HEADER


	! calculate energy distributions from c - s, e, nu
	! gnu decay in flight

	mc = 1.87d0 ! mass of charm quark
	ms = 0.50d0 ! mass of strange quark
	me = 0.000511d0!mass of electron
	mg = 0.00d0!mass of neutrino
	dz = 23.0d0 ! Initial momentum of incident charm in GeV/c
	n  = 100000  ! iterations in monte carlo

	! Write the values on screen

	write(*,15) mc, ms, dz, n
15 	format('mc, ms, dz, n = ', 2f7.4,f8.4,i10) !formatted to floating point and integer respectively
	
	! initialise the arrays, set all values to zero
	do  i = 1,25
		do  j = 1,3
			err(i,j) = 0.0d0
	 		bin(i,j) = 0.00d0
		end do
	end do
	

	gamma = 0.0d0 !decay width
	var   = 0.0d0 !Variance

	do  j = 1,n
		call cdec (mc, ms, me, mg,w, se, sx, sy,sz,ee, ex, ey, ez, ge, gx, gy, gz)
		gamma = gamma + w/n
		var = var + w**2/n
	end do
		var = var - gamma**2.0d0
		sd = sqrt(var/n)
		write(*,20) gamma, sd
20 		format('gamma, sd = ' , 2e12.4,'GeV')
		
! implement event generation

do  i = 1,n
	 
	call cdec  (mc, ms, 0.0d0, 0.0d0, w, se, sx, sy, sz, ee, ex, ey, ez, ge, gx, gy, gz)
	call boost (mc, 0.0d0, 0.0d0, dz, ee, ex, ey, ez)	

	call boost (mc, 0.0d0, 0.0d0, dz, se, sx, sy, sz)
	
	call boost (mc, 0.0d0, 0.0d0, dz, ge, gx, gy, gz)

! possible acceptance cuts can be inserted here
! sort each event into appropriate bin numbers nb(i) ,

	nb(1) = se
	nb(2) = ee
	nb(3) = ge

	do j=1,3,1
	bin(nb(j),j) = bin(nb(j),j) + w/real(n)
	err(nb(j),j) = err(nb(j),j)+ (w**2.0d0)/n
	end do	
	gamma = gamma + w/real(n)
	
end do
	write(*,*)"The decay width is:", gamma
	do  i = 1,25
		do  j = 1,3
		err(i,j) = err(i,j) - bin(i,j)**2
		err(i,j) = sqrt(err(i,j)/n)/gamma
	        bin(i,j) = (bin(i,j)/gamma)
		end do
 		write(20,25) i, (bin(i,j), j = 1,3)
25	        format (i4, 6e12.4)
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
write(30,*)"set output 'distribution.eps' "
write(30,*)"set title 'ENERGY DISTRIBUTION PLOT' "
write(30,*)"set xlabel '2*ENERGY/m_{c} ' "
write(30,*)"set xlabel 'ENERGY in GeV' "
write(30,*)"set xtics 0,1,25 "
write(30,*)"set xrange [0:25] "
write(30,*)"show mxtics"
write(30,*)"set ylabel '1/ {/Symbol G} d {/Symbol G} /dx_{i}'  "
write(30,*)"set ylabel 'weight distribution w (Barger and Phillips)'  "
!write(30,*)"plot for [col=2:4] 'bin.txt' using 1:col smooth bezier title columnheader "
write(30,*)"plot for [col=2:4] 'bin.txt' using 1:col w p title columnheader "

	rewind(30)
	close(30)
!close file

END PROGRAM charmdecay

      subroutine cdec (mc, ms, me, mg,w, se, sx, sy,sz,ee, ex, ey, ez, ge, gx, gy, gz)
	implicit none
! generates c to s + e + gnu (denoted g) in c rest frame.
! weight w
	real*8, parameter  :: Gf = 1.166*0.1**5.0d0 !Fermi Constant
	real*8, parameter  :: pie = 3.14159 !Surface area of your favorite pizza divided by it's radius squared
	real*8, parameter  :: mw = 81.0d0
	real*8, parameter  :: gw = 2.6d0
	real*8, intent(in) :: mc,ms,me,mg
	real*8, intent(out)::w,se,sx,sy,sz,ee,ex,ey,ez,ge,gx,gy,gz
        real*8 ::pmax2,p1,p1sq,e1,m23,pesq,pe,costh,sinth,phi,pcsq,ce,cz,g,bg,b,ct,st,cp,sp,msq,num
	integer::x
	real*8,external:: alam
	!INTEGER, DIMENSION (1) :: seed = (/2*64-1/) 
	CALL RANDOM_SEED (SIZE=x)

! start in c rest frame
	pmax2 = alam(mc**2.0d0, ms**2.0d0, (me+mg)**2.0d0)/(4.0d0*mc**2.0d0)
	CALL RANDOM_NUMBER(num)
	p1sq = pmax2*num
	p1 = sqrt(p1sq)
	e1 = sqrt(p1sq + ms**2)
	m23 = sqrt(mc**2.0d0 + ms**2.0d0 - 2.0d0*mc*e1)
! now work in eg rest frame; assume boosted along z-axis.
	pesq = alam(m23**2.0d0, me**2.0d0, mg**2.0d0)/(4.0d0*m23**2.0d0)
	pe = sqrt(pesq)
	ee = (m23**2.0d0 + me**2.0d0 - mg**2.0d0)/(2.0d0*m23)
	CALL RANDOM_NUMBER(num)
	costh = 1.0d0 - 2.0d0*num
	sinth = sqrt(1.0d0 - costh**2.0d0)
	CALL RANDOM_NUMBER(num)
	phi = 2.0d0*pie*num
	ex = pe*sinth*sin(phi)
	ey = pe*sinth*cos(phi)
	ez = pe*costh
	ge = m23 - ee
	gx = - ex
	gy = - ey
	gz = - ez
	pcsq = alam(m23**2.0d0, mc**2.0d0, ms**2.0d0)/(4.0d0*m23**2.0d0)
	ce = (mc**2.0d0 - ms**2.0d0 + m23**2.0d0)/(2.0d0*m23)
	cz = sqrt(pcsq)
	se = ce - m23
	sz = cz
! matrix element squared and averaged
! Gf in units GeV**(-2)
	msq = 64.0d0*Gf**2.0d0*(ce*ee - cz*ez)*(se*ge - sz*gz)
! for very heavy me, include w-propagator, too
	msq = msq*mw**4.0d0/((m23**2.0d0 - mw**2.0d0)**2.0d0 + (mw*gw)**2.0d0)
	w = msq*p1*pmax2*pe/(64.0d0*(pie**3.0d0)*mc*e1*m23)
! boost back to c rest frame.
! b,g, bg = beta, gamma, beta*gamma
	g  = ce/mc
	bg = - cz/mc
	b  = bg/g
	ee = ee*g + ez*bg
	ge = ge*g + gz*bg
	ez = ee*b + ez/g
	gz = ge*b + gz/g
! randomize orientation in c-restframe
	CALL RANDOM_NUMBER(num)
	ct = 1.0d0 - 2.0d0*num
	st = sqrt(1.0d0 - ct**2.0d0)
	CALL RANDOM_NUMBER(num)
	phi = 2.0d0*pie*num
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
end
subroutine boost (mc, dx, dy, dz, ee, ex, ey, ez)
	implicit none
	real*8,intent(in)   :: mc,dx,dy,dz
	real*8,intent(inout):: ee,ex,ey,ez
	real*8 		    :: bg,g,ct,cps,sps,rz,rx,ry
! boost e-mom from d-rest to lab; pd=dx.dy.dz
! rotate z-axis to d-direction. boost. rotate back
! define boost and angles; avoid theta = 0 ambiguities;
! ct = cos(theta). cps = cos(phi)sin(theta). etc.
	bg = sqrt(dx**2.0d0 + dy**2.0d0 + dz**2.0d0)/mc
	g = sqrt(1.0d0 + bg**2.0d0)
	ct = dz/(mc*bg)
	cps = dy/(mc*bg)
	sps = dx/(mc*bg)
	rz = ex*sps + ey*cps + ez*ct
	ee = ee*g + rz*bg
	rz = ee*bg/g + rz/g
	rx = ex*(1.0d0-sps**2.0d0)-ey*sps*cps-ez*sps*ct+rz*sps
	ry =-ex*sps*cps+ey*(1.0d0-cps**2.0d0)-ez*cps*ct+rz*cps
	ez =-ex*sps*ct-ey*cps*ct+ez*(1.0d0-ct**2.0d0)+rz*ct
	ey = ry
	ex = rx
	return
end


! Kallen Function
	real*8 function alam(a,b,c)
	implicit none
	real*8, intent(in)::a,b,c
	alam = a**2 + b**2 + c**2 - 2*a*b - 2*b*c - 2*c*a
	return
end


