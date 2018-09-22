/**SOUMYANANDA GOSWAMI This program calculates the invariant mass distribution in a  1--->3 type decay in particle physics*/
/**Uses std=c++11 and g++ version 4.9.2 (not default package provided) g++ 4.9.2 linux mint 17.2 std=c++11*/
//REFER BARGER & PHILLIPS COLLLIDER PHYSICS,REVISED ED, CHAPTER 11-MC SIMULATIONS --ALGO COPYRIGHTED

//HEADER DEFINITIONS

#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>
//USING NAMESPACE STANDARD 
using namespace std;

//FUNCTION PROTOTYPE DECLARATIONS

long double alam (long double,long double,long double) ;//kallen function
double random_uniform_0_1();//random number generator

//charm decay function that calculates decay width and four momenta 
long double* cdec   (long double const,long double const ,long double const,long double const,long double const,long double const,long double const,long double const);

//boost function boosts to lab frame
long double* boost (long double const ,long double ,long double ,long double ,long double ,long double ,long double ,long double,long double,long double,long double,long double,long double,long double,long double,long double );

//main function body, uses four functions declared above

int main(void){
			long double const mc=1.87,ms=0.50,mw=81.0,gw=2.5,Gf=1.166e-5,vcs=1.02;
/** mc=charm mass, ms =strange mass, me=electron mass, mg=neutrino mass, Gf= fermi constant, mw=w boson mass, gw=weak coupling constant, 
gamm=decay width,vcs=matrix element of V_cs from Cabibo quark matrix. All mass units in GeV */

			long double dz,w1=0.0, gamma=0.0,se1=0.0, sx1=0.0, sy1=0.0, sz1=0.0, ee1=0.0, ex1=0.0, ey1=0.0, ez1=0.0, ge1=0.0, gx1=0.0, gy1=0.0, gz1=0.0;

//proportionality constant

			long double const K=(pow(vcs,2)*pow(Gf,2)*pow(mw,4))/(mc*pow(M_PI,3));

					
			int  n;//no of energy divisions
			cout<<"Enter the value of the monte carlo shots:  "<<flush;
			cin >>n;

//declare and initialise other variables
			cout<<"Enter the initial momentum of the charm quark within 20 Gev: "<<flush;
			cin>>dz;
			cout<<endl;
			int dimension;

	
	int count =0;
	if((dz<1.0)&&(dz>0)){
		int digit=0;
		do{
    			digit=int(dz*10);
			count+=1;
		}while(digit==0);

		count+=3;
		dimension=int((dz+5.0)*pow(10,count));
	}
	
	if((dz>=1.0)&&(dz>0)){
	count=ceil(log10(fabs(dz)+1));
	dimension=int((dz+5.0)*1e+3);
	
	}
	if(dz==0){
	dimension=int((dz+5.0)*1e+3);
	}
	const int dim=dimension;
	
//allocating and initialising the energy vs gamma data array				
			long double   bin[dim][3];//this array bins decay widths of the three particles   
			long double   err[dim][3];//this array stores standard deviations
				int   nb[3]={0};//this array stores the random energy value out of cdec for each particle in a single loop

			cout<<endl;
	for(int l=0;l<dim;l++){
		for(int il=0;il<3;il++){
		bin[l][il]=0.0;
		}
	}			





//			long double euplim=dz+5;

//Monte Carlo starts
		
	for(int j=0;j<n;j++){

//   assign four momenta to the pointer slots after cdec operation and boosting respectively *t is for decay kinematics pointer; 
// * b is for boost pointer after basic kinematic operations
		
	long double * t=cdec(mc,ms,0.0,0.0,Gf,mw,gw,K);//pointer to cdec function

	w1=(*t), se1=*(t+1), sx1=*(t+2),sy1=*(t+3), sz1=*(t+4),ee1=*(t+5), ex1=*(t+6), ey1=*(t+7), ez1=*(t+8),ge1=*(t+9), gx1=*(t+10), gy1=*(t+11), gz1=*(t+12);
	
	if(dz>0){
	long double * b=boost (mc, 0.0, 0.0, dz, se1, sx1, sy1, sz1,ee1, ex1, ey1, ez1,ge1, gx1, gy1, gz1);//pointer to boost function

	se1=(*b), sx1=*(b+1),sy1=*(b+2), sz1=*(b+3),ee1=*(b+4), ex1=*(b+5), ey1=*(b+6), ez1=*(b+7),ge1=*(b+8), gx1=*(b+9), gy1=*(b+10), gz1=*(b+11);
	}

//captures energy in the form of greatest integer function for strange quark, positron and its corresponding antineutrino
	if((dz>=1.0)||(dz==0)){
	nb[0]=int((se1)*1e3);
	nb[1]=int((ee1)*1e3);
	nb[2]=int((ge1)*1e3);
	
	}

	if(dz>0){
	if(dz<1.0){
	nb[0]=int((se1)*count);
	nb[1]=int((ee1)*count);
	nb[2]=int((ge1)*count);
	
	}
	}
//loop over the decay widths and allocate them to respective bins
		for(int i=0;i<3;i++){
			bin[nb[i]][i] +=w1/double(n);//the data goes to that bin whose energy is in that range 
			err[nb[i]][i] +=pow(w1,2)/double(n);// similarly for standard deviation
		}
			
		gamma +=w1/n;//the integrated averaged decay width

	}

cout<<"Averaged decay width=......"<<gamma<<endl;

//loop over the matrices to normalise bins and assign the standard deviation matrix

	for (int i = 0; i<dim; i++) {

		    for (int j = 0; j<3; j++) {
				err[i][j] -=pow(bin[i][j],2.0);
				err[i][j]  =sqrt(err[i][j]/double(n))/gamma;
				bin[i][j] /=((2.0/mc)*gamma);
    		    }
	}

//open output file for writing distribution data
		ofstream outputFile;
		outputFile.open("dgammade.txt");

//label values on top so that gnuplot can plot

outputFile<<"energy		s      		e			g"<<endl;
			for(int i = 0;i<dim;i++){
			if((dz>=1.0)||(dz==0)){
				outputFile<<setprecision(12)<<double((i*2.0)/(1000.0*mc))<<scientific<<'\t';
			}
			if((dz>0)&&(dz<1.0)){
			
				    outputFile<<setprecision(12)<<double((i*2.0)/(pow(10,count)*mc))<<scientific<<'\t';


			}
				for(int j =0;j<3;j++){
					outputFile<<setprecision(12)<<bin[i][j]<<scientific<<"\t";//scientific formatting
				}
			  outputFile<<endl;
			}
	outputFile.close();	

// GNUPLOT PLOT COMMAND FILE

		ofstream outputFile1;
		outputFile1.open("plt.txt");

outputFile1<<"set terminal postscript eps size 12,10 enhanced color font 'Helvetica,20' linewidth 2"<<endl;
outputFile1<<"set output 'dist_curve.eps' "<<endl;
outputFile1<<"set title 'INVARIANT MASS DISTRIBUTION PLOT '   "<<endl;
outputFile1<<"set xlabel '(2*ENERGY(GeV)/charm mass)'  "<<endl;
outputFile1<<"set xtics 0.1 "<<endl;
if(dz==0){
	outputFile1<<"set xrange [0:2]"<<endl;
}
//outputFile1<<"set yrange [0:1e-2]"<<endl;

outputFile1<<" show mxtics   "<<endl;
outputFile1<<"set ylabel '1/{/Symbol G}d{/Symbol G}/dx_{i} =Const*f(x_{i})g_{i} ((X_{i})^{2})   INVARIANT MASS DISTRIBUTION proportional to GeV^{4}' "<<endl;

outputFile1<<"plot 'dgammade.txt' using 1:2 w yerrorbars smooth bezier title columnheader, '' using 1:3 w yerrorbars smooth bezier title columnheader, ''using 1:4 w yerrorbars smooth bezier title columnheader "<<endl;

//outputFile1<<"plot 'dgammade.txt' using 1:2 w p yerrorbars title columnheader, '' using 1:3 w p yerrorbars title columnheader, ''using 1:4 w p yerrorbars title columnheader"<<endl;

//outputFile1<<"plot 'dgammade.txt' using 1:2 smooth bezier , '' using 1:3 smooth bezier , ''using 1:4 smooth bezier  "<<endl;

outputFile1.close();	

return 0;

}

//end main

//begin decay kinematics
long double * cdec(long double const mc, long double const ms,long double const me,long double const mg,long double const Gf,long double const mw,long double const gw,long double const K)
{

	long double array[13]={0};
	// generates c to s + e + gnu (denoted g) in c rest frame, weight w
	
	// Start in c rest frame
	long double a=pow(mc,2.0);
	long double b=pow(ms,2.0);
	long double c=me+mg;
	long double d=pow(c,2.0);

	long double ran=random_uniform_0_1();

	long double pmax2 = alam(a,b,d)/(4.0*a);
	long double p1sq = pmax2*ran;
	long double p1 = sqrt(p1sq);
	long double e1 = sqrt(p1sq + b);
	long double m23 = sqrt(a + b - 2*mc*e1);

	//Then work in eg rest frame; assume boost along z-axis
	a=pow(m23,2.0);
	b=pow(me,2.0);
	c=pow(mg,2.0);
	d=pow(mc,2.0);

	long double e=pow(ms,2.0);
	long double pesq=alam(a,b,c)/(4.0*a);
	long double pe=sqrt(pesq);

	//ee,ex,ey,ez are four momentas of electron.. similarly for s quark se,sy,sx,sz and for neutrino.. ge,gx,gy,gz

	long double ee=(a+b-c)/(2.0*m23);

	ran=random_uniform_0_1();

	long double  costh=1.0-2.0*ran;
	long double  sinth=sqrt(1.0-pow(costh,2.0));
	long double  phi=2.0*M_PI*random_uniform_0_1();

	long double ex=pe*sinth*sin(phi);
	long double ey=pe*sinth*cos(phi);
	long double ez=pe*costh;

	long double ge=m23-ee;
	long double gx=-ex;
	long double gy=-ey;
	long double gz=-ez;


	long double pcsq=alam(a,d,e)/(4.0*a);
	long double ce=(d-e+a)/(2.0*m23);
	long double cz=sqrt(pcsq);
	long double se=ce-m23;
	long double sz=cz;

	//Matrix element squared and averaged
	//Gf in units of Gev^-2

	long double msq=(ce*ee-cz*ez)*(se*ge-sz*gz);

	//for very heavy mc, include w propagator too
	msq=(msq*K)/((pow(m23,2.0)-pow(mw,2.0))*(pow(m23,2.0)-pow(mw,2.0))+ pow(mw,2.0)*pow(gw,2.0));

	long double w=(msq*p1*pmax2*pe)/(e1*m23);

	// boost back to c rest frame
	//initialise b as beta and g as gamma here; bg=product : beta*gamma

	long double g=ce/mc;
	long double bg=-cz/mc;
	b=bg/g;

	ee=ee*g+ez*bg;
	ge=ge*g+gz*bg;
	ez=ee*b+ez/g;
	gz=ge*b+gz/g;

	//randomize orientation in c-restframe

	long double ct=1.0-2.0*random_uniform_0_1();
	long double st=sqrt(1-pow(ct,2.0));

	phi=2*M_PI*random_uniform_0_1();
	long double cp=cos(phi);
	long double sp=sin(phi);

	ez=ez*ct-ey*st;
	gz=gz*ct-gy*st;

	ey=(ez*st+ey)/ct;
	gy = (gz*st + gy)/ct;

	ex = ex*cp - ey*sp;
	gx = gx*cp - gy*sp;

	ey = (ex*sp + ey)/cp;
	gy = (gx*sp + gy)/cp;

	se = mc - ee - ge;
	long double sx = -ex - gx;
	long double sy = -ey - gy;
		   sz = -ez - gz;

//store four vectors
	array[0]=w;
	array[1]=se;
	array[2]=sx;
	array[3]=sy;
	array[4]=sz;
	array[5]=ee;
	array[6]=ex;
	array[7]=ey;
	array[8]=ez;
	array[9]=ge;
	array[10]=gx;
	array[11]=gy;
	array[12]=gz;

long double * cdecstr = array;
	return cdecstr;	
		
}

//begin boost to lab frame

long double * boost (long double const mc, long double dx, long double dy,long double dz,long double se,long double sx,long double sy,long double sz,long double ee,long double ex,long double ey,long double ez,long double ge,long double gx,long double gy,long double gz)
{
	long double array[12]={0};
	
// boost e-mom from d-rest to lab; pd=dx,dy,dz
// rotate z-axis to d-direction, boost, rotate back
// define boost and angles; avoid theta = 0 ambiguities;
// ct = cos(theta), cps = cos(phi)sin(theta) etc.

	long double bg = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0))/mc;
	long double   g = sqrt(1 + pow(bg,2.0));

	long double ct = dz/(mc*bg);
	long double cps = dy/(mc*bg);
	long double sps = dx/(mc*bg);

//charm scope
	long double rz = sx*sps + sy*cps + sz*ct;
	 se = se*g + rz*bg;
	rz = se*(bg/g) + (rz/g);

	long double rx = sx*(1-pow(sps,2.0))-(sy*sps*cps)-(sz*sps*ct)+(rz*sps);
	long double ry =-(sx*sps*cps)+(sy*(1-pow(cps,2.0)))-(sz*cps*ct)+(rz*cps);

	sz =-(sx*sps*ct)-(sy*cps*ct)+(sz*(1-pow(ct,2.0)))+(rz*ct);
	sy = ry;
	sx = rx;

//store four vector components into an array
	array[0]=se;
	array[1]=sx;
	array[2]=sy;
	array[3]=sz;

//positron scope
	rz = ex*sps + ey*cps + ez*ct;
	ee = ee*g + rz*bg;
	rz = ee*(bg/g) + (rz/g);

	rx = ex*(1-pow(sps,2.0))-(ey*sps*cps)-(ez*sps*ct)+(rz*sps);
	ry =-(ex*sps*cps)+(ey*(1-pow(cps,2.0)))-(ez*cps*ct)+(rz*cps);
	ez =-(ex*sps*ct)-(ey*cps*ct)+(ez*(1-pow(ct,2.0)))+(rz*ct);

	ey = ry;
	ex = rx;
	array[4]=ee;
	array[5]=ex;
	array[6]=ey;
	array[7]=ez;

//gnu scope
	rz = gx*sps + gy*cps + gz*ct;
	ge = ge*g + rz*bg;
	rz = ge*(bg/g) + (rz/g);
	rx = gx*(1-pow(sps,2.0))-(gy*sps*cps)-(gz*sps*ct)+(rz*sps);
	ry =-(gx*sps*cps)+(gy*(1-pow(cps,2.0)))-(gz*cps*ct)+(rz*cps);
	gz =-(gx*sps*ct)-(gy*cps*ct)+(gz*(1-pow(ct,2.0)))+(rz*ct);
	gy = ry;
	gx = rx;

	array[8]=ge;
	array[9]=gx;
	array[10]=gy;
	array[11]=gz;

	long double * booststr = array;
return booststr;	

}
//end boost

//Kallen Formula
long double alam(long double a,long double b,long double c)
{
	double ans=(pow(a,2.0) + pow(b,2.0) + pow(c,2.0) - 2*a*b - 2*b*c - 2*c*a);
	return ans;
}

//random number generator
double random_uniform_0_1()
{
	return  double(rand())/ double (RAND_MAX);
}

//Compilation : time g++ -std=c++11 -Wall -o  final final.cpp
//Run	      :	./final
//Plot	      : gnuplot 'plt.txt'
