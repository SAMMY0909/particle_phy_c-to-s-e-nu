/** NAME: SOUMYANANDA GOSWAMI, HRI. NUMERICAL PROJECT, 2ND SEMESTER CODING
*
*     This program performs an integral over two dimensional phase space to find out the distribution of a daughter product in a one to three type 
*decay process. Charm =========> strange + positron+ electron neutrino. The integral is done over that phase space of two particles whose invariant mass 
*distributions are not to be found out. The output is invariant mass distribution data points as a function of a mandelstam type variable in the three-body 
*decay. That mandelstam type variable is bound within cut-offs and that will be the range over which the loop will be run to get the data points. 
*The pathological delta functions are integrated out and the energy delta function therein can be reduced to a Heaviside function which is positive 
*definite therein. The only integral left is an angular integral, a much simpler integral to be evaluated by monte carlo. No azimuthal dependence. 
*The forms of the integrated expression are Lorentz invariant forms. Matrix element square is assumed from Feynman diagram calculations using 
*gamma matrices and trace mechanism. The process is mediated by W+ boson and hence we shall include its effect. Rest of the program is self explanatory.
* Refer Barger and Phillips, Collider Physics
*/
/**Uses std=c++11 and g++ version 4.9.2 (not package provided) g++ 4.9.2 linux mint 17.2 std=c++11*/

//HEADER DEFINITIONS

#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <fstream>
#include <cstdlib>
//USING NAMESPACE STANDARD 

using namespace std;
long double random_uniform_0_1();//generates random number between 0 and 1
long double integsolid(int);// generates the value of monte carlo integration over a solid angle 

//begin main
int main(void){
			//all constants used are defined below
		long double const mc=1.87L,ms=0.50L,mcsq=3.4969L,mssq=0.25L,me=0.0L,mg=0.0L,Gf=1.166e-5L,mw=81.0L,gw=2.5L,gamma=3.096e-13L,vcs=1.02;
/** mc=charm mass, ms =starnge mass, me=electron mass, mg=neutrino mass, Gf= fermi constant, mw=w boson mass, gw=weak coupling constant, 
gamm=decay width,vcs=matrix element of V_cs from Cabibo qquark matrix. All mass units in GeV */
			long double const K=(pow(vcs,2)*pow(Gf,2)*pow(mw,4))/(96*mc*pow(M_PI,5)*50);
//values of the invariant mass distrubution expressions
			long double vals=0.0L,vale=0.0L,valnu=0.0L;
			long double valsx=0.0L,valex=0.0L,valnux=0.0L;
// this array will store the ranges and their corresponding values for the three particles
			long double distn[100][6]; 
			int n=0;
			
			cout<<"ENTER THE NUMBER OF ITERATIONS for monte carlo=  "<<endl;
			cin>>n;
//initialise the array			
	for(int i = 0;i<100;i++){
		for(int j = 0;j<6;j++){
			distn[i][j]= 0.0L;
		}
	}
// define cutoff mandelstam type variable for s quark distribution
	long double slimu=mc-ms;

//define different types of combinations with masses
	long double diff=mcsq-mssq;
	long double sum=mcsq+mssq;
	long double prod=mcsq*mssq;
//enter distribution data into the array for s quark		
	for(int i = 0;i<100;i++){
			valsx=0.0L+(slimu/100)*(i+1);
			long double temp=pow(valsx,2)-pow(mw,2);
			long double temp1=((sum-pow(valsx,2))/4.0*mcsq)-mssq;
			vals=4.0*M_PI*(sqrt(temp1))*K*((pow(diff,2)+pow(valsx,2)*sum-2*pow(valsx,4))/(pow(temp,2)+pow(mw*gw,2)))*integsolid(n);
			distn[i][0]=valsx;
			distn[i][1]=vals/gamma;
		
	}
// define cutoff mandelstam type variable for positron distribution
	long double eliml=ms;
	long double elimu=mc;
//enter distribution data into the array for positron
			
	for(int i = 0;i<100;i++){
			valex=eliml+((elimu-eliml)/100)*(i+1); //break variable limit into sections 
			long double temp=pow(valex,2)-pow(mw,2);
			long double temp1=pow(valex,2)-mssq;
			long double temp2=((mcsq-pow(valex,2))/4.0*mcsq);
// inv mass distn value for some particular value of variable
			vale=4.0*M_PI*(sqrt(temp2))*3.0*K*((pow(temp1,2)*(mcsq-pow(valex,2)))/(pow(valex,2)*(pow(temp,2)+pow(mw*gw,2))))*integsolid(n);
			distn[i][2]=valex;
			distn[i][3]= vale/gamma;
		
	}

// define cutoff mandelstam type variable for neutrino distribution
	long double nuliml=ms;
	long double nulimu=mc;
//enter the distribution data into the array for neutrino			
	for(int i = 0;i<100;i++){
			valnux=nuliml+((nulimu-nuliml)/100)*(i+1);// break variable limit into sections
			long double temp=pow(valnux,2)-pow(mw,2);
			long double temp1=pow(valnux,2)-mssq;
			long double temp2=((mcsq-pow(valnux,2))/4.0*mcsq);
// inv mass distn value for some particular value of variable
valnu=4.0*M_PI*(sqrt(temp2))*K*((pow(temp1,2)*(mcsq-pow(valnux,2))*(2*pow(valnux,4)+pow(valnux,2)*sum+2*prod))/(pow(valnux,6)*(pow(temp,2)+pow(mw*gw,2))))*integsolid(n);
			distn[i][4]=valnux;
			distn[i][5]= valnu/gamma;
		
	}

//open output file for writing distribution data
ofstream outputFile;
outputFile.open("distrib.txt");
//label values on top so that gnuplot can plot
outputFile<<"	      		gs						ge						gnu"<<endl;
for(int i = 0;i<100;i++){
	for(int j =0 ;j<6;j++){
	outputFile<<setprecision(10)<<distn[i][j]<<scientific<<"\t";
	}
outputFile << endl;
}
outputFile.close();	

// GNUPLOT PLOT COMMAND FILE
ofstream outputFile1;

outputFile1.open("gplot.txt");
outputFile1<<"set terminal postscript eps size 12,10 enhanced color font 'Helvetica,20' linewidth 2"<<endl;
outputFile1<<"set output 'distnimage.eps' "<<endl;
outputFile1<<"set title 'INVARIANT MASS DISTRIBUTION PLOT'   "<<endl;
outputFile1<<"set xlabel 'ENERGY(GeV) in terms of X_{i}s'    "<<endl;
outputFile1<<"set xtics 0,0.1,2.5 "<<endl;
outputFile1<<"set xrange [0:2.5]"<<endl;
outputFile1<<"unset yrange"<<endl;
outputFile1<<" show mxtics   "<<endl;
outputFile1<<"set ylabel 'E_i d{/Symbol G} /d^{3} P_{i} =g_{i} ((X_{i})^{2})   INVARIANT MASS DISTRIBUTION proportional to GeV^{4}' "<<endl;
outputFile1<<"plot 'distrib.txt' using 1:2 smooth bezier title columnheader, '' using 3:4 smooth bezier title columnheader, ''using 5:6 smooth bezier title columnheader "<<endl;

outputFile1.close();	

return 0;

}
//end main


//Monte carlo integrator
long double integsolid(int p){
			long double intgval=0.0L,x;
	for(int i=1;i<p;i++){
					x=random_uniform_0_1();
					 intgval += 2*M_PI*sin(x);
	}
	
	intgval=intgval/double (p);

return intgval;

}


//random number generator
long double random_uniform_0_1()
{
return  double(rand())/double (RAND_MAX);
}
//Compilation : g++ -std=c++11 -Wall -o  distri distri.cpp	
//run: ./distri
//for plot: gnuplot 'gplot.txt'
//see distnimage.eps image file in the same directory
