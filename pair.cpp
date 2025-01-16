#include "nr3.h"
#include "ran.h"

#include "GNM.h"

#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <time.h>
#include <sstream>
using namespace std;

 //The pair energy calculation
long double pair_energy(int A,int B,Box& bx,pair_E pairlj[][PAIR],Atom at[],input& dt){
	
	long double dxx=abs(at[A].x-at[B].x);
	long double dyy=abs(at[A].y-at[B].y);
	long double dzz=abs(at[A].z-at[B].z);
	if(dxx> bx.LXX/2.0){dxx = dxx-bx.LXX;}
	if(dyy> bx.LYY/2.0){dyy = dyy-bx.LYY;}
	if(dzz> bx.LZZ/2.0){dzz = dzz-bx.LZZ;}
	long double potent=0.0;
	long double ljeps=0.0;
	long double ljsigma=0.0;
	long double ljRcut=0.0;
	
	ljeps=pairlj[at[A].type][at[B].type].eps;
	ljsigma=pairlj[at[A].type][at[B].type].sigma;
	ljRcut=pairlj[at[A].type][at[B].type].Rcut;

	long double rs=dxx*dxx+dyy*dyy+dzz*dzz;
	long double rcs=ljRcut*ljRcut;
	if (rs<rcs){
		long double rcs6=rcs*rcs*rcs;
		long double rcs12=rcs6*rcs6;
		long double rs6=rs*rs*rs;
		long double rs12=rs6*rs6;
		long double sigma6=ljsigma*ljsigma*ljsigma*ljsigma*ljsigma*ljsigma;
		long double sigma12=sigma6*sigma6;
		potent = 4.0*ljeps*((sigma12/rs12-sigma6/rs6)-(sigma12/rcs12-sigma6/rcs6));
	}
	else{
		potent=0.0;
	}
	return(potent);
}