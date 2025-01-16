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
#include <sstream>
using namespace std;

 //The total energy calculation
long double energy(int N,Box& bx,input& dt,pair_E pairlj[][PAIR],Atom at[]){
	
	long double Et0=0.0;
	for(int i=0;i<N-1;i++){
		for(int j=i+1;j<N;j++){
			Et0+=pair_energy(i,j,bx,pairlj,at,dt);
		}
	}

	return (Et0);
}

 //The change of energy calculation
long double delta_energy(int Ni,int N,int Nextra,neigh** cell,Box& bx,input& dt,pair_E pairlj[][PAIR],Atom at[]){
	
	long double delta=0.0;
	int *p;
	p=neigh_array(Ni,cell,bx,dt,at);
	if(Nextra<0){//Calculate the atom/particle will be moved/inserted/deleted
		for(int i=0;i<dt.neigh_arr_len;i++){
			if(*(p+i)!=Ni&&*(p+i)!=N){
				delta+=pair_energy(Ni,*(p+i),bx,pairlj,at,dt);
			}
		}
	}
	if(Nextra>=0){//Calculate the atom/particle has already moved
		for(int i=0;i<dt.neigh_arr_len;i++){
			if(*(p+i)!=Nextra&&*(p+i)!=Ni){
				delta+=pair_energy(Ni,*(p+i),bx,pairlj,at,dt);
			}
		}
	}

	return (delta);
}