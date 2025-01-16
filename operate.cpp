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

 //Insert a atom/particle in box randomly
int insert_atom(int N,Ran& ran,neigh** cell,Box& bx,input& dt,Atom at[]){
	
	long double lr_x=ran.doub();
	long double lr_y=ran.doub();
	long double lr_z=ran.doub();
	int type_new=int(ran.doub()*dt.type_num);
	at[N].x=lr_x*bx.LXX;
	at[N].y=lr_y*bx.LYY;
	at[N].z=lr_z*bx.LZZ;
	at[N].type=dt.Ntype[dt.type_name[type_new]];
	insert_cell(N,cell,bx,dt,at);
	
	return(dt.type_name[type_new]);
}

 //Delete a atom/particle in box randomly
int delete_atom(Ran& ran,input& dt,Atom at[]){
	
	int type_choose=int(ran.doub()*dt.type_num);
	int delete_atom_choose=0;
	if(dt.Natom_type[dt.type_name[type_choose]]==0){
		delete_atom_choose=-1;
	}else{
		int i=0;
		while(i<dt.Natoms){
			delete_atom_choose=ran.doub()*dt.Natoms;
			if(at[delete_atom_choose].type==dt.Ntype[dt.type_name[type_choose]]){
				break;
			}
			i++;
		}
	}
	
	return(delete_atom_choose);
}

 //Move a atom/particle in box randomly
void displace_atom(int N,long double r,int N2,Ran& ran,neigh** cell,Box& bx,input& dt,Atom at[]){
	
	long double ld_x=ran.doub();
	long double ld_y=ran.doub();
	long double ld_z=ran.doub();
	long double xx=(ld_x-0.5)*r+at[N].x;
	long double yy=(ld_y-0.5)*r+at[N].y;
	long double zz=(ld_z-0.5)*r+at[N].z;
	if(xx<0){xx+=bx.LXX;}
	if(xx>bx.LXX){xx-=bx.LXX;}
	if(yy<0){yy+=bx.LYY;}
	if(yy>bx.LYY){yy-=bx.LYY;}
	if(zz<0){zz+=bx.LZZ;}
	if(zz>bx.LXX){zz-=bx.LZZ;}
	at[N2].x=xx;
	at[N2].y=yy;
	at[N2].z=zz;
	at[N2].type=at[N].type;
	insert_cell(N2,cell,bx,dt,at);
}

 //The probobality of insertion
long double insert_prob(long double E_ins,Box& bx,input& dt,int Ni,Atom at[]){
	
	int N=dt.Natom_type[at[Ni].type];
	long double act=dt.Activity[at[Ni].type];
	
	long double pins=exp(-dt.Beta*E_ins)*act*(bx.vol)/(N+1);
	
#ifdef DEBUG	
	cout << "act " << act <<  " V/(N+1) " << (bx.vol)/(N+1) << " Boltz " << exp(-dt.Beta*E_ins) << " pins " << pins << endl;
#endif

	return (pins);
}

 //The probobality of deletion
long double delete_prob(long double E_del,Box& bx,input& dt,int Ni,Atom at[]){
	
	int N=dt.Natom_type[at[Ni].type];
	long double act=dt.Activity[at[Ni].type];
			
	long double pdel=exp(-dt.Beta*E_del)*N/(bx.vol*act);
	
#ifdef DEBUG	
	cout << "act " << act <<  " N/V " << N/bx.vol << " Boltz " << exp(-dt.Beta*E_del) << " pdel " << pdel << endl;
#endif

	return (pdel);
}

 //The probobality of displacement
long double displace_prob(long double E,input& dt){
	long double pdis=exp(-dt.Beta*E);
	return (pdis);
}

 //Update the atom/particle list after deletion
void del_update_atoms(int N,int Nall,neigh** cell,Box& bx,input& dt,Atom at[]){

	if(Nall!=N){
		delete_cell(Nall,cell,bx,dt,at);
		delete_cell(N,cell,bx,dt,at);
		at[N].x=at[Nall].x;
		at[N].y=at[Nall].y;
		at[N].z=at[Nall].z;  
		dt.Natom_type[at[N].type]--;
		at[N].type=at[Nall].type;
		insert_cell(N,cell,bx,dt,at);
	} else{//Delete the last one 
		delete_cell(Nall,cell,bx,dt,at);
		dt.Natom_type[at[N].type]--;
	}
}

 //Update the atom/particle list after displacement
void dis_update_atoms(int N,int Nall,neigh** cell,Box& bx,input& dt,Atom at[]){

	delete_cell(Nall,cell,bx,dt,at);
	delete_cell(N,cell,bx,dt,at);
	at[N].x=at[Nall].x;
	at[N].y=at[Nall].y;
	at[N].z=at[Nall].z; 
	insert_cell(N,cell,bx,dt,at);
}