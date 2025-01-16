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

 //Find which cell the particle in (the cell is start from 0)
int search_cell(int N,Box& bx,input& dt,Atom at[]){

	int nx=int(at[N].x/dt.lc);
	int ny=int(at[N].y/dt.lc);
	int nz=int(at[N].z/dt.lc);
	int num=nx+ny*bx.nlx+nz*bx.nlx*bx.nly;
	
	return(num);
}

 //Build a neighbor-cell list
void neigh_cell(neigh** cell,Box& bx){

	for(int i=0;i<bx.nlz;i++){
		for(int j=0;j<bx.nly;j++){
			for(int k=0;k<bx.nlx;k++){
				int nc=k+j*bx.nlx+i*bx.nlx*bx.nly;
				neigh* n=new neigh;
				n->list=-1;//Initial cell will store -1
				n->next=NULL;
				*(cell+nc)=n;
			}
		}
	}		
}

 //Insert the particle into neighbor-cell list
void insert_cell(int N,neigh** cell,Box& bx,input& dt,Atom at[]){

	int nc=search_cell(N,bx,dt,at);	
	neigh* p =new neigh;
	p->list=N; //create a new neigh pointer and store the particle into it
	p->next=*(cell+nc);//put this new pointer on the top of cell list
	*(cell+nc)=p;
}

 //Delete the particle which was in the neighbor-cell
void delete_cell(int N,neigh** cell,Box& bx,input& dt,Atom at[]){

	int nc=search_cell(N,bx,dt,at);
	neigh* r = *(cell+nc);	

	if(r->list !=N){//if the delete particle is not the first particle of the cell
		while (r->next->list !=N){//search one by one until the previous one
			r=r->next;
		}
		neigh* del=r->next;
		r->next=r->next->next;//skip this particle
		delete del;
	}else{
		*(cell+nc)=r->next;//skip this particle
		delete r;
	}
}

 //Build the array of the particle which in the neighbor-cell
int* neigh_array(int N,neigh** cell,Box& bx,input& dt,Atom at[]){
	
	static int array[NEIatom];
	int nc=search_cell(N,bx,dt,at);
	int nc_near;
	int l=0;
	int x0=int(nc%bx.nlx);
	int linshi=int(nc/bx.nlx);
	int y0=int(linshi%bx.nly);
	int z0=nc/(bx.nlx*bx.nly);
	int xc,yc,zc;
	for(int c1=-1;c1<2;c1++){
		for(int c2=-1;c2<2;c2++){
			for(int c3=-1;c3<2;c3++){
				
				if(l>=(NEIatom-1)){
					cout << "ERROR: Too many atoms/particles "<<endl;
					exit(0);
				}
				
				xc=x0+c1;
				yc=y0+c2;
				zc=z0+c3;
				if(xc<0){xc+=bx.nlx;}
				if(xc>=bx.nlx){xc-=bx.nlx;}
				if(yc<0){yc+=bx.nly;}
				if(yc>=bx.nly){yc-=bx.nly;}
				if(zc<0){zc+=bx.nlz;}
				if(zc>=bx.nlz){zc-=bx.nlz;}
				nc_near=xc*1+yc*bx.nlx+zc*bx.nlx*bx.nly;//find the near cell even at the periodic boundaries	
				neigh* t=*(cell+nc_near);
				while(t->list!=-1){
					array[l]=t->list;//store the particle in these near cells into a array
					t=t->next;
					l+=1;
				}
			}
		}
	}
	dt.neigh_arr_len=l;
	
	return array;
}

 //Check neighbor-cell list
void check_cell(int n,ofstream& outn,neigh** cell,Box& bx){

	outn << "Step: " <<n<<endl;
	outn << "Cell Memebers" <<endl;
	int num=0;
	for(int i=0;i<bx.nlz;i++){
		for(int j=0;j<bx.nly;j++){
			for(int k=0;k<bx.nlx;k++){
				int nc=k+j*bx.nlx+i*bx.nlx*bx.nly;
				outn << nc << " --> ";
				neigh* n=*(cell+nc);
				while(n->next!=NULL){
					outn << n->list<<" ";
					n=n->next;
					num++;
				}
				outn << endl;
			}
		}
	}
	outn << "All memebers are " <<num<<endl;
	outn << endl;
}