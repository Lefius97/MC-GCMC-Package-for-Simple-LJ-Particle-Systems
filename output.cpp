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

 //Output settting
void Output(int n,int N,Box& bx,ofstream& out,Atom at[]){

	out << "#Number of atoms: "<<N<<endl;
	out << "#MC running step is: "<<n<<endl;
	out << "xhi xlo "<<"0.0 " << bx.LXX<<endl;
	out << "yhi ylo "<<"0.0 " << bx.LYY<<endl;
	out << "zhi zlo "<<"0.0 " << bx.LZZ<<endl;
	for(int i=0;i<N;++i){
		out << i+1 <<" "<<at[i].type<<" "<<at[i].x <<" "<<at[i].y<<" "<<at[i].z<<endl;
	}
	out << endl;
}

 //Output settting
void Output_E(int n,int N,long double E,Box& bx,ofstream& out1,input& dt){

	out1 <<setw(8)<< n <<setw(16) <<E<<setw(16)<<double(N)/(bx.vol);
	if(dt.type_num>1){
		for(int k=0;k<dt.type_num;++k){
			out1 <<setw(16)<<double(dt.Natom_type[dt.type_name[k]])/(bx.vol);				
		}
	}
	out1 <<endl;
}