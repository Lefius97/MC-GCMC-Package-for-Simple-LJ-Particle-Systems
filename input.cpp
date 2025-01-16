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

 //Print the user guide
void Show_usage(){

	cout <<endl;
	cout <<"***"<<"USAGE: "<<endl;
	cout <<"***"<<"\t GNM {-input} {-output}"<<endl;
	cout <<"***"<<"\t --{-input} -in/-i input_file_name"<<endl;
	cout <<"***"<<"\t --{-output} -out/-o output_file_name"<<endl;
	cout <<"***"<<"\t FOR EXAMPLE: ./GNM -i ./mc/in.dat -o ./mc/out"<<endl;
	cout <<"***"<<"\t NOTE: only the input part is required, and others will default as:"<<endl;
	cout <<"***"<<"\t       MC simulation; the output is in the current path and the file name"<<endl;
	cout <<"***"<<"\t       starts with 'output'; the configuration won't output."<<endl;
	cout <<endl;
}

 //Print the input file writing guide
void Show_Inputfile(){

	cout <<"***"<<"Input file writing guide: "<<endl;
	cout <<"***"<<"\t Ntype {type} {atoms} * {chemical_potential}"<<endl;
	cout <<"***"<<"\t Pair {pair1} {pair2} {LJ_eps} {LJ_sigma} {LJ_Rcut} *"<<endl;
	cout <<"***"<<"\t Initial_density *"<<endl;
	cout <<"***"<<"\t Boxlength *"<<endl;
	cout <<"***"<<"\t Temp *"<<endl;
	cout <<"***"<<"\t MCsteps *"<<endl;
	cout <<"***"<<"\t Outputsteps_atoms #default -1 means won't output the configuration"<<endl;
	cout <<"***"<<"\t Outputfreq_atoms #default 1 means output the configuration every steps"<<endl;
	cout <<"***"<<"\t Outputsteps_E #default 0 means won't output the energy"<<endl;
	cout <<"***"<<"\t Outputfreq_E #default 1 means output the energy every steps"<<endl;
	cout <<"***"<<"\t Maxdisplace #default is 0.5*boxlength"<<endl;
	cout <<"***"<<"\t Cellsize #default is 1.0, means the 1.0*Rcut"<<endl;
	cout <<"***"<<"\t Randseed #default is 32"<<endl;
	cout <<"***"<<"\t Fraction_GCMC #default is -1.0"<<endl;
	cout <<"***"<<endl;
	cout <<"***"<<"\t NOTE: You should write down these keywords line by line, and then write "<<endl;
	cout <<"***"<<"\t       down your settings after each keyword, The keywords and setting"<<endl;
	cout <<"***"<<"\t       values can be separated by these symbols:' :/%#'. The keywords"<<endl;
	cout <<"***"<<"\t       marked with * are mandatory. If you don't set the rest keywords, "<<endl;
	cout <<"***"<<"\t       the default value will be used. It should be noted that the "<<endl;
	cout <<"***"<<"\t       'Initial_density' and 'Boxlength' can't exist at the same time."<<endl;
	cout <<"***"<<"\t       And it is better not write down the output keywords if you also"<<endl;
	cout <<"***"<<"\t       set it in command, because if these two conflict, they won't output"<<endl;
	cout <<"***"<<"\t       anything. But if you decide to output, the step to start output must"<<endl;
	cout <<"***"<<"\t       smaller than the simulation steps, and the output interval of"<<endl;
	cout <<"***"<<"\t       configuration must be an integer multiple of the output interval of"<<endl;
	cout <<"***"<<"\t       energy. Finally, there is no order requirement for the placement"<<endl;
	cout <<"***"<<"\t       of allkeywords, and there is a complete example file for reference."<<endl;
	cout <<endl;
}

 //Read the input file
void Input(const string& infile,input& dt,pair_E pairlj[][PAIR],Box& bx){

	const char delim[]=" :/%#";

	ifstream in((infile).c_str());
	if(!in.is_open()){cout << "Error: Input file is empty"<<endl;exit(0);}
	
	//Some default value setting
	dt.o_steps_atoms=-1;
	dt.freq_atoms=1;
	dt.o_steps_E=-1;
	dt.freq_E=1;
	dt.seed=32;
	dt.Maxdr=0.0;
	dt.cell_size=1.0;
	dt.frac_gcmc=-1;
	
	bx.LXX=0.0;
	
	//Read the file line by line and crop to char array
	int ww,str_lg;
	int ee=0;
	char inner[40*7][20];
	char* inreal[40][7];
	string str;
	while(getline(in,str))
	{
		ww=1;
		str_lg=str.length()+1;
		char str_array[str_lg];
		strcpy(str_array,str.c_str());
		char *zh=strtok(str_array,delim);
		while(zh!=NULL){
			sprintf(inner[ee*7+ww],"%s",zh);
			inreal[ee][ww]=&inner[ee*7+ww][0];
			ww+=1;
		zh=strtok(NULL,delim);}
		if(str_array!=NULL){ee+=1;}
	}
	dt.type_num=0;
	dt.pair_num=0;
	long double MaxRcut=2.5;
	long double chemical_potential[TYPEatom];
	int cp=0;
	int pair1=0,pair2=0;
	for(int i=0;i<ee;i++){
		if(inreal[i][1]!=NULL){
			if(strcmp(inreal[i][1],"Ntype")==0){
				dt.type_name[dt.type_num]=atoi(inreal[i][2]);
				dt.Ntype[dt.type_name[dt.type_num]]=dt.type_name[dt.type_num];
				dt.Natom_type[dt.type_name[dt.type_num]]=atoi(inreal[i][3]);
				dt.Natoms+=dt.Natom_type[dt.type_name[dt.type_num]];
				dt.acceptence_ins[dt.type_name[dt.type_num]]=0;
				dt.acceptence_del[dt.type_name[dt.type_num]]=0;
				dt.acceptence_dis[dt.type_name[dt.type_num]]=0;
				if(inreal[i][4]!=NULL){
					chemical_potential[dt.type_name[dt.type_num]]=atof(inreal[i][4]);
					cp++;
				}
				dt.type_num++;
			}
			if(strcmp(inreal[i][1],"Pair")==0){
				pair1=atoi(inreal[i][2]);
				pair2=atoi(inreal[i][3]);
				pairlj[pair1][pair2].eps=atof(inreal[i][4]);
				pairlj[pair1][pair2].sigma=atof(inreal[i][5]);
				pairlj[pair1][pair2].Rcut=atof(inreal[i][6]);
				if(pairlj[pair1][pair2].Rcut>MaxRcut){MaxRcut=pairlj[pair1][pair2].Rcut;}
				if(pair1!=pair2){
					pairlj[pair2][pair1].eps=pairlj[pair1][pair2].eps;
					pairlj[pair2][pair1].sigma=pairlj[pair1][pair2].sigma;
					pairlj[pair2][pair1].Rcut=pairlj[pair1][pair2].Rcut;
				}
				dt.pair_num++;
			}
			if(strcmp(inreal[i][1],"Initial_density")==0){dt.Initial_density=atof(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Boxlength")==0){bx.LXX=atof(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Temp")==0){dt.Temp=atof(inreal[i][2]);}
			if(strcmp(inreal[i][1],"MCsteps")==0){dt.MCsteps=atoi(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Outputsteps_atoms")==0){dt.o_steps_atoms=atoi(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Outputsteps_E")==0){dt.o_steps_E=atoi(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Outputfreq_atoms")==0){dt.freq_atoms=atoi(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Outputfreq_E")==0){dt.freq_E=atoi(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Maxdisplace")==0){dt.Maxdr=atof(inreal[i][2]);}//The max distance of atom/particle dispalcement
			if(strcmp(inreal[i][1],"Fraction_GCMC")==0){dt.frac_gcmc=atof(inreal[i][2]);}
			if(strcmp(inreal[i][1],"Cellsize")==0){dt.cell_size=atof(inreal[i][2]);}//The seed value using for random generation
			if(strcmp(inreal[i][1],"Randseed")==0){dt.seed=atoi(inreal[i][2]);}//The seed value using for random generation
		}
	}
	if(cp!=dt.type_num && dt.frac_gcmc>0){cout<<"ERROR: Incomplete setting of chemical potential"<<endl;exit(0);}
	if(dt.pair_num==0){cout<<"ERROR: No setting of LJ potential"<<endl;exit(0);}
	if(dt.pair_num < ((dt.type_num+1)*dt.type_num/2) ){cout<<"ERROR: Incomplete setting of LJ potential"<<endl;exit(0);}
	if(dt.Natoms==0){cout<<"ERROR: There are no atoms/particles"<<endl;exit(0);}
	if(dt.Natoms>MAXatom){cout<<"ERROR: There are tooo many atoms/particles"<<endl;exit(0);}
	if(dt.Initial_density==0 && bx.LXX<0.1){cout<<"ERROR: The box can't build by density and box size do not exist at the same time"<<endl;exit(0);}
	if(dt.Initial_density!=0 && bx.LXX>=0.1){cout<<"ERROR: The box can't build by density and box size setting at the same time"<<endl;exit(0);}
	if(dt.Temp==0){cout<<"ERROR: Temperature is 0.0"<<endl;exit(0);}
	if(dt.MCsteps==0){cout<<"ERROR: The number of MC loop is 0"<<endl;exit(0);}
	if(dt.o_steps_E==-1){dt.o_steps_E=0;}
	if(dt.o_steps_atoms==-1){dt.o_steps_atoms=dt.MCsteps;}
	if(dt.freq_E<=0){cout<<"ERROR: The output interval of energy can't smaller than 1"<<endl;exit(0);}
	if(dt.o_steps_E>dt.MCsteps || dt.o_steps_atoms>dt.MCsteps){cout<<"ERROR: The output can't start with the step larger than simulation steps"<<endl;exit(0);}
	if(dt.o_steps_atoms>0 && dt.o_steps_atoms<dt.MCsteps && dt.freq_atoms%dt.freq_E!=0){cout<<"ERROR: The output interval of configuration must be an integer multiple of the output interval of energy"<<endl;exit(0);}
	if(dt.cell_size<1.0){cout << "ERROR: Wrong cell size setting, should set the cell size >=1.0"<<endl;exit(0);}
	
	cout <<"\t"<<"   >> Number of all atoms/particles is "<<dt.Natoms<<endl;
	cout <<"\t"<<"   >> Number of atoms/particle types is "<<dt.type_num<<endl;
	
	if(dt.Initial_density!=0 && bx.LXX<0.1){
		bx.LXX=pow(dt.Natoms/dt.Initial_density,1.0/3);
	}
	
	dt.lc=dt.cell_size * MaxRcut;
	bx.nlx=int(bx.LXX/dt.lc);
	dt.lc=bx.LXX/double(bx.nlx);
	if(bx.nlx<3){
		cout<<"ERROR: The system you create is too small to use the version with neighbor list"<<endl;
		exit(0);
	}
	bx.nly=bx.nlx;
	bx.nlz=bx.nlx;
	bx.neighsize=bx.nlx*bx.nly*bx.nlz;
	
	bx.LYY=bx.LXX;
	bx.LZZ=bx.LXX;
	bx.vol=bx.LXX*bx.LYY*bx.LZZ;
	if(dt.Maxdr==0){dt.Maxdr=0.5*bx.LXX;}
	long double lamda=1.0;//Lamda in LJ system is 1.0
	dt.Beta=1.0/dt.Temp;
	for(int i=0;i<cp;++i){
		dt.Activity[dt.type_name[i]]=exp(dt.Beta*chemical_potential[dt.type_name[i]])/(lamda*lamda*lamda);
		cout <<"\t"<<"   >> Chemical_potential of atom/particle type"<<dt.Ntype[dt.type_name[i]]<<" is "<<chemical_potential[dt.type_name[i]]<<endl;
	}
	
	in.close();

 }
 
 //Build a initial configuration
void Config(const string& outfile,Ran& ran,neigh** cell,Box& bx,input& dt,Atom at[]){

	ofstream outc((outfile+"_initial_config.dat").c_str());

	neigh_cell(cell,bx);

	outc << "#Number of atoms: "<<dt.Natoms<<endl;
	outc << "#MC running step is: 0"<<endl;
	outc << "xhi xlo "<<"0.0 " << bx.LXX<<endl;
	outc << "yhi ylo "<<"0.0 " << bx.LYY<<endl;
	outc << "zhi zlo "<<"0.0 " << bx.LZZ<<endl;
	int i=0;
	for(int k=0;k<dt.type_num;++k){
		for(int j=0;j<dt.Natom_type[dt.type_name[k]];++j){
			long double nn_x=ran.doub();
			long double nn_y=ran.doub();
			long double nn_z=ran.doub();
			at[i].x=nn_x*bx.LXX;
			at[i].y=nn_y*bx.LYY;
			at[i].z=nn_z*bx.LZZ;
			at[i].type=dt.Ntype[dt.type_name[k]];
			outc << i+1 <<" "<<at[i].type<<" "<<at[i].x <<" "<<at[i].y<<" "<<at[i].z<<endl;
			insert_cell(i,cell,bx,dt,at);
			i++;
		}
	}
	
	outc.close();
}