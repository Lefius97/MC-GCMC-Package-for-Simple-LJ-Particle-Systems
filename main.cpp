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

int main(int argc, char *argv[])
{

	struct Atom atom[MAXatom];
	struct Box box;
	struct input data;
	struct pair_E LJ[PAIR][PAIR];

	time_t begin,end;
	long double ret;
	begin=clock();
	
	string input_file_name;
	string output_file_name;
	int com = 1;

	//Commands identification
	if(argc <= 1){
		Show_usage();
		Show_Inputfile();
		exit(0);
	}
	else{ 
		while(com < argc){
			if( (strcmp(argv[com],"-i") == 0) || (strcmp(argv[com],"-in") == 0) ){
				com++;
				input_file_name = argv[com];
			}else if( (strcmp(argv[com],"-o") == 0) || (strcmp(argv[com],"-out") == 0) ){
				com++;
				output_file_name = argv[com];
			}else {
				cout << "ERROR: Unrecognized option -- "<<argv[com]<<endl;
				exit(0);
			}
			com++;
		}
	}
	if (input_file_name.empty()) {
		cout << "ERROR: Wrong input file name setting "<<endl;
		exit(0);
	}

	if (output_file_name.empty()) {
		output_file_name = "output";
	}

	//Start to read the input file and get all the setting
	cout <<endl;
	cout <<"\t"<<"Simulation start: "<<endl;
	cout <<"\t"<<">> Read input file"<<endl;
	Input(input_file_name,data,LJ,box);
	
	//Define the neighbor list cell
	neigh** cell;
	cell =new neigh*[box.neighsize];
	
	//Random number
	Ran myran(data.seed);
	
	//Biuld the initial configuration
	cout <<"\t"<<">> Build initial configuration"<<endl;
	Config(output_file_name,myran,cell,box,data,atom);
	
	cout <<"\t"<<"   >> Volume of the box is "<<box.LXX<<"*"<<box.LYY<<"*"<<box.LZZ<<"="<<box.vol<<endl;
	
	//Calculate the total energy
	data.energy_t=energy(data.Natoms,box,data,LJ,atom);
	
	//Open the output file
	ofstream out;
	out.open((output_file_name+"_configration.dat").c_str());
	ofstream out1;
	out1.open((output_file_name+"_energy.dat").c_str());
	out1 <<"#Steps Energy Density";
	if(data.type_num>1){
		for(int k=0;k<data.type_num;++k){
			out1 <<" Density_type:"<<data.Ntype[data.type_name[k]];				
		}
	}
	out1 <<endl;
	//ofstream outn;
	//outn.open((output_file_name+"_check_nlist.dat").c_str());
	
	int screen=data.MCsteps/4;
	int accept_dis=0;
	int accept_ins=0;
	int accept_del=0;
	int freq_min=data.freq_E;
	int loops=data.MCsteps/data.freq_E;
	int n=0;
	for(int m=0;m<loops;++m){
		if(data.Natoms>=(MAXatom-1)){
			cout << "ERROR: Too many atoms/particles "<<endl;
			exit(0);
		}
		if(n>=data.o_steps_atoms && n%data.freq_atoms==0){
			Output(n,data.Natoms,box,out,atom);
		}
		if(n>=data.o_steps_E){
			Output_E(n,data.Natoms,data.energy_t,box,out1,data);
			//check_cell(n,outn,cell,box);
		}
		for(int l=0;l<freq_min;++l){
		
			if(n==0){cout <<"\t"<<">> Simulation progress :"<<endl;}	
			if(n==screen){cout <<"\t"<<"   >> "<< 25<<"%"<<endl;}
			if(n==2*screen){cout <<"\t"<<"   >> "<< 50<<"%"<<endl;}
			if(n==3*screen){cout <<"\t"<<"   >> "<< 75<<"%"<<endl;}
			
			long double rand_gcmc=myran.doub();	
#ifdef DEBUG				
				double e_before = energy(data.Natoms,box,data,LJ,atom);
				cout << "\nBefore " << setprecision(20) << e_before << endl;
#endif
			if(rand_gcmc<data.frac_gcmc){
				long double choose_ins_del=myran.doub();
				// Insertion of atom/particle
				if(choose_ins_del>=0.5){
					int Ins_atom=data.Natoms;//Number of the insert atom/particle is the last one
					int newtype=insert_atom(Ins_atom,myran,cell,box,data,atom);
					long double dE=delta_energy(Ins_atom,data.Natoms,-1,cell,box,data,LJ,atom);
					long double rand_ins=myran.doub();
					long double prob_ins=insert_prob(dE,box,data,Ins_atom,atom);
					if(rand_ins<=prob_ins){// Accepteance of insertion
						data.Natoms=data.Natoms+1;
						data.Natom_type[newtype]++;
						data.energy_t+=dE;
						accept_ins++;	
						data.acceptence_ins[atom[Ins_atom].type]++;
#ifdef DEBUG							
						cout << "I Accept, diff " << dE << endl;
						double e_after = energy(data.Natoms,box,data,LJ,atom);
						cout << "After " << e_after << " cmp " << e_after-e_before << endl;
						if(fabs(e_after-e_before-dE)>1e-7){exit(0);}
#endif
					}else{
						delete_cell(Ins_atom,cell,box,data,atom);
//cout << "Reject"  << endl;							
#ifdef DEBUG							
						cout << "I Reject, diff " << dE << endl;
						double e_after = energy(data.Natoms,box,data,LJ,atom);
						cout << "After " << e_after << " cmp " << e_after-e_before << endl;
						if(fabs(e_after-e_before)>1e-7){exit(0);}
#endif
					}
				}
				// Deletion of atom/particle
				else if(choose_ins_del<0.5){
					if(data.Natoms>0){
						int Del_atom=delete_atom(myran,data,atom);//Choose the number of the atom/particle to delete randomly
						if(Del_atom!=-1){
							long double dE=-delta_energy(Del_atom,data.Natoms,-1,cell,box,data,LJ,atom);
							long double rand_del=myran.doub();
							long double prob_del=delete_prob(dE,box,data,Del_atom,atom);
							if(rand_del<=prob_del){// Accepteance of deletion
								data.Natoms=data.Natoms-1;		
								del_update_atoms(Del_atom,data.Natoms,cell,box,data,atom);
								data.energy_t+=dE;
								accept_del++;
								data.acceptence_del[atom[Del_atom].type]++;								
#ifdef DEBUG							 
								cout << "D Accept, diff " << dE << endl;
								cout << "New E " << setprecision(16) << data.energy_t << endl;
								double e_after = energy(data.Natoms,box,data,LJ,atom);
								cout << "After " << setprecision(16) << e_after << " cmp " << e_after-e_before <<setprecision(20) << endl;
								if(fabs(e_after-e_before-dE)>1e-7){exit(0);}
#endif		
							}
						}
					}
				}
			}
			
			else{
				// Motion of atom/particle
				int Dis_atom=int(myran.doub()*data.Natoms);//Choose the number of the atom/particle to move randomly	
				int New_Dis_atom=data.Natoms;//Make a virtual atom/particle in new position
				displace_atom(Dis_atom,data.Maxdr,New_Dis_atom,myran,cell,box,data,atom);
				long double dE1=delta_energy(Dis_atom,data.Natoms,-1,cell,box,data,LJ,atom);
				long double dE2=delta_energy(New_Dis_atom,data.Natoms,Dis_atom,cell,box,data,LJ,atom);
				long double dE=-dE1+dE2;
				long double rand_dis=myran.doub();
				long double prob_dis=displace_prob(dE,data);
				if(rand_dis<=prob_dis){// Accepteance of dispalcement
					dis_update_atoms(Dis_atom,data.Natoms,cell,box,data,atom);
					data.energy_t+=dE;
					accept_dis++;
					data.acceptence_dis[atom[Dis_atom].type]++;
				}else{
					delete_cell(New_Dis_atom,cell,box,data,atom);
				}
			}
			n++;
		}
	}
	cout <<"\t"<<"   >> "<< "100%"<<endl;
	cout <<"\t"<<">> Simulation completed"<<endl;

	long double pa_dis=100*double(accept_dis)/double(data.MCsteps);
	long double pa_ins=100*double(accept_ins)/double(data.MCsteps);
	long double pa_del=100*double(accept_del)/double(data.MCsteps);
	if(data.frac_gcmc!=1.0){
		cout <<"\t"<<">> Acceptance percentage of displacement: "<<pa_dis<<"%"<<endl;
		for(int k=0;k<data.type_num;++k){
			long double pa_dis_type=100*double(data.acceptence_dis[data.type_name[k]])/double(data.MCsteps);
			cout <<"\t"<<">> Acceptance percentage of type"<<data.type_name[k]<<" displacement: "<<pa_dis_type<<"%"<<endl;
		}
	}
	if(data.frac_gcmc>0.0 && data.frac_gcmc<=1.0){
		cout <<"\t"<<">> Acceptance percentage of insertion: "<<pa_ins<<"%"<<endl;
		cout <<"\t"<<">> Acceptance percentage of deletion: "<<pa_del<<"%"<<endl;
		for(int k=0;k<data.type_num;++k){
			long double pa_del_type=100*double(data.acceptence_del[data.type_name[k]])/double(data.MCsteps);
			long double pa_ins_type=100*double(data.acceptence_ins[data.type_name[k]])/double(data.MCsteps);
			cout <<"\t"<<">> Acceptance percentage of type"<<data.type_name[k]<<" insertion: "<<pa_ins_type<<"%"<<endl;
			cout <<"\t"<<">> Acceptance percentage of type"<<data.type_name[k]<<" deletion: "<<pa_del_type<<"%"<<endl;
		}
	}
	
	out.close();
	out1.close();
	if(data.o_steps_atoms==data.MCsteps){
		remove((output_file_name+"_configration.dat").c_str());
	}
	//outn.close();

	end=clock();
	ret=double(end-begin)/CLOCKS_PER_SEC;
	int hour = int(ret/3600);
	int minute = int(ret/60); 
	int s=int(ret)%60;
	cout <<"\t"<<">> Runtime : "<<hour<<":"<<minute<<":"<<s<<endl;
	cout<<endl;
}