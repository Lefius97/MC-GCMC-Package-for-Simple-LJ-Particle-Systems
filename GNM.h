#include <stdio.h>

#define MAXatom 1000000000
#define NEIatom 1000000
#define TYPEatom 500
#define PAIR 12750
struct Atom
{
	long double x,y,z; 
	int type;
};
extern struct Atom atom[MAXatom];

struct Box
{
	long double LXX,LYY,LZZ,vol;
	int nlx,nly,nlz,neighsize;
};
extern struct Box box;

struct input
{
	int Ntype[TYPEatom];
	int Natom_type[TYPEatom];
	long double Activity[TYPEatom];
	int type_num,pair_num;
	int type_name[TYPEatom];
	int Natoms,MCsteps,o_steps_atoms,freq_atoms,o_steps_E,freq_E,seed;
	long double Initial_density,Temp,Beta,energy_t,Maxdr,frac_gcmc;
	long double lc,neigh_arr_len,cell_size;
	int acceptence_ins[TYPEatom];
	int acceptence_dis[TYPEatom];
	int acceptence_del[TYPEatom];
};
extern struct input data;

struct pair_E
{
	long double Rcut;
	long double eps;
	long double sigma;
};
extern struct pair_E LJ[PAIR][PAIR];

struct neigh
{
	int list;
	struct neigh *next;
};

void Show_usage();
void Show_Inputfile();
void Input(const string& infile,input&,pair_E pairlj[][PAIR],Box&);
void Config(const string& outfile,Ran& ran,neigh** cell,Box&,input&,Atom at[]);
int insert_atom(int N,Ran& ran,neigh** cell,Box&,input& dt,Atom at[]);
int delete_atom(Ran& ran,input& dt,Atom at[]);
void displace_atom(int N,long double r,int N2,Ran& ran,neigh** cell,Box&,input&,Atom at[]);
long double insert_prob(long double E_ins,Box&,input&,int Ni,Atom at[]);
long double delete_prob(long double E_del,Box&,input&,int Ni,Atom at[]);
long double displace_prob(long double E,input&);
long double energy(int N,Box&,input&,pair_E pairlj[][PAIR],Atom at[]);
long double delta_energy(int Ni,int N,int Nextra,neigh** cell,Box&,input&,pair_E pairlj[][PAIR],Atom at[]);
long double pair_energy(int A,int B,Box&,pair_E pairlj[][PAIR],Atom at[],input&);
void del_update_atoms(int N,int Nall,neigh** cell,Box&,input&,Atom at[]);
void dis_update_atoms(int N,int Nall,neigh** cell,Box&,input&,Atom at[]);
int search_cell(int N,Box&,input&,Atom at[]);
void neigh_cell(neigh** cell,Box&);
void insert_cell(int N,neigh** cell,Box&,input&,Atom at[]);
void delete_cell(int N,neigh** cell,Box&,input&,Atom at[]);
int* neigh_array(int N,neigh** cell,Box&,input&,Atom at[]);
void Output(int n,int N,Box& bx,ofstream& out,Atom at[]);
void Output_E(int n,int N,long double E,Box&,ofstream& out1,input&);
void check_cell(int n,ofstream& outn,neigh** cell,Box&);

