//
//  read_information.hpp
//  connect_bond


#ifndef read_information_hpp
#define read_information_hpp
//template <typename Anytype>
#include <stdio.h>
#include <string>
using namespace std;

class atom_aa
{
	static int total;
public:	
	int xuhao;//元素序号
	double rp[3];//真实坐标
	int flag ;//是否被用到了
	int layer ;//标记的O原子的层数
	int org_xuhao;
	atom_aa();
	atom_aa(int org_xuhao,int xuhao, double* rp);
	static int get_total()
	{
		return total;
	}
	
};

//把三角形封装为链表
class triangle
{
public:
	int org_xuhao[3];//三个顶点的原始序号
	double **rp;//三个顶点的真实坐标 3*3
	int sure;//确定每个三角形是什么，0表示四面体空位，1表示四面体，2表示八面体空位，3表示八面体，4表示三棱柱
	triangle* next;
	triangle(triangle* p = NULL) //默认构造函数
	{
		this->org_xuhao[0] = this->org_xuhao[1] = this->org_xuhao[2] = -1;
		this->next = p;
		this->sure = -1;
	}
};

class triangle_list
{
public:
	triangle* head;
	triangle* tail;
	int count;//计数放了多少个元素
	triangle_list() 
	{
		head = new triangle;
		tail = head->next;
		count = 0;				
	}
	~triangle_list()
	{
		delete[]head;
		delete[]head;
	}
	void insert(int *org);//根据原始的序号插入一个三角形


};
struct angle_save
{
	double angle;
	int xuhao[2];
};

class element{
public:
	element();
    int atomic_num;//atomic number
    char name[3];   // element name
    double vdw_radius_min;
    double vdw_radius_max;
    //van der waals radius,min and max
    double cov_radius;//covalence radius
    int num_metal_radius;  // number of metallic radius
    double * metal_radius; // the radius of metallic radius
    int num_common_val;  // number of common valence
    int * common_val;  // common valence
    int num_unusual_val; // number of unusual valence
    int * unusual_val;  // unusual valence
    double electron_negativity;   //electronic negativity
	double first_ionization_energy;    
};

class ionic_radius{
public:
    int atomic_num;    // ionic atomic number
    int coor_num;     // ionic coordination number
    int val_state;    // ionic valence state
    
    double radius;     // ionic radius
    int structure_type;    // ionic structure type, 0=common, 1=square, 2=pyrometric cone
    int spin_stat;     // ionic spin state, 0=common, 1=high spin, 2=low spin
    
};
class cell
{
private:
	char atom_name[200][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
public:
	cell(char *jiegou_name);
	int num;
	double **letice;
	double **p;
	double ***p_real;
	double ***real_position;
	int *type;
	~cell();
};

void read_radius(ionic_radius *****ir,string& name1,string &name2);
void read_element(element *e,string& file1,string& file2,string& file3,string& file4);
int analyse_pure_metal(string& file_name, string & path, string & out_path,element* e);//input is file_name,data path and element class information
double three_jie_chaji(double **a);
double vector_angle(double**a);
int inone_face(double **a);
int stick_label(atom_aa* aa, int num);
int get_triangel(triangle_list& list, atom_aa* aa, int &num, int layer=1 );
void buble_paixu_plus(double *data, int num,int *xuhao);
void buble_struct(struct angle_save* save, int num);
void buble_paixu_plus_forint(int *data, int num, int *xuhao);
#endif /* read_information_hpp */
