//
//  read_atom.hpp
//  WKM_module
//
//  Created by 翁谋毅 on 2018/12/4.
//  Copyright © 2018 翁谋毅. All rights reserved.
//

#ifndef read_atom_hpp
#define read_atom_hpp

#include <stdio.h>
#include <stdlib.h>


class atom {
private:
	char atom_name[200][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
public:
	enum {cengshu=3,yanshen=27};
	int num_atom;
	double **d;  //晶格常数
	double **rd;
	double **p;//原子的分量坐标
	double ***rp;//原子真实坐标
	int *atom_type;//原子序数的一维数组
	void read_atom(char *name);
	void print_name(int i, FILE *in);
	int type_num;//元素序号种类个数
	int type_name[30];
	int flag_read;
	void freeatom();
};




class cluster{
public:
    int nx, ny, nz;    // jingbao weizhi
    int nn;            // na yi hao yuanzi
    int level;        //di ji ceng
    int num;          //xuhao in total
    double rp[3];
    int ntype;
    int num_cnnct;     // lianjie shu
    int flag_find;
    int flag2;
    int n_type_xuhao;
    
    cluster **cnnct;
    cluster *next;
    cluster *last;
    cluster *mysf;
    void make_cluster(int n, atom aa, double **dist);
    void make_cnnct(double **dist);
    int find_last(int n1, int n2, int n3, int natom);
};

#endif /* read_atom_hpp */
