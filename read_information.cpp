//
//  read_information.cpp
//  connect_bond
//
//  Created by 翁谋毅 on 2019/9/7.
//  Copyright © 2019 翁谋毅. All rights reserved.
//

#include "read_information.hpp"
#include "tools.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include<cmath>
#pragma warning(disable : 4996)
const int cengshu = 3;
const int yanshen = cengshu * cengshu*cengshu;
const double LAYER_RULE = 0.5;
using namespace std;




atom_aa::atom_aa()
{
	this->xuhao = 0;
	flag = 0;
	for (int i = 0; i < 3; i++)
	{
		this->rp[i] = 0;
	}
}
atom_aa::atom_aa(int org_xuhao, int xuhao, double* rp)
{
	total = total + 1;
	flag = 0;
	this->org_xuhao = org_xuhao;
	this->xuhao = xuhao;
	for (int i = 0; i < 3; i++)
	{
		this->rp[i] = rp[i];
	}
}

element::element()
{
	atomic_num = 0;
	name[0] = '\0';
	num_common_val = 0;
	num_metal_radius = 0;
	num_unusual_val = 0;
}
cell::cell(char *name)
{
	int i, j, k;
	//cout << "expand the :" << cengshu << "layer" << endl;
	char temp[300];
	double x_pian = 0.0;
	double y_pian = 0.0;
	double z_pian = 0.0;
	//strcpy(wenjian, "atom1.config");
	FILE *in;
	in = fopen(name, "rt");
	//system("pause");
	if (in == NULL)
	{
		printf("error of rading atom.config!\n");
		printf("the filename is :%s\n", name);
		//cin.get();
		return;
	}
	fscanf(in, "%d", &num);

	type = (int *)malloc(num * sizeof(int));
	letice = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		letice[i] = (double *)malloc(3 * sizeof(double));
	}
	p = (double **)malloc(num * sizeof(double *));
	for (i = 0; i < num; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	real_position = (double ***)malloc(yanshen * sizeof(double **));
	for (i = 0; i < yanshen; i++)
	{
		real_position[i] = (double **)malloc(num * sizeof(double *));
		for (k = 0; k < num; k++)
			real_position[i][k] = (double *)malloc(3 * sizeof(double));
	}

	p_real = new double **[yanshen];
	for (i = 0; i < yanshen; i++)
	{
		p_real[i] = new double *[num];
		for (k = 0; k < num; k++)
		{
			p_real[i][k] = new double[3];
		}
	}
	while (fgets(temp, 300, in) != NULL)
	{
		if (strstr(temp, "vector") != NULL || strstr(temp, "LATTICE") != NULL)
			break;
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &letice[i][j]);
		}
	}

	fgets(temp, 300, in);
	fgets(temp, 300, in);
	for (i = 0; i < num; i++)
	{

		fscanf(in, "%d", &type[i]);
		fscanf(in, "%lf", &p[i][0]);
		fscanf(in, "%lf", &p[i][1]);
		fscanf(in, "%lf", &p[i][2]);
		fgets(temp, 300, in);
	}
	//int x_xishu = 0;
	//int y_xishu = 0;
	//int z_zishu = 0;

	for (i = 0; i < yanshen; i++)
	{

		for (j = 0; j < num; j++)
		{
			//x_xishu = i/3;
			//y_xishu = i / 3;
			//z_zishu = (i % 9) / 3;

			real_position[i][j][0] = letice[0][0] * p[j][0] + letice[1][0] * p[j][1] + letice[2][0] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][0] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][0] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][0];
			real_position[i][j][1] = letice[0][1] * p[j][0] + letice[1][1] * p[j][1] + letice[2][1] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][1] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][1] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][1];
			real_position[i][j][2] = letice[0][2] * p[j][0] + letice[1][2] * p[j][1] + letice[2][2] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][2] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][2] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][2];
		}
	}
	for (i = 0; i < yanshen; i++)
	{
		for (j = 0; j < num; j++)
		{
			p_real[i][j][0] = (i % cengshu - ((cengshu - 1) / 2)) + p[j][0];
			p_real[i][j][1] = (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) + p[j][1];
			p_real[i][j][2] = -(i / (cengshu * cengshu) - ((cengshu - 1) / 2)) + p[j][2];
		}
	}

	fclose(in);
}

cell::~cell()
{
	int i = 0, j = 0;
	if (num == 0)
		return;
	else
	{
		free(type);
		for (i = 0; i < num; i++)
		{
			free(p[i]);
		}
		free(p);
		for (i = 0; i < 3; i++)
		{
			free(letice[i]);
		}
		free(letice);
		for (i = 0; i < yanshen; i++)
		{
			for (j = 0; j < num; j++)
			{
				free(real_position[i][j]);
				free(p_real[i][j]);
			}
		}
		free(p_real);
		free(real_position);
	}
}



void triangle_list::insert(int *org)//先不动，来一个放一个
{
	int i, j;
	int *temp = new int[3];
	
	buble_paixu_plus_forint(org, 3, temp);

	if (head->next == NULL)
	{
		triangle* p = new triangle;
		for (i = 0; i < 3; i++)
		{
			p->org_xuhao[i] = org[i];
		}
		head->next = p;
		p->next = tail;
		count = 1;
		delete[]temp;
	}
	else
	{
		int falg = 0;//表示是不是之前有过了
		//先检查之前是不是已经有了
		triangle* tt;
		tt = head->next;
		while (tt)
		{
			if (tt->org_xuhao[0] == org[0] && tt->org_xuhao[1] == org[1] && tt->org_xuhao[2] == org[2])
			{
				falg = 1;
				break;
			}
			tt = tt->next;
		}
		if (falg == 1)
			return;


		triangle* p = new triangle;
		for (i = 0; i < 3; i++)
		{
			p->org_xuhao[i] = org[i];
		}
		tail->next = p;
		tail = p;
		tail->next = NULL;
		count++;
		delete[]temp;
	}

}
void read_radius(ionic_radius *****ir,string& file_name1,string& file_name2)
{

    int i,j,k,l,m;
    int line_in;
    int coor_in;
    int val_in, structure_in, spin_state_in;
    double radius_in;
    char c;
    char str[5];
    char str_temp[500];
    FILE *in;
    double max=-9999;
    char atom[120][3]={"D","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"};
    
    for (i=0;i<120;i++)
    {
        for (j=0;j<20;j++)
        {
            for (k=0;k<15;k++)
            {
				for (l = 0; l < 3; l++)
				{
					for (m = 0; m < 3; m++)
					{
						ir[i][j][k][l][m].atomic_num = i;
						ir[i][j][k][l][m].coor_num = k;
						ir[i][j][k][l][m].val_state = j - 8;
						ir[i][j][k][l][m].radius = -9999;
						ir[i][j][k][l][m].structure_type = -9999;
						ir[i][j][k][l][m].spin_stat = -9999;
					}
				}				
            }
        }
    }
    in=fopen(file_name1.c_str(), "r");
    if (in==NULL)
    {
        printf("input_radius does not exist!\n");
        return;
    }
    c='0';
    for (line_in=0;line_in<493;line_in++)
    {
        fscanf(in, "%s", str);
        for (i=0;i<120;i++)
        {
            if (str[0]==atom[i][0] && str[1]==atom[i][1] && str[1]=='\0')
            {
                break;
            }
            if (str[0]==atom[i][0] && str[1]==atom[i][1] && str[2]==atom[i][2] && str[2]=='\0')
            {
                break;
            }
        }
        if (i==120)
        {
            printf("ERROR :%s\n", str);
			cout << line_in << endl;
			cin.get();
            return;
        }
        fscanf(in, "%d", &val_in);
        fscanf(in, "%d", &coor_in);
        fscanf(in, "%d", &structure_in);
        fscanf(in, "%d", &spin_state_in);
        fscanf(in, "%lf", &radius_in);
        ir[i][val_in+8][coor_in][structure_in][spin_state_in].atomic_num=i;
        ir[i][val_in+8][coor_in][structure_in][spin_state_in].val_state=val_in;
        ir[i][val_in+8][coor_in][structure_in][spin_state_in].coor_num=coor_in;
        ir[i][val_in+8][coor_in][structure_in][spin_state_in].spin_stat=spin_state_in;
        ir[i][val_in+8][coor_in][structure_in][spin_state_in].radius=radius_in;
        ir[i][val_in+8][coor_in][structure_in][spin_state_in].structure_type=structure_in;
        /*for (i=0;i<15;i++)
        {
            if (max<ir[i][val_in+8][i].radius)
            {
                max=ir[i][val_in+8][i].radius;
            }
        }
        ir[i][val_in+8][0].radius=max*1.1;*/
    }    
    fclose(in);
    //继续读取我补充的信息
	ifstream fin;
	fin.open(file_name2, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_name2 << endl;
		cin.get();
	}	
	while (fin.peek() != EOF && fin.good())
	{
		fin >> str;
		//cout << str;
		for (i = 0; i < 120; i++)
		{
			if (str[0] == atom[i][0] && str[1] == atom[i][1] && str[1] == '\0')
			{
				break;
			}
			if (str[0] == atom[i][0] && str[1] == atom[i][1] && str[2] == atom[i][2] && str[2] == '\0')
			{
				break;
			}
		}
		if (i == 120)
		{
			//printf("ERROR :%s\n", str);
			//cout << file_name2 << endl;
			//cin.get();
			return;
		}
		fin >> val_in;
		fin >> coor_in;
		fin >> radius_in;
		ir[i][val_in + 8][coor_in][0][0].atomic_num = i;
		ir[i][val_in + 8][coor_in][0][0].coor_num= coor_in;
		ir[i][val_in + 8][coor_in][0][0].val_state = val_in;
		ir[i][val_in + 8][coor_in][0][0].spin_stat =0;
		ir[i][val_in + 8][coor_in][0][0].structure_type = 0;
		ir[i][val_in + 8][coor_in][0][0].radius = radius_in;

	}
	fin.close();


	//做完这个需要补充一点事情，就是取出一个最大值用来预估配位数，我们把不同元素和价态的最大半径放到一个位置，以后方便取用
	for (i = 0; i < 120; i++)
	{
		for (j = 0; j < 20; j++)
		{
			for (k = 0; k < 15; k++)
			{
				for (l = 0; l < 3; l++)
				{
					for (m = 0; m < 3; m++)
					{
						if (max < ir[i][j][k][l][m].radius)
						{
							max = ir[i][j][k][l][m].radius;
						}
					}
				}
			}
			//每种情况取出来的最大值在这里赋值,放到对应的元素和价态，0配位，0结构，0自旋的半径下
			max = ir[i][j][0][0][0].radius;
			max = -9999;
		}
	}
	return;
}

void read_element(element *e, string& file_element_r, string& file_colvance, string& file_electronic_negativity, string& file_first_ionization_energy)
{
	int i, j, k;
	char atom_name[120][3] = { "D","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

	//先初始化信息
	for (i = 0; i < 120; i++)
	{
		strcpy(e[i].name, atom_name[i]);
		e[i].atomic_num = i;
		e[i].num_common_val = 0;
		e[i].num_metal_radius = 0;
		e[i].num_unusual_val = 0;
	}
	i = 0, j = 0, k = 0;
    int line_in;
	ifstream fin;
    int atom_num_in;
    double vdw_r_min_in, vdw_r_max_in;
    double cov_r;
    int num_metal_r;
    double *temp;
    char str[5];
    char str_temp[500];
    char atom[120][3]={"D","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"};
	int num_common_in;
	int num_unsual_in;
	//读取共价键范德华键和金属键

	fin.open(file_element_r.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_element_r << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		fin >> vdw_r_min_in;
		fin >> vdw_r_max_in;
		fin >> cov_r;
		fin >> num_metal_r;
		
		e[atom_num_in].metal_radius = new double[num_metal_r];
		double *temp = new double[num_metal_r];
		for (i = 0; i < num_metal_r; i++)
		{
			fin >> temp[i];
		}
		e[atom_num_in].atomic_num = atom_num_in;
		e[atom_num_in].vdw_radius_max = vdw_r_max_in;
		e[atom_num_in].vdw_radius_min = vdw_r_min_in;
		e[atom_num_in].cov_radius = cov_r;
		e[atom_num_in].num_metal_radius = num_metal_r;
		for (i = 0; i < num_metal_r; i++)
		{
			e[atom_num_in].metal_radius[i] = temp[i];
		}
		free(temp);
	}
	fin.close();
    
    
    //然后开始读取价态相关信息
	fin.clear();
	fin.open(file_colvance.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_colvance.c_str() << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
        fin>>atom_num_in;
        //printf("%d\n", atom_num_in);
		//fin >> num_common_in;
		fin >> e[atom_num_in].num_common_val;        
		e[atom_num_in].common_val = new int[e[atom_num_in].num_common_val];
        for (i=0;i<e[atom_num_in].num_common_val;i++)
        {
			fin >> e[atom_num_in].common_val[i];
        }
		fin >> e[atom_num_in].num_unusual_val;
		e[atom_num_in].unusual_val = new int[e[atom_num_in].num_unusual_val];
        for (i=0;i<e[atom_num_in].num_unusual_val;i++)
        {
			fin >> e[atom_num_in].unusual_val[i];
        }
		//delete[]e[atom_num_in].common_val;
    }
	fin.close();
   
	//开始读取电负性相关信息
	string temp_test;
	fin.clear();
    fin.open(file_electronic_negativity.c_str(), ios::in);	
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_electronic_negativity << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		//cout << atom_num_in << endl;
		fin >> str;
		//fin >> temp_test;
		//cout <<temp_test << endl;
		/*if (str[0] == '\0')
		{
			break;
		}*/
		fin >> e[atom_num_in].electron_negativity;   
		//cout << e[atom_num_in].electron_negativity << endl;
        if (str[0]!=atom[atom_num_in][0] || str[1]!= atom[atom_num_in][1])
        {
            printf("ERROR atom name %s %s %d\n", str, atom[atom_num_in], atom_num_in);
			cout << file_electronic_negativity << endl;
			cin.get();
		}
    } 
	fin.close();
    
	//最后开始读取第一电离能的相关信息
	fin.clear();
	fin.open(file_first_ionization_energy.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_first_ionization_energy << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		//cout << atom_num_in << endl;
		fin >> str;
		/*if (str[0] = '\0')
		{
			break;
		}*/
		fin >> e[atom_num_in].first_ionization_energy;
		//cout << e[atom_num_in].first_ionization_energy << endl;
		if (str[0] != atom[atom_num_in][0] || str[1] != atom[atom_num_in][1])
		{
			printf("ERROR atom name %s %s %d\n", str, atom[atom_num_in], atom_num_in);
			cout << file_first_ionization_energy << endl;
			cin.get();
		}
	}
	fin.close();
    return;
}

int analyse_pure_metal(string& file_name, string & path, string & out_path,element* e)
{
	char atom_name[120][3] = { "D","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	int falg = 0;
	int i = 0, j = 0, k = 0,ii=0;
	int type_num = 0;
	int num_line;
	string temp_name;
	ifstream fin;
	fin.open(file_name, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file :" << file_name << endl;
		cin.get();
		return 0;
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> temp_name;
		if (temp_name.empty())
		{
			break;
		}
		i++;
	}
	cout << "i got " << i << " lines" << endl;
	num_line = i;
	fin.close();
	fin.clear();
	fin.open(file_name, ios::in);
	fin.seekg(ios::beg);	

	//建立stringname的空间和对象
	string *name = new string[i];
	i = 0;
	while (fin.peek() != EOF && fin.good())
	{
		fin >> name[i];
		//name[i] = path + name[i];
		//cout << name[i] << endl;
		i++;		
	}
	fin.close();
	fin.clear();
	/*for (i = 0; i < num_line; i++)
	{
		cout << name[i] << endl;
	}*/
	
	//已经读完了的文件名字下面开始统计金属键和真正键长的相关信息
	i = 0;
	
	double diss = 0;
	double dis_rule;
	double max_dis = 7;//这里我们统计的是7A以内的
	double max_rule[3];//找的时候还有一个规则就是，小于晶胞常数
	cout << num_line << endl;
	string path_data = path;
	int sure_flag[120] = { 0 };
	ofstream ftt;
	ftt.open("sb_name_15.txt", ios::out);
	for (i = 0; i < num_line; i++)
	{		
		ofstream fout;
		//cout << name[i] << endl;
		cell cell_a(const_cast<char*>((path_data+name[i]).c_str()));
		for (ii = 0; ii < 120; ii++)
		{
			sure_flag[ii] = 0;
		}
		type_num = 0;
		for (ii = 0; ii < cell_a.num; ii++)
		{
			sure_flag[cell_a.type[ii]] = 1;
		}
		ii = 0;
		for (ii = 0; ii < 120; ii++)
		{
			type_num += sure_flag[ii];
		}
		//如果原子数过多，我们不管，这种合金不常见
		if (type_num != 2 )
			continue;
		for (ii = 0; ii < 3; ii++)
		{
			max_rule[ii] = pow(pow(cell_a.letice[ii][0], 2) + pow(cell_a.letice[ii][1], 2) + pow(cell_a.letice[ii][2], 2), 0.5)*0.9;
		}
		fout.open(out_path+name[i]+"_result.txt",ios::out);
		for (j = 0; j < cell_a.num; j++)
		{
			for (k = 0; k < yanshen*cell_a.num; k++)
			{
				diss = dis(cell_a.real_position[(yanshen - 1) / 2][j], cell_a.real_position[k / cell_a.num][k%cell_a.num]);
				if (diss > max_dis || diss> max_rule[0]||diss> max_rule[1]||diss> max_rule[2]||abs(diss)<1e-2)
				{
					continue;
				}
				else 
				{
					//先写真实距离再写半径之和
					if (e[cell_a.type[j]].num_metal_radius > 0 && e[cell_a.type[k%cell_a.num]].num_metal_radius > 0)
					{
						
						dis_rule = (e[cell_a.type[j]].metal_radius[0] + e[cell_a.type[k%cell_a.num]].metal_radius[0]) / 100.0;
						if (abs(diss - dis_rule) > 1.5)
							falg = 1;
						fout << abs(e[cell_a.type[j]].electron_negativity - e[cell_a.type[k%cell_a.num]].electron_negativity) << "\t" << diss/dis_rule << "\t"<<atom_name[cell_a.type[j]]<<"\t"<< atom_name[cell_a.type[k%cell_a.num]]<<endl;
					}

				}
			}
			
		}
		fout.close();		
		//cout << "has generated the metal result :" << name[i] << endl;
		if (i % 100 == 0&&i>0)
		{
			cout << "has generate:" << i << endl;
			//cin.get();
		}

		if (falg == 1)
		{
			ftt << name[i] << endl;
		}

		falg = 0;


		
	}

	ftt.close();
	delete[]name;
	return falg;

}

//double dis(double*a, double*b)
//{
//	return(pow(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2), 0.5));
//}
inline double  three_jie_chaji(double **a)
{
	int i = 0, j = 0;
	double temp = 0;
	temp = pow(pow((a[1][1] * a[2][2] - a[1][2] * a[2][1]), 2) + pow((a[1][2] * a[2][0] - a[1][0] * a[2][2]), 2) + pow((a[1][0] * a[2][1] - a[1][1] * a[2][0]), 2), 0.5);
	return temp;
}

inline int inone_face(double **a) //���淵��1�������淵��0
{
	int i = 0, j = 0;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	x1 = a[0][0] - a[3][0];
	x2 = a[1][0] - a[3][0];
	x3 = a[2][0] - a[3][0];

	y1 = a[0][1] - a[3][1];
	y2 = a[1][1] - a[3][1];
	y3 = a[2][1] - a[3][1];

	z1 = a[0][2] - a[3][2];
	z2 = a[1][2] - a[3][2];
	z3 = a[2][2] - a[3][2];
	double k = 0;
	k = (x1 * y2 * z3) + (x2 * y3 * z1) + (x3 * y1 * z2) - (x3 * y2 * z1) - (y3 * z2 * x1) - (z3 * x2 * y1);

	if (-1e-1 < k && k < 1e-1)
		return 1;
	else
		return 0;
}

inline double vector_angle(double**a)
{
	//用来求两个向量的夹角，结果以角度制返回

	double a_mu = pow(pow(a[0][0], 2) + pow(a[0][1], 2) + pow(a[0][2], 2), 0.5);
	double b_mu = pow(pow(a[1][0], 2) + pow(a[1][1], 2) + pow(a[1][2], 2), 0.5);
	double diancheng = a[0][0] * a[1][0] + a[0][1] * a[1][1] + a[0][2] * a[1][2];
	double jungel_orig = acos(diancheng / (a_mu*b_mu)) * 180 / 3.1415926;
	//cout << "the judge angle is:" << jungel_orig << endl;
	return jungel_orig;


}


int  stick_label(atom_aa* aa, int num)
{
	//这个函数输入是一堆骨架原子，结果根据纵坐标的高低排出层的高低
	double max = -100;
	int layer_sure = 0;
	int i, j;
	int *flagg = new int[num];
	for (i = 0; i < num; i++)
	{
		flagg[i] = 0;
	}
	
	while (num > 0)
	{
		for (i = 0; i < num; i++)
		{
			if (aa[i].rp[2] > max && flagg[i] == 0)
			{
				max = aa[i].rp[2];
			}
		}
		//先找到最大值
		for (i = 0; i < num; i++)
		{
			if (aa[i].rp[2] > max - LAYER_RULE && flagg[i] == 0)
			{
				aa[i].layer = layer_sure;
				flagg[i] = 1;
				num--;
			}
		}
		layer_sure++;
		max = -100;
	}
	delete[]flagg;
	return layer_sure;

}


int get_triangel(triangle_list &list, atom_aa* aa,int &num,int layer)
{
	//输入是一个三角形的链表和一堆原子类，作用是在list后面补充链表节点
	int all_count = 0;
	struct angle_save  *angelaa=new struct angle_save[15];

	const double DIS_RULE = 1.0;
	const double ANGLE_RULE = 70.0;
	int neighbor[6];
	for (int i = 0; i < 6; i++)
	{
		neighbor[i] = -1;
	}
	int i, j,k,m;
	int ii, jj;
	double temp_dis;
	//各种需要储存的东西
	double* dis_save = new double[num];
	int *xuhao_save = new int[num];
	double *angle_save = new double[15];
	double **temp_angle = new double*[2];
	for (i = 0; i < 2; i++)
	{
		temp_angle[i] = new double[3];
	}
	int *insert = new int[3];
	
	for (i = 0; i < num; i++) //从每个任意点出发找到最近距离的氧原子
	{
		for (j = 0; j < num; j++)
		{
			temp_dis = dis(aa[i].rp, aa[j].rp);
			dis_save[j] = temp_dis;
			
		}

		buble_paixu_plus(dis_save, num, xuhao_save);
		//1：是选取其中六个最小的
		//先测试排的对不对
		cout << "min is :" << dis_save[0] << endl;
		cin.get();
		for (k = 0; k < 6; k++)
		{
			if (dis_save[k + 1] < (dis_save[1] + DIS_RULE))
			{
				neighbor[k] = xuhao_save[k + 1];
			}

		}
		//:2循环得到角度，并进行排序
		int angele_count = 0;
		for (k = 0; k < 6; k++)
		{
			for (m = 0; m < k; m++)
			{
				//先向temp_angel中储存东西
				for (ii = 0; ii < 3; ii++)
				{
					temp_angle[0][ii] = aa[i].rp[ii] - aa[neighbor[k]].rp[ii];
					temp_angle[1][ii] = aa[i].rp[ii] - aa[neighbor[m]].rp[ii];
				}	
				angle_save[angele_count] = vector_angle(temp_angle);
				angelaa[angele_count].angle = angle_save[angele_count];
				angelaa[angele_count].xuhao[0] = neighbor[k];
				angelaa[angele_count].xuhao[1] = neighbor[m];
				angele_count++;
			}
		}
		buble_struct(angelaa, 15);



		//3：最后是根据得到的结果插入三角形的信息
		for (ii = 0; ii < 15; ii++)
		{
			if (angelaa[ii].angle > ANGLE_RULE)
				break;
			else
			{
				insert[0] = aa[i].org_xuhao;
				insert[1] = aa[angelaa[ii].xuhao[0]].org_xuhao;
				insert[2] = aa[angelaa[ii].xuhao[1]].org_xuhao;
				list.insert(insert);
				all_count++;
			}
		}

	}
	delete[]dis_save;
	delete[]xuhao_save;
	delete[]angle_save;
	for (i = 0; i < 2; i++)
	{
		delete[]angle_save;
	}
	delete[]angle_save;
	delete[]angelaa;
	delete[]insert;
	return all_count;
}



void buble_struct(struct angle_save* save,int num)
{
	int i, j;
	struct angle_save temp_save;
	for (i = 0; i < num-1; i++)
	{
		for (j = 0; j < num - i - 1; j++)
		{
			if (save[j].angle > save[i].angle)
			{
				temp_save = save[j];
				save[j] = save[j + 1];
				save[j + 1] = temp_save;
			}
		}
	}
	return;
}

void buble_paixu_plus(double *data, int num, int *xuhao)
{
	int i, j;
	double temp;
	int temp_xuhao;
	for (i = 0; i < num; i++)
	{
		xuhao[i] = i;
	}
	for (i = 0; i < num - 1; i++)
	{
		for (j = 0; j < num - 1 - i; j++)
		{
			if (data[j] > data[j + 1])
			{
				temp = data[j];
				data[j] = data[j + 1];
				data[j + 1] = temp;


				//同时要记录下他的序号
				temp_xuhao = xuhao[j];
				xuhao[j] = xuhao[j + 1];
				xuhao[j + 1] = temp_xuhao;
			}
		}
	}
	return;
}

void buble_paixu_plus_forint(int *data, int num, int *xuhao)
{
	int i, j;
	int temp;
	int temp_xuhao;
	for (i = 0; i < num; i++)
	{
		xuhao[i] = i;
	}
	for (i = 0; i < num - 1; i++)
	{
		for (j = 0; j < num - 1 - i; j++)
		{
			if (data[j] > data[j + 1])
			{
				temp = data[j];
				data[j] = data[j + 1];
				data[j + 1] = temp;


				//同时要记录下他的序号
				temp_xuhao = xuhao[j];
				xuhao[j] = xuhao[j + 1];
				xuhao[j + 1] = temp_xuhao;
			}
		}
	}
	return;
}