//
//  main.cpp
//  connect_bond
//
//  Created by wangzhi 2019.9.10
//  Copyright © 2019 王志. All rights reserved.
//

#include "read_information.hpp"
#include "read_atom.hpp"
#include <iostream>
#include "tools.hpp"
#include <stdlib.h>
#include <fstream>
#include <cstring>
#pragma warning(disable : 4996)
using namespace std;
const double ratio = 1.2;
const int NUM_POSITIVE = 5;
const int NUM_NEGATIVE = 5;
char atom_name[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

const int cengshu = 3;
const int yanshen = cengshu * cengshu*cengshu;
int* find_coor(atom& cell_a,  double** dist);//输入是该结构，以及120个元素之间的距离信息，返回每个原子的配位数
int generate_new_atom(string path,element * e, ionic_radius***** ir);
int main(int argc, const char * argv[]) {
    // insert code here...
    int i,j,k,l;
	string supporting_data_path = "E:\\ionc_ridus_series\\ionc_ridus_series\\data\\";
	string file_name1 = supporting_data_path+"input_ionic";
	string file_name2 = supporting_data_path+"input_ionic_plus";
	string file_element_r = supporting_data_path+"ridus";
	string file_colvance = supporting_data_path+"colvance";
	string file_file_nagetivity = supporting_data_path+"negativity.txt";
	string file_first_ionization_energy = supporting_data_path+"first_ionazation_energy.txt";
	//建立离子半径的存储结构	
	ionic_radius***** ir;
	ir= new ionic_radius****[120];
	for (i = 0; i < 120; i++)
	{
		ir[i] = new ionic_radius***[20];
		for (j = 0; j < 20; j++)
		{
			ir[i][j] = new ionic_radius**[15];
			for (k = 0; k < 15; k++)
			{
				ir[i][j][k] = new ionic_radius*[3];
				for (l = 0;l < 3; l++)
				{
					ir[i][j][k][l] = new ionic_radius[3];
				}

			}
		}
	}
	//建立元素的	
	element *e;
	e = new element[120];
	//开始依次读取相关信息
    read_radius(ir,file_name1,file_name2);
    read_element(e,file_element_r,file_colvance,file_file_nagetivity,file_first_ionization_energy);
	//测试一下是不是完整的录入了信息
	//cout << e[6].common_val[1] << endl;
	
	//建立对金属的分析结果
	/*string metal_name = "E:\\ionc_ridus_series\\ionc_ridus_series\\analyse_pure_metal\\metal_name";
	string metal_path = "E:\\ionc_ridus_series\\ionc_ridus_series\\analyse_pure_metal\\metal\\";
	string out_path = "E:\\ionc_ridus_series\\ionc_ridus_series\\analyse_pure_metal\\result\\";
	analyse_pure_metal(metal_name, metal_path,out_path, e);
	cout << "finished!" << endl;*/


	//然后产生一个离子半径的结果
	/*generate
	cin.get();*/



	////产生电负性的表格
	//string file_name[3] = { "1T","2H","BI3" };
	//for (int i = 0; i < 3; i++)
	//{
	//	generate_negativity_indor(file_name[i], e);
	//}
	
	//建立新的批量的结构文件	
	string dir_gne_path = "C:\\Users\\王志\\Desktop\\新二维材料基元2019.10.9\\new_str\\C_like\\config\\";
	if (generate_new_atom(dir_gne_path,e, ir) == 1)
	{
		cout << "generate done!" << endl;
		cin.get();
	}

	//find_matel(ir, e);
	delete[]e;
	for (i = 0; i < 120; i++)
	{
		delete[]ir[i];
	}
	delete[]ir;
    
    
    return 0;
}

int* find_coor(cell& cell_a, double** dist,int * coor)//输入该结构，该原子在结构的次序，统计的信息（离子半径和元素信息），返回它和周围的配位数
{
	int i = 0, j = 0,k=0;
	int coor_num = 0;
	double dis_temp;
	double dis_rule;
	//这个地方需要注意，还要知道价态信息，需要明天补充
	for (i=0; i < cell_a.num; i++)
	{
		//对每个原子寻找其配位数
		for (j = 0; j < yanshen*cell_a.num; j++)
		{
			dis_temp = dis(cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[j/cell_a.num][j%cell_a.num]);
			if (abs(dis_temp) > 1e-4  && dis_temp < dist[cell_a.type[i]][cell_a.type[j%cell_a.num]])
			{
				coor[i]++;
			}	
		}
	}	
	return coor;
}

int find_matel(ionic_radius ***ir, element *e)
{
    char str[500];
    char str1[500];
    FILE *in;
    atom a;
    int i,j;
    in=fopen("name", "r");
    int temp;
    int temp2;
    int metal_xuhao[85] = { 3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111 };
    int flag[120];
    for (i=0;i<120;i++)
    {
        flag[i]=1;   //nonmetal
    }
    for (i=0;i<85;i++)
    {
        flag[metal_xuhao[i]]=-1;   // metal
    }
    //flag[15]=1;
    //flag[16]=1;
    FILE *out;
    int atom_a, atom_b;
    int flag1;
    int valence_atom=-1;
    atom_a=-1;
    atom_b=-1;
    while (fgets(str, 500, in)!=NULL) 
	{
        str[strlen(str)-1]='\0';
        sprintf(str1, "E:\\ionc_ridus_series\\icsd_data\\%s\\atom.config", str);
        a.read_atom(str1);
        temp=1;
        temp2=0;
        for (i=0;i<a.num_atom;i++)
        {
            temp=temp*flag[a.atom_type[i]];
            temp2=temp2+flag[a.atom_type[i]];
        }
        //printf("%s\n", str1);
        if (temp2 ==-1*a.num_atom && a.type_num>1)
        {
            for (i=0;i<2;i++)
            {
                if (flag[a.type_name[i]]==1)
                {
                    atom_a=a.type_name[i];
                }
                if(flag[a.type_name[i]]==-1)
                {
                    atom_b=a.type_name[i];
                }
            }
            flag1=0;
            for (i=0;i<e[atom_a].num_common_val;i++)
            {
                if (e[atom_a].common_val[i]<0)
                {
                    for (j=0;j<e[atom_b].num_common_val;j++)
                    {
                        if (e[atom_a].common_val[i]+ e[atom_b].common_val[j]==0)
                        {
                            valence_atom=e[atom_a].common_val[i];
                            flag1++;
                            break;
                        }
                    }
                }
            }
            if (flag1==1)
            {
                sprintf(str1,"/Users/wengmouyi/Desktop/icsd/ionic_two/%s.config", str);
                out=fopen(str1, "wb");
                fprintf(out, "    %d\n", a.num_atom);
                fprintf(out, "Lattice vector\n");
                for (i=0;i<3;i++)
                {
                    for (j=0;j<3;j++)
                    {
                        fprintf(out, "    %.10lf ", a.d[i][j]);
                    }
                    fprintf(out, "\n");
                }
                fprintf(out, "Position\n");
                for (i=0;i<a.num_atom;i++)
                {
                    
                    if (a.atom_type[i]==atom_a)
                    {
                        fprintf(out, "    %d    %.10lf    %.10lf    %.10lf    1  1  1  !  0   %lf\n",a.atom_type[i], a.p[i][0], a.p[i][1], a.p[i][2],  ir[atom_a][8+valence_atom][0].radius);
                    }
                    if (a.atom_type[i]==atom_b)
                    {
                        fprintf(out, "    %d    %.10lf    %.10lf    %.10lf    1  1  1  !  0   %lf\n",a.atom_type[i], a.p[i][0], a.p[i][1], a.p[i][2],  ir[atom_a][8-valence_atom][0].radius);
                    }
                    

                }
                
                
                printf("%s\n", str);
                fclose(out);
                //sprintf(str1, "mv /Users/wengmouyi/Desktop/icsd/icsd_files_2018/%s/ /Users/wengmouyi/Desktop/icsd/finished/", str);
                //system(str1);
            }
            
            
        }
        a.freeatom();
    }
    //a.read_atom(str1);
    fclose(in);
    return 0;
}



int generate_new_atom(string path,element * e,ionic_radius***** ir)
{
	//这个函数的目的是输入几个文件名字，批量产生对应的结构
	FILE* fout_a=NULL;
	double *positive = new double[120];
	double *negative = new double[120];
	int i = 0, j = 0, k = 0,l=0,m=0;

	//对应的正负情况
	//+3int p[] = {21,39,71, 24,26,44,27,45,77,79,5,13,31,49,81,5,33,51,83};
	//B 这一族int p[] = {5,13,31,49,81};
	//+2 int p[] = { 4,12,20,38,56,88,21,39,22,40,72,23,41,73,24,42,74,25,43,75,26,44,76,27,45,77,28,46,78,29,47,79,30,48,80,31,49,50,82};
	//int p[] = {};
	int p[] = { 5,13,31,49,81 };
	int n[NUM_NEGATIVE] = { 7,15,33,51,83};
	//第一步是找到半径的最大值
	
	double max = -9999;
	string file_out_name;	
	
	for (i = 0; i < 120; i++)
	{
		for (j = 0; j < 8; j++)
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
			//每种情况取出来的最大值在这里赋值,放到对应的元素和价态，0配位，0结构，0自旋的半径下negative情况			
		}
		negative[i] = max;
		max = -9999;
		for (j = 8; j < 20; j++)
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
			//每种情况取出来的最大值在这里赋值,放到对应的元素和价态，0配位，0结构，0自旋的半径下positive情况			
		}
		positive[i] = max;
		max = -9999;
	}	
	const double moren = 0.8;
	for (i = 0; i < 120; i++)
	{
		if (negative[i] > 10 || negative[i] <= 0)
		{
			negative[i] = moren;
		}
		if (positive[i] > 10 || positive[i] <= 0)
		{
			positive[i] = moren;
		}
	}




	char file_name_real[300];
	for (i = 0; i < NUM_POSITIVE; i++)
	{
		for (j = 0; j < NUM_NEGATIVE; j++)
		{
			
			file_out_name = string(atom_name[p[i]]) + "" + string(atom_name[n[j]])+"";
			strcpy(file_name_real, path.c_str());
			strcat(file_name_real,file_out_name.c_str());
			printf(file_name_real);
			//fout.open(file_out_name, ios::out);
			cout << "has generate the file :" << file_out_name << endl;
			//下面开始正式写结构文件了
			fout_a = fopen(file_name_real, "wb");
			if (NULL == fout_a )
			{
				cout << "wrong" << endl;
				cin.get();
			}
			fprintf(fout_a, "  2\n");
			fprintf(fout_a, "Lattice Vector\n");
			fprintf(fout_a, "   %lf    %lf    %lf\n", 2.51/ 1.73*(positive[p[i]] + negative[n[j]]), 0.0, 0.0);
			fprintf(fout_a, "   %lf    %lf    %lf\n", -1.25 / 1.73*(positive[p[i]] + negative[n[j]]),2.17/1.73*(positive[p[i]] + negative[n[j]]), 0.0);
			fprintf(fout_a, "   %lf    %lf    %lf\n", 0.0, 0.0, 25.0);
			fprintf(fout_a, "Position\n");

			//c_like
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 1.0/3, 2.0/3, 0.5);
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 2.0/3, 1.0/3, 0.5);


			//纯4
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", p[i], 0.25, 0.25, 0.0);
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", p[i], 0.75,0.75, 0.0);
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.25,0.75, 0.918);
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.75, 0.25,0.082);

			//缺陷4
		/*	fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", p[i], 0.0, 0.0, 0.5);
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", p[i], 0.5,0.5, 0.5);
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.209, 0.283, 0.472);
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.291, 0.783, 0.472);
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.709, 0.217, 0.528);
			fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.791, 0.717, 0.528);
*/


			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.709, 0.217, 0.528);
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.791, 0.717, 0.528);
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.0, 0.75, 0.329);
			//fprintf(fout_a, "%d    %lf    %lf    %lf  1  1  1\n", n[j], 0.5, 0.25, 0.671);
			
			fclose(fout_a);						/*fout<<"  5\n";
			fout<<"Lattice Vector\n";
			fout <<4.157 / 4.15*(positive[p[i]] + negative[n[j]])<<"\t"<<0.0<<"\t"<< 0.0;
			fout << endl;
			fout<<0.0<<"\t"<< 4.157 / 4.15*(positive[p[i]] + negative[n[j]])<<"\t"<< 0.0;
			fout << endl;
			fout<<0.0<<"\t"<< 0.0<<"\t"<< 25.0;
			fout << endl;
			fout<<"Position\n";
			fout << p[i]<<"\t"<< 0.33<<"\t"<< 0.66<<"\t"<< 0.427;
			fout << "1  1  1\n" << endl;
			
			fout<<p[i]<<"\t"<< 0.0<<"\t"<< 0.0<<"\t"<< 0.57;
			fout <<"\t"<< "1  1  1\n" << endl;
			fout<<n[j]<<"\t"<< 0.33<<"\t"<< 0.66<<"\t"<< 0.631;
			fout << "\t" << "1  1  1\n" << endl;

			fout<< n[j]<<"\t"<< 0.66<<"\t"<< 0.33<<"\t"<< 0.5;
			fout << "\t" << "1  1  1\n" << endl;
			fout<< n[j]<<"\t"<< 0.0<<"\t"<< 0.0<<"\t"<< 0.37;	
			fout << "\t" << "1  1  1\n" << endl;
			fout.close();*/
		}
	}

	delete[]negative;
	delete[]positive;

	return 1;

}