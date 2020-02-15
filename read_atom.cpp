//
//  read_atom.cpp
//  WKM_module
//
//  Created by 翁谋毅 on 2018/12/4.
//  Copyright © 2018 翁谋毅. All rights reserved.
//

#define BOHR 0.529177
#define PI  3.1415926535897932

#include "read_atom.hpp"
#include "tools.hpp"
#include <iostream>
#include <cstdlib>
const int yanshen=3;
const int cengshu=cengshu*cengshu*cengshu;
#pragma warning(disable : 4996)




void atom::read_atom(char *name)
{
	FILE *in;
	char temp[300];
	int i, j;
	char str[500];
	double det;
	int flag;
	in = fopen(name, "r");
	if (in == NULL)
	{
		flag_read = 0;
		printf("Do not have file %s\n", name);
		return;
	}
	flag_read = 1;
	fscanf(in, "%d", &num_atom);
	d = (double **)malloc(3 * sizeof(double *));
	rd = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		d[i] = (double *)malloc(3 * sizeof(double));
		rd[i] = (double *)malloc(3 * sizeof(double));
	}
	atom_type = (int *)malloc(num_atom * sizeof(int));
	p = (double **)malloc(num_atom * sizeof(double *));
	for (i = 0; i < num_atom; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	rp = (double***)malloc(yanshen * sizeof(double**));
	for (j = 0; j < yanshen; j++)
	{
		rp[j] = (double **)malloc(num_atom * sizeof(double *));
		for (i = 0; i < num_atom; i++)
		{
			rp[j][i] = (double *)malloc(3 * sizeof(double));
		}
	}
	
	while (fgets(temp, 300, in) != NULL)
	{
		if (strstr(temp, "vector") != NULL)
			break;
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &d[i][j]);
		}
		fgets(str, 500, in);
	}
	
	fgets(str, 500, in);
	//fgets(str, 500, in);
	for (i = 0; i < num_atom; i++)
	{
		fscanf(in, "%d", &atom_type[i]);
		//printf("%d\n", atom_type[i]);
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &p[i][j]);
		}
		fgets(str, 500, in);
	}
	type_num = 0;


	for (i = 0; i < 30; i++)
	{
		type_name[i] = -1;
	}

	for (i = 0; i < num_atom; i++)
	{
		flag = 0;
		for (j = 0; j < type_num; j++)
		{
			if (atom_type[i] == type_name[j])
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			type_name[type_num] = atom_type[i];
			type_num++;
		}
	}

	for (i = 0; i < yanshen; i++)
	{

		for (j = 0; j < num_atom; j++)
		{
			//x_xishu = i/3;
			//y_xishu = i / 3;
			//z_zishu = (i % 9) / 3;
			rp[i][j][0] = d[0][0] * p[j][0] + d[1][0] * p[j][1] + d[2][0] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * d[0][0] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * d[1][0] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * d[2][0];
			rp[i][j][1] = d[0][1] * p[j][0] + d[1][1] * p[j][1] + d[2][1] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * d[0][1] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * d[1][1] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * d[2][1];
			rp[i][j][2] = d[0][2] * p[j][0] + d[1][2] * p[j][1] + d[2][2] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * d[0][2] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * d[1][2] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * d[2][2];
		}
	}
	det = d[0][0] * d[1][1] * d[2][2] + d[0][1] * d[1][2] * d[2][0] + d[0][2] * d[1][0] * d[2][1] - d[0][0] * d[1][2] * d[2][1] - d[0][1] * d[1][0] * d[2][2] - d[0][2] * d[1][1] * d[2][0];
	rd[0][0] = (d[1][1] * d[2][2] - d[1][2] * d[2][1]) / det * 2 * PI*BOHR;
	rd[0][1] = (d[1][2] * d[2][0] - d[1][0] * d[2][2]) / det * 2 * PI*BOHR;
	rd[0][2] = (d[1][0] * d[2][1] - d[2][0] * d[1][1]) / det * 2 * PI*BOHR;

	rd[1][0] = (d[2][1] * d[0][2] - d[2][2] * d[0][1]) / det * 2 * PI*BOHR;
	rd[1][1] = (d[2][2] * d[0][0] - d[2][0] * d[0][2]) / det * 2 * PI*BOHR;
	rd[1][2] = (d[2][0] * d[0][1] - d[2][1] * d[0][0]) / det * 2 * PI*BOHR;

	rd[2][0] = (d[0][1] * d[1][2] - d[1][1] * d[0][2]) / det * 2 * PI*BOHR;
	rd[2][1] = (d[0][2] * d[1][0] - d[1][2] * d[0][0]) / det * 2 * PI*BOHR;
	rd[2][2] = (d[0][0] * d[1][1] - d[1][0] * d[0][1]) / det * 2 * PI*BOHR;
	fclose(in);
	return;
}

void atom::print_name(int i, FILE *in)
{
	fprintf(in, "%s", atom_name[i]);
	return;
}

void atom::freeatom() {
	int i;
	if (flag_read == 1)
	{
		for (i = 0; i < 3; i++)
		{
			free(d[i]);
			free(rd[i]);
		}
		free(d);
		free(rd);
		free(atom_type);
		for (i = 0; i < num_atom; i++)
		{
			free(p[i]);
		}
		free(p);
		for (int j = 0; j < yanshen; j++)
		{
			for (i = 0; i < num_atom; i++)
			{
				free(rp[j][i]);
			}
			free(rp[j]);
		}		
		free(rp);
	}
}

void cluster::make_cluster(int n, atom aa, double **dist)
{
    int i,j,k;
    int natom;
    int ii;
    double rp1[3];
    cluster *pointer1, *pointer2;
    num_cnnct=0;
    int m=0;
    cluster *temp, *temp1;;
    flag_find=1;
    for (i=nx-2;i<nx+3;i++)
    {
        for (j=ny-2;j<ny+3;j++)
        {
            for (k=nz-2;k<nz+3;k++)
            {
                for(natom=0;natom<aa.num_atom;natom++)
                {
                    rp1[0]=aa.rp[(aa.yanshen-1)/2][natom][0]+i*aa.d[0][0]+j*aa.d[1][0]+k*aa.d[2][0];
                    rp1[1]=aa.rp[(aa.yanshen - 1) / 2][natom][1]+i*aa.d[0][1]+j*aa.d[1][1]+k*aa.d[2][1];
                    rp1[2]=aa.rp[(aa.yanshen - 1) / 2][natom][2]+i*aa.d[0][2]+j*aa.d[1][2]+k*aa.d[2][2];
                    if (dis(rp1, rp)<dist[ntype][aa.atom_type[natom]])
                    {
                        num_cnnct++;
                        if (find_last(i, j, k, natom)==0)
                        {
                            temp=mysf;
                            //temp1=next;
                            while (temp->next!=NULL) {
                                temp=temp->next;
                            }
                            temp->next=(cluster *)malloc(sizeof(cluster));
                            temp->next->last=temp;
                            temp->next->next=NULL;
                            temp->next->mysf=temp->next;
                            temp=temp->next;
                            if (m==0)
                            {
                                //next=temp;
                                //temp->last=mysf;
                                //temp->next=NULL;
                                temp->nx=i;
                                temp->ny=j;
                                temp->nz=k;
                                temp->nn=natom;
                                temp->level=level+1;
                                temp->num=num+1;
                                temp->ntype=aa.atom_type[natom];
                                temp->rp[0]=rp1[0];
                                temp->rp[1]=rp1[1];
                                temp->rp[2]=rp1[2];
                                temp->flag_find=0;
                            }
                            else
                            {
                                //temp1->next=temp;
                                //temp->last=temp1;
                                //temp->next=NULL;
                                temp->nx=i;
                                temp->ny=j;
                                temp->nz=k;
                                temp->nn=natom;
                                temp->level=level+1;
                                temp->num=temp->last->num+1;
                                temp->ntype=aa.atom_type[natom];
                                temp->rp[0]=rp1[0];
                                temp->rp[1]=rp1[1];
                                temp->rp[2]=rp1[2];
                                temp->flag_find=0;
                            }
                            m++;
                        }
                    }
                }
            }
        }
    }
    cnnct=(cluster **)malloc(num_cnnct*sizeof(cluster *));
    pointer1=last;
    pointer2=next;
    m=0;
    while (pointer1!=NULL)
    {
        if (dis(pointer1->rp, rp)<dist[ntype][pointer1->ntype])
        {
            cnnct[m]=pointer1;
            m++;
        }
        pointer1=pointer1->last;
    }
    while (pointer2!=NULL)
    {
        if (dis(pointer2->rp, rp)<dist[ntype][pointer2->ntype])
        {
            cnnct[m]=pointer2;
            m++;
        }
        pointer2=pointer2->next;
    }
    num_cnnct=num_cnnct-1;
    pointer2=next;
    while (pointer2!=NULL)
    {
        if (pointer2->flag_find==0 && pointer2->level < n)
        {
            pointer2->make_cluster(n, aa, dist);
        }
        pointer2=pointer2->next;
    }
    return;
}

void cluster::make_cnnct(double **dist)
{
    cluster *p, *q;
    p=mysf;
    int count=0;
    while (p->last!=NULL) {
        p=p->last;
    }
    while (p!=NULL) {
        p->num=count;
        count++;
        p->num_cnnct=0;
        p=p->next;
    }
    p=mysf;
    while (p->last!=NULL) {
        p=p->last;
    }
    
    while (p!=NULL)
    {
        q=mysf;
        while (q->last!=NULL) {
            q=q->last;
        }
        while (q!=NULL)
        {
            if (dis(p->rp, q->rp)<dist[p->ntype][q->ntype])
            {
                p->num_cnnct=p->num_cnnct+1;
            }
            q=q->next;
        }
        p=p->next;
    }
    p=mysf;
    while (p->last!=NULL) {
        p=p->last;
    }
    while (p!=NULL) {
        p->cnnct=(cluster **)malloc(p->num_cnnct*sizeof(cluster *));
        p=p->next;
    }
    
    p=mysf;
    while (p->last!=NULL) {
        p=p->last;
    }
    while (p!=NULL)
    {
        count=0;
        q=mysf;
        while (q->last!=NULL)
        {
            q=q->last;
        }
        while (q!=NULL)
        {
            if (dis(p->rp, q->rp)<dist[p->ntype][q->ntype])
            {
                p->cnnct[count]=q;
                count++;
            }
            q=q->next;
        }
        if (count!=p->num_cnnct)
        {
            printf("ERROR !\n");
        }
        p=p->next;
    }
    return;
}

int cluster::find_last(int n1, int n2, int n3, int natom)
{
    cluster *pointer1, *pointer2;
    pointer1=last;
    pointer2=next;
    if (nx==n1 && ny==n2 && nz==n3 && nn == natom)
    {
        return 1;
    }
    while (pointer1!=NULL)
    {
        if (pointer1->nx==n1 &&pointer1->ny==n2 && pointer1->nz==n3 && pointer1->nn == natom)
        {
            if (pointer1->level > level+1)
            {
                pointer1->level = level+1;
                pointer1->flag_find=0;
            }
            return 1;
        }
        pointer1=pointer1->last;
    }
    while (pointer2!=NULL)
    {
        if (pointer2->nx==n1 &&pointer2->ny==n2 && pointer2->nz==n3 && pointer2->nn == natom)
        {
            if (pointer2->level > level+1)
            {
                pointer2->level = level+1;
                pointer2->flag_find=0;
            }
            return 1;
        }
        pointer2=pointer2->next;
    }
    return 0;
}


