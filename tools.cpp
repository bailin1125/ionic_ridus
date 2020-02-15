//
//  tools.cpp
//  wannier_band
//
//  Created by 翁谋毅 on 2018/12/1.
//  Copyright © 2018 翁谋毅. All rights reserved.
//

#include "tools.hpp"
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
double det(double **a, int n)
{
    int i,j,k;
    int ii,jj;
    int a_i,a_j;
    a_i=-1;
    a_j=-1;
    double m;
    m=0;
    double **b;
    int f;
    b=(double **)malloc ((n-1)*sizeof(double*));
    for (i=0;i<n-1;i++)
    {
        b[i]=(double *)malloc((n-1)*sizeof(double));
    }
    if (n==1)
    {
        m=a[0][0];
    }
    else
    {
        j=0;
        for(i=0;i<n;i++)
        {
            f=1;
            for (k=0;k<i+j;k++)
            {
                f=f*(-1);
            }
            for (ii=0;ii<n;ii++)
            {
                for (jj=0;jj<n;jj++)
                {
                    if (ii<i)
                    {
                        a_i=ii;
                    }
                    if (ii>i)
                    {
                        a_i=ii-1;
                    }
                    if (jj<j)
                    {
                        a_j=jj;
                    }
                    if (jj>j)
                    {
                        a_j=jj-1;
                    }
                    if (ii!=i && jj!=j)
                    {
                        b[a_i][a_j]=a[ii][jj];
                    }
                }
            }
            m=m+a[i][j]*det(b, n-1)*f;
        }
    }
    for (i=0;i<n-1;i++)
    {
        free(b[i]);
    }
    free(b);
    return m;
}
void inv(int n, double **a, double **b)
{
    int i,j,k;
    int ii,jj;
    int a_i,a_j;
    a_i=-1;
    a_j=-1;
    double det_a;
    double **c;
    int f;
    
    c=(double **)malloc ((n-1)*sizeof(double*));
    for (i=0;i<n-1;i++)
    {
        c[i]=(double *)malloc((n-1)*sizeof(double));
    }
    
    det_a=det(a, n);
    if (fabs(det_a)<0.000000001)
    {
        for (i=0;i<n;i++)
        {
            for (j=0;j<n;j++)
            {
                b[i][j]=0;
            }
        }
    }
    else
    {
        for (i=0;i<n;i++)
        {
            for (j=0;j<n;j++)
            {
                f=1;
                for (k=0;k<i+j;k++)
                {
                    f=f*(-1);
                }
                for (ii=0;ii<n;ii++)
                {
                    for (jj=0;jj<n;jj++)
                    {
                        if (ii<i)
                        {
                            a_i=ii;
                        }
                        if (ii>i)
                        {
                            a_i=ii-1;
                        }
                        if (jj<j)
                        {
                            a_j=jj;
                        }
                        if (jj>j)
                        {
                            a_j=jj-1;
                        }
                        if (ii!=i && jj!=j)
                        {
                            c[a_i][a_j]=a[ii][jj];
                        }
                    }
                }
                b[j][i]=det(c, n-1)/det_a*f;
            }
        }
    }
    for(i=0;i<n-1;i++)
    {
        free(c[i]);
    }
    free(c);
    return;
}

void trans(int ax, int ay, double **a, double **b)
{
    int i,j;
    for (i=0;i<ay;i++)
    {
        for (j=0;j<ax;j++)
        {
            b[j][i]=a[i][j];
        }
    }
    return;
}

void cal_rp(double **d, double x, double y, double z, double *rp)
{
    rp[0]=x*d[0][0]+y*d[1][0]+z*d[2][0];
    rp[1]=x*d[0][1]+y*d[1][1]+z*d[2][1];
    rp[2]=x*d[0][2]+y*d[1][2]+z*d[2][2];
    return ;
}

double dot_product(double ***a, double ***b, int nx, int ny, int nz)
{
    int i,j,k;
    double sum=0;
    for (i=0;i<nx;i++)
    {
        for (j=0;j<ny;j++)
        {
            for (k=0;k<nz;k++)
            {
                sum=sum+a[i][j][k]*b[i][j][k];
            }
        }
    }
    return sum;
}
double dot_product(double ***a_r, double ***a_i, double ***b_r, double ***b_i, int nx, int ny, int nz)
{
    int i,j,k;
    double sum=0;
    for (i=0;i<nx;i++)
    {
        for (j=0;j<ny;j++)
        {
            for (k=0;k<nz;k++)
            {
                sum=sum+a_r[i][j][k]*b_r[i][j][k]+a_i[i][j][k]*b_i[i][j][k];
            }
        }
    }
    return sum;
}

double dis(double *a, double *b)
{
    return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]));
}

