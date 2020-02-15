//
//  tools.hpp
//  wannier_band
//
//  Created by 翁谋毅 on 2018/12/1.
//  Copyright © 2018 翁谋毅. All rights reserved.
//

#ifndef tools_hpp
#define tools_hpp

#include <stdio.h>

double det(double **a, int n);
void inv(int n, double **a, double **b);
void trans(int ax, int ay, double **a, double **b);
void cal_rp(double **d, double x, double y, double z, double *rp);
double dot_product(double ***a, double ***b, int nx, int ny, int nz);
double dot_product(double ***a_r, double ***a_i, double ***b_r, double ***b_i, int nx, int ny, int nz);
double dis(double *a, double *b);

#endif /* tools_hpp */
