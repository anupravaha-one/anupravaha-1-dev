#ifndef __APKuser_defined_funcs2_h
#define __APKuser_defined_funcs2_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include"APKData_Structures.h"

double userfunc_Energysource(double X,double Y,double Z,double time,Zone *zone);
double user_func_14(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc);
double user_func_15(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc);
double user_func_16(double X,double Y,double Z,double time,Zone *zone,int i,int j,int k,ZoneBoundary *zonebc);

double CalcLagrangeCoeff(double X, double XC, double B, double C);

#endif

