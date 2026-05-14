#ifndef __APKuserDefinedVarCoefficients_h
#define __APKuserDefinedVarCoefficients_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

double userEnterconstant(Zone *z,int cellindex,Base *base);
double userfuncDensity(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncViscosity(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncThermalConductivity(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncSpecificHeat(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff1(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff2(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff3(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff4(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff5(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff6(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff7(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff8(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff9(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userfuncScalCoeff10(Zone *z,int cellindex,double X,double Y,double Z,Base *base);

double userfuncChargeDiffusion(Zone *z,int cellindex,double X,double Y,double Z,Base *base);

double userfuncCPdrop(Zone *z,int cellindex,double X,double Y,double Z,Base *base);
double userDefinedCD(double Re);
double userDefinedCDDisperse(double Re);
double userDefinedDropPfYf(double Tsurf,double P);
double userDefinedDropCpa(double Tref);
double userDefinedDropCpv(double Tref);
double userDefinedDropKa(double Tref);
double userDefinedDropKv(double Tref,double Cpv);
double userDefinedDropL(double Tsurf);
double userDefinedDropCpf(double Tsurf);
double userDefinedDropRhoF(double Tsurf);
#endif

