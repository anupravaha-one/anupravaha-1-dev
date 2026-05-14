#ifndef __APKBC_userDefinedVarCoefficients_h
#define __APKBC_userDefinedVarCoefficients_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include"APKData_Structures.h"

double BCuserfuncDensity(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncViscosity(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncThermalConductivity(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncSpecificHeat(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff1(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff2(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff3(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff4(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff5(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff6(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff7(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff8(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff9(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);
double BCuserfuncScalCoeff10(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);

double BCuserfuncChargeDiffusion(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);

double BCuserfuncCPdrop(ZoneBoundary *bcPtr,int faceCenterindex,double X,double Y,double Z,Base *base);



#endif

