#include "APKCallCalcCoefficients.h"

/**
This function calls CalcCoefficients2d,CalcCoefficients,CalcCoefficientsCorner2d,CalcCoefficientsCorner which calculate
the coefficients of pressure Poisson equation which when multiplied by diffusion coeffient (e.g., kinematic viscosity),
can be used in convection-diffusion equation as diffusion term coefficients but only if it is a constant diffusion
coefficients case.
*/
double CallCalcCoefficients(Base *base,Zone *zone,int g)
{
	int NI,NJ,NK,i,j,k,index[3],indexP,*ptrSize,flagrho;
	double *diffdia,*CoeffN[18],*slabVolume[6];
	double *beta,*volumeDomain;
	double Volg,Vol[18],rho;
	
//	FILE *coef;
//	coef = fopen("coefficients_fine.dat","w");
	
	diffdia = (zone->grid+g)->DiffusionCoeff_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid+g)->DiffusionCoeff_nb[i];
	}
	ptrSize = zone->Size ;
	beta=zone->ZGeom->Beta;
	volumeDomain = zone->ZGeom->Volume;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	flagrho = zone->Flowsol->flagRho;
	
	if(base->BaseParameters->flag2d == 1)
	{
//--------------------------------------------- Boundary cells ----------------------------------(1)---------
		k = 0; j = 0; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
		Vol[7] = 0.0;
		Vol[8] = 0.0;
		Vol[9] = 0.0;
		Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
		Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
		Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
		Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
		Vol[15] = 0.0;
		Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
		Vol[17] = 0.0;
		
		*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
		
//		fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
//---------------------------------------------------------------------------------------------(2)---------
		k = 0; j = 0;
		
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
			Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
			Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
			Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
			Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = 0.0;
			Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
			Vol[17] = 0.0;
			
			*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
			
//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
		}
//-------------------------------------------------------------------------------------------------(3)---------
		k = 0; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
			Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
			Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
			Vol[14] = 0.0;
			Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
			Vol[16] = 0.0;
			Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
			
//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
		}
//---------------------------------------------------------------------------------------------------(4)-------
		k = 0; j = NJ-2; i = NI-2; 
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
		Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
		Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
		Vol[10] = 0.0;
		Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
		Vol[12] = 0.0;
		Vol[13] = 0.0;
		Vol[14] = 0.0;
		Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
		Vol[16] = 0.0;
		Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
		
//		fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
//----------------------------------------------------------------------------------------------------(5)------
		k = 0; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
		Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
		Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
		Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
		Vol[11] = 0.0;
		Vol[12] = 0.0;
		Vol[13] = 0.0;
		Vol[14] = *(slabVolume[4]+(NI-1)*(j-1)+i);
		Vol[15] = 0.0;
		Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
		Vol[17] = 0.0;
		
		*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
		
//		fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
//------------------------------------------------------------------------------------------------------(6)-------
		k = 0; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
			Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
			Vol[8] = 0.0;
			Vol[9] = 0.0;
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
			Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
			Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
			
//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
		}
//---------------------------------------------------------------------------------------------------------(7)-----
		k = 0; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
			Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
			Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
			Vol[12] = 0.0;
			Vol[13] = 0.0;
			Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
			
//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
		}
//-------------------------------------------------------------------------------------------------------(8)---
		k = 0; i = 0; j = NJ-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = 0.0;
		Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
		Vol[8] = 0.0;
		Vol[9] = 0.0;
		Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
		Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
		Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
		Vol[14] = 0.0;
		Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
		Vol[16] = 0.0;
		Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner2d(index,Volg,Vol,ptrSize,beta,CoeffN);
		
//		fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
//--------------------------------------------Interior Cells------------------------------------------------------------(9)------
		k = 0;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
				Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
				Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
				Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
				Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
				Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
				Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
				Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
				Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
				Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
				Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
				Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
				Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
				Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients2d(index,Volg,Vol,ptrSize,beta,CoeffN);
				
//				fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
		}
//---------------------------------------End of Interior Cells------------------------------------------
	}
	else if(base->BaseParameters->flag2d == 0)
	{
//--------------------------------------------- Boundary cells ----------------------------------(1)---------
		k = 0; j = 0; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
		Vol[7] = 0.0;
		Vol[8] = *(slabVolume[1]+(NJ-1)*(k+1)+j);
		Vol[9] = 0.0;
		Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
		Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
		Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		Vol[15] = *(slabVolume[3]+(NI-1)*(k+1)+i);
		Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
		Vol[17] = 0.0;
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//---------------------------------------------------------------------------------------------(2)---------
		k = 0; j = 0;
		
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
			Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
			Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[3]+(NI-1)*(k+1)+i);
			Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
			Vol[17] = 0.0;
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//-------------------------------------------------------------------------------------------------(3)---------
		k = 0; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
			Vol[14] = *(slabVolume[2]+(NI-1)*(k+1)+i);
			Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			Vol[16] = 0.0;
			Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//---------------------------------------------------------------------------------------------------(4)-------
		j = NJ-2;
		for(k=1;k<NK-2;k++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
				Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
				Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
				Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
				Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
				Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
				Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
				Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
				Vol[14] = *(slabVolume[2]+(NI-1)*(k+1)+i);
				Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
				Vol[16] = *(slabVolume[2]+(NI-1)*(k-1)+i);
				Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			}
		}
//--------------------------------------------------------------------------------------------------(5)--------
		k = 0; j = NJ-2; i = NI-2; 
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
		Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
		Vol[10] = 0.0;
		Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
		Vol[12] = *(slabVolume[0]+(NJ-1)*(k+1)+j);
		Vol[13] = 0.0;
		Vol[14] = *(slabVolume[2]+(NI-1)*(k+1)+i);
		Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		Vol[16] = 0.0;
		Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//----------------------------------------------------------------------------------------------------(6)------
		j = NJ-2; i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			Vol[10] = 0.0;
			Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
			Vol[12] = *(slabVolume[0]+(NJ-1)*(k+1)+j);
			Vol[13] = *(slabVolume[0]+(NJ-1)*(k-1)+j);
			Vol[14] = *(slabVolume[2]+(NI-1)*(k+1)+i);
			Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[2]+(NI-1)*(k-1)+i);
			Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//------------------------------------------------------------------------------------------------------(7)-------
		k = NK-2; j = 0;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
			Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
			Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
			Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
			Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
			Vol[15] = 0.0;
			Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[3]+(NI-1)*(k-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//------------------------------------------------------------------------------------------------------(8)-------
		k = 0; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
		Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
		Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
		Vol[11] = 0.0;
		Vol[12] = *(slabVolume[0]+(NJ-1)*(k+1)+j);
		Vol[13] = 0.0;
		Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		Vol[15] = *(slabVolume[3]+(NI-1)*(k+1)+i);
		Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
		Vol[17] = 0.0;
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);

//------------------------------------------------------------------------------------------------------(9)-------
		j = 0; i = 0;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
			Vol[7] = 0.0;
			Vol[8] = *(slabVolume[1]+(NJ-1)*(k+1)+j);
			Vol[9] = *(slabVolume[1]+(NJ-1)*(k-1)+j);
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
			Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[3]+(NI-1)*(k+1)+i);
			Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[3]+(NI-1)*(k-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//-----------------------------------------------------------------------------------------------------(10)--------
		j = 0; i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
			Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
			Vol[11] = 0.0;
			Vol[12] = *(slabVolume[0]+(NJ-1)*(k+1)+j);
			Vol[13] = *(slabVolume[0]+(NJ-1)*(k-1)+j);
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[3]+(NI-1)*(k+1)+i);
			Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[3]+(NI-1)*(k-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//-----------------------------------------------------------------------------------------------------(11)-----
		j = NJ-2; i = 0;
		for(k=1;k<NK-2;k++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = 0.0;
			Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
			Vol[8] = *(slabVolume[1]+(NJ-1)*(k+1)+j);
			Vol[9] = *(slabVolume[1]+(NJ-1)*(k-1)+j);
			Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			Vol[14] = *(slabVolume[2]+(NI-1)*(k+1)+i);
			Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[2]+(NI-1)*(k-1)+i);
			Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//------------------------------------------------------------------------------------------------------(12)-------
		j = 0; i = 0; k = NK-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
		Vol[7] = 0.0;
		Vol[8] = 0.0;
		Vol[9] = *(slabVolume[1]+(NJ-1)*(k-1)+j);
		Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
		Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
		Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
		Vol[15] = 0.0;
		Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		Vol[17] = *(slabVolume[3]+(NI-1)*(k-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//------------------------------------------------------------------------------------------------------(13)----
		k = 0; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
			Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
			Vol[8] = *(slabVolume[1]+(NJ-1)*(k+1)+j);
			Vol[9] = 0.0;
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//---------------------------------------------------------------------------------------------------------(14)-----
		k = NK-2; i = 0;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
			Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
			Vol[8] = 0.0;
			Vol[9] = *(slabVolume[1]+(NJ-1)*(k-1)+j);
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
			Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
			Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//---------------------------------------------------------------------------------------------------------(15)-----
		k = 0; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
			Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
			Vol[12] = *(slabVolume[0]+(NJ-1)*(k+1)+j);
			Vol[13] = 0.0;
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
			Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//-------------------------------------------------------------------------------------------------------(16)---
		k = 0; i = 0; j = NJ-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
		Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		Vol[6] = 0.0;
		Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
		Vol[8] = *(slabVolume[1]+(NJ-1)*(k+1)+j);
		Vol[9] = 0.0;
		Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
		Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
		Vol[14] = *(slabVolume[2]+(NI-1)*(k+1)+i);
		Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		Vol[16] = 0.0;
		Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//--------------------------------------------------------------------------------------------------------(17)------
		j = 0;
		for(k=1;k<NK-2;k++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
				Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
				Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
				Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
				Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
				Vol[11] = *(slabVolume[3]+(NI-1)*k+i+1);
				Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
				Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
				Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
				Vol[15] = *(slabVolume[3]+(NI-1)*(k+1)+i);
				Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
				Vol[17] = *(slabVolume[3]+(NI-1)*(k-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			}
		}
//------------------------------------------------------------------------------------------------------(18)------
		k = 0;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(slabVolume[5]+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
				Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
				Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
				Vol[9] = *(slabVolume[5]+(NI-1)*j+i-1);
				Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
				Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
				Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
				Vol[13] = *(slabVolume[5]+(NI-1)*j+i+1);
				Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
				Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
				Vol[16] = *(slabVolume[5]+(NI-1)*(j+1)+i);
				Vol[17] = *(slabVolume[5]+(NI-1)*(j-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			}
		}
//-------------------------------------------------------------------------------------------------------------(19)---
		i=0;
		for(k=1;k<NK-2;k++)
		{
			for(j=1;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				Vol[6] = *(slabVolume[1]+(NJ-1)*k+j+1);
				Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
				Vol[8] = *(slabVolume[1]+(NJ-1)*(k+1)+j);
				Vol[9] = *(slabVolume[1]+(NJ-1)*(k-1)+j);
				Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
				Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
				Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
				Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
				Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
				Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
				Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
				Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			}
		}
//----------------------------------------------------------------------------------------------------------(20)-----
		k = NK-2; j = NJ-2; i = 0;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(slabVolume[1]+(NJ-1)*k+j);
		Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = 0.0;
		Vol[7] = *(slabVolume[1]+(NJ-1)*k+j-1);
		Vol[8] = 0.0;
		Vol[9] = *(slabVolume[1]+(NJ-1)*(k-1)+j);
		Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
		Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
		Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		Vol[14] = 0.0;
		Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
		Vol[16] = *(slabVolume[2]+(NI-1)*(k-1)+i);
		Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//--------------------------------------------------------------------------------------------------------(21)-----
		k = NK-2; j = NJ-2;
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
			Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[2]+(NI-1)*k+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
			Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			Vol[14] = 0.0;
			Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
			Vol[16] = *(slabVolume[2]+(NI-1)*(k-1)+i);
			Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//-------------------------------------------------------------------------------------------------------(22)------
		k = NK-2;
		for(j=1;j<NJ-2;j++)
		{
			for(i=1;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
				Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
				Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
				Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
				Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
				Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
				Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
				Vol[12] = *(slabVolume[4]+(NI-1)*j+i+1);
				Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
				Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
				Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
				Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
				Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			}
		}
//-------------------------------------------------------------------------------------------------------(23)------
		k = NK-2; j = 0; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(slabVolume[3]+(NI-1)*k+i);
		Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		Vol[7] = *(slabVolume[3]+(NI-1)*k+i-1);
		Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
		Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
		Vol[11] = 0.0;
		Vol[12] = 0.0;
		Vol[13] = *(slabVolume[0]+(NJ-1)*(k-1)+j);
		Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
		Vol[15] = 0.0;
		Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		Vol[17] = *(slabVolume[3]+(NI-1)*(k-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//---------------------------------------------------------------------------------------------------------(24)-----
		k = NK-2; i = NI-2;
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
			Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
			Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
			Vol[12] = 0.0;
			Vol[13] = *(slabVolume[0]+(NJ-1)*(k-1)+j);
			Vol[14] = *(slabVolume[4]+(NI-1)*(j+1)+i);
			Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
			Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
		}
//---------------------------------------------------------------------------------------------------------(25)-----
		i = NI-2;
		for(k=1;k<NK-2;k++)
		{
			for(j=1;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
				Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
				Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
				Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
				Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
				Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
				Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
				Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
				Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
				Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
				Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
				Vol[10] = *(slabVolume[0]+(NJ-1)*k+j+1);
				Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
				Vol[12] = *(slabVolume[0]+(NJ-1)*(k+1)+j);
				Vol[13] = *(slabVolume[0]+(NJ-1)*(k-1)+j);
				Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
				Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
				Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
				Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
				
				*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			}
		}
//---------------------------------------------------------------------------------------------------------(26)-------
		k = NK-2; j = NJ-2; i = NI-2;
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		Vol[1] = *(slabVolume[0]+(NJ-1)*k+j);
		Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		Vol[3] = *(slabVolume[2]+(NI-1)*k+i);
		Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
		Vol[6] = *(slabVolume[2]+(NI-1)*k+i-1);
		Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		Vol[8] = *(slabVolume[4]+(NI-1)*j+i-1);
		Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		Vol[10] = 0.0;
		Vol[11] = *(slabVolume[0]+(NJ-1)*k+j-1);
		Vol[12] = 0.0;
		Vol[13] = *(slabVolume[0]+(NJ-1)*(k-1)+j);
		Vol[14] = 0.0;
		Vol[15] = *(slabVolume[4]+(NI-1)*(j-1)+i);
		Vol[16] = *(slabVolume[2]+(NI-1)*(k-1)+i);
		Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		*(diffdia+indexP)  = CalcCoefficientsCorner(index,Volg,Vol,ptrSize,beta,CoeffN);
//------------------------------------------- End of Boundary cells ----------------------------------------
	
//------------------------------------------- Interior cells ----------------------------------------------(27)------
		for(k=1;k<NK-2;k++)
		for(j=1;j<NJ-2;j++)
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Volg = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			Vol[0] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			Vol[1] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			Vol[2] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			Vol[3] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			Vol[4] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			Vol[5] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			Vol[6] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			Vol[7] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			Vol[8] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			Vol[9] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			Vol[10] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			Vol[11] = *(volumeDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			Vol[12] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			Vol[13] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			Vol[14] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			Vol[15] = *(volumeDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			Vol[16] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			Vol[17] = *(volumeDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			*(diffdia+indexP)  = CalcCoefficients(index,Volg,Vol,ptrSize,beta,CoeffN);
			
		}
//---------------------------------------End of Interior Cells------------------------------------------
	}
	
//	fclose(coef);
	return(0);
}

