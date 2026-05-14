#include "APKCoarseGridCoefficients.h"

/**
This function is used to calculate coarse grid coefficients of the convection-diffusion equation.
*/

double CoarseGridCoefficients(Base *base,Zone *zone,int g)
{
	int NI,NJ,NK,i,j,k,NIF,NJF,NKF,iF,jF,kF,indexP,*ptrSize,*ptrSizeF;
	double *diffdia,*CoeffN[18],*diffdiaF,*CoeffNF[18];
	
//	FILE *coef;
//	coef = fopen("coefficients.dat","w");

	diffdia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
	}
	ptrSize = (zone->grid+g)->GSize;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	diffdiaF = (zone->grid+g-1)->ConvDiffusionCoeffPhi_P;
	for(i=0;i<18;i++)
	{
		CoeffNF[i] = (zone->grid+g-1)->ConvDiffusionCoeffPhi_nb[i];
	}
	ptrSizeF = (zone->grid+g-1)->GSize;
	
	NIF=*(ptrSizeF+0);
	NJF=*(ptrSizeF+1);
	NKF=*(ptrSizeF+2);
	
	if(base->BaseParameters->flag2d == 1)//2D case
	{
		if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
		{
			k = 0;i = NI-2;
			for(j=0;j<NJ-1;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
							
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
		{
			k = 0;j = NJ-2;
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
							
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
							
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
				
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
							
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
		{
			k = 0;i = NI-2;
			for(j=0;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
							
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
			
			k = 0;j = NJ-2;
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
							
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
							
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
				
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
							
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
			
			k = 0;i = NI-2;j = NJ-2;
			
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
			iF = 2*i;jF = 2*j;kF = k;
			
			*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
			*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
			*(CoeffN[1]+indexP)=	0.0;
				
			*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
						
			*(CoeffN[3]+indexP)=	0.0;
				
			*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
			
			*(CoeffN[5]+indexP)=	0.0;
			
			*(CoeffN[6]+indexP)=	0.0;
			
			*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
			
			*(CoeffN[8]+indexP)=	0.0;
			
			*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
			
			*(CoeffN[10]+indexP)=	0.0;
			
			*(CoeffN[11]+indexP)=	0.0;
			
			*(CoeffN[12]+indexP)=	0.0;
			
			*(CoeffN[13]+indexP)=	0.0;
			
			*(CoeffN[14]+indexP)=	0.0;
			
			*(CoeffN[15]+indexP)=	0.0;
			
			*(CoeffN[16]+indexP)=	0.0;
			
			*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
			
	//		fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
	//			fprintf(coef,"g=%d	i=%d	j=%d	k=%d	%lf\n",g,i,j,k,*(diffdia+indexP));
			}
		}
	}
	else//3D case
	{
		if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
		{
			i = NI-2;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
		{
			j = NJ-2;
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
		{
			k = NK-2;
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-2;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
		{
			k = NK-2;
			j = NJ-2;
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			
			j = NJ-2;
			for(k=0;k<NK-2;k++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			k = NK-2;
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-2;k++)
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
		{
			k = NK-2;
			i = NI-2;
			for(j=0;j<NJ-1;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			i = NI-2;
			for(k=0;k<NK-2;k++)
			for(j=0;j<NJ-1;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			k = NK-2;
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-2;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
		{
			j = NJ-2;
			i = NI-2;
			for(k=0;k<NK-1;k++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			i = NI-2;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			j = NJ-2;
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
		{
			k = NK-2;
			j = NJ-1;
			i = NI-1;
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
			iF = 2*i;jF = 2*j;kF = 2*k;
			
			*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
			*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
			*(CoeffN[1]+indexP)=	0.0;
					
			*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
			*(CoeffN[3]+indexP)=	0.0;
					
			*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
			*(CoeffN[5]+indexP)=	0.0;
				
			*(CoeffN[6]+indexP)=	0.0;
				
			*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
			*(CoeffN[8]+indexP)=	0.0;
				
			*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
			*(CoeffN[10]+indexP)=	0.0;
				
			*(CoeffN[11]+indexP)=	0.0;
				
			*(CoeffN[12]+indexP)=	0.0;
				
			*(CoeffN[13]+indexP)=	0.0;
				
			*(CoeffN[14]+indexP)=	0.0;
				
			*(CoeffN[15]+indexP)=	0.0;
				
			*(CoeffN[16]+indexP)=	0.0;
				
			*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
//			printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			k = NK-2;i = NI-2;
			for(j=0;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
								+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
								+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
								+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			k = NK-2;j = NJ-2;
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
								+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
								+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
								+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			i = NI-2;j = NJ-2;
			for(k=0;k<NK-2;k++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			k = NK-2;
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	0.0;
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	0.0;
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	0.0;
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			j = NJ-2;
			for(k=0;k<NK-2;k++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	0.0;
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[6]+indexP)=	0.0;
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
				
				*(CoeffN[14]+indexP)=	0.0;
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	0.0;
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			i = NI-2;
			for(k=0;k<NK-2;k++)
			for(j=0;j<NJ-2;j++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	0.0;
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	0.0;
				
				*(CoeffN[11]+indexP)=	0.0;
				
				*(CoeffN[12]+indexP)=	0.0;
				
				*(CoeffN[13]+indexP)=	0.0;
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
			for(k=0;k<NK-2;k++)
			for(j=0;j<NJ-2;j++)
			for(i=0;i<NI-2;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				iF = 2*i;jF = 2*j;kF = 2*k;
				
				*(diffdia+indexP) =	*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[5]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(diffdiaF+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[4]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[0]+indexP)=	*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[0]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
					
				*(CoeffN[1]+indexP)=	*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[1]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[2]+indexP)=	*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[2]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
					
				*(CoeffN[3]+indexP)=	*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[3]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
					
				*(CoeffN[4]+indexP)=	*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[4]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[5]+indexP)=	*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF)
							+*(CoeffNF[5]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[6]+indexP)=	*(CoeffNF[6]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[6]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[7]+indexP)=	*(CoeffNF[7]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[7]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF);
				
				*(CoeffN[8]+indexP)=	*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[8]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[9]+indexP)=	*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[9]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[10]+indexP)=	*(CoeffNF[10]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[10]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[11]+indexP)=	*(CoeffNF[11]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[11]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[12]+indexP)=	*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1)+*(CoeffNF[12]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[13]+indexP)=	*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1)+*(CoeffNF[13]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[14]+indexP)=	*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF+1)+*(CoeffNF[14]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*(jF+1)+iF);
				
				*(CoeffN[15]+indexP)=	*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF)+*(CoeffNF[15]+(NIF-1)*(NJF-1)*(kF+1)+(NIF-1)*jF+iF+1);
				
				*(CoeffN[16]+indexP)=	*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF)+*(CoeffNF[16]+(NIF-1)*(NJF-1)*kF+(NIF-1)*(jF+1)+iF+1);
				
				*(CoeffN[17]+indexP)=	*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF)+*(CoeffNF[17]+(NIF-1)*(NJF-1)*kF+(NIF-1)*jF+iF+1);
		//		printf("%d	%lf	%lf\n",g,*(CoeffN[0]+indexP),*(CoeffN[1]+indexP));
			}
		}
	}
//	fclose(coef);
	return(0);
}
