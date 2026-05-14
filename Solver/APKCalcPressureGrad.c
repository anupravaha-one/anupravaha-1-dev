#include "APKCalcPressureGrad.h"
/**
The function CalcPressureGrad() is used to compute the pressure gradient term in the Navier Stokes equation.
This function takes the individual velocity component's variable tag as input and calculates the coresponding
component of Pressure gradient. This function is also used for calculating the pressure gradient term while
solving the massvelocities also.
Note: This function ALSO includes the "negative sign" involved in the pressure gradient term occuring in the 
Navier Stokes equations.
*/
double CalcPressureGrad(double *variable,double *slabVariable[6],double *volume,double *slabVolume[6],double *SurfArea,int index[3],int *ptrSize,int variableTag)
{
	
	int k,j,i,comp,NK,NJ,NI,indexP,indexE,indexW,indexN,indexS,indexB,indexT;
	double sP,sW,sE,sS,sN,sB,sT;
	double volP,volW,volE,volS,volN,volB,volT;
	double s_PE,s_PW,s_PS,s_PN,s_PB,s_PT;
	double source[6][3],sumSource;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	k=index[0],j=index[1],i=index[2];
	indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
	indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
	indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
	indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
	indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
	indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
	
	sP=*(variable+indexP);
	volP=*(volume+indexP);
	
	if(i==(NI-2))
	{
		sE=*(slabVariable[0]+(NJ-1)*k+j);
		volE=*(slabVolume[0]+(NJ-1)*k+j);
	}
	else
	{
		sE=*(variable+indexE);
		volE=*(volume+indexE);
	}
	
	if(i==0)
	{
		sW =*(slabVariable[1]+(NJ-1)*k+j);
		volW =*(slabVolume[1]+(NJ-1)*k+j);
	
	}
	else
	{
		sW=*(variable+indexW);
		volW=*(volume+indexW);
	}
	
	if(j==(NJ-2))
	{
		sN=*(slabVariable[2]+(NI-1)*k+i);
		volN=*(slabVolume[2]+(NI-1)*k+i);
	}
	else
	{
		sN=*(variable+indexN);
		volN=*(volume+indexN);
	}
	
	if(j==0)
	{
		sS=*(slabVariable[3]+(NI-1)*k+i);
		volS=*(slabVolume[3]+(NI-1)*k+i);
	}
	else
	{
		sS=*(variable+indexS);
		volS=*(volume+indexS);
	}
	
	if(k==(NK-2))
	{
		sT=*(slabVariable[4]+(NI-1)*j+i);
		volT=*(slabVolume[4]+(NI-1)*j+i);
	}
	else
	{
		sT=*(variable+indexT);
		volT=*(volume+indexT);
	}
	
	if(k==0)
	{
		sB=*(slabVariable[5]+(NI-1)*j+i);
		volB=*(slabVolume[5]+(NI-1)*j+i);
	}
	else
	{
		sB=*(variable+indexB);
		volB=*(volume+indexB);
	}
	
	s_PE = (sP*volE + sE*volP )/( volP + volE);
	s_PW = (sP*volW + sW*volP )/( volP + volW);
	s_PN = (sP*volN + sN*volP )/( volP + volN);
	s_PS = (sP*volS + sS*volP )/( volP + volS);
	s_PT = (sP*volT + sT*volP )/( volP + volT);
	s_PB = (sP*volB + sB*volP )/( volP + volB);

	
	switch(variableTag)
	{
		case 1:
			
			comp=0;
			source[0][comp]=s_PW*(*(SurfArea+indexP*6*3+3*0+comp));//West direction, x-component
			source[1][comp]=s_PE*(*(SurfArea+indexP*6*3+3*1+comp));//East direction, x-component
			source[2][comp]=s_PS*(*(SurfArea+indexP*6*3+3*2+comp));//South direction, x-component
			source[3][comp]=s_PN*(*(SurfArea+indexP*6*3+3*3+comp));//North direction, x-component
			source[4][comp]=s_PB*(*(SurfArea+indexP*6*3+3*4+comp));//Bottom direction, x-component
			source[5][comp]=s_PT*(*(SurfArea+indexP*6*3+3*5+comp));//Top direction, x-component

			sumSource=(-1.0)*(source[0][comp]+source[1][comp]+source[2][comp]+source[3][comp]+source[4][comp]+source[5][comp]);
		break;
		case 100:
			comp=0;
			source[0][comp]=s_PW*(*(SurfArea+indexP*6*3+3*0+comp));//West direction, x-component
			source[1][comp]=s_PE*(*(SurfArea+indexP*6*3+3*1+comp));//East direction, x-component
			source[2][comp]=s_PS*(*(SurfArea+indexP*6*3+3*2+comp));//South direction, x-component
			source[3][comp]=s_PN*(*(SurfArea+indexP*6*3+3*3+comp));//North direction, x-component
			source[4][comp]=s_PB*(*(SurfArea+indexP*6*3+3*4+comp));//Bottom direction, x-component
			source[5][comp]=s_PT*(*(SurfArea+indexP*6*3+3*5+comp));//Top direction, x-component

			sumSource=(-1.0)*(source[0][comp]+source[1][comp]+source[2][comp]+source[3][comp]+source[4][comp]+source[5][comp]);
		break;
		case 2:

			comp=1;
			source[0][comp]=s_PW*(*(SurfArea+indexP*6*3+3*0+comp));//West direction, y-component
			source[1][comp]=s_PE*(*(SurfArea+indexP*6*3+3*1+comp));//East direction, y-component
			source[2][comp]=s_PS*(*(SurfArea+indexP*6*3+3*2+comp));//South direction, y-component
			source[3][comp]=s_PN*(*(SurfArea+indexP*6*3+3*3+comp));//North direction, y-component
			source[4][comp]=s_PB*(*(SurfArea+indexP*6*3+3*4+comp));//Bottom direction, y-component
			source[5][comp]=s_PT*(*(SurfArea+indexP*6*3+3*5+comp));//Top direction, y-component

			sumSource=(-1.0)*(source[0][comp]+source[1][comp]+source[2][comp]+source[3][comp]+source[4][comp]+source[5][comp]);

		break;
		case 200:
			comp=1;
			source[0][comp]=s_PW*(*(SurfArea+indexP*6*3+3*0+comp));//West direction, y-component
			source[1][comp]=s_PE*(*(SurfArea+indexP*6*3+3*1+comp));//East direction, y-component
			source[2][comp]=s_PS*(*(SurfArea+indexP*6*3+3*2+comp));//South direction, y-component
			source[3][comp]=s_PN*(*(SurfArea+indexP*6*3+3*3+comp));//North direction, y-component
			source[4][comp]=s_PB*(*(SurfArea+indexP*6*3+3*4+comp));//Bottom direction, y-component
			source[5][comp]=s_PT*(*(SurfArea+indexP*6*3+3*5+comp));//Top direction, y-component

			sumSource=(-1.0)*(source[0][comp]+source[1][comp]+source[2][comp]+source[3][comp]+source[4][comp]+source[5][comp]);

		break;
		case 3:

			comp=2;
			source[0][comp]=s_PW*(*(SurfArea+indexP*6*3+3*0+comp));//West direction, z-component
			source[1][comp]=s_PE*(*(SurfArea+indexP*6*3+3*1+comp));//East direction, z-component
			source[2][comp]=s_PS*(*(SurfArea+indexP*6*3+3*2+comp));//South direction, z-component
			source[3][comp]=s_PN*(*(SurfArea+indexP*6*3+3*3+comp));//North direction, z-component
			source[4][comp]=s_PB*(*(SurfArea+indexP*6*3+3*4+comp));//Bottom direction, z-component
			source[5][comp]=s_PT*(*(SurfArea+indexP*6*3+3*5+comp));//Top direction, z-component

			sumSource=(-1.0)*(source[0][comp]+source[1][comp]+source[2][comp]+source[3][comp]+source[4][comp]+source[5][comp]);

		break;
		case 300:
			comp=2;
			source[0][comp]=s_PW*(*(SurfArea+indexP*6*3+3*0+comp));//West direction, z-component
			source[1][comp]=s_PE*(*(SurfArea+indexP*6*3+3*1+comp));//East direction, z-component
			source[2][comp]=s_PS*(*(SurfArea+indexP*6*3+3*2+comp));//South direction, z-component
			source[3][comp]=s_PN*(*(SurfArea+indexP*6*3+3*3+comp));//North direction, z-component
			source[4][comp]=s_PB*(*(SurfArea+indexP*6*3+3*4+comp));//Bottom direction, z-component
			source[5][comp]=s_PT*(*(SurfArea+indexP*6*3+3*5+comp));//Top direction, z-component

			sumSource=(-1.0)*(source[0][comp]+source[1][comp]+source[2][comp]+source[3][comp]+source[4][comp]+source[5][comp]);

		break;
	}
	return(sumSource);
}
