#include "APKCalcPressureGrad_Common.h"
/**
The function CalcPressureGrad_Common() is used to compute the pressure gradient term in the Navier Stokes equation. This function takes the individual velocity component's variable tag as input and calculates the coresponding
component of Pressure gradient. Note: This function ALSO includes the "negative sign" involved in the pressure gradient term occuring in the Navier Stokes equations.
*/
int CalcPressureGrad_Common(Base * base,double *variable,double *slabVariable[6],double *volume,double *slabVolume[6],double *SurfArea,int index[3],int *ptrSize,double grad[3])
{
	
	int k,j,i,NK,NJ,NI,indexP,indexE,indexW,indexN,indexS,indexB,indexT;
	double sP,sW,sE,sS,sN,sB,sT;
	double volP,volW,volE,volS,volN,volB,volT;
	double s_PE,s_PW,s_PS,s_PN,s_PB,s_PT;
	
	double *Surf_indexP;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	k=index[0],j=index[1],i=index[2];
	indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
	indexE=indexP+1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
	indexW=indexP-1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
	indexN=indexP+(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
	indexS=indexP-(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
	indexT=indexP+(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
	indexB=indexP-(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
	
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

	Surf_indexP = (SurfArea+indexP*6*3);
	//West direction, x-component+East direction, x-component+South direction, x-component+North direction,x-component+Bottom direction, x-component+Top direction, x-component
	grad[0]=s_PW*(*(Surf_indexP))+s_PE*(*(Surf_indexP+3))+s_PS*(*(Surf_indexP+6))+s_PN*(*(Surf_indexP+9))+
			s_PB*(*(Surf_indexP+12))+s_PT*(*(Surf_indexP+15));

	grad[0]=(-1.0)*grad[0];
	
	//West direction, y-component//East direction, y-component//South direction, y-component//North direction, y-component//Bottom direction, y-component//Top direction, y-component
	grad[1]=s_PW*(*(Surf_indexP+1))+s_PE*(*(Surf_indexP+4))+s_PS*(*(Surf_indexP+7))+s_PN*(*(Surf_indexP+10))+s_PB*(*(Surf_indexP+13))+s_PT*(*(Surf_indexP+16));
	
	grad[1]=(-1.0)*grad[1];
	
	if(base->BaseParameters->flag2d==0)
	{
		//West direction, z-component+East direction, z-component+South direction, z-component+North direction, z-component+Bottom direction, z-component+Top direction, z-component
		grad[2]=s_PW*(*(Surf_indexP+2))+s_PE*(*(Surf_indexP+5))+s_PS*(*(Surf_indexP+8))+s_PN*(*(Surf_indexP+11))+s_PB*(*(Surf_indexP+14))+s_PT*(*(Surf_indexP+17));
		grad[2] = -1.0 *grad[2];
	}
	return(0);
}
