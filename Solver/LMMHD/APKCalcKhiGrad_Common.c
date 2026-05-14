#include "APKCalcKhiGrad_Common.h"

int CalcKhiGrad_Common(Base * base,double *variable1,double *slabVariable1[6],double *variable2,double *slabVariable2[6],double *variable3,double *slabVariable3[6],double *volume,double *slabVolume[6],double *SurfArea,int index[3],int *ptrSize,double grad[3])
{
	
	int k,j,i,NK,NJ,NI,indexP,indexE,indexW,indexN,indexS,indexB,indexT;
	double sP1,sW1,sE1,sS1,sN1,sB1,sT1;
	double sP2,sW2,sE2,sS2,sN2,sB2,sT2;
	double sP3,sW3,sE3,sS3,sN3,sB3,sT3;

	double KhiP,KhiW,KhiE,KhiS,KhiN,KhiB,KhiT;

	double volP,volW,volE,volS,volN,volB,volT;
// 	double s_PE1,s_PW1,s_PS1,s_PN1,s_PB1,s_PT1;
// 	double s_PE2,s_PW2,s_PS2,s_PN2,s_PB2,s_PT2;
// 	double s_PE3,s_PW3,s_PS3,s_PN3,s_PB3,s_PT3;

	double Khi_PE,Khi_PW,Khi_PS,Khi_PN,Khi_PB,Khi_PT;
	
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
	
	sP1=*(variable1+indexP);
	sP2=*(variable2+indexP);
	sP3=*(variable3+indexP);
	volP=*(volume+indexP);
	
	if(i==(NI-2))
	{
		sE1=*(slabVariable1[0]+(NJ-1)*k+j);
		sE2=*(slabVariable2[0]+(NJ-1)*k+j);
		sE3=*(slabVariable3[0]+(NJ-1)*k+j);
		volE=*(slabVolume[0]+(NJ-1)*k+j);
	}
	else
	{
		sE1=*(variable1+indexE);
		sE2=*(variable2+indexE);
		sE3=*(variable3+indexE);
		volE=*(volume+indexE);
	}
	
	if(i==0)
	{
		sW1 =*(slabVariable1[1]+(NJ-1)*k+j);
		sW2 =*(slabVariable2[1]+(NJ-1)*k+j);
		sW3 =*(slabVariable3[1]+(NJ-1)*k+j);
		volW =*(slabVolume[1]+(NJ-1)*k+j);
	
	}
	else
	{
		sW1=*(variable1+indexW);
		sW2=*(variable2+indexW);
		sW3=*(variable3+indexW);
		volW=*(volume+indexW);
	}
	
	if(j==(NJ-2))
	{
		sN1=*(slabVariable1[2]+(NI-1)*k+i);
		sN2=*(slabVariable2[2]+(NI-1)*k+i);
		sN3=*(slabVariable3[2]+(NI-1)*k+i);
		volN=*(slabVolume[2]+(NI-1)*k+i);
	}
	else
	{
		sN1=*(variable1+indexN);
		sN2=*(variable2+indexN);
		sN3=*(variable3+indexN);
		volN=*(volume+indexN);
	}
	
	if(j==0)
	{
		sS1=*(slabVariable1[3]+(NI-1)*k+i);
		sS2=*(slabVariable2[3]+(NI-1)*k+i);
		sS3=*(slabVariable3[3]+(NI-1)*k+i);
		volS=*(slabVolume[3]+(NI-1)*k+i);
	}
	else
	{
		sS1=*(variable1+indexS);
		sS2=*(variable2+indexS);
		sS3=*(variable3+indexS);
		volS=*(volume+indexS);
	}
	
	if(k==(NK-2))
	{
		sT1=*(slabVariable1[4]+(NI-1)*j+i);
		sT2=*(slabVariable2[4]+(NI-1)*j+i);
		sT3=*(slabVariable3[4]+(NI-1)*j+i);
		volT=*(slabVolume[4]+(NI-1)*j+i);
	}
	else
	{
		sT1=*(variable1+indexT);
		sT2=*(variable2+indexT);
		sT3=*(variable3+indexT);
		volT=*(volume+indexT);
	}
	
	if(k==0)
	{
		sB1=*(slabVariable1[5]+(NI-1)*j+i);
		sB2=*(slabVariable2[5]+(NI-1)*j+i);
		sB3=*(slabVariable3[5]+(NI-1)*j+i);
		volB=*(slabVolume[5]+(NI-1)*j+i);
	}
	else
	{
		sB1=*(variable1+indexB);
		sB2=*(variable2+indexB);
		sB3=*(variable3+indexB);
		volB=*(volume+indexB);
	}

	KhiP=sP1*sP1+sP2*sP2+sP3*sP3;
	KhiE=sE1*sE1+sE2*sE2+sE3*sE3;
	KhiW=sW1*sW1+sW2*sW2+sW3*sW3;
	KhiN=sN1*sN1+sN2*sN2+sN3*sN3;
	KhiS=sS1*sS1+sS2*sS2+sS3*sS3;
	KhiT=sT1*sT1+sT2*sT2+sT3*sT3;
	KhiB=sB1*sB1+sB2*sB2+sB3*sB3;

	
	
	Khi_PE = (KhiP*volE + KhiE*volP )/( volP + volE);
	Khi_PW = (KhiP*volW + KhiW*volP )/( volP + volW);
	Khi_PN = (KhiP*volN + KhiN*volP )/( volP + volN);
	Khi_PS = (KhiP*volS + KhiS*volP )/( volP + volS);
	Khi_PT = (KhiP*volT + KhiT*volP )/( volP + volT);
	Khi_PB = (KhiP*volB + KhiB*volP )/( volP + volB);

	Surf_indexP = (SurfArea+indexP*6*3);
	//West direction, x-component+East direction, x-component+South direction, x-component+North direction,x-component+Bottom direction, x-component+Top direction, x-component
	grad[0]=Khi_PW*(*(Surf_indexP))+Khi_PE*(*(Surf_indexP+3))+Khi_PS*(*(Surf_indexP+6))+Khi_PN*(*(Surf_indexP+9))+
			Khi_PB*(*(Surf_indexP+12))+Khi_PT*(*(Surf_indexP+15));

	//grad[0]=(-1.0)*grad[0];
	
	//West direction, y-component//East direction, y-component//South direction, y-component//North direction, y-component//Bottom direction, y-component//Top direction, y-component
	grad[1]=Khi_PW*(*(Surf_indexP+1))+Khi_PE*(*(Surf_indexP+4))+Khi_PS*(*(Surf_indexP+7))+Khi_PN*(*(Surf_indexP+10))+Khi_PB*(*(Surf_indexP+13))+Khi_PT*(*(Surf_indexP+16));
	
	//grad[1]=(-1.0)*grad[1];
	
	if(base->BaseParameters->flag2d==0)
	{
		//West direction, z-component+East direction, z-component+South direction, z-component+North direction, z-component+Bottom direction, z-component+Top direction, z-component
		grad[2]=Khi_PW*(*(Surf_indexP+2))+Khi_PE*(*(Surf_indexP+5))+Khi_PS*(*(Surf_indexP+8))+Khi_PN*(*(Surf_indexP+11))+Khi_PB*(*(Surf_indexP+14))+Khi_PT*(*(Surf_indexP+17));
		//grad[2] = -1.0 *grad[2];
	}
	return(0);
}
