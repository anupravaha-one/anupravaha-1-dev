#include "APKTURBgrad3D.h"
/**
This function calculates the product of the velocity gradients which appear in the production function of the k-epsilon 
equations. (This function is used for 3D gemoetries).
*/
double grad3D(Base *base,Zone *zone,double *variable,double *slabVariable[6],double *variable1,double *slabVariable1[6],double *variable2,double *slabVariable2[6],int INDEX[3])
{
	
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize;
	double UE,UW,UN,US,UT,UB,UP,VE,VW,VN,VS,VT,VB,VP,WE,WW,WN,WS,WT,WB,WP;
	double DUx,DUy,DUz,DVx,DVy,DVz,DWx,DWy,DWz,ephsilon,pi[3];
	double e[3][3],*deltainv;
	double *rho,*slabrho[6];
	
	ptrsize = zone->Size;
	NI=*(ptrsize+0);
	NJ=*(ptrsize+1);
	NK=*(ptrsize+2);
	
	deltainv = zone->ZGeom->DeltaInverse; // Delta inverse matrix
	rho=MakeDomainArray(zone,710);
	MakeSlabPtrArray(zone,slabrho,710);
	
	k=INDEX[0];
	j=INDEX[1];
	i=INDEX[2];
	
	P=(NI-1)*(NJ-1)*k+(NI-1)*j+i; 
	E=(NI-1)*(NJ-1)*k+(NI-1)*j+(i+1);
	W=(NI-1)*(NJ-1)*k+(NI-1)*j+(i-1);
	N=(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i;
	S=(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i;
	T=(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i;
	B=(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i;
	
	if(base->flagMPFV==0)
	{
		UP=*(variable+P);
		VP=*(variable1+P);
		WP=*(variable2+P);
		
		if(i==(NI-2))
		{
			UE=*(slabVariable[0]+(NJ-1)*k+j);
			VE=*(slabVariable1[0]+(NJ-1)*k+j);
			WE=*(slabVariable2[0]+(NJ-1)*k+j);
		}
		else
		{
			UE=*(variable+E);
			VE=*(variable1+E);
			WE=*(variable2+E);
		}
		if(i==0)
		{
			UW=*(slabVariable[1]+(NJ-1)*k+j);
			VW=*(slabVariable1[1]+(NJ-1)*k+j);
			WW=*(slabVariable2[1]+(NJ-1)*k+j);
		}
		else
		{
			UW=*(variable+W);
			VW=*(variable1+W);
			WW=*(variable2+W);
		}
		if(j==(NJ-2))
		{
			UN=*(slabVariable[2]+(NI-1)*k+i);
			VN=*(slabVariable1[2]+(NI-1)*k+i);
			WN=*(slabVariable2[2]+(NI-1)*k+i);
		}
		else
		{
			UN=*(variable+N);
			VN=*(variable1+N);
			WN=*(variable2+N);
		}
		if(j==0)
		{
			US=*(slabVariable[3]+(NI-1)*k+i);
			VS=*(slabVariable1[3]+(NI-1)*k+i);
			WS=*(slabVariable2[3]+(NI-1)*k+i);
		}
		else
		{
			US=*(variable+S);
			VS=*(variable1+S);
			WS=*(variable2+S);
		}
		if(k==0)
		{
			UB=*(slabVariable[5]+(NI-1)*j+i);
			VB=*(slabVariable1[5]+(NI-1)*j+i);
			WB=*(slabVariable2[5]+(NI-1)*j+i);
		}
		else
		{
			UB=*(variable+B);
			VB=*(variable1+B);
			WB=*(variable2+B);
		}
		if(k==(NK-2))
		{
			UT=*(slabVariable[4]+(NI-1)*j+i);
			VT=*(slabVariable1[4]+(NI-1)*j+i);
			WT=*(slabVariable2[4]+(NI-1)*j+i);
		}
		else
		{
			UT=*(variable+T);
			VT=*(variable1+T);
			WT=*(variable2+T);
		}
	}
	else if(base->flagMPFV==1)
	{
		UP=*(variable+P)/(*(rho+P));
		VP=*(variable1+P)/(*(rho+P));
		WP=*(variable2+P)/(*(rho+P));
		if(i==(NI-2))
		{
			UE=*(slabVariable[0]+(NJ-1)*k+j)/(*(slabrho[0]+(NJ-1)*k+j));
			VE=*(slabVariable1[0]+(NJ-1)*k+j)/(*(slabrho[0]+(NJ-1)*k+j));
			WE=*(slabVariable2[0]+(NJ-1)*k+j)/(*(slabrho[0]+(NJ-1)*k+j));
		}
		else
		{
			UE=*(variable+E)/(*(rho+E));
			VE=*(variable1+E)/(*(rho+E));
			WE=*(variable2+E)/(*(rho+E));
		}
		if(i==0)
		{
			UW=*(slabVariable[1]+(NJ-1)*k+j)/(*(slabrho[1]+(NJ-1)*k+j));
			VW=*(slabVariable1[1]+(NJ-1)*k+j)/(*(slabrho[1]+(NJ-1)*k+j));
			WW=*(slabVariable2[1]+(NJ-1)*k+j)/(*(slabrho[1]+(NJ-1)*k+j));
		}
		else
		{
			UW=*(variable+W)/(*(rho+W));
			VW=*(variable1+W)/(*(rho+W));
			WW=*(variable2+W)/(*(rho+W));
		}
		if(j==(NJ-2))
		{
			UN=*(slabVariable[2]+(NI-1)*k+i)/(*(slabrho[2]+(NI-1)*k+i));
			VN=*(slabVariable1[2]+(NI-1)*k+i)/(*(slabrho[2]+(NI-1)*k+i));
			WN=*(slabVariable2[2]+(NI-1)*k+i)/(*(slabrho[2]+(NI-1)*k+i));
		}
		else
		{
			UN=*(variable+N)/(*(rho+N));
			VN=*(variable1+N)/(*(rho+N));
			WN=*(variable2+N)/(*(rho+N));
		}
		if(j==0)
		{
			US=*(slabVariable[3]+(NI-1)*k+i)/(*(slabrho[3]+(NI-1)*k+i));
			VS=*(slabVariable1[3]+(NI-1)*k+i)/(*(slabrho[3]+(NI-1)*k+i));
			WS=*(slabVariable2[3]+(NI-1)*k+i)/(*(slabrho[3]+(NI-1)*k+i));
		}
		else
		{
			US=*(variable+S)/(*(rho+S));
			VS=*(variable1+S)/(*(rho+S));
			WS=*(variable2+S)/(*(rho+S));
		}
		if(k==0)
		{
			UB=*(slabVariable[5]+(NI-1)*j+i)/(*(slabrho[5]+(NI-1)*j+i));
			VB=*(slabVariable1[5]+(NI-1)*j+i)/(*(slabrho[5]+(NI-1)*j+i));
			WB=*(slabVariable2[5]+(NI-1)*j+i)/(*(slabrho[5]+(NI-1)*j+i));
		}
		else
		{
			UB=*(variable+B)/(*(rho+B));
			VB=*(variable1+B)/(*(rho+B));
			WB=*(variable2+B)/(*(rho+B));
		}
		if(k==(NK-2))
		{
			UT=*(slabVariable[4]+(NI-1)*j+i)/(*(slabrho[4]+(NI-1)*j+i));
			VT=*(slabVariable1[4]+(NI-1)*j+i)/(*(slabrho[4]+(NI-1)*j+i));
			WT=*(slabVariable2[4]+(NI-1)*j+i)/(*(slabrho[4]+(NI-1)*j+i));
		}
		else
		{
			UT=*(variable+T)/(*(rho+T));
			VT=*(variable1+T)/(*(rho+T));
			WT=*(variable2+T)/(*(rho+T));
		}
	}
	//printf("\nUE=%lf$UW=%lf$UN=%lf$US=%lf$UT=%lf$UB=%lf",UE,UW,UN,US,UT,UB);
	//printf("\nVE=%lf$VW=%lf$VN=%lf$VS=%lf$VT=%lf$VB=%lf",VE,VW,VN,VS,VT,VB);
	//printf("\nWE=%lf$WW=%lf$WN=%lf$WS=%lf$WT=%lf$WB=%lf",VE,VW,VN,VS,VT,VB);
	pi[0]=UE-UW;
	pi[1]=UN-US;
	pi[2]=UT-UB;
	
	DUx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
	DUy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
	DUz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
	
	pi[0]=VE-VW;
	pi[1]=VN-VS;
	pi[2]=VT-VB;
	
	DVx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
	DVy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
	DVz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
	
	pi[0]=WE-WW;
	pi[1]=WN-WS;
	pi[2]=WT-WB;
	
	DWx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
	DWy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
	DWz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
	
	e[0][0]=DUx;e[0][1]=DUy;e[0][2]=DUz;
	e[1][0]=DVx;e[1][1]=DVy;e[1][2]=DVz;
	e[2][0]=DWx;e[2][1]=DWy;e[2][2]=DWz;
	
	ephsilon=0;
	for(j=0;j<3;j++)
	for(i=0;i<3;i++)
		ephsilon=ephsilon+(e[i][j]+e[j][i])*(e[i][j]+e[j][i]);
	
	ephsilon=0.5*ephsilon;
	
	return (ephsilon);
}

