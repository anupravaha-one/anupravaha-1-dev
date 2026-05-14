#include "APKTURBgrad2D.h"
/**
This function calculates the product of the velocity gradients which appear in the production function of the k-epsilon 
equations. (This function is used for 2D gemoetries).
 */
double grad2D(Base *base,Zone *zone,double *variable,double *slabVariable[6],double *variable1,double *slabVariable1[6],int INDEX[3])
{              // variable is notation for U and vraiable1 is notation for V
	
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize;
	double UE,UW,UN,US,UT,UB,UP,VE,VW,VN,VS,VT,VB,VP;
	double DUx,DUy,DVx,DVy,ephsilon,pi[2],*deltainv;
	double e[2][2];
	double *rho,*slabrho[6];
	
	ptrsize = zone->Size;
	deltainv = zone->ZGeom->DeltaInverse; // Delta inverse matrix
	rho=MakeDomainArray(zone,710);
	MakeSlabPtrArray(zone,slabrho,710);
	
	NI=*(ptrsize+0);
	NJ=*(ptrsize+1);
	NK=*(ptrsize+2);
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
		
		if(i==(NI-2))
		{
			UE=*(slabVariable[0]+(NJ-1)*k+j);
			VE=*(slabVariable1[0]+(NJ-1)*k+j);
		}
		else
		{
			UE=*(variable+E);
			VE=*(variable1+E);
		}
		if(i==0)
		{
			UW=*(slabVariable[1]+(NJ-1)*k+j);
			VW=*(slabVariable1[1]+(NJ-1)*k+j);
		}
		else
		{
			UW=*(variable+W);
			VW=*(variable1+W);
		}
		if(j==(NJ-2))
		{
			UN=*(slabVariable[2]+(NI-1)*k+i);
			VN=*(slabVariable1[2]+(NI-1)*k+i);
		}
		else
		{
			UN=*(variable+N);
			VN=*(variable1+N);
		}
		if(j==0)
		{
			US=*(slabVariable[3]+(NI-1)*k+i);
			VS=*(slabVariable1[3]+(NI-1)*k+i);
		}
		else
		{
			US=*(variable+S);
			VS=*(variable1+S);
		}
		if(k==0)
		{
			UB=*(slabVariable[5]+(NI-1)*j+i);
			VB=*(slabVariable1[5]+(NI-1)*j+i);
		}
		else
		{
			UB=*(variable+B);
			VB=*(variable1+B);
		}
		if(k==(NK-2))
		{
			UT=*(slabVariable[4]+(NI-1)*j+i);
			VT=*(slabVariable1[4]+(NI-1)*j+i);
		}
		else
		{
			UT=*(variable+T);
			VT=*(variable1+T);
		}
	}
	else if(base->flagMPFV==1)
	{
		UP=*(variable+P)/(*(rho+P));
		VP=*(variable1+P)/(*(rho+P));
		
		if(i==(NI-2))
		{
			UE=*(slabVariable[0]+(NJ-1)*k+j)/(*(slabrho[0]+(NJ-1)*k+j));
			VE=*(slabVariable1[0]+(NJ-1)*k+j)/(*(slabrho[0]+(NJ-1)*k+j));
		}
		else
		{
			UE=*(variable+E)/(*(rho+E));
			VE=*(variable1+E)/(*(rho+E));
		}
		if(i==0)
		{
			UW=*(slabVariable[1]+(NJ-1)*k+j)/(*(slabrho[1]+(NJ-1)*k+j));
			VW=*(slabVariable1[1]+(NJ-1)*k+j)/(*(slabrho[1]+(NJ-1)*k+j));
		}
		else
		{
			UW=*(variable+W)/(*(rho+W));
			VW=*(variable1+W)/(*(rho+W));
		}
		if(j==(NJ-2))
		{
			UN=*(slabVariable[2]+(NI-1)*k+i)/(*(slabrho[2]+(NI-1)*k+i));
			VN=*(slabVariable1[2]+(NI-1)*k+i)/(*(slabrho[2]+(NI-1)*k+i));
		}
		else
		{
			UN=*(variable+N)/(*(rho+N));
			VN=*(variable1+N)/(*(rho+N));
		}
		if(j==0)
		{
			US=*(slabVariable[3]+(NI-1)*k+i)/(*(slabrho[3]+(NI-1)*k+i));
			VS=*(slabVariable1[3]+(NI-1)*k+i)/(*(slabrho[3]+(NI-1)*k+i));
		}
		else
		{
			US=*(variable+S)/(*(rho+S));
			VS=*(variable1+S)/(*(rho+S));
		}
		if(k==0)
		{
			UB=*(slabVariable[5]+(NI-1)*j+i)/(*(slabrho[5]+(NI-1)*j+i));
			VB=*(slabVariable1[5]+(NI-1)*j+i)/(*(slabrho[5]+(NI-1)*j+i));
		}
		else
		{
			UB=*(variable+B)/(*(rho+B));
			VB=*(variable1+B)/(*(rho+B));
		}
		if(k==(NK-2))
		{
			UT=*(slabVariable[4]+(NI-1)*j+i)/(*(slabrho[4]+(NI-1)*j+i));
			VT=*(slabVariable1[4]+(NI-1)*j+i)/(*(slabrho[4]+(NI-1)*j+i));
		}
		else
		{
			UT=*(variable+T)/(*(rho+T));
			VT=*(variable1+T)/(*(rho+T));
		}
	}
	//printf("\nUE=%lf$UW=%lf$UN=%lfUS=%lfUT=%lfUB=%lf",UE,UW,UN,US,UT,UB);
	//printf("\nVE=%lf$VW=%lf$VN=%lfVS=%lfVT=%lfVB=%lf",VE,VW,VN,VS,VT,VB);
	pi[0]=UE-UW;
	pi[1]=UN-US;
	DUx=*(deltainv+P*4) * pi[0] + *((deltainv+P*4)+1) * pi[1];
	DUy=*((deltainv+P*4)+2) * pi[0] + *((deltainv+P*4)+3) * pi[1];
	
	pi[0]=VE-VW;
	pi[1]=VN-VS;
	DVx=*(deltainv+P*4) * pi[0] + *((deltainv+P*4)+1) * pi[1];
	DVy=*((deltainv+P*4)+2) * pi[0] + *((deltainv+P*4)+3) * pi[1];
	
	e[0][0]=DUx;e[0][1]=DUy;
	e[1][0]=DVx;e[1][1]=DVy;
	
	ephsilon=0.0;
	
	for(j=0;j<2;j++)
		for(i=0;i<2;i++)
			ephsilon=ephsilon+(e[i][j]+e[j][i])*(e[i][j]+e[j][i]);
		
	ephsilon=0.5*ephsilon;
	return (ephsilon);
}

