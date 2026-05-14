#include "APKCalcGradVar_Disperse.h"
/**
The function CalcGradVar_Disperse() calculates and returns the derivatives of a particle velocities in x,y,z directions
using Finite Difference approximation. This function is used in calulating the source terms in momentum
equations that arises from varaible viscosity of particle.
*/
double CalcGradVar_Disperse(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3],double *U_g,double *slabU_g[6])
{
	
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize,flagRhod;
	double UE,UW,UN,US,UT,UB,UP;
	double pi[2];
	
	//Turbulent Gas-Particle
	double *SlabRho[6],*SlabVartheta_Disperse[6];
	double *rho,*Vartheta_Disperse;
	
	//Volume fraction of the disperse
	Vartheta_Disperse=MakeDomainArray(zone,650);
	MakeSlabPtrArray(zone,SlabVartheta_Disperse,650);
	//Density of the disperse phase
	rho=MakeDomainArray(zone,66);
	MakeSlabPtrArray(zone,SlabRho,66);
	flagRhod=base->BaseParameters->flagRhoDisperse;
	//end
	double EPSILONU;
	EPSILONU = base->BaseParameters->EPSILONDrop;
	ptrsize = zone->Size;
	
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
	
	if(*(Vartheta_Disperse+P)>EPSILONU)
		UP=*(variable+P)/(*(Vartheta_Disperse+P) * *(rho+P*flagRhod));
	else
		UP=*(U_g+P);
	
	if(i==(NI-2))
	{
		if(*(SlabVartheta_Disperse[0]+(NJ-1)*k+j)>EPSILONU)
			UE=*(slabVariable[0]+(NJ-1)*k+j)/(*(SlabVartheta_Disperse[0]+(NJ-1)*k+j) * *(SlabRho[0]+((NJ-1)*k+j)*flagRhod));
		else
			UE=*(slabU_g[0]+(NJ-1)*k+j);
	}
	else
	{
		if(*(Vartheta_Disperse+E)>EPSILONU)
			UE=*(variable+E)/(*(Vartheta_Disperse+E) * *(rho+E*flagRhod));
		else
			UE=*(U_g+E);
	}
	if(i==0)
	{
		if(*(SlabVartheta_Disperse[1]+(NJ-1)*k+j)>EPSILONU)
			UW=*(slabVariable[1]+(NJ-1)*k+j)/(*(SlabVartheta_Disperse[1]+(NJ-1)*k+j) * *(SlabRho[1]+((NJ-1)*k+j)*flagRhod));
		else
			UW=*(slabU_g[1]+(NJ-1)*k+j);
	}
	else
	{	
		if(*(Vartheta_Disperse+W)>EPSILONU)
			UW=*(variable+W)/(*(Vartheta_Disperse+W) * *(rho+W*flagRhod));
		else
			UW=*(U_g+W);
					
	}
	if(j==(NJ-2))
	{
		if(*(SlabVartheta_Disperse[2]+(NI-1)*k+i)>EPSILONU)
			UN=*(slabVariable[2]+(NI-1)*k+i)/(*(SlabVartheta_Disperse[2]+(NI-1)*k+i) * *(SlabRho[2]+((NI-1)*k+i)*flagRhod));
		else
			UN=*(slabU_g[2]+(NI-1)*k+i);
	}
	else
	{
		if(*(Vartheta_Disperse+N)>EPSILONU)
			UN=*(variable+N)/(*(Vartheta_Disperse+N) * *(rho+N*flagRhod));
		else
			UN=*(U_g+N);
	}
	if(j==0)
	{
		if(*(SlabVartheta_Disperse[3]+(NI-1)*k+i)>EPSILONU)
			US=*(slabVariable[3]+(NI-1)*k+i)/(*(SlabVartheta_Disperse[3]+(NI-1)*k+i) * *(SlabRho[3]+((NI-1)*k+i)*flagRhod));
		else
			US=*(slabU_g[3]+(NI-1)*k+i);
	}
	else
	{
		if(*(Vartheta_Disperse+S)>EPSILONU)
			US=*(variable+S)/(*(Vartheta_Disperse+S) * *(rho+S*flagRhod));
		else
			US=*(U_g+S);
			
	}
	if(k==0)
	{
		if(*(SlabVartheta_Disperse[5]+(NI-1)*j+i)>EPSILONU)
			UB=*(slabVariable[5]+(NI-1)*j+i)/(*(SlabVartheta_Disperse[5]+(NI-1)*j+i) * *(SlabRho[5]+((NI-1)*j+i)*flagRhod));
		else
			UB=*(slabU_g[5]+(NI-1)*j+i);
	}
	else
	{
		if(*(Vartheta_Disperse+B)>EPSILONU)
			UB=*(variable+B)/(*(Vartheta_Disperse+B) * *(rho+B*flagRhod));
		else
			UB=*(U_g+B);
	}
	if(k==(NK-2))
	{
		if(*(SlabVartheta_Disperse[4]+(NI-1)*j+i)>EPSILONU)
			UT=*(slabVariable[4]+(NI-1)*j+i)/(*(SlabVartheta_Disperse[4]+(NI-1)*j+i) * *(SlabRho[4]+((NI-1)*j+i)*flagRhod));
		else
			UT=*(slabU_g[4]+(NI-1)*j+i);
	}
	else
	{
		if(*(Vartheta_Disperse+T)>EPSILONU)
			UT=*(variable+T)/(*(Vartheta_Disperse+T) * *(rho+T*flagRhod));
		else
			UT=*(U_g+T);
	}
	
	
	//printf("\nUP=%e,UE=%e,UW=%e,UN=%e,US=%e,UT=%e,UB=%e--[Zone %s]",UP,UE,UW,UN,US,UT,UB,zone->zonename);
	if(base->BaseParameters->flag2d == 1)
	{
		pi[0]=UE-UW;
		pi[1]=UN-US;
		D[0]=*((zone->ZGeom->DeltaInverse+P*4)+0) * pi[0] + *((zone->ZGeom->DeltaInverse+P*4)+1) * pi[1];
		D[1]=*((zone->ZGeom->DeltaInverse+P*4)+2) * pi[0] + *((zone->ZGeom->DeltaInverse+P*4)+3) * pi[1];
		D[2]=0.0;
		//printf("\n[Zone %s]$$[k=%d][j=%d][i=%d]$D[1]=%lf$D[2]=%lf$D[3]=%lf",zone->zonename,k,j,i,D[0],D[1],D[2]);	
	}
	else
	{
		pi[0]=UE-UW;
		pi[1]=UN-US;
		pi[2]=UT-UB;
		D[0]=*((zone->ZGeom->DeltaInverse+P*9)+0) * pi[0] + *((zone->ZGeom->DeltaInverse+P*9)+1) * pi[1] + *((zone->ZGeom->DeltaInverse+P*9)+2) * pi[2];
		D[1]=*((zone->ZGeom->DeltaInverse+P*9)+3) * pi[0] + *((zone->ZGeom->DeltaInverse+P*9)+4) * pi[1] + *((zone->ZGeom->DeltaInverse+P*9)+5) * pi[2];
		D[2]=*((zone->ZGeom->DeltaInverse+P*9)+6) * pi[0] + *((zone->ZGeom->DeltaInverse+P*9)+7) * pi[1] + *((zone->ZGeom->DeltaInverse+P*9)+8) * pi[2];
		//printf("\n[Zone %s]$$[k=%d][j=%d][i=%d]$D[1]=%lf$D[2]=%lf$D[3]=%lf",zone->zonename,k,j,i,D[0],D[1],D[2]);
	}
	return 0.0;
}