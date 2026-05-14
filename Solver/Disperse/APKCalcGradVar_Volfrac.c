#include "APKCalcGradVar_Volfrac.h"
/**
The function CalcGradVar_Volfrac() calculates and returns the derivatives of a variable(Velocity multiplied by the Density of disperse phase and Diffusion coefficient of disperse phase) in x,y,z directions
using Finite Difference approximation. This function is used in calulating the source terms in momentum
equations that arises from turbulent dispersion term in particle momentum equations.
*/
double CalcGradVar_Volfrac(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3],int variableTag)
{
	
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize,flagRhod,flagRhoc;
	double UE,UW,UN,US,UT,UB,UP;
	double pi[2];
	
	//Turbulent Gas-Particle
	double *SlabRho[6],*SlabVartheta_Disperse[6],*SlabD_Disperse[6];
	double *rho,*Vartheta_Disperse,*D_Disperse,Cons;
	double *U_g,*SlabU_g[6];
	//Volume fraction of the disperse
	Vartheta_Disperse=MakeDomainArray(zone,650);
	MakeSlabPtrArray(zone,SlabVartheta_Disperse,650);
	//Density of the disperse phase
	rho=MakeDomainArray(zone,66);
	MakeSlabPtrArray(zone,SlabRho,66);
	flagRhod=base->BaseParameters->flagRhoDisperse;
	flagRhoc=zone->Flowsol->flagRho;
	//Turbulent Viscosity of the gas phase
	D_Disperse=MakeDomainArray(zone,727);
	MakeSlabPtrArray(zone,SlabD_Disperse,727);
	
	double EPSILONU;
	EPSILONU = base->BaseParameters->EPSILONDrop;
	
	if(variableTag==62)
	{
		U_g=zone->PrevFlowsol->u;
		MakeSlabPtrArray(zone,SlabU_g,10);
	}
	else if(variableTag==63)
	{
		U_g=zone->PrevFlowsol->v;
		MakeSlabPtrArray(zone,SlabU_g,20);
	
	}
	else if(variableTag==64)
	{
		U_g=zone->PrevFlowsol->w;
		MakeSlabPtrArray(zone,SlabU_g,30);
	
	}
	
	
	
	
	
	Cons=1.0;
	//printf("\n Cons=%e",Cons);
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
	//printf("\n *(D_Disperse+P)=%e",*(D_Disperse+P));
	if(*(Vartheta_Disperse+P)>EPSILONU)
		UP=*(variable+P)/(*(Vartheta_Disperse+P) * *(rho+P*flagRhod));
	else
		UP=*(U_g+P);
				
	UP=UP* *(rho+P*flagRhod)* *(D_Disperse+P) *Cons;
	
	if(i==(NI-2))
	{
		if(*(SlabVartheta_Disperse[0]+(NJ-1)*k+j)>EPSILONU)
			UE=*(slabVariable[0]+(NJ-1)*k+j)/(*(SlabVartheta_Disperse[0]+(NJ-1)*k+j) * *(SlabRho[0]+((NJ-1)*k+j)*flagRhod));
		else
			UE=*(SlabU_g[0]+(NJ-1)*k+j);
		
		UE=UE* *(SlabRho[0]+((NJ-1)*k+j)*flagRhod)* *(SlabD_Disperse[0]+(NJ-1)*k+j) *Cons;
	}
	else
	{
		if(*(Vartheta_Disperse+E)>EPSILONU)
			UE=*(variable+E)/(*(Vartheta_Disperse+E) * *(rho+E*flagRhod));
		else
			UE=*(U_g+E);
		
		UE=UE* *(rho+E*flagRhod)* *(D_Disperse+E) *Cons;
	}
	if(i==0)
	{
		if(*(SlabVartheta_Disperse[1]+(NJ-1)*k+j)>EPSILONU)
			UW=*(slabVariable[1]+(NJ-1)*k+j)/(*(SlabVartheta_Disperse[1]+(NJ-1)*k+j) * *(SlabRho[1]+((NJ-1)*k+j)*flagRhod));
		else
			UW=*(SlabU_g[1]+(NJ-1)*k+j);
		
		UW=UW* *(SlabRho[1]+((NJ-1)*k+j)*flagRhod)* *(SlabD_Disperse[1]+(NJ-1)*k+j) *Cons;
	}
	else
	{
		if(*(Vartheta_Disperse+W)>EPSILONU)
			UW=*(variable+W)/(*(Vartheta_Disperse+W) * *(rho+W*flagRhod));
		else
			UW=*(U_g+W);
			
		UW=UW* *(rho+W*flagRhod) * *(D_Disperse+W)*Cons;
	}
	if(j==(NJ-2))
	{
		if(*(SlabVartheta_Disperse[2]+(NI-1)*k+i)>EPSILONU)
			UN=*(slabVariable[2]+(NI-1)*k+i)/(*(SlabVartheta_Disperse[2]+(NI-1)*k+i) * *(SlabRho[2]+((NI-1)*k+i)*flagRhod));
		else
			UN=*(SlabU_g[2]+(NI-1)*k+i);
			
		UN=UN* *(SlabRho[2]+((NI-1)*k+i)*flagRhod)* *(SlabD_Disperse[2]+(NI-1)*k+i) *Cons;
	}
	else
	{
		if(*(Vartheta_Disperse+N)>EPSILONU)
			UN=*(variable+N)/(*(Vartheta_Disperse+N) * *(rho+N*flagRhod));
		else
			UN=*(U_g+N);
		
		UN=UN* *(rho+N*flagRhod)* *(D_Disperse+N)*Cons;
	}
	if(j==0)
	{
		if(*(SlabVartheta_Disperse[3]+(NI-1)*k+i)>EPSILONU)
			US=*(slabVariable[3]+(NI-1)*k+i)/(*(SlabVartheta_Disperse[3]+(NI-1)*k+i) * *(SlabRho[3]+((NI-1)*k+i)*flagRhod));
		else
			US=*(SlabU_g[3]+(NI-1)*k+i);
					
		US=US* *(SlabRho[3]+((NI-1)*k+i)*flagRhod)* *(SlabD_Disperse[3]+(NI-1)*k+i) *Cons;
	}
	else
	{
		if(*(Vartheta_Disperse+S)>EPSILONU)
			US=*(variable+S)/(*(Vartheta_Disperse+S) * *(rho+S*flagRhod));
		else
			US=*(U_g+S);
		
		US=US* *(rho+S*flagRhod)* *(D_Disperse+S) *Cons;
	}
	if(k==0)
	{
		if(*(SlabVartheta_Disperse[5]+(NI-1)*j+i)>EPSILONU)
			UB=*(slabVariable[5]+(NI-1)*j+i)/(*(SlabVartheta_Disperse[5]+(NI-1)*j+i) * *(SlabRho[5]+((NI-1)*j+i)*flagRhod));
		else
			UB=*(SlabU_g[5]+(NI-1)*j+i);
		
		UB=UB* *(SlabRho[5]+((NI-1)*j+i)*flagRhod)* *(SlabD_Disperse[5]+(NI-1)*j+i) *Cons;
	}
	else
	{	
		if(*(Vartheta_Disperse+B)>EPSILONU)
			UB=*(variable+B)/(*(Vartheta_Disperse+B) * *(rho+B*flagRhod));
		else
			UB=*(U_g+B);
		
		UB=UB* *(rho+B*flagRhod)* *(D_Disperse+B) *Cons;
	}
	if(k==(NK-2))
	{
		if(*(SlabVartheta_Disperse[4]+(NI-1)*j+i)>EPSILONU)
			UT=*(slabVariable[4]+(NI-1)*j+i)/(*(SlabVartheta_Disperse[4]+(NI-1)*j+i) * *(SlabRho[4]+((NI-1)*j+i)*flagRhod));
		else
			UT=*(SlabU_g[4]+(NI-1)*j+i);
		
		UT=UT* *(SlabRho[4]+((NI-1)*j+i)*flagRhod)* *(SlabD_Disperse[4]+(NI-1)*j+i) *Cons;
	}
	else
	{
		if(*(Vartheta_Disperse+T)>EPSILONU)
			UT=*(variable+T)/(*(Vartheta_Disperse+T) * *(rho+T*flagRhod));
		else
			UT=*(U_g+T);
		
		UT=UT* *(rho+T*flagRhod)* *(D_Disperse+T) *Cons;
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