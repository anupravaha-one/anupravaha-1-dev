#include "APKCalcGradVar.h"
/**
The function CalcGradVar() calculates and returns the derivatives of a variable in x,y,z directions
using Finite Difference approximation. This function is used in calulating the source terms in momentum
equations that arises from varaible viscosity.
 */
double CalcGradVar(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3])
{
	
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize;
	double UE,UW,UN,US,UT,UB,UP;
	double pi[2];
	
	//Turbulent MPFV//chandrasekhar
	double *SlabRho[6];
	double *rho;
	rho=MakeDomainArray(zone,710);
	MakeSlabPtrArray(zone,SlabRho,710);
	//end
	
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
	
	if(base->flagMPFV==0)
	{
		UP=*(variable+P);
		
		if(i==(NI-2))
			UE=*(slabVariable[0]+(NJ-1)*k+j);
		
		else
			UE=*(variable+E);
		
		if(i==0)
			UW=*(slabVariable[1]+(NJ-1)*k+j);
		
		else
			UW=*(variable+W);
		
		if(j==(NJ-2))
			UN=*(slabVariable[2]+(NI-1)*k+i);
		
		else
			UN=*(variable+N);
		
		if(j==0)
			US=*(slabVariable[3]+(NI-1)*k+i);
		else
			US=*(variable+S);
		
		if(k==0)
			UB=*(slabVariable[5]+(NI-1)*j+i);
		
		else
			UB=*(variable+B);
		
		if(k==(NK-2))
			UT=*(slabVariable[4]+(NI-1)*j+i);
		
		else
			UT=*(variable+T);
	
	}
	//As momentum stored in the variable and the slabVariable , we need to divide all these variables with rho at that 
	//particlular location (cell or slab)-->chandrasekhar
	else if(base->flagMPFV==1)
	{
		UP=*(variable+P)/(*(rho+P));
		
		if(i==(NI-2))
			UE=*(slabVariable[0]+(NJ-1)*k+j)/(*(SlabRho[0]+(NJ-1)*k+j));
		
		else
			UE=*(variable+E)/(*(rho+E));
		
		if(i==0)
			UW=*(slabVariable[1]+(NJ-1)*k+j)/(*(SlabRho[1]+(NJ-1)*k+j));
		
		else
			UW=*(variable+W)/(*(rho+W));
		
		if(j==(NJ-2))
			UN=*(slabVariable[2]+(NI-1)*k+i)/(*(SlabRho[2]+(NI-1)*k+i));
		
		else
			UN=*(variable+N)/(*(rho+N));
		
		if(j==0)
			US=*(slabVariable[3]+(NI-1)*k+i)/(*(SlabRho[3]+(NI-1)*k+i));
		else
			US=*(variable+S)/(*(rho+S));
		
		if(k==0)
			UB=*(slabVariable[5]+(NI-1)*j+i)/(*(SlabRho[5]+(NI-1)*j+i));
		
		else
			UB=*(variable+B)/(*(rho+B));
		
		if(k==(NK-2))
			UT=*(slabVariable[4]+(NI-1)*j+i)/(*(SlabRho[4]+(NI-1)*j+i));
		
		else
			UT=*(variable+T)/(*(rho+T));
	}//end--chandrasekhar
	//printf("\nUP=%lf,UE=%lf,UW=%lf,UN=%lf,US=%lf,UT=%lf,UB=%lf--[Zone %s]",UP,UE,UW,UN,US,UT,UB,zone->zonename);
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
/**
The function CalcGradViscoeff() calculates and returns the derivatives of a viscosity in x,y,z directions
using Finite Difference approximation. This function is used in calulating the source terms in momentum
equations that arises from varaible viscosity.
 */
double CalcGradViscoeff(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3],double D[3])
{
	
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize,flagViscoeff;
	double UE,UW,UN,US,UT,UB,UP;
	double pi[2];
	
	flagViscoeff=zone->Flowsol->flagViscoeff;
	
	if(base->flagTurb==1&&base->flagParticle==1)
	{
	flagViscoeff=1;
	}
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
	
	UP=*(variable+P*flagViscoeff);
	
	if(i==(NI-2))
		UE=*(slabVariable[0]+((NJ-1)*k+j)*flagViscoeff);
	
	else
		UE=*(variable+E*flagViscoeff);
	
	if(i==0)
		UW=*(slabVariable[1]+((NJ-1)*k+j)*flagViscoeff);
	
	else
		UW=*(variable+W*flagViscoeff);
	
	if(j==(NJ-2))
		UN=*(slabVariable[2]+((NI-1)*k+i)*flagViscoeff);
	
	else
		UN=*(variable+N*flagViscoeff);
	
	if(j==0)
		US=*(slabVariable[3]+((NI-1)*k+i)*flagViscoeff);
	else
		US=*(variable+S*flagViscoeff);
	
	if(k==0)
		UB=*(slabVariable[5]+((NI-1)*j+i)*flagViscoeff);
	
	else
		UB=*(variable+B*flagViscoeff);
	
	if(k==(NK-2))
		UT=*(slabVariable[4]+((NI-1)*j+i)*flagViscoeff);
	
	else
		UT=*(variable+T*flagViscoeff);
	
	//printf("\nUP=%e,UE=%e,UW=%e,UN=%e,US=%e,UT=%e,UB=%e--[Zone %s]",UP,UE,UW,UN,US,UT,UB,zone->zonename);
	if(base->BaseParameters->flag2d == 1)
	{
		pi[0]=UE-UW;
		pi[1]=UN-US;
		D[0]=*((zone->ZGeom->DeltaInverse+P*4)+0) * pi[0] + *((zone->ZGeom->DeltaInverse+P*4)+1) * pi[1];
		D[1]=*((zone->ZGeom->DeltaInverse+P*4)+2) * pi[0] + *((zone->ZGeom->DeltaInverse+P*4)+3) * pi[1];
		D[2]=0.0;
		//printf("\n[Zone %s]$$[k=%d][j=%d][i=%d]$D[1]=%lf$D[2]=%lf$D[3]=%lf\n",zone->zonename,k,j,i,D[0],D[1],D[2]);
	}
	else
	{
		pi[0]=UE-UW;
		pi[1]=UN-US;
		pi[2]=UT-UB;
		D[0]=*((zone->ZGeom->DeltaInverse+P*9)+0) * pi[0] + *((zone->ZGeom->DeltaInverse+P*9)+1) * pi[1] + *((zone->ZGeom->DeltaInverse+P*9)+2) * pi[2];
		D[1]=*((zone->ZGeom->DeltaInverse+P*9)+3) * pi[0] + *((zone->ZGeom->DeltaInverse+P*9)+4) * pi[1] + *((zone->ZGeom->DeltaInverse+P*9)+5) * pi[2];
		D[2]=*((zone->ZGeom->DeltaInverse+P*9)+6) * pi[0] + *((zone->ZGeom->DeltaInverse+P*9)+7) * pi[1] + *((zone->ZGeom->DeltaInverse+P*9)+8) * pi[2];
		//printf("\n[Zone %s]$$[k=%d][j=%d][i=%d]$D[1]=%lf$D[2]=%lf$D[3]=%lf\n",zone->zonename,k,j,i,D[0],D[1],D[2]);
	}
	return 0.0;
}


