#include "APKUpdatePressure2d.h"

/**
The function UpdatePressure2d() does the Point Gauss Seidel Iterations for the Pressure Poisson equation for 2D domain.
*/

double UpdatePressure2d(Base *base,Zone *zone)
{

	double residZone=0.0,*pressureDomain,*diffdia;
	double *LHS,*RHS,resid,rho;
	int NI,NJ,NK,i,j,indexP=0,*ptrSize,flagrho,cntDir;
	double delta_t,*CoeffN[18],V[18],Vg,omega=1.0,*slabPressure[6];
	
	//int k;
	int J_Increment;
	
	LHS=(double *)malloc(sizeof(double));
	RHS=(double *)malloc(sizeof(double));
	
	omega = *(base->BaseParameters->OMEGA+1);

	diffdia = (zone->grid)->DiffusionCoeff_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
	}
	delta_t = base->BaseParameters->TimeStep;
	pressureDomain = zone->Flowsol->P;
	ptrSize = zone->Size ;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	J_Increment = NI-1;
	
	for(cntDir=0;cntDir<6;cntDir++)
	{
		slabPressure[cntDir] = zone->Flowsol->SlabPress[cntDir];
	}
	
	flagrho = zone->Flowsol->flagRho;
	indexP=0;
//--------------------------------------------- Boundary cells -----------------------------------(1)--------
	//for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		if(i!=0)
			V[0]=*(pressureDomain+indexP-1);
		else
			V[0] =*(slabPressure[1]+j);
		
		if(i!=(NI-2))
			V[1]=*(pressureDomain+indexP+1);
		else
			V[1]= *(slabPressure[0]+j);
		
		if(j!=0)
			V[2]=*(pressureDomain+indexP-J_Increment);
		else
			V[2]=*(slabPressure[3]+i);
		
		if(j!=(NJ-2))
			V[3]=*(pressureDomain+indexP+J_Increment);
		else
			V[3]=*(slabPressure[2]+i);
		
			
		if(!(i==0 && j==(NJ-2)))
		{
			if(i==0)
				V[6]=*(slabPressure[1]+(j+1));//from western slab
			else
			{
				if(j==(NJ-2))
					V[6]=*(slabPressure[2]+(i-1));//from northern slab
				else
					V[6]=*(pressureDomain+indexP+J_Increment-1);
			}
		}
		else
			V[6]=0.0;
		
		if(!(i==0 && j==0))
		{
			if(i==0)
				V[7]=*(slabPressure[1]+(j-1));//from west slab
			else
			{
				if(j==0)
					V[7]=*(slabPressure[3]+(i-1));//from southern slab
				else
					V[7]=*(pressureDomain+indexP-J_Increment-1);
			}
		}
		else
			V[7]=0.0;
		
			
		if(!(i==(NI-2) && j==(NJ-2)))
		{
			if(i==(NI-2))
				V[10]=*(slabPressure[0]+(j+1));
			else
			{
				if(j==(NJ-2))
					V[10]=*(slabPressure[2]+(i+1));
				else
					V[10]=*(pressureDomain+indexP+J_Increment+1);
			}
		}
		else
			V[10]=0.0;
		
		if(!(i==(NI-2) && j==0))
		{
			if(i==(NI-2))
				V[11]=*(slabPressure[0]+(j-1));
			else
			{	
				if(j==0)
					V[11]=*(slabPressure[3]+(i+1));
				else
					V[11]=*(pressureDomain+indexP-J_Increment+1);
			}
		}
		else
			V[11]=0.0;
		
		if(base->flagMPFV==1)
		{
			rho =1.0;
		}
		if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
		{
			rho=1.0;
		}
		else
		{
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
		}

		Vg = *(pressureDomain+indexP);
		
		//MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		*LHS = *(diffdia+indexP) * Vg;
		*RHS = *(CoeffN[0]+indexP) *V[0] + *(CoeffN[1]+indexP) *V[1] + *(CoeffN[2]+indexP) *V[2] +*(CoeffN[3]+indexP) *V[3] + *(CoeffN[6]+indexP) *V[6] + *(CoeffN[7]+indexP) *V[7] + *(CoeffN[10]+indexP) *V[10]+*(CoeffN[11]+indexP) *V[11];
		//printf("\n LHS")	
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		//printf("\n *(pressureDomain+indexP)=%e",*(pressureDomain+indexP));
		residZone=residZone+resid*resid;
		
		
		indexP++;
	}
	//residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	
	free(LHS);
	free(RHS);
	return(residZone);
}

