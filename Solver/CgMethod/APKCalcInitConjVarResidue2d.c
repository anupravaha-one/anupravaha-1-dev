/**
It is similar to the function "APKCalcInitPressureResidue.c" except for it is used to provide the initial residue for the potential flow case.
This function is used to find the initial residual from the initial guess. This provides the residue vector for both BI-CGSTAB and Opti. GS algo.s
This function is for 2-D cases.
 **/



#include "APKCalcInitConjVarResidue2d.h"


int CalcInitConjVarResidue2d(Base *base,int vartag)
{
	int NI,NJ,NK,cntDir,iZone,numZones;
	
	int IMAX,JMAX,I_Increment,J_Increment;//,ab;
	
	int INDEX,k,j,i,flagSpheat,flagrho,flagheat_var;
	
	double *variable,*prevTstepVar,*slabVariable[6],*ptrVolume,*Energysource;

	double *residue,*R_init,*P_i;
	
	double rho,Diag,Sum_Nb,*Rho,source=0.0;
	
	double *diffdia,*CoeffN[18],V[18],Vg;
	
	
	Zone *zone;
	
	numZones=base->Number_of_zones;
	
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		
		NI=*(zone->Size+0);
		NJ=*(zone->Size+1);
		NK=*(zone->Size+2);
		
		ptrVolume= zone->ZGeom->Volume;
		
		residue = zone->Residue;
		
		if(base->BiCgflag==1)
		{
			R_init = zone->Residue_init;
			
			P_i = zone->Pi;
		}
		
		diffdia = (zone->grid)->DiffusionCoeff_P;
		for(i=0;i<18;i++)
		{
			CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
		}
		
		Rho = zone->Flowsol->Rho;
		
		flagrho = zone->Flowsol->flagRho;
		
		if(vartag==6)
		{
			variable = zone->Flowsol->T;
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabT[cntDir];
			}
			
			prevTstepVar = zone->PrevFlowsol->T;
		}
		else if(vartag==990)
		{		
			variable = zone->Flowsol->Flow_Potential;
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Flow_Potential[cntDir];
			}
			
			prevTstepVar = zone->PrevFlowsol->Flow_Potential;
		}
		else
		{
			printf("No case found for the Vartag(%d), exiting from CalcInitConjVarResidue2d\n",vartag);
			exit(0);
		}	
		
		
		IMAX = NI-2;
		JMAX = NJ-2;
		
		I_Increment =1;
		J_Increment =NI-1;
		
		Energysource = (zone->Flowsol)->Energysource;
		
		flagheat_var = zone->flagHeat_con_var;
		
		flagSpheat = zone->Flowsol->flagSpecificHeat;
		
		INDEX=0;
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
	
			if(i!=0)
				V[0]=*(variable+INDEX-1);
			else
				V[0] =*(slabVariable[1]+j);
			
			if(i!=IMAX)
				V[1]=*(variable+INDEX+1);
			else
				V[1]= *(slabVariable[0]+j);
			
			if(j!=0)
				V[2]=*(variable+INDEX-J_Increment);
			else
				V[2]=*(slabVariable[3]+i);
			
			if(j!=JMAX)
				V[3]=*(variable+INDEX+J_Increment);
			else
				V[3]=*(slabVariable[2]+i);
			
			if(!(i==0 && j==JMAX))
			{
				if(!(i==0 || j==JMAX))
					V[6]=*(variable+INDEX+J_Increment-1);
				else if(i==0)
					V[6]=*(slabVariable[1]+(j+1));//from western slab
				else if(j==JMAX)
					V[6]=*(slabVariable[2]+(i-1));//from northern slab
					
			}
			else
				V[6]=0.0;
			
			if(!(i==0 && j==0))
			{
				if(!(i==0 || j==0))
					V[7]=*(variable+INDEX-J_Increment-1);
				else if(i==0)
					V[7]=*(slabVariable[1]+(j-1));//from west slab
				else if(j==0)
					V[7]=*(slabVariable[3]+(i-1));//from southern slab
			}
			else
				V[7]=0.0;
			
			if(!(i==IMAX && j==JMAX))
			{
				if(!(i==IMAX || j==JMAX))
					V[10]=*(variable+INDEX+J_Increment+1);
				else if(i==IMAX)
					V[10]=*(slabVariable[0]+(j+1));
				else if(j==JMAX)
					V[10]=*(slabVariable[2]+(i+1));
	
			}
			else
				V[10]=0.0;
			
			if(!(i==IMAX && j==0))
			{
				if(!(i==IMAX || j==0))
					V[11]=*(variable+INDEX-J_Increment+1);
				else if(i==IMAX)
					V[11]=*(slabVariable[0]+(j-1));
				else if(j==0)
					V[11]=*(slabVariable[3]+(i+1));
			}
			else
				V[11]=0.0;
			
			rho = *(Rho+INDEX*flagrho);
			//rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			
			Vg = *(variable+INDEX);
			
			Diag = *(diffdia+INDEX) * Vg;
			
			Sum_Nb = *(CoeffN[0]+INDEX) *V[0] + *(CoeffN[1]+INDEX) *V[1] + *(CoeffN[2]+INDEX) *V[2] +*(CoeffN[3]+INDEX) *V[3] + *(CoeffN[6]+INDEX) *V[6] + *(CoeffN[7]+INDEX) *V[7] + *(CoeffN[10]+INDEX) *V[10]+*(CoeffN[11]+INDEX) *V[11];
			
			if(zone->flagHeatGeneration==1)
				source=(1.0/rho)*(*(Energysource+INDEX*flagheat_var))*(*(ptrVolume+INDEX));
	
			*(residue+INDEX) = source - (Sum_Nb+Diag)/rho;//*(-(1.0/rho));
			
			if(base->BiCgflag==1)
				*(R_init+INDEX) = *(residue+INDEX);
			
			///*(P_i+INDEX) = *(residue+INDEX);//
			
			///rms += *(residue+INDEX) * *(residue+INDEX);//
			
			INDEX++;
		}
	}
	return(0);
}


