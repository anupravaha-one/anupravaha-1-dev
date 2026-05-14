/**
It is similar to the function "APKCalcInitPressureResidue.c" except for it is used to provide the initial residue for the potential flow case.
This function is used to find the initial residual from the initial guess. This provides the residue vector for both BI-CGSTAB and Opti. GS algo.s
This function is for 3-D cases.
 **/

#include "APKCalcInitConjVarResidue.h"


int CalcInitConjVarResidue(Base *base,int vartag)
{
	
	int NI,NJ,NK,cntDir,iZone,numZones;
	
	int K_Increment,J_Increment;//,ab;
	
	int INDEX,k,j,i,flagSpheat,flagrho,flagheat_var;
	
	double *variable,*slabVariable[6],*ptrVolume,*Energysource;//,*prevTstepVar;

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
		
		R_init = zone->Residue_init;
		
		P_i = zone->Pi;
		
		diffdia = (zone->grid)->DiffusionCoeff_P;
		for(i=0;i<18;i++)
		{
			CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
		}
		
		Rho = zone->PrevFlowsol->Rho;
		
		flagrho = zone->Flowsol->flagRho;
		
		if(vartag==6)
		{
			variable = zone->Flowsol->T;
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabT[cntDir];
			}
			
			//prevTstepVar = zone->PrevFlowsol->T;
		}
		else if(vartag==990)
		{		
			variable = zone->Flowsol->Flow_Potential;
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Flow_Potential[cntDir];
			}
			
			//prevTstepVar = zone->PrevFlowsol->Flow_Potential;
		}
		else
		{
			printf("No case found for the Vartag(%d), exiting from CalcInitConjVarResidue2d\n",vartag);
			exit(0);
		}	
		
		K_Increment = (NI-1)*(NJ-1);
		J_Increment = NI-1;
		
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
				V[0] =*(slabVariable[1]+(NJ-1)*k+j);
				
			if(i!=(NI-2))
				V[1]=*(variable+INDEX+1);
			else
				V[1]= *(slabVariable[0]+(NJ-1)*k+j);
			
			if(j!=0)
				V[2]=*(variable+INDEX-J_Increment);
			else
				V[2]=*(slabVariable[3]+(NI-1)*k+i);
			
			if(j!=(NJ-2))
				V[3]=*(variable+INDEX+J_Increment);
			else
				V[3]=*(slabVariable[2]+(NI-1)*k+i);
			
			if(k!=0)
				V[4]=*(variable+INDEX-K_Increment);
			else
				V[4]=*(slabVariable[5]+(NI-1)*j+i);
			
			if(k!=(NK-2))
				V[5]=*(variable+INDEX+K_Increment);
			else
				V[5]=*(slabVariable[4]+(NI-1)*j+i);
			
				
			if(!(i==0 && j==(NJ-2)))
			{
				if(i==0)
					V[6]=*(slabVariable[1]+(NJ-1)*k+(j+1));//from western slab
				else
				{
					if(j==(NJ-2))
						V[6]=*(slabVariable[2]+(NI-1)*k+(i-1));//from northern slab
					else
						V[6]=*(variable+INDEX+J_Increment-1);
				}
			}
			else
				V[6]=0.0;
			
			if(!(i==0 && j==0))
			{
				if(i==0)
					V[7]=*(slabVariable[1]+(NJ-1)*k+(j-1));//from west slab
				else
				{
					if(j==0)
						V[7]=*(slabVariable[3]+(NI-1)*k+(i-1));//from southern slab
					else
						V[7]=*(variable+INDEX-J_Increment-1);
				}
			}
			else
				V[7]=0.0;
			
			if(!(i==0 && k==(NK-2)))
			{
				if(i==0)
					V[8]=*(slabVariable[1]+(NJ-1)*(k+1)+j);
				else
				{	
					if(k==(NK-2))
						V[8]=*(slabVariable[4]+(NI-1)*j+(i-1));
					else
						V[8]=*(variable+INDEX+K_Increment-1);
				}
				
			}
			else
				V[8]=0.0;
			
			if(!(i==0 && k==0))
			{
				if(i==0)
					V[9]=*(slabVariable[1]+(NJ-1)*(k-1)+j);
				else
				{
					if(k==0)
						V[9]=*(slabVariable[5]+(NI-1)*j+(i-1));
					else
						V[9]=*(variable+INDEX-K_Increment-1);
				}
				
			}
			else
				V[9]=0.0;
				
			if(!(i==(NI-2) && j==(NJ-2)))
			{
				if(i==(NI-2))
					V[10]=*(slabVariable[0]+(NJ-1)*k+(j+1));
				else
				{
					if(j==(NJ-2))
						V[10]=*(slabVariable[2]+(NI-1)*k+(i+1));
					else
						V[10]=*(variable+INDEX+J_Increment+1);
				}
			}
			else
				V[10]=0.0;
			
			if(!(i==(NI-2) && j==0))
			{
				if(i==(NI-2))
					V[11]=*(slabVariable[0]+(NJ-1)*k+(j-1));
				else
				{	
					if(j==0)
						V[11]=*(slabVariable[3]+(NI-1)*k+(i+1));
					else
						V[11]=*(variable+INDEX-J_Increment+1);
				}
			}
			else
				V[11]=0.0;
			
			if(!(i==(NI-2) && k==(NK-2)))
			{
				if(i==(NI-2))
					V[12]=*(slabVariable[0]+(NJ-1)*(k+1)+j);
				else
				{
					if(k==(NK-2))
						V[12]=*(slabVariable[4]+(NI-1)*j+(i+1));
					else
						V[12]=*(variable+INDEX+K_Increment+1);
				}
			}
			else
				V[12]=0.0;	
			
			if(!(i==(NI-2) && k==0))
			{
				if(i==(NI-2))
					V[13]=*(slabVariable[0]+(NJ-1)*(k-1)+j);
				else
				{
					if(k==0)
						V[13]=*(slabVariable[5]+(NI-1)*j+(i+1));
					else
						V[13]=*(variable+INDEX-K_Increment+1);
				}
			}
			else
				V[13]=0.0;
			
			if(!(k==(NK-2) && j==(NJ-2)))
			{
				if(k==(NK-2))
					V[14]=*(slabVariable[4]+(NI-1)*(j+1)+i);
				else
				{
					if(j==(NJ-2))
						V[14]=*(slabVariable[2]+(NI-1)*(k+1)+i);
					else
						V[14]=*(variable+INDEX+J_Increment+K_Increment);
				}
			}
			else
				V[14]=0.0;
			
			if(!(k==(NK-2) && j==0))
			{
				if(k==(NK-2))
					V[15]=*(slabVariable[4]+(NI-1)*(j-1)+i);
				else
				{
					if(j==0)
						V[15]=*(slabVariable[3]+(NI-1)*(k+1)+i);
					else
						V[15]=*(variable+INDEX-J_Increment+K_Increment);
				}
			}
			else
				V[15]=0.0;
			
			if(!(k==0 && j==(NJ-2)))
			{
				if(k==0)
					V[16]=*(slabVariable[5]+(NI-1)*(j+1)+i);
				else
				{
					if(j==(NJ-2))
						V[16]=*(slabVariable[2]+(NI-1)*(k-1)+i);
					else
						V[16]=*(variable+INDEX+J_Increment-K_Increment);
				}
			}
			else
				V[16]=0.0;
			
			if(!(k==0 && j==0))
			{
				if(k==0)
					V[17]=*(slabVariable[5]+(NI-1)*(j-1)+i);
				else
				{
					if(j==0)
						V[17]=*(slabVariable[3]+(NI-1)*(k-1)+i);
					else
						V[17]=*(variable+INDEX-J_Increment-K_Increment);
				}
			}
			else
				V[17]=0.0;
			
			Vg = *(variable+INDEX);
			rho = *(Rho+INDEX*flagrho);
		
			Diag = *(diffdia+INDEX) * Vg;
		
			Sum_Nb = 
				*(CoeffN[0 ]+INDEX)*V[0 ]+*(CoeffN[1 ]+INDEX)*V[1 ]+*(CoeffN[2 ]+INDEX)*V[2 ]+
				*(CoeffN[3 ]+INDEX)*V[3 ]+*(CoeffN[4 ]+INDEX)*V[4 ]+*(CoeffN[5 ]+INDEX)*V[5 ]+
				*(CoeffN[6 ]+INDEX)*V[6 ]+*(CoeffN[7 ]+INDEX)*V[7 ]+*(CoeffN[8 ]+INDEX)*V[8 ]+
				*(CoeffN[9 ]+INDEX)*V[9 ]+*(CoeffN[10]+INDEX)*V[10]+*(CoeffN[11]+INDEX)*V[11]+
				*(CoeffN[12]+INDEX)*V[12]+*(CoeffN[13]+INDEX)*V[13]+*(CoeffN[14]+INDEX)*V[14]+
				*(CoeffN[15]+INDEX)*V[15]+*(CoeffN[16]+INDEX)*V[16]+*(CoeffN[17]+INDEX)*V[17];
	
			if(zone->flagHeatGeneration==1)
				source=(1.0/rho)*(*(Energysource+INDEX*flagheat_var))*(*(ptrVolume+INDEX));
	
			*(residue+INDEX) = source - (Sum_Nb+Diag)/rho;
			
			*(R_init+INDEX) = *(residue+INDEX);
			
			INDEX++;
		}
	}
	
	return(0);
}


