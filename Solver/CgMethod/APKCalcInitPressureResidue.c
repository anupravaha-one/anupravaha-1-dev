/**
This function is used to find the initial residual from the initial guess. This provides the residue vector for both BI-CGSTAB and Opti. GS algo.s
This function is for 3-D cases.

 **/


#include "APKCalcInitPressureResidue.h"


double CalcInitPressureResidue(Base *base)
{

	int NI,NJ,NK,i,j,k,numZones,iZone;
	
	int indexP,flagrho,cntDir;
	
	int K_Increment,J_Increment,I_Increment;
	
	double *Rho,*RHS;
	
	double *slabPressure[6],rho,Diag,Sum_Nb;
	
	double delta_t,*CoeffN[18],V[18],Vg;
	
	double *pressureDomain,*diffdia,*residue,*R_init;
	
	Zone *zone;
	
	numZones=base->Number_of_zones;
	delta_t = base->BaseParameters->TimeStep;
		
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		
		
		NI=*(zone->Size+0);
		NJ=*(zone->Size+1);
		NK=*(zone->Size+2);
		
		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;
		
		RHS = (zone->grid)->GRHS;
	
		diffdia = (zone->grid)->DiffusionCoeff_P;
		for(i=0;i<18;i++)
		{
			CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
		}
		
		Rho = zone->Flowsol->Rho;
		
		flagrho = zone->Flowsol->flagRho;
		
		pressureDomain = zone->Flowsol->P;
		
		for(cntDir=0;cntDir<6;cntDir++)
		{
			slabPressure[cntDir] = zone->Flowsol->SlabPress[cntDir];
		}
		
		residue = zone->Residue;
		
		R_init = zone->Residue_init;
		
		indexP =0;
	//--------------------------------------------- Boundary cells ------------------------------(1)--------
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			if(i!=0)
				V[0]=*(pressureDomain+indexP-I_Increment);
			else
				V[0] =*(slabPressure[1]+(NJ-1)*k+j);
			
			if(i!=(NI-2))
				V[1]=*(pressureDomain+indexP+I_Increment);
			else
				V[1]= *(slabPressure[0]+(NJ-1)*k+j);
			
			if(j!=0)
				V[2]=*(pressureDomain+indexP-J_Increment);
			else
				V[2]=*(slabPressure[3]+(NI-1)*k+i);
			
			if(j!=(NJ-2))
				V[3]=*(pressureDomain+indexP+J_Increment);
			else
				V[3]=*(slabPressure[2]+(NI-1)*k+i);
			
			if(k!=0)
				V[4]=*(pressureDomain+indexP-K_Increment);
			else
				V[4]=*(slabPressure[5]+(NI-1)*j+i);
		
			if(k!=(NK-2))
				V[5]=*(pressureDomain+indexP+K_Increment);
			else
				V[5]=*(slabPressure[4]+(NI-1)*j+i);
			
				
			if(!(i==0 && j==(NJ-2)))
			{
				if(i==0)
					V[6]=*(slabPressure[1]+(NJ-1)*k+(j+1));//from western slab
				else
				{
					if(j==(NJ-2))
						V[6]=*(slabPressure[2]+(NI-1)*k+(i-1));//from northern slab
					else
						V[6]=*(pressureDomain+indexP+J_Increment-I_Increment);
				}
			}
			else
				V[6]=0.0;
			
			if(!(i==0 && j==0))
			{
				if(i==0)
					V[7]=*(slabPressure[1]+(NJ-1)*k+(j-1));//from west slab
				else
				{
					if(j==0)
						V[7]=*(slabPressure[3]+(NI-1)*k+(i-1));//from southern slab
					else
						V[7]=*(pressureDomain+indexP-J_Increment-I_Increment);
				}
			}
			else
				V[7]=0.0;
			
			if(!(i==0 && k==(NK-2)))
			{
				if(i==0)
					V[8]=*(slabPressure[1]+(NJ-1)*(k+1)+j);
				else
				{	
					if(k==(NK-2))
						V[8]=*(slabPressure[4]+(NI-1)*j+(i-1));
					else
						V[8]=*(pressureDomain+indexP+K_Increment-I_Increment);
				}
				
			}
			else
				V[8]=0.0;
			
			if(!(i==0 && k==0))
			{
				if(i==0)
					V[9]=*(slabPressure[1]+(NJ-1)*(k-1)+j);
				else
				{
					if(k==0)
						V[9]=*(slabPressure[5]+(NI-1)*j+(i-1));
					else
						V[9]=*(pressureDomain+indexP-K_Increment-I_Increment);
				}
				
			}
			else
				V[9]=0.0;
				
			if(!(i==(NI-2) && j==(NJ-2)))
			{
				if(i==(NI-2))
					V[10]=*(slabPressure[0]+(NJ-1)*k+(j+1));
				else
				{
					if(j==(NJ-2))
						V[10]=*(slabPressure[2]+(NI-1)*k+(i+1));
					else
						V[10]=*(pressureDomain+indexP+J_Increment+I_Increment);
				}
			}
			else
				V[10]=0.0;
			
			if(!(i==(NI-2) && j==0))
			{
				if(i==(NI-2))
					V[11]=*(slabPressure[0]+(NJ-1)*k+(j-1));
				else
				{	
					if(j==0)
						V[11]=*(slabPressure[3]+(NI-1)*k+(i+1));
					else
						V[11]=*(pressureDomain+indexP-J_Increment+I_Increment);
				}
			}
			else
				V[11]=0.0;
			
			if(!(i==(NI-2) && k==(NK-2)))
			{
				if(i==(NI-2))
					V[12]=*(slabPressure[0]+(NJ-1)*(k+1)+j);
				else
				{
					if(k==(NK-2))
						V[12]=*(slabPressure[4]+(NI-1)*j+(i+1));
					else
						V[12]=*(pressureDomain+indexP+K_Increment+I_Increment);
				}
			}
			else
				V[12]=0.0;	
			
			if(!(i==(NI-2) && k==0))
			{
				if(i==(NI-2))
					V[13]=*(slabPressure[0]+(NJ-1)*(k-1)+j);
				else
				{
					if(k==0)
						V[13]=*(slabPressure[5]+(NI-1)*j+(i+1));
					else
						V[13]=*(pressureDomain+indexP-K_Increment+I_Increment);
				}
			}
			else
				V[13]=0.0;
			
			if(!(k==(NK-2) && j==(NJ-2)))
			{
				if(k==(NK-2))
					V[14]=*(slabPressure[4]+(NI-1)*(j+1)+i);
				else
				{
					if(j==(NJ-2))
						V[14]=*(slabPressure[2]+(NI-1)*(k+1)+i);
					else
						V[14]=*(pressureDomain+indexP+K_Increment+J_Increment);
				}
			}
			else
				V[14]=0.0;
			
			if(!(k==(NK-2) && j==0))
			{
				if(k==(NK-2))
					V[15]=*(slabPressure[4]+(NI-1)*(j-1)+i);
				else
				{
					if(j==0)
						V[15]=*(slabPressure[3]+(NI-1)*(k+1)+i);
					else
						V[15]=*(pressureDomain+indexP+K_Increment-J_Increment);
				}
			}
			else
				V[15]=0.0;
			
			if(!(k==0 && j==(NJ-2)))
			{
				if(k==0)
					V[16]=*(slabPressure[5]+(NI-1)*(j+1)+i);
				else
				{
					if(j==(NJ-2))
						V[16]=*(slabPressure[2]+(NI-1)*(k-1)+i);
					else
						V[16]=*(pressureDomain+indexP-K_Increment+J_Increment);
				}
			}
			else
				V[16]=0.0;
			
			if(!(k==0 && j==0))
			{
				if(k==0)
					V[17]=*(slabPressure[5]+(NI-1)*(j-1)+i);
				else
				{
					if(j==0)
						V[17]=*(slabPressure[3]+(NI-1)*(k-1)+i);
					else
						V[17]=*(pressureDomain+indexP-K_Increment-J_Increment);
				}
			}
			else
				V[17]=0.0;
			
			rho = *(Rho+indexP*flagrho);
			
			Vg = *(pressureDomain+indexP);
			
			Diag = *(diffdia+indexP) * Vg;
			
			Sum_Nb = *(CoeffN[0 ]+indexP)*V[0 ]+*(CoeffN[1 ]+indexP)*V[1 ]+*(CoeffN[2 ]+indexP)*V[2 ]+
				 *(CoeffN[3 ]+indexP)*V[3 ]+*(CoeffN[4 ]+indexP)*V[4 ]+*(CoeffN[5 ]+indexP)*V[5 ]+
				 *(CoeffN[6 ]+indexP)*V[6 ]+*(CoeffN[7 ]+indexP)*V[7 ]+*(CoeffN[8 ]+indexP)*V[8 ]+
				 *(CoeffN[9 ]+indexP)*V[9 ]+*(CoeffN[10]+indexP)*V[10]+*(CoeffN[11]+indexP)*V[11]+
				 *(CoeffN[12]+indexP)*V[12]+*(CoeffN[13]+indexP)*V[13]+*(CoeffN[14]+indexP)*V[14]+
				 *(CoeffN[15]+indexP)*V[15]+*(CoeffN[16]+indexP)*V[16]+*(CoeffN[17]+indexP)*V[17];
			
			
			*(residue+indexP) = *(RHS+indexP) - (delta_t/rho)*(Sum_Nb+Diag);
		
			if(base->BiCgflag==1)
				*(R_init+indexP) = *(residue+indexP);
			
// 			printf(" %d Resid %.18e\n",indexP,*(residue+indexP));
			//if(indexP<20)
			//	printf("Updated Pressure=%e\t%e\t%e\t%e\n",*(pressureDomain+indexP),*((zone->grid)->GRHS+indexP),(delta_t/rho)*(*LHS),(delta_t/rho)*(*RHS));
			
			indexP++;
		}
	}
	
	return(0);
}

