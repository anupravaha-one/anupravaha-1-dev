/**
This function is used to find the initial residual from the initial guess. This provides the residue vector for both BI-CGSTAB and Opti. GS algo.s
This function is for 2-D cases.

 **/



#include "APKCalcInitPressureResidue2d.h"

double CalcInitPressureResidue2d(Base *base)
{

	double *residue,*pressureDomain,*diffdia,*R_init,*P_i;
	
	double rho,Diag,Sum_Nb,delta_t,*Rho,*RHS;//,rms=0.0;
	
	double *CoeffN[18],V[18],Vg,*slabPressure[6];
	
	int NI,NJ,NK,i,j,*ptrSize,flagrho,iZone,numZones;
	
	int indexP,J_Increment,TCell;//,ab;//,k;
	
	Zone *zone;
	
	numZones=base->Number_of_zones;
	
	delta_t = base->BaseParameters->TimeStep;
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		
		ptrSize = zone->Size ;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		TCell=(NI-1)*(NJ-1);
		
		residue = zone->Residue;
		
		R_init = zone->Residue_init;
		
		P_i = zone->Pi;
			
		pressureDomain = zone->Flowsol->P;
		
		diffdia = (zone->grid)->DiffusionCoeff_P;
		
		for(i=0;i<18;i++)
		{
			CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
		}
		
		for(i=0;i<6;i++)
		{
			slabPressure[i] = zone->Flowsol->SlabPress[i];
			
		}
		
		RHS=(zone->grid)->GRHS;
		
		Rho = zone->Flowsol->Rho;
		
		flagrho = zone->Flowsol->flagRho;
		
		J_Increment = NI-1;
		
		indexP=0; 
		
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
			
			
			rho = *(Rho+indexP*flagrho);
			
			Vg = *(pressureDomain+indexP);
			
			Diag = *(diffdia+indexP) * Vg;
			
			Sum_Nb = *(CoeffN[0]+indexP) *V[0] + *(CoeffN[1]+indexP) *V[1] + *(CoeffN[2]+indexP) *V[2] +*(CoeffN[3]+indexP) *V[3] + *(CoeffN[6]+indexP) *V[6] + *(CoeffN[7]+indexP) *V[7] + *(CoeffN[10]+indexP) *V[10]+*(CoeffN[11]+indexP) *V[11];
	
			
			*(residue+indexP) = *(RHS+indexP)-(delta_t/rho)*(Sum_Nb+Diag);
			
			if(base->BiCgflag==1)
				*(R_init+indexP) = *(residue+indexP);
			
			indexP++;
		}
	}		
			//*(P_i+indexP) = *(residue+indexP);
			
// 			rms += *(residue+indexP) * *(residue+indexP);
			
// 			printf(" %d Resid %.18e\n",indexP,*(residue+indexP));
				
/*
		if(indexP==16)
		{	
			printf(" %d Resid %.18e\n",indexP,*(residue+indexP));
				
			printf("Diff %.18e Vg %.18e\n",*(diffdia+indexP),Vg);
			for(ab=0;ab<4;ab++)
			{
				printf("Coeff %.18e Vg %.18e\n",*(CoeffN[ab]+indexP),V[ab]);
			}
			exit(0);
		}
*/			
			//printf("%d %lf %lf %lf %lf %lf \n",j*(NI-1)+i,*(diffdia+indexP),*(CoeffN[0]+indexP),*(CoeffN[1]+indexP),*(CoeffN[2]+indexP),*(CoeffN[3]+indexP));
			
			
// 			indexP++;
// 		}
// 		exit(0);
// 	}
	//exit(0);
	
	return(0);
}
