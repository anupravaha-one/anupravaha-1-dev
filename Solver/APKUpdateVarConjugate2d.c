#include "APKUpdateVarConjugate2d.h"

/**
This is exactly similar to "UpdateVar" but is only used for Conjugate case for solving temperature.
This function solves for temperature using Gauss-Seidal point method and returns the
residual of temperature for that particular iteration. This is implemented for 2D
case having either constant or variable diffusion coefficient.
*/


double UpdateVarConjugate2d(Base *base,Zone *zone,int vartag)
{
	int NI,NJ,NK,*ptrSize,cntDir;
	double AX,RS,*CoeffN[18],*diffdia;
	int k,j,i,flagSpheat,flagrho;
	double residue,res=0.0,source=0.0,omega=1.8,Vg,V[18],*LHS,*RHS;
	double *variable,*prevTstepVar,*slabVariable[6],*ptrVolume;

	//int indexP,index[3];
	omega = *(base->BaseParameters->OMEGA+3);

	//double DiagSource=0.0;
	double rho;
	
	source=0.0;
	ptrSize= zone->Size;
	ptrVolume= zone->ZGeom->Volume;

	LHS = malloc(sizeof(double));
	RHS = malloc(sizeof(double));
	switch(vartag)
	{
		case 6:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabT[cntDir];
			}
			diffdia = (zone->grid)->DiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid)->DiffusionCoeffPhi_nb[i];
				V[i] = 0.0;
			}
			variable = zone->Flowsol->T;
			prevTstepVar = zone->PrevFlowsol->T;
		break;
		case 990:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Flow_Potential[cntDir];
			}
			diffdia = (zone->grid)->DiffusionCoeff_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
			}
			variable = zone->Flowsol->Flow_Potential;
			prevTstepVar = zone->PrevFlowsol->Flow_Potential;
		break;
		case 59://rajendra
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Psi_MHD[cntDir];
			}
			diffdia = (zone->grid)->DiffusionCoeff_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
			}
			variable = zone->Flowsol->Psi;
			prevTstepVar = zone->DeferredCorrCoeffDummy_P;//rajendra
			/**Don't get confused.  This is just given to assure that the vlaues in prevTstepVar are zero**/
			omega=0.8;
		break;
		default:
			printf("Wrong Tag supplied into UpdateVarConjugate2d.c...\n");
			exit(1);
		break;
	}
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	int IMAX,JMAX,I_Increment,J_Increment;
	
 	IMAX = NI-2;
 	JMAX = NJ-2;
	I_Increment =1;
	J_Increment =NI-1;
	
	
	flagSpheat = zone->Flowsol->flagSpecificHeat;
	flagrho = zone->Flowsol->flagRho;
	
	int INDEX=0;
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
		
		rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
		Vg = *(variable+INDEX);
		
		*LHS = *(diffdia+INDEX) * Vg;
		
 		*RHS = *(CoeffN[0]+INDEX) *V[0] + *(CoeffN[1]+INDEX) *V[1] + *(CoeffN[2]+INDEX) *V[2] +*(CoeffN[3]+INDEX) *V[3] + *(CoeffN[6]+INDEX) *V[6] + *(CoeffN[7]+INDEX) *V[7] + *(CoeffN[10]+INDEX) *V[10]+*(CoeffN[11]+INDEX) *V[11];
		
		if(zone->flagHeatGeneration==1)
		source=(1.0/rho)*(*((zone->Flowsol)->Energysource+INDEX*(zone->flagHeat_con_var)))*(*(ptrVolume+INDEX));

		AX = -(1.0/rho)*(*LHS+*RHS);
		
		RS = (1.0*source);
 		
		if(vartag==59)//rajendra
		{
		source=*((zone->grid)->GRHS+INDEX);
		
		AX = (*LHS+*RHS);
 		RS = (1.0*source);
		}
		residue=RS-AX;

		
		if(vartag==59)
		{
		*(variable+INDEX) = *(variable+INDEX) + omega*residue/(*(diffdia+INDEX));
		
		}
		else
		*(variable+INDEX) = *(variable+INDEX) + omega*residue/((-1.0/rho)*(*(diffdia+INDEX)));
		res=res+(residue*residue);
			
		
		INDEX++;
	}
	return(res);
}


