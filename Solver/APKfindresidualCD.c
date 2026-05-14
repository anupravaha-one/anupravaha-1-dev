#include "APKfindresidualCD.h"

/**
This function solves for the mass velocities only using Gauss-Seidal point method and returns the
residual of the mass velocities for that particular iteration. It is called only for 3D case in which central difference
conevction scheme is used and having either constant or variable diffusion coefficient.
*/

double findresidualCD(Base *base,Zone *zone,int vartag)
{

	int NI,NJ,NK,index[3],INDEX,numBC,k,j,i,*ptrSize,g,cntDir;
	double vol,velP,AX,RS,residue,res=0.0,DiagSource,source;
	double delta_t,omega=0.6,*slabVelocity[6],diacoeff;
	double Vg,V[18],*LHS,*RHS,DefCor,*CoeffN[18],*dia,*DefCoeffN[12],*Defdia,defer[6];
	double *massVelDomain,*prevTstepVel,*prevIterMassVel;
	
	
	int K_Increment,J_Increment,I_Increment;
	
	//ZoneBoundary *zonebc;
	//int nBc;
	omega = *(base->BaseParameters->OMEGA);
	ptrSize= zone->Size;
	
	g = 0;
	switch(vartag)
	{
		case 1:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffPhi_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabUmass[cntDir];
			}
			massVelDomain = zone->Flowsol->umass;
			prevTstepVel = zone->PrevFlowsol->u;
		break;
		case 2:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffPhi_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabVmass[cntDir];
			}
			massVelDomain = zone->Flowsol->vmass;
			prevTstepVel = zone->PrevFlowsol->v;
		break;
		case 3:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffPhi_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabWmass[cntDir];
			}
			massVelDomain = zone->Flowsol->wmass;
			prevTstepVel = zone->PrevFlowsol->w;
		break;
		//Added by sushil for multiphase flow problem
		case 31:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffPhi_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabVxstar[cntDir];
			}
			massVelDomain = zone->Flowsol->Momentumxstar;
			prevTstepVel = zone->PrevFlowsol->Momentumx;
			break;
		case 32:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffPhi_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabVystar[cntDir];
			}
			massVelDomain = zone->Flowsol->Momentumystar;
			prevTstepVel = zone->PrevFlowsol->Momentumy;
			break;
		
		case 33:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffPhi_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffPhi_nb[i];
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabVzstar[cntDir];
			}
			massVelDomain = zone->Flowsol->Momentumzstar;
			prevTstepVel = zone->PrevFlowsol->Momentumz;
//for dilute gas-droplet flows with evaporation
		case 62:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] =zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabUDmass[cntDir];
			}
			massVelDomain = zone->Flowsol->U_Dmass;
			prevTstepVel = zone->PrevFlowsol->U_Disperse;
		break;
		case 63:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabVDmass[cntDir];
			}
			massVelDomain = zone->Flowsol->V_Dmass;
			prevTstepVel = zone->PrevFlowsol->V_Disperse;
		break;
		case 64:
			dia = (zone->grid+g)->ConvDiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid+g)->ConvDiffusionCoeffPhi_nb[i];
			}
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVelocity[cntDir] = zone->Flowsol->SlabWDmass[cntDir];
			}
			massVelDomain = zone->Flowsol->W_Dmass;
			prevTstepVel = zone->PrevFlowsol->W_Disperse;
		break;
		//end of multiphase flow 
		
	}
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	K_Increment = (NI-1)*(NJ-1);
	J_Increment = (NI-1);
	I_Increment = 1;
	
	prevIterMassVel = zone->Flowsol->prevIterVariable;
	
	INDEX=0;
	
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		*(prevIterMassVel+INDEX) = *(massVelDomain+INDEX);
		INDEX++;
	}
	
	delta_t = base->BaseParameters->TimeStep;
	
	LHS=(double *)malloc(sizeof(double));
	RHS=(double *)malloc(sizeof(double));
	
	numBC= zone->Number_of_BC_Patches;
	double MULTIPLE=1.0;
	
	if(base->BaseParameters->TimeStepingSchemeTag == 1)
		MULTIPLE = 0.5;
	
	INDEX=0;
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		if(i!=0)
		{
			V[0]=*(massVelDomain+INDEX-1);
			defer[0]=*(prevIterMassVel+INDEX-1);
		}
		else
		{
			V[0] =*(slabVelocity[1]+(NJ-1)*k+j);
			defer[0] = V[0];
		}
		if(i!=(NI-2))
		{
			V[1]=*(massVelDomain+INDEX+1);
			defer[1]=*(prevIterMassVel+INDEX+1);
		}
		else
		{
			V[1]= *(slabVelocity[0]+(NJ-1)*k+j);
			defer[1]= V[1];
		}
	
		if(j!=0)
		{
			V[2]=*(massVelDomain+INDEX-J_Increment);
			defer[2] = *(prevIterMassVel+INDEX-J_Increment);
		}
		else
		{
			V[2]=*(slabVelocity[3]+(NI-1)*k+i);
			defer[2]=V[2];
		}
	
		if(j!=(NJ-2))
		{
			V[3]=*(massVelDomain+INDEX+J_Increment);
			defer[3] =*(prevIterMassVel+INDEX+J_Increment);
		}
		else
		{
			V[3]=*(slabVelocity[2]+(NI-1)*k+i);
			defer[3] = V[3];
		}
	
		if(k!=0)
		{
			V[4]=*(massVelDomain+INDEX-K_Increment);
			defer[4] = *(prevIterMassVel+INDEX-K_Increment);
		}
		else
		{
			V[4]=*(slabVelocity[5]+(NI-1)*j+i);
			defer[4] = V[4];
		}
	
		if(k!=(NK-2))
		{
			V[5]=*(massVelDomain+INDEX+K_Increment);
			defer[5]=*(prevIterMassVel+INDEX+K_Increment);
		}
		else
		{
			V[5]=*(slabVelocity[4]+(NI-1)*j+i);
			defer[5] = V[5];
		}
		
			
		if(!(i==0 && j==(NJ-2)))
		{
			if(i==0)
				V[6]=*(slabVelocity[1]+(NJ-1)*k+(j+1));//from western slab
			else
			{
				if(j==(NJ-2))
					V[6]=*(slabVelocity[2]+(NI-1)*k+(i-1));//from northern slab
				else
					V[6]=*(massVelDomain+INDEX+J_Increment-1);
			}
		}
		else
			V[6]=0.0;
		
		if(!(i==0 && j==0))
		{
			if(i==0)
				V[7]=*(slabVelocity[1]+(NJ-1)*k+(j-1));//from west slab
			else
			{
				if(j==0)
					V[7]=*(slabVelocity[3]+(NI-1)*k+(i-1));//from southern slab
				else
					V[7]=*(massVelDomain+INDEX-J_Increment-1);
			}
		}
		else
			V[7]=0.0;
		
		if(!(i==0 && k==(NK-2)))
		{
			if(i==0)
				V[8]=*(slabVelocity[1]+(NJ-1)*(k+1)+j);
			else
			{	
				if(k==(NK-2))
					V[8]=*(slabVelocity[4]+(NI-1)*j+(i-1));
				else
					V[8]=*(massVelDomain+INDEX+K_Increment-1);
			}
			
		}
		else
			V[8]=0.0;
		
		if(!(i==0 && k==0))
		{
			if(i==0)
				V[9]=*(slabVelocity[1]+(NJ-1)*(k-1)+j);
			else
			{
				if(k==0)
					V[9]=*(slabVelocity[5]+(NI-1)*j+(i-1));
				else
					V[9]=*(massVelDomain+INDEX-K_Increment-1);
			}
			
		}
		else
			V[9]=0.0;
			
		if(!(i==(NI-2) && j==(NJ-2)))
		{
			if(i==(NI-2))
				V[10]=*(slabVelocity[0]+(NJ-1)*k+(j+1));
			else
			{
				if(j==(NJ-2))
					V[10]=*(slabVelocity[2]+(NI-1)*k+(i+1));
				else
					V[10]=*(massVelDomain+INDEX+J_Increment+1);
			}
		}
		else
			V[10]=0.0;
		
		if(!(i==(NI-2) && j==0))
		{
			if(i==(NI-2))
				V[11]=*(slabVelocity[0]+(NJ-1)*k+(j-1));
			else
			{	
				if(j==0)
					V[11]=*(slabVelocity[3]+(NI-1)*k+(i+1));
				else
					V[11]=*(massVelDomain+INDEX-J_Increment+1);
			}
		}
		else
			V[11]=0.0;
		
		if(!(i==(NI-2) && k==(NK-2)))
		{
			if(i==(NI-2))
				V[12]=*(slabVelocity[0]+(NJ-1)*(k+1)+j);
			else
			{
				if(k==(NK-2))
					V[12]=*(slabVelocity[4]+(NI-1)*j+(i+1));
				else
					V[12]=*(massVelDomain+INDEX+K_Increment+1);
			}
		}
		else
			V[12]=0.0;	
		
		if(!(i==(NI-2) && k==0))
		{
			if(i==(NI-2))
				V[13]=*(slabVelocity[0]+(NJ-1)*(k-1)+j);
			else
			{
				if(k==0)
					V[13]=*(slabVelocity[5]+(NI-1)*j+(i+1));
				else
					V[13]=*(massVelDomain+INDEX-K_Increment+1);
			}
		}
		else
			V[13]=0.0;
		
		if(!(k==(NK-2) && j==(NJ-2)))
		{
			if(k==(NK-2))
				V[14]=*(slabVelocity[4]+(NI-1)*(j+1)+i);
			else
			{
				if(j==(NJ-2))
					V[14]=*(slabVelocity[2]+(NI-1)*(k+1)+i);
				else
					V[14]=*(massVelDomain+INDEX+J_Increment+K_Increment);
			}
		}
		else
			V[14]=0.0;
		
		if(!(k==(NK-2) && j==0))
		{
			if(k==(NK-2))
				V[15]=*(slabVelocity[4]+(NI-1)*(j-1)+i);
			else
			{
				if(j==0)
					V[15]=*(slabVelocity[3]+(NI-1)*(k+1)+i);
				else
					V[15]=*(massVelDomain+INDEX-J_Increment+K_Increment);
			}
		}
		else
			V[15]=0.0;
		
		if(!(k==0 && j==(NJ-2)))
		{
			if(k==0)
				V[16]=*(slabVelocity[5]+(NI-1)*(j+1)+i);
			else
			{
				if(j==(NJ-2))
					V[16]=*(slabVelocity[2]+(NI-1)*(k-1)+i);
				else
					V[16]=*(massVelDomain+INDEX+J_Increment-K_Increment);
			}
		}
		else
			V[16]=0.0;
		
		if(!(k==0 && j==0))
		{
			if(k==0)
				V[17]=*(slabVelocity[5]+(NI-1)*(j-1)+i);
			else
			{
				if(j==0)
					V[17]=*(slabVelocity[3]+(NI-1)*(k-1)+i);
				else
					V[17]=*(massVelDomain+INDEX-J_Increment-K_Increment);
			}
		}
		else
			V[17]=0.0;
		
		index[0]=k,index[1]=j,index[2]=i;
	
		vol = *(zone->ZGeom->Volume+INDEX);
		velP=*(prevTstepVel+INDEX);
	
		Vg = *(massVelDomain+INDEX);
		
		source = *(zone->Flowsol->Source_term+INDEX);
		DiagSource = *(zone->Flowsol->Diag_Source+INDEX);
		source = source - DiagSource*Vg;
		
		*LHS = *(dia+INDEX) * Vg;
		
		*RHS = 
		*(CoeffN[0 ]+INDEX)*V[0 ]+*(CoeffN[1 ]+INDEX)*V[1 ]+*(CoeffN[2 ]+INDEX)*V[2 ]+
		*(CoeffN[3 ]+INDEX)*V[3 ]+*(CoeffN[4 ]+INDEX)*V[4 ]+*(CoeffN[5 ]+INDEX)*V[5 ]+
		*(CoeffN[6 ]+INDEX)*V[6 ]+*(CoeffN[7 ]+INDEX)*V[7 ]+*(CoeffN[8 ]+INDEX)*V[8 ]+
		*(CoeffN[9 ]+INDEX)*V[9 ]+*(CoeffN[10]+INDEX)*V[10]+*(CoeffN[11]+INDEX)*V[11]+
		*(CoeffN[12]+INDEX)*V[12]+*(CoeffN[13]+INDEX)*V[13]+*(CoeffN[14]+INDEX)*V[14]+
		*(CoeffN[15]+INDEX)*V[15]+*(CoeffN[16]+INDEX)*V[16]+*(CoeffN[17]+INDEX)*V[17];
				
		AX = *LHS + *RHS;
		
		DefCor = (*(Defdia+INDEX))*(*(prevIterMassVel+INDEX))+(*(DefCoeffN[0]+INDEX))*defer[0]+(*(DefCoeffN[1]+INDEX))*defer[1]+(*(DefCoeffN[2]+INDEX))*defer[2]+(*(DefCoeffN[3]+INDEX))*defer[3]+(*(DefCoeffN[4]+INDEX))*defer[4]+(*(DefCoeffN[5]+INDEX))*defer[5];
		
		RS = delta_t*source + vol*velP - DefCor*MULTIPLE;
		residue = RS-AX;
		diacoeff = *(dia + INDEX) + delta_t*DiagSource;
		
		*(massVelDomain+INDEX)= *(massVelDomain+INDEX)+(omega*residue)/diacoeff;
		res=res+(residue*residue);
		/*
		if(INDEX==976)
		printf("RS =%e\tAX=%e\t%d\t%d\t%d\t%d\t%e\tdia=%e\n",RS,AX,INDEX,i,j,k,*(massVelDomain+INDEX),diacoeff);*/
		INDEX++;
		
	}

// 	ExtractMassVelocities(vartag*100,base,zone);
// 	for(nBc=0;nBc<numBC;nBc++)
// 	{
// 		zonebc = (zone->ZoneBc+nBc);	
// 		velocitymassBc(base,vartag,zone,zonebc);/* Ensure that the right tag is supplied for */
// 		res=res+(zonebc->NueRms); /*Ensure that the right tag is supplied for */
// 	}
	
	free(LHS);
	free(RHS);
	return(res);
}
