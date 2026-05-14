#include "APKfindresidualQVK2d.h"

/**
This function solves for the mass velocities only using Gauss-Seidal point method and returns the
residual of the mass velocities for that particular iteration. It is called only for 2D case in which QUICK
conevction scheme is used and having either constant or variable diffusion coefficient.
*/

double findresidualQVK2d(Base *base,Zone *zone,int vartag)
{
	int NI,NJ,NK,index[3],INDEX,numBC,k,j,i,*ptrSize,g,cntDir;
	double vol,velP,AX,RS,residue,res=0.0,DiagSource,source;
	double delta_t,omega=0.6,*slabVelocity[6],diacoeff;
	double Vg,V[18],*LHS,*RHS,DefCor,*CoeffN[18],*dia,*DefCoeffN[12],*Defdia;
	//double Vww,Vee,Vss,Vnn;
	double *massVelDomain,*prevTstepVel,*prevIterMassVel,defer[12];
	
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
		
		
		//end of multiphase flow 
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
	}
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	int J_Increment;
	
	J_Increment = NI-1;
	
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
	
//	MakeSlabPtrArray(zone,slabVelocity,vartag*100);/*Called for making slabs of u */
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
			if(i!=1)
				defer[6] = *(prevIterMassVel+INDEX-2);
			else
				defer[6] = *(slabVelocity[1]+(NJ-1)*k+j);
		}
		else
		{
			V[0] =*(slabVelocity[1]+(NJ-1)*k+j);
			defer[0] = V[0];
			defer[6] = 0.0;
		}
		
		if(i!=(NI-2))
		{
			V[1]=*(massVelDomain+INDEX+1);
			defer[1]=*(prevIterMassVel+INDEX+1);
			
			if(i!=NI-3)
				defer[7] = *(prevIterMassVel+INDEX+2);
			else
				defer[7] = *(slabVelocity[0]+(NJ-1)*k+j);
		}
		else
		{
			V[1]= *(slabVelocity[0]+(NJ-1)*k+j);
			defer[1]= V[1];
			defer[7]= 0.0;
		}
	
		if(j!=0)
		{
			V[2]=*(massVelDomain+INDEX-J_Increment);
			defer[2] = *(prevIterMassVel+INDEX-J_Increment);
			
			if(j!=1)
				defer[8] = *(prevIterMassVel+INDEX-J_Increment-J_Increment);
			else
				defer[8] = *(slabVelocity[3]+(NI-1)*k+i);
		}
		else
		{
			V[2]=*(slabVelocity[3]+(NI-1)*k+i);
			defer[2]=V[2];
			defer[8] = 0.0;
		}
	
		if(j!=(NJ-2))
		{
			V[3]=*(massVelDomain+INDEX+J_Increment);
			defer[3] =*(prevIterMassVel+INDEX+J_Increment);
			
			if(j!=NJ-3)
				defer[9] = *(prevIterMassVel+INDEX+J_Increment+J_Increment);
			else
				defer[9] = *(slabVelocity[2]+(NI-1)*k+i);
		}
		else
		{
			V[3]=*(slabVelocity[2]+(NI-1)*k+i);
			defer[3] = V[3];
			defer[9] = 0.0;
		}
	
	
		if(!(i==0 && j==(NJ-2)))
		{
			if(i==0)
				V[6]=*(slabVelocity[1]+(j+1));//from western slab
			else
			{
				if(j==(NJ-2))
					V[6]=*(slabVelocity[2]+(i-1));//from northern slab
				else
					V[6]=*(massVelDomain+INDEX+(NI-1)-1);
			}
		}
		else
			V[6]=0.0;
		
		if(!(i==0 && j==0))
		{
			if(i==0)
				V[7]=*(slabVelocity[1]+(j-1));//from west slab
			else
			{
				if(j==0)
					V[7]=*(slabVelocity[3]+(i-1));//from southern slab
				else
					V[7]=*(massVelDomain+INDEX-(NI-1)-1);
					//V[7]=*(massVelDomain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i-1));
			}
		}
		else
			V[7]=0.0;
		
		if(!(i==(NI-2) && j==(NJ-2)))
		{
			if(i==(NI-2))
				V[10]=*(slabVelocity[0]+(j+1));
			else
			{
				if(j==(NJ-2))
					V[10]=*(slabVelocity[2]+(i+1));
				else
					V[10]=*(massVelDomain+INDEX+(NI-1)+1);
			}
		}
		else
			V[10]=0.0;
		
		if(!(i==(NI-2) && j==0))
		{
			if(i==(NI-2))
				V[11]=*(slabVelocity[0]+(j-1));
			else
			{	
				if(j==0)
					V[11]=*(slabVelocity[3]+(i+1));
				else
					V[11]=*(massVelDomain+INDEX-(NI-1)+1);
					//V[11]=*(massVelDomain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i+1));
			}
		}
		else
			V[11]=0.0;
		
		index[0]=k,index[1]=j,index[2]=i;
	
		vol = *(zone->ZGeom->Volume+INDEX);
		velP=*(prevTstepVel+INDEX);
	
		Vg = *(massVelDomain+INDEX);
		
		source = *(zone->Flowsol->Source_term+INDEX);
		DiagSource = *(zone->Flowsol->Diag_Source+INDEX);
		source = source - DiagSource*Vg;
		
		*LHS = *(dia+INDEX) * Vg;
		*RHS = *(CoeffN[0]+INDEX) *V[0] + *(CoeffN[1]+INDEX) *V[1] + *(CoeffN[2]+INDEX) *V[2] +*(CoeffN[3]+INDEX) *V[3] + *(CoeffN[6]+INDEX) *V[6] + *(CoeffN[7]+INDEX) *V[7] + *(CoeffN[10]+INDEX) *V[10]+*(CoeffN[11]+INDEX) *V[11];

		AX = *LHS + *RHS;
		
		DefCor = 	(*(Defdia+INDEX))*(*(prevIterMassVel+INDEX)) +
				(*(DefCoeffN[0]+INDEX))*defer[0]+ (*(DefCoeffN[1]+INDEX))*defer[1]+ (*(DefCoeffN[2]+INDEX))*defer[2]+ (*(DefCoeffN[3]+INDEX))*defer[3]+
				(*(DefCoeffN[6]+INDEX))*defer[6]+ (*(DefCoeffN[7]+INDEX))*defer[7]+
				(*(DefCoeffN[8]+INDEX))*defer[8]+ (*(DefCoeffN[9]+INDEX))*defer[9];
		
		RS = delta_t*source + vol*velP - DefCor*MULTIPLE;
		
		residue = RS-AX;
		diacoeff = *(dia + INDEX) + delta_t*DiagSource;
		
		*(massVelDomain+INDEX)= *(massVelDomain+INDEX)+(omega*residue)/diacoeff;
		res=res+(residue*residue);

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
