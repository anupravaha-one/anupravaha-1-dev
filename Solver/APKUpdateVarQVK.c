#include "APKUpdateVarQVK.h"

/**
This is a general function unlike "findresidual" (which only solves for mass velocities).
This function solves for any variable using Gauss-Seidal point method and returns the
residual of that variable for that particular iteration. It is called only for 3D case in which QUICK
conevction scheme is used and having either constant or variable diffusion coefficient.
*/

double UpdateVarQVK(Base *base,Zone *zone,int vartag)
{
	int NI,NJ,NK,index[3],INDEX,k,j,i,*ptrSize,g,cntDir;
	double vol,velP,AX,RS,diacoeff,residue,res=0.0,omega=0.6;
	double delta_t,*slabVariable[6],varCoefficientPrev,source,DiagSource;
	double Vg,V[18],*LHS,*RHS,DefCor,*CoeffN[18],*dia,*DefCoeffN[12],*Defdia;
	double *prevIterVar,*variable,*prevTstepVar,defer[12];
	
	omega = *(base->BaseParameters->OMEGA);
	delta_t = base->BaseParameters->TimeStep;
	
	ptrSize= zone->Size;
	
	g = 0;
	
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
	
	switch(vartag)
	{
		case 1:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabU[cntDir];
			}
			variable = zone->Flowsol->u;
			prevTstepVar = zone->PrevFlowsol->u;
		break;
		case 2:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabV[cntDir];
			}
			variable = zone->Flowsol->v;
			prevTstepVar = zone->PrevFlowsol->v;
		break;
		case 3:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabW[cntDir];
			}
			variable = zone->Flowsol->w;
			prevTstepVar = zone->PrevFlowsol->w;
		break;
		case 6:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabT[cntDir];
			}
			variable = zone->Flowsol->T;
			prevTstepVar = zone->PrevFlowsol->T;
		break;
		case 7:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Enthalpy[cntDir];
			}
			variable = zone->Flowsol->enthalpy;
			prevTstepVar = zone->PrevFlowsol->enthalpy;
		break;
		case 9:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_RhoIon[cntDir];
			}
			variable = zone->Flowsol->RhoIon;
			prevTstepVar = zone->PrevFlowsol->RhoIon;
			//printf("Entered for Ion Density\n");
			omega=0.8;
		break;
		
		case 14:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabMixFrac[cntDir];
			}
			variable = zone->Flowsol->MixFrac;
			prevTstepVar = zone->PrevFlowsol->MixFrac;
			//printf("Entered for Ion Density\n");
			//omega=0.8;
			break;
			
		case 15:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabFuelFrac[cntDir];
			}
			variable = zone->Flowsol->FuelMassFrac;
			prevTstepVar = zone->PrevFlowsol->FuelMassFrac;
			//printf("Entered for Ion Density\n");
			//omega=0.8;
		break;
		case 600:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[0][cntDir];
			}
			variable = zone->Flowsol->Scal[0];
			prevTstepVar = zone->PrevFlowsol->Scal[0];
			
			omega = *(base->BaseParameters->OMEGA+4);
			delta_t = base->BaseParameters->Turb->DualTimeStep;//delta_t =delta_t/10.0;
			break;
		case 601:
			
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[1][cntDir];
			}
			variable = zone->Flowsol->Scal[1];
			prevTstepVar = zone->PrevFlowsol->Scal[1];
			
			omega = *(base->BaseParameters->OMEGA+5);
			delta_t = base->BaseParameters->Turb->DualTimeStep;
			break;
		case 602:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[2][cntDir];
			}
			variable = zone->Flowsol->Scal[2];
			prevTstepVar = zone->PrevFlowsol->Scal[2];
		break;
		case 603:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[3][cntDir];
			}
			variable = zone->Flowsol->Scal[3];
			prevTstepVar = zone->PrevFlowsol->Scal[3];
		break;
		case 604:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[4][cntDir];
			}
			variable = zone->Flowsol->Scal[4];
			prevTstepVar = zone->PrevFlowsol->Scal[4];
		break;
		case 605:
			
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[5][cntDir];
			}
			variable = zone->Flowsol->Scal[5];
			prevTstepVar = zone->PrevFlowsol->Scal[5];
		break;
		case 606:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[6][cntDir];
			}
			variable = zone->Flowsol->Scal[6];
			prevTstepVar = zone->PrevFlowsol->Scal[6];
		break;
		case 607:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[7][cntDir];
			}
			variable = zone->Flowsol->Scal[7];
			prevTstepVar = zone->PrevFlowsol->Scal[7];
		break;
		case 608:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[8][cntDir];
			}
			variable = zone->Flowsol->Scal[8];
			prevTstepVar = zone->PrevFlowsol->Scal[8];
		break;
		case 609:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabScalar[9][cntDir];
			}
			variable = zone->Flowsol->Scal[9];
			prevTstepVar = zone->PrevFlowsol->Scal[9];
		break;
		
		//added for multiphase flow problem.
		case 31:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVx[cntDir];
			}
			variable = zone->Flowsol->Momentumx;
			prevTstepVar = zone->PrevFlowsol->Momentumx;
			break;
		case 32:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVy[cntDir];
			}
			variable = zone->Flowsol->Momentumy;
			prevTstepVar = zone->PrevFlowsol->Momentumy;
			break;
		case 33:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVz[cntDir];
			}
			variable = zone->Flowsol->Momentumz;
			prevTstepVar = zone->PrevFlowsol->Momentumz;
			break;
		case 41:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVxPhase[0][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasex[0];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasex[0];
			
			break;
		case 44:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVxPhase[1][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasex[1];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasex[1];
			
			break;
		case 47:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVxPhase[2][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasex[2];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasex[2];
			
			break;
		case 50:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVxPhase[3][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasex[3];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasex[3];
			
			break;
		case 53:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVxPhase[4][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasex[4];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasex[4];
			
			break;
		case 43:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVzPhase[0][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasez[0];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasez[0];
			
			break;
		case 46:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVzPhase[1][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasez[1];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasez[1];
			
			break;
		case 49:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVzPhase[2][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasez[2];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasez[2];
			
			break;
		case 52:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVzPhase[3][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasez[3];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasez[3];
			
			break;
		case 55:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVzPhase[4][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasez[4];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasez[4];
			
			break;
		case 42:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVyPhase[0][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasey[0];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasey[0];
			
			break;
		case 45:	
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVyPhase[1][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasey[1];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasey[1];
			
			break;
		case 48:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVyPhase[2][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasey[2];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasey[2];
			
			break;
		case 51:	
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVyPhase[3][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasey[3];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasey[3];
			
			break;
		case 54:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVyPhase[4][cntDir];
			}
			variable = zone->Flowsol->Momentumofphasey[4];
			prevTstepVar = zone->PrevFlowsol->Momentumofphasey[4];
			
			break;
		case 34://for volume fraction
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVartheta[0][cntDir];
			}
			variable = zone->Flowsol->VarthetaPhase[0];
			prevTstepVar = zone->PrevFlowsol->VarthetaPhase[0];
			break;
		case 35:
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVartheta[1][cntDir];
			}
			variable = zone->Flowsol->VarthetaPhase[1];
			prevTstepVar = zone->PrevFlowsol->VarthetaPhase[1];
			break;
		case 36:
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVartheta[2][cntDir];
			}
			variable = zone->Flowsol->VarthetaPhase[2];
			prevTstepVar = zone->PrevFlowsol->VarthetaPhase[2];
			break;
		case 37:
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVartheta[3][cntDir];
			}
			variable = zone->Flowsol->VarthetaPhase[3];
			prevTstepVar = zone->PrevFlowsol->VarthetaPhase[3];
			break;
		case 38:
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVartheta[4][cntDir];
			}
			variable = zone->Flowsol->VarthetaPhase[4];
			prevTstepVar = zone->PrevFlowsol->VarthetaPhase[4];
			break;
		//added by adhi+shash for disperse phase flow//changed by mrunalini
		case 62:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVxDisperse[cntDir];
			}
			
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			variable = zone->Flowsol->U_Disperse;
			prevTstepVar = zone->PrevFlowsol->U_Disperse;
			omega=1.0;
			break;
		case 63:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVyDisperse[cntDir];
			}
			variable = zone->Flowsol->V_Disperse;
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			prevTstepVar = zone->PrevFlowsol->V_Disperse;
			omega=1.0;
			break;
		case 64:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVzDisperse[cntDir];
			}
			variable = zone->Flowsol->W_Disperse;
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			prevTstepVar = zone->PrevFlowsol->W_Disperse;
			omega=1.0;
			break;
		case 65:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVarthetaDisperse[cntDir];
			}
			variable = zone->Flowsol->VarthetaDisperse;
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			prevTstepVar = zone->PrevFlowsol->VarthetaDisperse;
			omega=1.0;
			break;
	/*	case 67:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabNumDensity[cntDir];
			}
			variable = zone->Flowsol->NumDensity;
			prevTstepVar = zone->PrevFlowsol->NumDensity;
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			omega=1.0;
			break;
*/
		case 68:
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
			   slabVariable[cntDir] = zone->Flowsol->SlabDropFuelMassFrac[cntDir];
			}
			variable = zone->Flowsol->DropFuelMassFrac;
			prevTstepVar = zone->PrevFlowsol->DropFuelMassFrac;
			break;
		case 69:
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			for(cntDir=0;cntDir<6;cntDir++)
			{
			   slabVariable[cntDir] = zone->Flowsol->SlabTdrop[cntDir];
			}
			variable = zone->Flowsol->Tdrop;
			prevTstepVar = zone->PrevFlowsol->Tdrop;
			break;
		case 710:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabDensity[cntDir];
			}
			variable = zone->Flowsol->Rho;
			prevTstepVar = zone->PrevFlowsol->Rho;
			Defdia = zone->DeferredCorrCoeffDummy_P;
			for(i=0;i<12;i++)
			{
				DefCoeffN[i] = zone->DeferredCorrCoeffDummy_P;
			}
			
			omega=1.0;
			break;
		case 56:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Total_MfieldX[cntDir];
			}
			variable = zone->Flowsol->TotalMfieldX;
			prevTstepVar = zone->PrevFlowsol->TotalMfieldX;
			break;

		case 57:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Total_MfieldY[cntDir];
			}
			variable = zone->Flowsol->TotalMfieldY;
			prevTstepVar = zone->PrevFlowsol->TotalMfieldY;
			break;
		case 58:
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->Slab_Total_MfieldZ[cntDir];
			}
			variable = zone->Flowsol->TotalMfieldZ;
			prevTstepVar = zone->PrevFlowsol->TotalMfieldZ;
			
			break;
		default:
			printf("This vartag %d is not found in UpdateVarQVK.c\n\n",vartag);
			exit(12);
			break;
		//end
	}
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	int J_Increment,K_Increment;
	
	K_Increment = (NI-1)*(NJ-1);
	J_Increment = (NI-1);
	
	prevIterVar = zone->Flowsol->prevIterVariable;
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(prevIterVar+INDEX) = *(variable+INDEX);
	}
	
	double MULTIPLE=1.0;
	
	if(base->BaseParameters->TimeStepingSchemeTag == 1)
		MULTIPLE = 0.5;
	
	LHS=(double *)malloc(sizeof(double));
	RHS=(double *)malloc(sizeof(double));
	
//	MakeSlabPtrArray(zone,slabVariable,vartag);//Called for making slabs of u
	
//---------------------------------------------------------------------------------------------------------(1)------------------
	INDEX=0;
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		if(i!=0)
		{
			V[0]=*(variable+INDEX-1);
			defer[0]=*(prevIterVar+INDEX-1);
			if(i!=1)
				defer[6] = *(prevIterVar+INDEX-2);
			else
				defer[6] = *(slabVariable[1]+(NJ-1)*k+j);
		}
		else
		{
			V[0] =*(slabVariable[1]+(NJ-1)*k+j);
			defer[0] = V[0];
			defer[6] = 0.0;
		}
		if(i!=(NI-2))
		{
			V[1]=*(variable+INDEX+1);
			defer[1]=*(prevIterVar+INDEX+1);
			
			if(i!=NI-3)
				defer[7] = *(prevIterVar+INDEX+2);
			else
				defer[7] = *(slabVariable[0]+(NJ-1)*k+j);
		}
		else
		{
			V[1]= *(slabVariable[0]+(NJ-1)*k+j);
			defer[1]= V[1];
			defer[7]= 0.0;
		}
	
		if(j!=0)
		{
			V[2]=*(variable+INDEX-J_Increment);
			defer[2] = *(prevIterVar+INDEX-J_Increment);
			
			if(j!=1)
				defer[8] = *(prevIterVar+INDEX-J_Increment-J_Increment);
			else
				defer[8] = *(slabVariable[3]+(NI-1)*k+i);
		}
		else
		{
			V[2]=*(slabVariable[3]+(NI-1)*k+i);
			defer[2]=V[2];
			defer[8] = 0.0;
		}
	
		if(j!=(NJ-2))
		{
			V[3]=*(variable+INDEX+J_Increment);
			defer[3] =*(prevIterVar+INDEX+J_Increment);
			
			if(j!=NJ-3)
				defer[9] = *(prevIterVar+INDEX+J_Increment+J_Increment);
			else
				defer[9] = *(slabVariable[2]+(NI-1)*k+i);
		}
		else
		{
			V[3]=*(slabVariable[2]+(NI-1)*k+i);
			defer[3] = V[3];
			defer[9] = 0.0;
		}
	
		if(k!=0)
		{
			V[4]=*(variable+INDEX-K_Increment);
			defer[4] = *(prevIterVar+INDEX-K_Increment);
			
			if(k!=1)
				defer[10] = *(prevIterVar+INDEX-K_Increment-K_Increment);
			else
				defer[10] = *(slabVariable[5]+(NI-1)*j+i);
		}
		else
		{
			V[4]=*(slabVariable[5]+(NI-1)*j+i);
			defer[4] = V[4];
			defer[10]= 0.0;
		}
	
		if(k!=(NK-2))
		{
			V[5]=*(variable+INDEX+K_Increment);
			defer[5]=*(prevIterVar+INDEX+K_Increment);
			
			if(k!=NK-3)
				defer[11] = *(prevIterVar+INDEX+K_Increment+K_Increment);
			else
				defer[11] = *(slabVariable[4]+(NI-1)*j+i);
		}
		else
		{
			V[5]=*(slabVariable[4]+(NI-1)*j+i);
			defer[5] = V[5];
			defer[11]= 0.0;
		}
		
			
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
		
		index[0]=k,index[1]=j,index[2]=i;
	
		vol = *(zone->ZGeom->Volume+INDEX);
		velP=*(prevTstepVar+INDEX);
	
		Vg = *(variable+INDEX);
		
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
		
		DefCor = (*(Defdia+INDEX))*(*(prevIterVar+INDEX))	+
		
		(*(DefCoeffN[0]+INDEX))*defer[0]+(*(DefCoeffN[1] +INDEX))*defer[1] +(*(DefCoeffN[2] +INDEX))*defer[2] +
		(*(DefCoeffN[3]+INDEX))*defer[3]+(*(DefCoeffN[4] +INDEX))*defer[4] +(*(DefCoeffN[5] +INDEX))*defer[5] +
		(*(DefCoeffN[6]+INDEX))*defer[6]+(*(DefCoeffN[7] +INDEX))*defer[7] +(*(DefCoeffN[8] +INDEX))*defer[8] +
		(*(DefCoeffN[9]+INDEX))*defer[9]+(*(DefCoeffN[10]+INDEX))*defer[10]+(*(DefCoeffN[11]+INDEX))*defer[11];
		
		
	
		varCoefficientPrev = CalcvarCoefficientPrev(base,zone,vartag,INDEX);
	
		RS = delta_t*source + varCoefficientPrev*vol*velP - DefCor*MULTIPLE*varCoefficientPrev;
	
		
		
		
		residue = RS-AX;
		diacoeff = *(dia + INDEX) + delta_t*DiagSource;
		
		*(variable+INDEX)= *(variable+INDEX)+(omega*residue)/diacoeff;
		res=res+(residue*residue);
		
		if(vartag==9 && ((*(variable+INDEX)) < 0.0))
			*(variable+INDEX) = 0.0;
		

		INDEX++;
	}

//------------------------------------------------------------------------------------------------------------------------
	
	free(LHS);
	free(RHS);
	return(res);
}



