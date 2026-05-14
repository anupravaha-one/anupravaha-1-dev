/**
This function solves for pressure in the current iteration after the mass velocities have been
converged. This function loops over all the zones and solves for pressure until the total residue
is less than the specified value.This function calls other functions for solving pressure
for different cases of density (constant or variable) and geometry (2D or 3D) in which the actual
Bi-CG iterations are executed. 

 */


#include "APKConj_PressureIter.h"

#define max(a,b) a>b?a:b
#define min(a,b) a>b?b:a

double Conj_PressureIter(Base *base)
{

	int iZone,nzones,i,j,k,tag1=0,tag2=0,TCell;//,*NCell;//,dir,zBc,numBc;
	
	int NI,NJ,NK,indexP,iter=0,niter;	
	
	double alpha=1.0,beta=0.0,omega=1.0,residBase=0.0,residZone,epsilon=1e-6;
	
	double Num=0.0,Dnum=0.0,SRms=0.0,rsum=0.0,Prev_Dp=1.0,Curr_Dp=0.0;
	
	double *massFlux;
	
	double *R_ini,*Res,*Vi,*Si,*Variable,*Pi,*Ti;
	
	FILE *pressConvcri;
	Zone *zone;
	
	nzones=base->Number_of_zones;
	
// 	NCell= (int *)malloc(nzones*sizeof(int));
	
	pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
	fscanf(pressConvcri,"%lf%d",&epsilon,&niter);
	fclose(pressConvcri);
	
	base->BaseParameters->PressLoopConv = epsilon;
	
	if((base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 0)||(base->BaseParameters->flag2d == 1 &&base->flagMPFV==1))
	{	
		
		for(iZone=0;iZone<nzones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			NI=*(zone->Size+0);
			NJ=*(zone->Size+1);
			NK=*(zone->Size+2);
			
			indexP=0;
				
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux=zone->Flowsol->Massflux+indexP*6;
					
				*((zone->grid)->GRHS+indexP) = *(massFlux+0)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3);
					
				indexP++;
			}
		}
			
		niter=51;
		do
		{
			for(iZone=0;iZone<nzones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					ExtractVariables(4,base,zone);
						
					residZone=UpdatePressure2d(base,zone);
						
					residBase=residBase+residZone;
				}
			}
				
			residBase=sqrt(residBase/(base->TotalCellFluid));
	
			iter++;
			
			if((residBase<epsilon)||(iter>51))
				tag1=1;
			
		}while(tag1==0);
		
		if(residBase<epsilon) 
		{
			//printf("Pressure Resid Converged\n");
			RefreshSurfaceCells(4,base);
			return(0);
		}
		
		base->BaseParameters->PressLoopConv = epsilon;
		
		if(zone->flagSolidFluid==0)
		{
			CalcInitPressureResidue2d(base);//calc.s both InitResids nd Resid^
		}
		
		iter=0;
		
		do
		{	
			iter++;
			residBase=0.0;
			Num=0.0;Dnum=0.0;
			SRms=0.0;
			Curr_Dp=0.0;	
				
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					
					TCell=(NI-1)*(NJ-1);
					
					Res = zone->Residue;
				
					R_ini = zone->Residue_init;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						Curr_Dp += *(R_ini+indexP) * *(Res+indexP);
					}
				}
			}
			
//  			if(Prev_Dp==0.0||Curr_Dp==0.0||omega==0.0) 
// 			{
				//printf("Method fails %e %e %e(Beta)\n",Prev_Dp,Curr_Dp,omega);
				//exit(0);
// 			}
			
			beta = (alpha/omega)*(Curr_Dp/Prev_Dp);
	
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					
					TCell=(NI-1)*(NJ-1);
				
					Res = zone->Residue;
				
					Vi = zone->Vi;
					
					Pi = zone->Pi;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(Pi+indexP) = *(Res+indexP) + beta * (*(Pi+indexP) - omega * *(Vi+indexP));
					}
				}
			}
			
			Dnum=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					CalcLHS(base,zone,1,4); ///Finds APi and returns Resid^ * Vi for all cells in a zone
					
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
						
					TCell=(NI-1)*(NJ-1);
					
					R_ini = zone->Residue_init;
					
					Vi = zone->Vi;
					
					for(indexP=0;indexP<TCell;indexP++)
					{	
						Dnum+= *(R_ini+indexP) * *(Vi+indexP);
						
					}
				}
			}
	// 		exit(0);
			Prev_Dp = Curr_Dp;/// since C_Dp is from the latest residual..no need to find Num again;
			
			if(Dnum==0.0) 
			{
				printf("Method fails %e (alpha)\n",Dnum);
				exit(0);
			}
			
			alpha = Curr_Dp/Dnum;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
						
					TCell=(NI-1)*(NJ-1);
						
					Si = zone->Si;
							
					Res = zone->Residue;
					
					Vi = zone->Vi;
						
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(Si+indexP)= *(Res+indexP) - alpha * *(Vi+indexP);
						
						SRms += *(Si+indexP) * *(Si+indexP);
					}
				}
			}
			
			residZone = sqrt(SRms/(base->TotalCellFluid));
			
			Dnum=0.0;Num=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					CalcLHS(base,zone,2,4);
						
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
						
					TCell=(NI-1)*(NJ-1);
					
					Ti = zone->Ti;
					
					Si = zone->Si;
						
					for(indexP=0;indexP<TCell;indexP++)
					{
						Num+= *(Ti+indexP) * *(Si+indexP);
							
						Dnum+= *(Ti+indexP) * *(Ti+indexP);
					}
				}
			}
			
			if(Dnum==0.0) 
			{
				printf("Method fails %e (omega)\n",Dnum);
				exit(0);
			}
			omega = Num/Dnum;
			
			SRms=0.0;
			
			Curr_Dp=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
						
					TCell=(NI-1)*(NJ-1);
						
					Si = zone->Si;
					
					Pi = zone->Pi;
					
					Ti = zone->Ti;
					
					Res = zone->Residue;
					
					R_ini = zone->Residue_init;
						
					Variable = zone->Flowsol->P;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						rsum = alpha * *(Pi+indexP) + omega * *(Si+indexP);
						
						*(Variable+indexP) = *(Variable+indexP) + rsum;
						
						SRms += rsum * rsum;
						
						*(Res+indexP) = *(Si+indexP) - omega * *(Ti+indexP);
						
					}
					
					ExtractVariables(4,base,zone);
				}
			}
			
			SRms=sqrt(SRms/(base->TotalCellFluid));
			
			residBase = max(residZone,SRms);
			
			if(iter%1000 == 0)
			{
				printf("CGIterCount=%d Pressure Residual= %12.8e\n",iter,residBase);
				pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
				fscanf(pressConvcri,"%lf",&epsilon);
				fclose(pressConvcri);
			}
			
		}while(residBase > epsilon);
		
		RefreshSurfaceCells(4,base);
	}
	
	else if((base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 0)||(base->BaseParameters->flag2d == 0 &&base->flagMPFV==1))
	{
		for(iZone=0;iZone<nzones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			NI=*(zone->Size+0);
			NJ=*(zone->Size+1);
			NK=*(zone->Size+2);

			indexP=0;
				
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux=zone->Flowsol->Massflux+indexP*6;
					
				*((zone->grid)->GRHS+indexP) = *(massFlux+0)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
					
				indexP++;
			}
		}
		niter=101;
		do
		{
			for(iZone=0;iZone<nzones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					residZone=UpdatePressure(base,zone);
						
					residBase=residBase+residZone;
					
					ExtractVariables(4,base,zone);
				}
			}
				
			residBase=sqrt(residBase/(base->TotalCellFluid));
	
			iter++;
			
			if((residBase<epsilon)||(iter>niter))
				tag1=1;
			
		}while(tag1==0);
		
		if(residBase<epsilon) 
		{
			//printf("Pressure Resid Converged %e\n",residBase);
			RefreshSurfaceCells(4,base);
			return(0);
		}
		
		base->BaseParameters->PressLoopConv = epsilon;
		
		if(zone->flagSolidFluid==0)
		{
			CalcInitPressureResidue(base);//calc.s both InitResids nd Resid^
		}
		
		iter=0;
		do
		{	
			iter++;
			residBase=0.0;
// 			neures=0.0;
			Num=0.0;Dnum=0.0;
			SRms=0.0;
			Curr_Dp=0.0;	
				
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
					
					TCell=(NI-1)*(NJ-1)*(NK-1);
					
					Res = zone->Residue;
				
					R_ini = zone->Residue_init;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						Curr_Dp += *(R_ini+indexP) * *(Res+indexP);
					}
				}
			}
			
			if(Prev_Dp==0.0||Curr_Dp==0.0||omega==0.0) 
			{
				printf("Method fails %e %e %e(Beta)\n",Prev_Dp,Curr_Dp,omega);
				//exit(0);
			}
			
			beta = (alpha/omega)*(Curr_Dp/Prev_Dp);
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
					
					TCell=(NI-1)*(NJ-1)*(NK-1);
				
					Res = zone->Residue;
				
					Vi = zone->Vi;
					
					Pi = zone->Pi;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(Pi+indexP) = *(Res+indexP) + beta * (*(Pi+indexP) - omega * *(Vi+indexP));
					}
				}
			}
			
			Dnum=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					CalcLHS3d(base,zone,1,4); 
					
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
						
					TCell=(NI-1)*(NJ-1)*(NK-1);
					
					R_ini = zone->Residue_init;
					
					Res = zone->Residue;
					
					Vi = zone->Vi;
					
					for(indexP=0;indexP<TCell;indexP++)
					{	
						Dnum+= *(R_ini+indexP) * *(Vi+indexP);
					}
				}
			}
	 		
			Prev_Dp = Curr_Dp;/// since C_Dp is from the latest residual..no need to find Num again;
			
			if(Dnum==0.0) 
			{
				printf("Method fails %e (alpha)\n",Dnum);
				//exit(0);
			}
			
			alpha = Curr_Dp/Dnum;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
						
					TCell=(NI-1)*(NJ-1)*(NK-1);
						
					Si = zone->Si;
							
					Res = zone->Residue;
					
					Vi = zone->Vi;
						
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(Si+indexP)= *(Res+indexP) - alpha * *(Vi+indexP);
						
						SRms += *(Si+indexP) * *(Si+indexP);
					}
				}
			}
			
			residZone = sqrt(SRms/(base->TotalCellFluid));
			
			Dnum=0.0;Num=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					CalcLHS3d(base,zone,2,4);
						
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
						
					TCell=(NI-1)*(NJ-1)*(NK-1);
					
					Ti = zone->Ti;
					
					Si = zone->Si;
						
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						Num+= *(Ti+indexP) * *(Si+indexP);
							
						Dnum+= *(Ti+indexP) * *(Ti+indexP);
					}
				}
			}
			if(Dnum==0.0) 
			{
				printf("Method fails %e (omega)\n",Dnum);
// 				exit(0);
			}
			omega = Num/Dnum;
	
			SRms=0.0;
	
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
						
					TCell=(NI-1)*(NJ-1)*(NK-1);
						
					Si = zone->Si;
					
					Pi = zone->Pi;
					
					Ti = zone->Ti;
					
					Res = zone->Residue;
					
					R_ini = zone->Residue_init;
						
					Variable = zone->Flowsol->P;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						rsum = alpha * *(Pi+indexP) + omega * *(Si+indexP);
						
						*(Variable+indexP) = *(Variable+indexP) + rsum;
						
						SRms += rsum * rsum;
						
						*(Res+indexP) = *(Si+indexP) - omega * *(Ti+indexP);
						
					}
					
					ExtractVariables(4,base,zone);
				}
			}
			
			SRms=sqrt(SRms/(base->TotalCellFluid));
			
			residBase = max(residZone,SRms);
			
			if(iter%1000 == 0)
			{
				printf("CGIterCount=%d Pressure residual= %e\n",iter,residBase);
				pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
				fscanf(pressConvcri,"%lf",&epsilon);
				fclose(pressConvcri);
			}
			
			tag1=0;tag2=0;
			
		}while(residBase > epsilon);
		
		RefreshSurfaceCells(4,base);
	}
	
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 1&&(base->flagMPFV==0))
	{
		if((base->flagTwoPhase == 1)||(base->flagTurbComb==1))
		{
				
			for(iZone=0;iZone<nzones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				
				NI=*(zone->Size+0);
				NJ=*(zone->Size+1);
				NK=*(zone->Size+2);
				
				zone->Flowsol->VarDensityCoeff=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				
				CalcExtraTermsOfVarDensityMassflux2d(base,zone);
			}
			
			do
			{
				for(iZone=0;iZone<nzones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						ExtractVariables(4,base,zone);
							
						residZone=UpdatePressure2dVarDensityMassFlux(base,zone);
							
						residBase=residBase+residZone;
					}
				}
					
				residBase=sqrt(residBase/(base->TotalCellFluid));
		
				iter++;
				
				if((residBase<epsilon)||(iter>51))
					tag1=1;
				
			}while(tag1==0);
		
			if(residBase<epsilon) 
			{
				printf("Pressure Resid Converged\n");
				
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					free(zone->Flowsol->VarDensityCoeff);
				}
				RefreshSurfaceCells(4,base);
				return(0);
			}
		
			if(zone->flagSolidFluid==0)
			{
				residBase=CalcInitPressureResidue2dVarDensityMassFlux(base);
			}
		
			iter=0;
		
			do
			{	
				iter++;
				residBase=0.0;
				Num=0.0;Dnum=0.0;
				SRms=0.0;
				Curr_Dp=0.0;	
					
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					
					if(zone->flagSolidFluid==0)
					{
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
						
						TCell=(NI-1)*(NJ-1);
						
						Res = zone->Residue;
					
						R_ini = zone->Residue_init;
						
						for(indexP=0;indexP<TCell;indexP++)
						{
							Curr_Dp += *(R_ini+indexP) * *(Res+indexP);
						}
					}
				}
			
// 			if(Prev_Dp==0.0||Curr_Dp==0.0||omega==0.0) 
// 			{
				//printf("Method fails %e %e %e(Beta)\n",Prev_Dp,Curr_Dp,omega);
				//exit(0);
// 			}
			
				beta = (alpha/omega)*(Curr_Dp/Prev_Dp);//cd/pdp;
	
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
						
						TCell=(NI-1)*(NJ-1);
					
						Res = zone->Residue;
					
						Vi = zone->Vi;
						
						Pi = zone->Pi;
						
						for(indexP=0;indexP<TCell;indexP++)
						{
							*(Pi+indexP) = *(Res+indexP) + beta * (*(Pi+indexP) - omega * *(Vi+indexP));
						}
					}
				}
			
				Dnum=0.0;
			
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						CalcLHS(base,zone,1,4); 
						
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
							
						TCell=(NI-1)*(NJ-1);
						
						R_ini = zone->Residue_init;
						
						Vi = zone->Vi;
						
						for(indexP=0;indexP<TCell;indexP++)
						{	
							Dnum+= *(R_ini+indexP) * *(Vi+indexP);
							
						}
					}
				}
				
				Prev_Dp = Curr_Dp;
			
				if(Dnum==0.0) 
				{
					printf("Method fails %e (alpha)\n",Dnum);
					exit(0);
				}
				
				alpha = Curr_Dp/Dnum;
				
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
							
						TCell=(NI-1)*(NJ-1);
							
						Si = zone->Si;
								
						Res = zone->Residue;
						
						Vi = zone->Vi;
							
						for(indexP=0;indexP<TCell;indexP++)
						{
							*(Si+indexP)= *(Res+indexP) - alpha * *(Vi+indexP);
							
							SRms += *(Si+indexP) * *(Si+indexP);
						}
					}
				}
				
				residZone = sqrt(SRms/(base->TotalCellFluid));
		
				Dnum=0.0;Num=0.0;
				
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						CalcLHS(base,zone,2,4);
							
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
							
						TCell=(NI-1)*(NJ-1);
						
						Ti = zone->Ti;
						
						Si = zone->Si;
							
						for(indexP=0;indexP<TCell;indexP++)
						{
							Num+= *(Ti+indexP) * *(Si+indexP);
								
							Dnum+= *(Ti+indexP) * *(Ti+indexP);
						}
					}
				}
				if(Dnum==0.0) 
				{
					printf("Method fails %e (omega)\n",Dnum);
					exit(0);
				}
				omega = Num/Dnum;
				
				SRms=0.0;
				
				Curr_Dp=0.0;
				
				for(iZone=0;iZone<nzones;iZone++)
				{	
					zone=base->CgnsZone+iZone;
					
					if(zone->flagSolidFluid==0)
					{
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
							
						TCell=(NI-1)*(NJ-1);
							
						Si = zone->Si;
						
						Pi = zone->Pi;
						
						Ti = zone->Ti;
						
						Res = zone->Residue;
						
						R_ini = zone->Residue_init;
							
						Variable = zone->Flowsol->P;
						
						for(indexP=0;indexP<TCell;indexP++)
						{
							rsum = alpha * *(Pi+indexP) + omega * *(Si+indexP);
							
							*(Variable+indexP) = *(Variable+indexP) + rsum;
							
							SRms += rsum * rsum;
							
							*(Res+indexP) = *(Si+indexP) - omega * *(Ti+indexP);
							
						}
						
						ExtractVariables(4,base,zone);
					}
				}
				
				SRms=sqrt(SRms/(base->TotalCellFluid));
				
				residBase = max(residZone,SRms);
				
				if(iter%1000 == 0)
				{
					printf("CGIterCount=%d Pressure Residual= %e\n",iter,residBase);
					pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
					fscanf(pressConvcri,"%lf",&epsilon);
					fclose(pressConvcri);
				}
				
			}while(residBase > epsilon);
			
			RefreshSurfaceCells(4,base);
		}
	}
	return(0);
}




