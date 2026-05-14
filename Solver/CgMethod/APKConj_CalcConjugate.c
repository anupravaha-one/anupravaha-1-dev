/**
The function is similar to Conj_PressureIter except for that it is used to find the potential flow solution i.e., it solves the Potential Flow equation. 
This can be called even when the user choses Opti GS as the main iteration scheme.
 */


#include "APKConj_CalcConjugate.h"

#define max(a,b) a>b?a:b

double Conj_CalcConjugate(Base *base,int vartag)
{

	int iZone,nzones,tag1=0,TCell;//,dir,zBc,numBc;
	
	int NI,NJ,indexP,iter=0;	
	
	double alpha=1.0,beta=0.0,omega=1.0,residBase=0.0,residZone,neures,epsilon=1e-6;
	
	double Num=0.0,Dnum=0.0,SRms=0.0,rsum=0.0,Prev_Dp=1.0,Curr_Dp=0.0;
	
	double *R_ini,*Res,*Vi,*Si,*Variable,*Pi,*Ti;
	
	FILE *Convcri;
	Zone *zone;
	
	nzones=base->Number_of_zones;
	
	Convcri=fopen("Settings/convergenceCriteria.dat","r");
	fscanf(Convcri,"%lf%lf",&epsilon,&epsilon);
	fclose(Convcri);
	
	*(base->BaseParameters->convergenceCriterion+1) = epsilon;
	
	EvaluateEnergysource(base);
	
	if(base->BiCgflag==2)
	{	
		for(iZone=0;iZone<nzones;iZone++)
		{	
			zone=base->CgnsZone+iZone;
			
			(zone->Residue_init)= (double *)initialize_variable(1,zone,0.0);	
			(zone->Pi)= (double *)initialize_variable(1,zone,0.0);
			(zone->Vi)= (double *)initialize_variable(1,zone,0.0);
			(zone->Si)= (double *)initialize_variable(1,zone,0.0);
			(zone->Ti)= (double *)initialize_variable(1,zone,0.0);
		}
	}	
	printf("Solving for Variable-%d:\n",vartag);
	
	if(base->BaseParameters->flag2d==1) //2D
	{
		CalcInitConjVarResidue2d(base,vartag);
		
		do
		{	
			iter++;
			residBase=0.0;
			neures=0.0;
			Num=0.0;Dnum=0.0;
			SRms=0.0;
			Curr_Dp=0.0;	
				
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				
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
			
			if(Prev_Dp==0.0||Curr_Dp==0.0||omega==0.0) 
			{
				printf("%dMethod fails %e %e %e(Beta)\n",iter,Prev_Dp,Curr_Dp,omega);
				exit(0);
			}
			
			beta = (alpha/omega)*(Curr_Dp/Prev_Dp);
			
	// 		printf("CDp is %.18e DNum is %.18e beta %.18e\n",alpha*Curr_Dp,omega*Prev_Dp,beta);
	// 		if(iter==2) exit(0);
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				
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
			
			Dnum=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
					
				CalcLHS(base,zone,1,vartag); ///Finds APi and returns Resid^ * Vi for all cells in a zone
				//printf("LHS1 Num is %e DNum is %e \n",Num,Dnum);
				
				NI=*(zone->Size+0);
				NJ=*(zone->Size+1);
					
				TCell=(NI-1)*(NJ-1);
				
				R_ini = zone->Residue_init;
				
				Res = zone->Residue;
				
				Vi = zone->Vi;
				//Dnum=0.0;	
				for(indexP=0;indexP<TCell;indexP++)
				{	
					//Num+= *(R_ini+indexP) * *(Res+indexP);
						
					Dnum+= *(R_ini+indexP) * *(Vi+indexP);
					
				}
				
				//printf("Num is %e DNum is %e \n",Num,Dnum);
			}
			
			Prev_Dp = Curr_Dp;/// since C_Dp is from the latest residual..no need to find Num again;
			
			if(Dnum==0.0) 
			{
				printf("Method fails %e (alpha)\n",Dnum);
				exit(0);
			}
			
			alpha = Curr_Dp/Dnum;
			
	// 		printf("CDp is %.18e DNum is %.18e Alpha %.18e\n",Curr_Dp,Dnum,alpha);
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
					
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
			
			residZone = sqrt(SRms/(base->TotalCellFluid));
			
			if(residZone<epsilon)
				tag1=1;
			
			Dnum=0.0;Num=0.0;
			
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
					
				CalcLHS(base,zone,2,vartag);
					
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
			
			if(Dnum==0.0) 
			{
				printf("Method fails %e (omega)\n",Dnum);
				exit(0);
			}
			
			omega = Num/Dnum;
	// 		printf("Num is %.18e DNum is %.18e omega %.18e\n",Num,Dnum,omega);//exit(0);
			SRms=0.0;
			Curr_Dp=0.0;
			for(iZone=0;iZone<nzones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
					
				NI=*(zone->Size+0);
				NJ=*(zone->Size+1);
					
				TCell=(NI-1)*(NJ-1);
					
				Si = zone->Si;
				
				Pi = zone->Pi;
				
				Ti = zone->Ti;
				
				Res = zone->Residue;
				
				R_ini = zone->Residue_init;
				
				if(vartag==6)
					Variable = zone->Flowsol->T;
				else if(vartag==990)
					Variable = zone->Flowsol->Flow_Potential;
				
				for(indexP=0;indexP<TCell;indexP++)
				{
					rsum = alpha * *(Pi+indexP) + omega * *(Si+indexP);
					
					*(Variable+indexP) = *(Variable+indexP) + rsum;
					
					SRms += rsum * rsum;
					
					*(Res+indexP) = *(Si+indexP) - omega * *(Ti+indexP);
					
					//SRms += *(Res+indexP) * *(Res+indexP);
					
				}
				
				ExtractVariables(vartag,base,zone);
			}
			
			SRms=sqrt(SRms/(base->TotalCellFluid));
			
			residBase = max(residZone,SRms);
			
			if(iter%1000 == 0)
			{
				Convcri=fopen("Settings/convergenceCriteria.dat","r");
				fscanf(Convcri,"%lf%lf",&epsilon,&epsilon);
				fclose(Convcri);
				
			}
		
		}while(residBase > epsilon);
		
		RefreshSurfaceCells(vartag,base);
	}
	else if(base->BaseParameters->flag2d==0) //3D
	{
		printf("exit\n");
		exit(0);
	}
	
	printf(" %d residzone %e srms %e\n",iter,residZone,SRms);
	
	if(base->BiCgflag==2)
	{	
		for(iZone=0;iZone<nzones;iZone++)
		{	
			zone=base->CgnsZone+iZone;
			
			free(zone->Residue_init);	
			free(zone->Pi);
			free(zone->Vi);
			free(zone->Si);
			free(zone->Ti);
		}
	}	
	
	
	return(0);
}




