/**
This function solves for pressure in the current iteration after the mass velocities have been
converged. This function loops over all the zones and solves for pressure until the total residue
is less than the specified value.This function calls other functions for solving pressure
for different cases of density (constant or variable) and geometry (2D or 3D) in which the actual
Optimized G-S iterations are executed. 

 */



#include "APKModiPressureIter.h"


int ModiPressureIter(Base *base)
{
	FILE *pressConvcri;
	Zone *zone;
	
	int i,j,k,NI,NJ,NK,*ptrSize,indexP;
	int totalBaseCells=0,iZone,numZones,bc_no;
	
	double residBase,residZone,neures,convcri=1e-6;
	double *massFlux;
	
	
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
	fscanf(pressConvcri,"%lf",&(convcri));
	fclose(pressConvcri);

	int TCell,iter=0;
	
	double *pressureDomain;
	double Alpha,RFactor=1.0,Num=0.0,Dnum=0.0;
	
	double *Ax_O,*delx_O,*Ax,*delx,*Residue;
	
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 0)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			ptrSize = (zone->grid)->GSize;
			
			NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
			
			indexP=0;
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				//indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				massFlux=zone->Flowsol->Massflux+indexP*6;
				
				*((zone->grid)->GRHS+indexP) = *(massFlux+0)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3);
				
				indexP++;
			}
		}
		
		for(iZone=0;iZone<numZones;iZone++)
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
		
		if(residBase<convcri)
		{
// 			printf("Pressure Resid Converged\n");
			RefreshSurfaceCells(4,base);
			return(0);
		}

		if(zone->flagSolidFluid==0)
		{
			CalcInitPressureResidue2d(base);
		}
		do
		{
			residBase=0.0;
			neures=0.0;
		/**	
			for(iZone=0;iZone<numZones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				
				delx = zone->Correction;
				
				NI=*(zone->Size+0);
				NJ=*(zone->Size+1);
				
				TCell=(NI-1)*(NJ-1);
			
				for(i=0;i<TCell;i++) 
					*(delx+i)=0.0;
			
			}***/
			
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				
				if(zone->flagSolidFluid==0)
				{
					residZone=CalcCorrection(base,zone,4);
				}
				residBase=residBase+residZone;
			}
		
			Num=0.0;Dnum=0.0;
			
			for(iZone=0;iZone<numZones;iZone++)
			{	
// 				printf("Zone Num %d\n",iZone);
				zone=base->CgnsZone+iZone;
				
				if(zone->flagSolidFluid==0)
				{
					CalcLHS(base,zone,0,4);//chk the syntax
				
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					
					TCell=(NI-1)*(NJ-1);
					
					Ax = zone->AP;
						
					Ax_O = zone->AP_Prev;
					
					if(iter!=0)
					for(indexP=0;indexP<TCell;indexP++)
					{
						Num+= *(Ax_O+indexP) * *(Ax+indexP);
						
						Dnum+= *(Ax_O+indexP) * *(Ax_O+indexP);
					}
				}
			}
			
//#################################################################################################### chk below
			
			if(iter!=0)
			{
				//Alpha=0.0;
				Alpha = Num/Dnum;
				//printf("%d Alpha %e %e %e\n",iter,Alpha,Num,Dnum);
				Num=0.0;Dnum=0.0;
				
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					
					TCell=(NI-1)*(NJ-1);
					
					delx = zone->Correction;
					
					Ax = zone->AP;
					
					delx_O = zone->Correction_Prev;
					
					Ax_O = zone->AP_Prev;
					
					Residue=zone->Residue;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(delx+indexP) = *(delx+indexP) - Alpha * *(delx_O+indexP);
				
						*(Ax+indexP) = *(Ax+indexP) - Alpha * *(Ax_O+indexP);
						
						Num+= *(Residue+indexP) * *(Ax+indexP);
					
						Dnum+= *(Ax+indexP) * *(Ax+indexP);
					}
				
				}
				
				RFactor = Num/Dnum;
			}
			
			
// 			RFactor=1.0;
// 			if(iter!=0)
// 			RFactor = Num/Dnum;
// 			printf("Omeg=%12.8e\n",RFactor);
// 			exit(0);
			if(RFactor<0.8) RFactor=1.0;
// 			if(RFactor>2.0) RFactor=2.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
					
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					
					TCell=(NI-1)*(NJ-1);
					
					Residue=zone->Residue;
					
					delx = zone->Correction;
					
					Ax = zone->AP;
						
					delx_O = zone->Correction_Prev;
						
					Ax_O = zone->AP_Prev;
					
					pressureDomain=zone->Flowsol->P;//domain=zone->Flowsol->Flow_Potential;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(pressureDomain+indexP) = *(pressureDomain+indexP)+ RFactor* *(delx+indexP); 
						
						*(Residue+indexP) = *(Residue+indexP)-RFactor * *(Ax+indexP);
						
						*(Ax_O+indexP) = *(Ax+indexP);
						
						*(delx_O+indexP) = *(delx+indexP);
					}
			
					ExtractVariables(4,base,zone);
	// 				printf("residzone=%12.8e\n",residBase);
					for(bc_no=0;bc_no<zone->Number_of_BC_Patches;bc_no++)
						neures=neures+(zone->ZoneBc+bc_no)->NueRms;
				}
			}
			
			residBase=sqrt(residBase/(base->TotalCellFluid))+neures;
			
			if(iter%1000 == 0)
			{
// 				printf("residual pressure=%12.8e\tneures=%e\n",residBase,RFactor);
				pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
				fscanf(pressConvcri,"%lf",&convcri);//(base->BaseParameters->PressLoopConv));
				fclose(pressConvcri);
			}
			base->BaseParameters->PressLoopConv = convcri;
			iter++;
		}
		while(residBase > convcri);//base->BaseParameters->PressLoopConv);
		printf("It took %d iterations to converge pressure %e\n\n",iter,residBase);
		
	}
	
	else if(base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 0)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			ptrSize = (zone->grid)->GSize;
			
			NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
			
			indexP=0;
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				//indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				massFlux=zone->Flowsol->Massflux+indexP*6;
				
				*((zone->grid)->GRHS+indexP) = *(massFlux+0)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
				
				indexP++;
			}
		}
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				ExtractVariables(4,base,zone);
				
				residZone=UpdatePressure(base,zone);
					
				residBase=residBase+residZone;
			}
		}
			
		residBase=sqrt(residBase/(base->TotalCellFluid));
		
		if(residBase<convcri) 
		{
			RefreshSurfaceCells(4,base);
			return(0);
		}

		if(zone->flagSolidFluid==0)
		{
			CalcInitPressureResidue(base);
		}
		do
		{
			residBase=0.0;
			neures=0.0;
			/**
			for(iZone=0;iZone<numZones;iZone++)
			{	
				zone=base->CgnsZone+iZone;
				
				delx = zone->Correction;
				
				NI=*(zone->Size+0);
				NJ=*(zone->Size+1);
				NK=*(zone->Size+2);
				
				TCell=(NI-1)*(NJ-1)*(NK-1);
			
				for(i=0;i<TCell;i++) 
					*(delx+i)=0.0;
			
			} **/
			
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				
				if(zone->flagSolidFluid==0)
				{
					residZone=CalcCorrection3d(base,zone);
				}
				residBase=residBase+residZone;
			}
		
			Num=0.0;Dnum=0.0;
			
			for(iZone=0;iZone<numZones;iZone++)
			{	
// 				printf("Zone Num %d\n",iZone);
				zone=base->CgnsZone+iZone;
				
				if(zone->flagSolidFluid==0)
				{
					CalcLHS3d(base,zone,0,4);
				
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
				
					TCell=(NI-1)*(NJ-1)*(NK-1);
					
					Ax = zone->AP;
						
					Ax_O = zone->AP_Prev;
					
					if(iter!=0)
					for(indexP=0;indexP<TCell;indexP++)
					{
						Num+= *(Ax_O+indexP) * *(Ax+indexP);
						
						Dnum+= *(Ax_O+indexP) * *(Ax_O+indexP);
					}
				}
			}
			
//#################################################################################################### chk below
			
			if(iter!=0)
			{
				//Alpha=0.0;
				Alpha = Num/Dnum;
				//printf("%d Alpha %e %e %e\n",iter,Alpha,Num,Dnum);
				Num=0.0;Dnum=0.0;
				
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					
					if(zone->flagSolidFluid==0)
					{
						NI=*(zone->Size+0);
						NJ=*(zone->Size+1);
						NK=*(zone->Size+2);
				
						TCell=(NI-1)*(NJ-1)*(NK-1);
					
						delx = zone->Correction;
						
						Ax = zone->AP;
						
						delx_O = zone->Correction_Prev;
						
						Ax_O = zone->AP_Prev;
						
						Residue=zone->Residue;
						
						for(indexP=0;indexP<TCell;indexP++)
						{
							*(delx+indexP) = *(delx+indexP) - Alpha * *(delx_O+indexP);
					
							*(Ax+indexP) = *(Ax+indexP) - Alpha * *(Ax_O+indexP);
							
							Num+= *(Residue+indexP) * *(Ax+indexP);
						
							Dnum+= *(Ax+indexP) * *(Ax+indexP);
						}
					}
				}
				
				RFactor = Num/Dnum;
			}
			
			
// 			RFactor=1.0;
// 			if(iter!=0)
// 			RFactor = Num/Dnum;
// 			printf("Omeg=%12.8e\n",RFactor);
// 			exit(0);
			if(RFactor<0.8) RFactor=1.0;
// 			if(RFactor>2.0) RFactor=2.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
					
				if(zone->flagSolidFluid==0)
				{
					NI=*(zone->Size+0);
					NJ=*(zone->Size+1);
					NK=*(zone->Size+2);
				
					TCell=(NI-1)*(NJ-1)*(NK-1);
					
					Residue=zone->Residue;
					
					delx = zone->Correction;
					
					Ax = zone->AP;
						
					delx_O = zone->Correction_Prev;
						
					Ax_O = zone->AP_Prev;
					
					pressureDomain=zone->Flowsol->P;//domain=zone->Flowsol->Flow_Potential;
					
					for(indexP=0;indexP<TCell;indexP++)
					{
						*(pressureDomain+indexP) = *(pressureDomain+indexP)+ RFactor* *(delx+indexP); 
						
						*(Residue+indexP) = *(Residue+indexP)-RFactor * *(Ax+indexP);
						
						*(Ax_O+indexP) = *(Ax+indexP);
						
						*(delx_O+indexP) = *(delx+indexP);
					}
			
					ExtractVariables(4,base,zone);
	// 				printf("residzone=%12.8e\n",residBase);
					for(bc_no=0;bc_no<zone->Number_of_BC_Patches;bc_no++)
						neures=neures+(zone->ZoneBc+bc_no)->NueRms;
				}
			}
			
			residBase=sqrt(residBase/(base->TotalCellFluid))+neures;
			
			if(iter%1000 == 0)
			{
				//printf("residual pressure=%12.8e\tneures=%e\n",residBase,RFactor);
				pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
				fscanf(pressConvcri,"%lf",&convcri);//(base->BaseParameters->PressLoopConv));
				fclose(pressConvcri);
			}
			base->BaseParameters->PressLoopConv=convcri;
			iter++;
		}
		while(residBase > convcri);//base->BaseParameters->PressLoopConv);
		printf("It took %d iterations to converge pressure %e\n\n",iter,residBase);
		
	}
	
	
	
	
	
	
	RefreshSurfaceCells(4,base);
// 	exit(0);
	
	return 0;
}




/*

			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					ExtractVariables(4,base,zone);
					
					for(bc_no=0;bc_no<zone->Number_of_BC_Patches;bc_no++)
						neures=neures+(zone->ZoneBc+bc_no)->NueRms;
					
					if(base->NewGS==0)
						residZone=UpdatePressure(base,zone);
					else if(base->NewGS==1)
						residZone=ConjGS_Pressure2d(base,zone);
					//residZone=UpdatePressure2d(base,zone);
					
					residBase=residBase+residZone;
				}
			}

*/



