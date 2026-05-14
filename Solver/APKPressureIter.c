#include "APKPressureIter.h"

 /**
This function solves for pressure in the current iteration after the mass velocities have been
converged. This function loops over all the zones and solves for pressure until the total residue
is less than the specified value.This function calls other functions for solving pressure
for different cases of density (constant or variable) and geometry (2D or 3D) in which the actual
G-S iterations are executed. 

*/
int PressureIter(Base *base)
{
	//printf("I am into pressureIter\n");
	
	FILE *pressConvcri;
	Zone *zone;
	int pressureIterCount=0;
	int g,i,j,k,NI,NJ,NK,*ptrSize,indexP,iFace;
	double residBase,residZone;
	double *massFlux,sumMassFlux;
	int totalBaseCells=0,iZone,numZones;
//Added by Mrunalini
	double *massFluxDisperse;
	double rhoDisperse;
	double *rho_f, *varthetaDisperse_f;
	double sumMassFluxDisperse;
	int flagRhoDisperse,index;

	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	#ifdef FILE_OPEN
	{
		pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
		fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
		fclose(pressConvcri);
	}
	#endif
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		g = 0;
		ptrSize = (zone->grid+g)->GSize;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		(zone->grid+g)->GRHS=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
	}
	
	if((base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 0)||(base->BaseParameters->flag2d == 1 &&base->flagMPFV==1))
	{
		if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
		{
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
			
				ptrSize = (zone->grid)->GSize;
				
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
					
				indexP=0;
				
				flagRhoDisperse = base->BaseParameters->flagRhoDisperse;

				varthetaDisperse_f = zone->PrevFlowsol->VarthetaD_face;
				//rho_f = zone->PrevFlowsol->RhoG_face;

				massFlux=zone->Flowsol->Massflux;
				massFluxDisperse=zone->PrevFlowsol->MassfluxDisperse;

				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					
					indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
					rhoDisperse = *(zone->PrevFlowsol->RhoDisperse+(indexP*flagRhoDisperse));
					sumMassFlux=0.0;
					sumMassFluxDisperse=0.0;

					for(iFace=0;iFace<4;iFace++)
					{
						index = (indexP*6)+iFace;	
					
						sumMassFlux=sumMassFlux + (*(massFlux+index))*(*(rho_f+index))*(1-(*(varthetaDisperse_f+index)*base->BaseParameters->InitialVolumeFraction)); //(*(rho_f+index))
						sumMassFluxDisperse=sumMassFluxDisperse + (*(massFluxDisperse+index))*(*(varthetaDisperse_f+index)*base->BaseParameters->InitialVolumeFraction);	
					}

					sumMassFluxDisperse=sumMassFluxDisperse*rhoDisperse;//This is done as density of droplet is considered a constant!
				/*	if(i==(int)(NI/2) && j==(int)(NJ/2))
					{
						printf("\nsum of gas phase mass flux = %lf",sumMassFlux);
						printf("\nsum of gas droplet mass flux = %lf\n",sumMassFluxDisperse);
					}	*/
					*((zone->grid)->GRHS+indexP) = sumMassFlux + sumMassFluxDisperse;
				}
			}


		}
		else
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
				
				sumMassFlux =  *(massFlux+0)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3);

				*((zone->grid)->GRHS+indexP) = sumMassFlux;
				
				indexP++;
				}
			}
		}
		
		do
		{
			residBase=0.0;
			
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
			//residBase=(residBase/numZones);
			residBase=sqrt(residBase/totalBaseCells);

			pressureIterCount++;
			if(pressureIterCount%10000 == 0)
			{
				
				printf("GSIterCount=%d Pressure residual= %12.8e\n",pressureIterCount,residBase);
				#ifdef FILE_OPEN
				{
					pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
					fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
					fclose(pressConvcri);
				}
				#endif
			}
			if(pressureIterCount > 10000000) residBase = 0.0;		//Terminating while loop if pressure iteration is greater than 1 crore
		}
		while(residBase > base->BaseParameters->PressLoopConv);
		//printf("It took %d iterations to converge\n\n",pressureIterCount);
	}
	if((base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 0)||(base->BaseParameters->flag2d == 0 &&base->flagMPFV==1))
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			ptrSize = (zone->grid+0)->GSize;
			NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				massFlux=zone->Flowsol->Massflux+indexP*6;
				sumMassFlux=0.0;
				for(iFace=0;iFace<6;iFace++)
				{
					sumMassFlux=sumMassFlux +(*(massFlux+iFace));
				}
				*((zone->grid)->GRHS+indexP) = sumMassFlux;
			}
		}
		do
		{
			residBase=0.0;
			
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
			//residBase=(residBase/numZones);
			residBase=sqrt(residBase/totalBaseCells);

			pressureIterCount++;
			if(pressureIterCount%1000 == 0)
			{
				printf("GSIterCount=%d Pressure residual= %12.8e\n",pressureIterCount,residBase);
				#ifdef FILE_OPEN
				{
					pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
					fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
					fclose(pressConvcri);
				}
				#endif
			}
			if(pressureIterCount > 10000000) residBase = 0.0;		//Terminating while loop if pressure iteration is greater than 1 crore
		}
		while(residBase > base->BaseParameters->PressLoopConv);
	}
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 1&&(base->flagMPFV==0))
	{
		if((base->flagTwoPhase == 1)||(base->flagTurbComb==1))
		{
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				ptrSize = (zone->grid+0)->GSize;
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
				zone->Flowsol->VarDensityCoeff=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				CalcExtraTermsOfVarDensityMassflux2d(base,zone);
			}
			do
			{
				residBase=0.0;
				
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						ExtractVariables(4,base,zone);
						
						residZone=UpdatePressure2dVarDensityMassFlux(base,zone);
						
						residBase=residBase+residZone;
					}
				}
				//residBase=(residBase/numZones);
				residBase=sqrt(residBase/totalBaseCells);

				pressureIterCount++;
				if(pressureIterCount%1000 == 0)
				{
					printf("GSIterCount=%d Pressure residual= %12.8e\n",pressureIterCount,residBase);
					#ifdef FILE_OPEN
					{
						pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
						fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
						fclose(pressConvcri);
					}
					#endif
				}
				if(pressureIterCount > 10000000) residBase = 0.0;		//Terminating while loop if pressure iteration is greater than 1 crore
			}
			while(residBase > base->BaseParameters->PressLoopConv);
		}
		else
		{
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				ptrSize = (zone->grid+0)->GSize;
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
				zone->Flowsol->VarDensityCoeff=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				CalcExtraTermsOfVarDensity2d(base,zone);
			}
			do
			{
				residBase=0.0;

				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						ExtractVariables(4,base,zone);
						
						residZone=UpdatePressure2dVarDensity(base,zone);
						
						residBase=residBase+residZone;
					}
				}
				//residBase=(residBase/numZones);
				residBase=sqrt(residBase/totalBaseCells);

				pressureIterCount++;
				if(pressureIterCount%1000 == 0)
				{
					printf("GSIterCount=%d Pressure residual= %12.8e\n",pressureIterCount,residBase);
					#ifdef FILE_OPEN
					{
						pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
						fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
						fclose(pressConvcri);
					}
					#endif
				}
				if(pressureIterCount > 10000000) residBase = 0.0;		//Terminating while loop if pressure iteration is greater than 1 crore
			}
			while(residBase > base->BaseParameters->PressLoopConv);
		}
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			free(zone->Flowsol->VarDensityCoeff);
		}
	}
	if(base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 1&&(base->flagMPFV==0))
	{
		if((base->flagTwoPhase == 1)||(base->flagTurbComb==1))
		{
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				ptrSize = (zone->grid+0)->GSize;
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
				zone->Flowsol->VarDensityCoeff=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				CalcExtraTermsOfVarDensityMassflux(base,zone);
			}
			do
			{
				residBase=0.0;
				
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						ExtractVariables(4,base,zone);
						
						residZone=UpdatePressureVarDensityMassFlux(base,zone);
						
						residBase=residBase+residZone;
					}
				}
				//residBase=(residBase/numZones);
				residBase=sqrt(residBase/totalBaseCells);

				pressureIterCount++;
				if(pressureIterCount%1000 == 0)
				{
					printf("GSIterCount=%d Pressure residual= %12.8e\n",pressureIterCount,residBase);
					#ifdef FILE_OPEN
					{
						pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
						fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
						fclose(pressConvcri);
					}
					#endif
				}
				if(pressureIterCount > 10000000) residBase = 0.0;		//Terminating while loop if pressure iteration is greater than 1 crore
				
			
			}
			while(residBase > base->BaseParameters->PressLoopConv);
		}
                else
		{
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				ptrSize = (zone->grid+0)->GSize;
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
				zone->Flowsol->VarDensityCoeff=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				///2d fuctions is called incase of 3d problem.  Please check this before uploading.
				///very important as it can halt other calculations.
				//CalcExtraTermsOfVarDensity2d(base,zone);
				CalcExtraTermsOfVarDensity(base,zone);
			}
			do
			{
				residBase=0.0;
				
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						ExtractVariables(4,base,zone);
						
						residZone=UpdatePressureVarDensity(base,zone);
						
						residBase=residBase+residZone;
					}
				}
				//residBase=(residBase/numZones);
				residBase=sqrt(residBase/totalBaseCells);

				pressureIterCount++;
				if(pressureIterCount%1000 == 0)
				{
					printf("GSIterCount=%d Pressure residual= %12.8e\n",pressureIterCount,residBase);
					#ifdef FILE_OPEN
					{
						pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
						fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
						fclose(pressConvcri);
					}
					#endif
				}
				if(pressureIterCount > 10000000) residBase = 0.0;		//Terminating while loop if pressure iteration is greater than 1 crore
			}
			while(residBase > base->BaseParameters->PressLoopConv);
		}
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			free(zone->Flowsol->VarDensityCoeff);
		}
	}
	RefreshSurfaceCells(4,base);
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		free((zone->grid)->GRHS);
	}
	return 0;
}
