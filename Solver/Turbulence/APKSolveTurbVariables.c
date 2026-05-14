#include "APKSolveTurbVariables.h"
/**
This function solves the governing equations for turbulence variables with a dual-time step approach
*/

int SolveTurbVariables(Root CgnsRoot,double variableResidual[11], int CounterTstep)
{
	int countdual,numDualTimeIter;
	Base *base;
	//double time,dt;
	//time = CgnsRoot.CgnsBase->BaseParameters->Time;
	//dt = CgnsRoot.CgnsBase->BaseParameters->TimeStep;
	//dt = CgnsRoot.CgnsBase->BaseParameters->Turb->DualTimeStep;
	
	base = CgnsRoot.CgnsBase;
	
	if(base->flagTurbModel == 8)	//Spalart-Allmaras one-equation turbulence model
	{
		CalcAbsVorticity(base);
		
		///Calculating the effective diffusion coefficient for the kinematic eddy viscosity (nu_t) equation.
		CalcNuDiffCoef(base);
		CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,600);
		CDcoeffCalculation(base,600);
		
		CalcSource_opti(CgnsRoot.CgnsBase,600);
		///Solving equation for kinematic eddy visocity (nu_t)
		variableResidual[1] = CalcVariable(base,600);
	}
	else if(base->flagTurb == 1 && base->flagTurbModel != 9)
	{
		if(base->flagTurbModel==4)
		{
			CalcKWallTerm(base);
			CalcEpsWallTerm(base);
		}
				
		if(base->flagTurbModel == 3)// || base->flagTurbModel == 4 || base->flagTurbModel == 5)
		{
			CalcC1epsi(base);
		}
	
		if(base->flagTurbModel == 6)
		{
			CalcAbsVorticity(base);
		}
	
		CalcTurbulentProduction(base);
		
		if(CounterTstep == 1 && base->BaseParameters->Turb->flagAutoInit == 1)
		{
			printf("*** Initializing TURBULENT KINETIC ENERGY (k) & DISSIPATION RATE (epsilon) using equillibrium assumption ***\n");
			initializeTurbVariables(base);
		}
		
		CalcMuEff(base,600);
		CalcMuEff(base,601);
// 		CDcoeffCalculation(base,600);
// 		CDcoeffCalculation(base,601);
		
		//What this function does??? <<< Narsimha 21st May, 2008
		Conv_Coeff_Calculation(base,0,600);///This function needs to be called again as the delta_t different for K and Epsilon
		
		CalcGamma(base);
		
		if(base->BaseParameters->Turb->TimeScaleMinLimitCount != 0 || base->BaseParameters->Turb->TimeScaleMaxLimitCount != 0)
			printf(">>> TimeScaleMinLimitCount = %d TimeScaleMaxLimitCount = %d <<<\n", base->BaseParameters->Turb->TimeScaleMinLimitCount,base->BaseParameters->Turb->TimeScaleMaxLimitCount);
		
		numDualTimeIter = base->BaseParameters->Turb->numDualTimeIter;
		countdual=0;
		do
		{
			countdual++;
			FillNewsolToPrevsolKepsilon(base);
			
			///Solving equation for turbulent kinetic energy (k)
			if(base->CgnsZone->Flowsol->Scal[0]!=NULL)
			{
				if(base->flagTurbModel==6)
				{
					CalcMuEff(base,600);
					//CDcoeffCalculation(base,600);//Check whether this is required here? <<< Narsimha.
				}
				CallCalcCoefficientsForVariables(base,0,600);
				CDcoeffCalculation(base,600);
				
				CalcSource_opti(base,600);
				variableResidual[1] = CalcVariable(base,600);
			}
			
			LimitVariable(base,600);
			
			///Calculating friction velocity which is used in the boundary condition for epsilon
			if(!(base->flagTurbModel==4 || base->flagTurbModel==5))
			{
				CalcUtau(base);
			}
	
			///Solving equation for (specific) dissipation rate of turbulent kinetic energy (epsilon or omega)
			if(base->CgnsZone->Flowsol->Scal[1]!=NULL)
			{
				if(base->flagTurbModel==4 || base->flagTurbModel==5)
				{
					CalcC1epsi(base);
				}
						
// 				if(base->flagTurbModel==7)
// 				{
// 					TimeDepend(base->BaseParameters->Time,601,base);
// 					RefreshSurfaceCells(601,base);
// 				}
				
				if(base->flagTurbModel==6) 
				{
					CalcMuEff(base,601);	//<-- Should be called once in each dual time loop
					//CDcoeffCalculation(base,601);
					CrossDiffusion(base);
				}
				
				CallCalcCoefficientsForVariables(base,0,601);
				CDcoeffCalculation(base,601);
				
				CalcSource_opti(base,601);
				variableResidual[2] = CalcVariable(base,601);
				
				if(base->flagTurbModel==6)
				{
					freeCrossDiff(base);
				}
			}
			
			LimitVariable(base,601);
		}while(countdual < numDualTimeIter);//10);
		
// 		LimitVariable(base,601);
// 		LimitVariable(base,600);
		if(base->BaseParameters->Turb->TProductionMaxLimitCount != 0)//This should be done only for k-epsilon model
			printf("######### Max. Turbulent Production is limited in %d Cell(s) ########\n", base->BaseParameters->Turb->TProductionMaxLimitCount);
		
		freeGamma(base);
		
		if(base->flagMPFV==0)
		Conv_Coeff_Calculation(base,0,1);
		else 
		Conv_Coeff_Calculation(base,0,31);
		//printf("Out of dual time stepping loop\n");
	}// End of else loop ie Turbulent case

	return 0;
}
