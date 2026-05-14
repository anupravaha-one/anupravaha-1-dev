#include "APKsolverCoupled2d.h"
/**
The function solverCoupled2d() is the solver for two dimensional geometries which solves all the variables
with a coupled approach. This function is a control function which calls different other functions in an order
as given in the coupled algorithm. See Project Report No.1, December 2005, Section 3.10.0.4 Page 49.
Brief about Algo:
This function calls the massvelocityoneIter() to calculate massvelocities and converges them simultaneously.
Calculates vilumeflux from the calculated value of massvelocities.
Finds out Pressure using the function PressureIter().
Repeats the above steps till Volume Flux converges.
Computes velocities by using the function CalcVariableOneIter() and converges them simultaneously.
Computes all other solvable Scalars like Temperature etc and converges simultaneusly.
Computes variables coefficients like Density, Viscosity etc. if they are userdefined functions,
using the function EvaluateVarCoefficients().
Saves the current timestep data if the user chooses so, or saves after a desired no of time steps as choosen by the user
using the function SaveTimeStepData()
Returns to massvelocity calculation till the total residual of all the solvable variables are less than their
individual convergence criteria.
If the user chooses to stop the computation at any point(by changing the value to 0 from 1 in the file "iwrite.dat" in the
directory "Settings"), then the solution at that time level is written using the function WriteCgnsFile().
*/
int solverCoupled2d(Root CgnsRoot)
{
	int cntVar,timetag,vartag,CounterTstep=0;
	int total,index_var,iScal;
	double time=0.0,dt;
	FILE *save,*odt,*Q,*fp;
	char filename2[] = "SOLUTIONS/restart.cgns";
	int timeStepToken,insideTimeStepConverToken,variableIterTag,velocityIterTag,cntVariable,massVelocityIterTag,stop;
	double massVelocityResidual[3],velocityResidual[3],variableResidual[11];
	double ErrorTimeStepCoupled[16],ErrorTimeStep[16];
	double residMassFlux;
	
	Q=fopen("SOLUTIONS/MassConservation.dat","w");
	dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
	printf("\nInside SolverCoupled2d\n");
	
	timetag = CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag;
	CgnsRoot.CgnsBase->BaseParameters->Time=0.0;
	total = CgnsRoot.CgnsBase->total_var;
	
	CDcoeffMemAlloc(CgnsRoot.CgnsBase);
	
	DiffusionCoeffCalculation(CgnsRoot.CgnsBase);
	
	DiffusionCoeffCalculationForVariables(CgnsRoot.CgnsBase);
	
	do/*time stepping loop*/
	{
		if (timetag==1 && time==0) /*Even the scheme is CN, the first timestep is done by BTCS scheme*/
			CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=0;
		else if (timetag==1 && time>0)
			CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=1;
		
			
		fp=fopen("SOLUTIONS/No_of_TimeSteps","w");
		fprintf(fp,"No of Time Steps: %d\nTime: %lf", CounterTstep,time);
		fclose(fp);
		
		time = time+dt;
		CounterTstep=CounterTstep+1;
		
		save=fopen("Settings/iwrite.dat","r");
		fscanf(save,"%d",&stop);
		fclose(save);
		
		if(stop==0)
		{
			printf("Writing the solution data ...\n");
			WriteCgnsFile(filename2,CgnsRoot);
                        WriteCgnsFile_vertex("SOLUTIONS/vertex.cgns",CgnsRoot);
			printsolver(CgnsRoot.CgnsBase);
			dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
			odt=fopen("Settings/old_dt.dat","w");
			fprintf(odt,"%12.10lf %12.10lf",dt,CgnsRoot.CgnsBase->BaseParameters->Time);
			fclose(odt);
			Write_Solution(CgnsRoot.CgnsBase,"SOLUTIONS/RESTART_DATA/solution.dat");
                        printf("Stopping AnuPravaha as per user request\n");
			exit(1);
		}
	
		CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;
	
		FillNewsolToPrevsol(CgnsRoot.CgnsBase);
	
		if(CgnsRoot.CgnsBase->SolverScheme==0)
		{
			
			do/*iterates till everything converges inside the time loop*/
			{
				for(cntVar=0;cntVar<16;cntVar++)
					ErrorTimeStepCoupled[cntVar]=ErrorTimeStep[cntVar];
	
				do/*mass flux loop*/
				{
					savevolumeflux(CgnsRoot.CgnsBase);
					
					CDcoeffCalculation(CgnsRoot.CgnsBase,1);
	
					do//Iterates mass velocities till convergence
					{
						massVelocityResidual[0]=massVelocityOneIter(CgnsRoot.CgnsBase, 1);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityResidual[1]=massVelocityOneIter(CgnsRoot.CgnsBase, 2);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityResidual[2]=0.00;//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,massVelocityResidual);
						//printf("massvelocity residuals=\n%e\n%e\n%e\n",massVelocityResidual[0],massVelocityResidual[1],massVelocityResidual[2]);
					}
					while(massVelocityIterTag==1);

					volumefluxstar2d(CgnsRoot.CgnsBase);

					PressureIter(CgnsRoot.CgnsBase);

					CalcNewVolumeFlux2d(CgnsRoot.CgnsBase);

					residMassFlux=FindVolumeFluxDifference(CgnsRoot.CgnsBase);

					//printf("mass flux residual =%e\n",residMassFlux);
				}
				while(residMassFlux > CgnsRoot.CgnsBase->BaseParameters->MassFluxConv);
				
				CDcoeffCalculation(CgnsRoot.CgnsBase,1);

				do//Iterates velocities till convergence
				{
					velocityResidual[0] = CalcVariableOneIter(CgnsRoot.CgnsBase,1);
	
					velocityResidual[1] = CalcVariableOneIter(CgnsRoot.CgnsBase,2);
	
					velocityResidual[2] = 0.0;
	
					velocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,velocityResidual);
	
					//printf("velocity residuals=\n%e\n%e\n%e\n",velocityResidual[0],velocityResidual[1],velocityResidual[2]);
				}
				while(velocityIterTag==1);
				
				if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
				{
					CDcoeffCalculation(CgnsRoot.CgnsBase,6);
				}
				for(iScal=0;iScal<10;iScal++)
				{
					if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
					{
						CDcoeffCalculation(CgnsRoot.CgnsBase,600+iScal);
					}
				}
				
				do//Iterates all the variables simulataneously till convergence
				{
					if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
					{
						if(CgnsRoot.CgnsBase->flagConjugateProblem==0 && CgnsRoot.CgnsBase->flagSteadyState==0)/*Conjugate Problem and Steady State*/
						{
							//printf("I am in CalcVariableOneIterConjugate\n");
							variableResidual[0] = CalcVariableOneIterConjugate(CgnsRoot.CgnsBase,6);
							//printf("Temperature Resid[0]=%e\n",variableResidual[0]);
							//ShearStress(CgnsRoot.CgnsBase,6);
						}
						else
						{
							//printf("i am in CalcVariableOneIter");
							variableResidual[0] = CalcVariableOneIter(CgnsRoot.CgnsBase,6);
							//ShearStress(CgnsRoot.CgnsBase,6);
							//printf("Temperature Residual = %e\n",variableResidual[0]);
						}
						//CalcNusseltNo(CgnsRoot.CgnsBase,6);/*This is not a generalised function*/
					}
					else
						variableResidual[0]=0.0;
					//printf("i am in solverCoupled\n");
					for(cntVariable=0;cntVariable<10;cntVariable++)
					{
						if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[cntVariable]!=NULL)
						{
							variableResidual[cntVariable+1]=CalcVariableOneIter(CgnsRoot.CgnsBase,600+cntVariable);
							//ShearStress(CgnsRoot.CgnsBase,600+cntVariable);
							//printf("residual scalar[%d]=%e\n",cntVariable,variableResidual[cntVariable+1]);
						}
						else
							variableResidual[cntVariable+1]=0.0;
					}
					
					variableIterTag = FindVariableIterToken(variableResidual,CgnsRoot.CgnsBase);
					//printf("I am before while variableIterTag==1 \n");
				}
				while(variableIterTag==1);
				
				EvaluateVarCoefficients(CgnsRoot.CgnsBase);// Calculates density,viscosity etc. if they are functions.
				if(CgnsRoot.CgnsBase->flagVariableDensityProblem == 1)
				{
					//EvaluateVarCoefficients(CgnsRoot.CgnsBase);
					for(cntVariable=700;cntVariable<=713;cntVariable++)
						RefreshVariableCoefficients(cntVariable,CgnsRoot.CgnsBase);
				}
				
				FindTimeStepErrorArray(CgnsRoot.CgnsBase,ErrorTimeStep);                              
				insideTimeStepConverToken = FindInsideTimeStepConverToken(ErrorTimeStep,ErrorTimeStepCoupled);
			}//End of inside-time-step convergence loop
			while(insideTimeStepConverToken==1);
		}
		else if(CgnsRoot.CgnsBase->SolverScheme==1)
		{
			do/*iterates till everything converges inside the time loop*/
			{
				for(cntVar=0;cntVar<16;cntVar++)
					ErrorTimeStepCoupled[cntVar]=ErrorTimeStep[cntVar];

				if(CounterTstep>1)
				{
					savevolumeflux(CgnsRoot.CgnsBase);
					
					CDcoeffCalculation(CgnsRoot.CgnsBase,1);
					
					do//Iterates mass velocities till convergence
					{
						massVelocityResidual[0]=massVelocityOneIter(CgnsRoot.CgnsBase, 1);//supply 1 for u, 2 for v and  3 for w and no other value
						
						massVelocityResidual[1]=massVelocityOneIter(CgnsRoot.CgnsBase, 2);//supply 1 for u, 2 for v and  3 for w and no other value
						
						massVelocityResidual[2]=0.00;//supply 1 for u, 2 for v and  3 for w and no other value
						
						massVelocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,massVelocityResidual);
					}
					while(massVelocityIterTag==1);

					volumefluxstar2d(CgnsRoot.CgnsBase);

					PressureIter(CgnsRoot.CgnsBase);

					CalcNewVolumeFlux2d(CgnsRoot.CgnsBase);
				}
				else
				{
					do
					{
						savevolumeflux(CgnsRoot.CgnsBase);
						
						CDcoeffCalculation(CgnsRoot.CgnsBase,1);
						
						do//Iterates mass velocities till convergence
						{
							massVelocityResidual[0]=massVelocityOneIter(CgnsRoot.CgnsBase, 1);//supply 1 for u, 2 for v and  3 for w and no other value
							
							massVelocityResidual[1]=massVelocityOneIter(CgnsRoot.CgnsBase, 2);//supply 1 for u, 2 for v and  3 for w and no other value
							
							massVelocityResidual[2]=0.00;//supply 1 for u, 2 for v and  3 for w and no other value
							
							massVelocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,massVelocityResidual);
						}
						while(massVelocityIterTag==1);
						
						volumefluxstar2d(CgnsRoot.CgnsBase);
						
						PressureIter(CgnsRoot.CgnsBase);
						
						CalcNewVolumeFlux2d(CgnsRoot.CgnsBase);
						
						residMassFlux=FindVolumeFluxDifference(CgnsRoot.CgnsBase);
						
						//printf("mass flux residual =%e\n",residMassFlux);
					}
					while(residMassFlux > CgnsRoot.CgnsBase->BaseParameters->MassFluxConv);
				}
				
				CDcoeffCalculation(CgnsRoot.CgnsBase,1);

				do//Iterates velocities till convergence
				{
					velocityResidual[0] = CalcVariableOneIter(CgnsRoot.CgnsBase,1);
					
					velocityResidual[1] = CalcVariableOneIter(CgnsRoot.CgnsBase,2);
					
					velocityResidual[2] = 0.0;
					
					velocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,velocityResidual);
					
					//printf("velocity residuals=\n%e\n%e\n%e\n",velocityResidual[0],velocityResidual[1],velocityResidual[2]);
				}
				while(velocityIterTag==1);
				
				if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
				{
					CDcoeffCalculation(CgnsRoot.CgnsBase,6);
				}
				for(iScal=0;iScal<10;iScal++)
				{
					if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
					{
						CDcoeffCalculation(CgnsRoot.CgnsBase,600+iScal);
					}
				}
				
				do//Iterates all the variables simulataneously till convergence
				{
					if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
					{
                                        	//printf("I am into solving Temp. Equation\n");
						EvaluateEnergysource(CgnsRoot.CgnsBase);
						if(CgnsRoot.CgnsBase->flagConjugateProblem==0 && CgnsRoot.CgnsBase->flagSteadyState==0)/*Conjugate Problem and Steady State*/
						{
							//printf("I am in CalcVariableOneIterConjugate\n");
							variableResidual[0] = CalcVariableOneIterConjugate(CgnsRoot.CgnsBase,6);
							//printf("variableResidual[0]=%lf\n",variableResidual[0]);
							//ShearStress(CgnsRoot.CgnsBase,6);
						}
						else
						{
							//printf("I am in CalcVariableOneIter\n");
							variableResidual[0] = CalcVariableOneIter(CgnsRoot.CgnsBase,6);
							//ShearStress(CgnsRoot.CgnsBase,6);
							//printf("Temperature Residual = %e\n",variableResidual[0]);
						}
						//CalcNusseltNo(CgnsRoot.CgnsBase,6);/*This is not a generalised function*/
					}
					else
						variableResidual[0]=0.0;
					//printf("i am in solverCoupled\n");
					for(cntVariable=0;cntVariable<10;cntVariable++)
					{
						if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[cntVariable]!=NULL)
						{
							variableResidual[cntVariable+1]=CalcVariableOneIter(CgnsRoot.CgnsBase,600+cntVariable);
							//ShearStress(CgnsRoot.CgnsBase,600+cntVariable);
							//printf("residual scalar[%d]=%e\n",cntVariable,variableResidual[cntVariable+1]);
						}
						else
							variableResidual[cntVariable+1]=0.0;
					}
					
					variableIterTag = FindVariableIterToken(variableResidual,CgnsRoot.CgnsBase);
				}
				while(variableIterTag==1);
				EvaluateVarCoefficients(CgnsRoot.CgnsBase);// Calculates density,viscosity etc. if they are functions.	
				if(CgnsRoot.CgnsBase->flagVariableDensityProblem == 1)
				{
					//EvaluateVarCoefficients(CgnsRoot.CgnsBase);
					for(cntVariable=700;cntVariable<=713;cntVariable++)
						RefreshVariableCoefficients(cntVariable,CgnsRoot.CgnsBase);
				}
				
				FindTimeStepErrorArray(CgnsRoot.CgnsBase,ErrorTimeStep);
				
				insideTimeStepConverToken = FindInsideTimeStepConverToken(ErrorTimeStep,ErrorTimeStepCoupled);
			}//End of inside-time-step convergence loop
			while(insideTimeStepConverToken==1);
		}//End of solverScheme==1
		
		timeStepToken = FindtimeStepToken(CgnsRoot.CgnsBase,ErrorTimeStep,CounterTstep);
		
		SaveTimeStepData(CgnsRoot,CounterTstep,CgnsRoot.CgnsBase->BaseParameters->Time);
		
		for(index_var=0;index_var<total;index_var++)
		{
			vartag = CgnsRoot.CgnsBase->VARTAG[index_var];
			TimeDepend(CgnsRoot.CgnsBase->BaseParameters->Time,vartag,CgnsRoot.CgnsBase);/*Alpha Beta Gamma values have to be updated here*/
			RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);
		}
	}//End of time stepping loop
	while(timeStepToken==1);
	
	/*CON = MassFluxConv(CgnsRoot.CgnsBase);
	fprintf(Q,"Massflux Conservation at all Boundary =%12.8e\n",CON);
	
	if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
	{
		//CON=VarMassFluxConv(CgnsRoot.CgnsBase,6);
		fprintf(Q,"Energy fluxConv at all Boundary =%12.8e\n",CON);
	}
	for(iScal=0;iScal<10;iScal++)
	{
		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
		{
			//CON=VarMassFluxConv(CgnsRoot.CgnsBase,600+iScal);
			fprintf(Q,"Scalar[%d] fluxConv at all Boundary =%12.8e\n",iScal,CON);
		}
	}*/
	fclose(Q);
	return 0;
}
