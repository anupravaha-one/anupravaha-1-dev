#include "APKCalcSteady_vel.h"


/**
The function solverSegregated() is the solver for two dimensional geometries which solves all the variables
in a segregated approach. This function is a control function which calls different other functions in an order
as given in the segregated algorithm. See Project Report No.1, December 2005, Section 3.10.0.2 Page 48.

Brief about Algo:
This function calls the massvelocity() to calculate massvelocities.
Calculates volumeflux from the calculated value of massvelocities.
Finds out Pressure using the function PressureIter().
Repeats the above steps till Volume Flux converges.
Computes velocities by using the function CalcVariabler().
Returns to massvelocity calculation till the total residual of all the velocity components are less than their
individual convergence criteria.
Computes all other solvable Scalars like Temperature etc seperately till they individually converge.
Saves the current timestep data if the user chooses so, or saves after a desired no of time steps as choosen by the user
using the function SaveTimeStepData()
If the user chooses to stop the computation at any point(by changing the value to 0 from 1 in the file "iwrite.dat" in the
directory "Settings"), then the solution at that time level is written using the function WriteCgnsFile().

Note: The part of calculation of Mass velocities, Pressure Iteration and Velocities computation is done by the
functions Solutions2d() or Solutions() depending on 2d or 3d domains.
*/
int CalcSteady_vel(Root CgnsRoot)
{
	double time=0.0,dt;
	
	int vartag;
	
	char filename2[] = "SOLUTIONS/restart.cgns";
	int timetag,counter=0;
	double maxerror;
	
	
	dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
	
	timetag = CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag;
	//CgnsRoot.CgnsBase->BaseParameters->Time=0.0;
	
	CDcoeffMemAlloc(CgnsRoot.CgnsBase);
	
	DiffusionCoeffCalculation(CgnsRoot.CgnsBase);
	
	DiffusionCoeffCalculationForVariables(CgnsRoot.CgnsBase);
	
	do
	{
		if (timetag==1 && time==0) /*Even the scheme is CN, the first timestep is done by BTCS scheme*/
			CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=0;
		else if (timetag==1 && time>0)
			CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=1;
		
		//EvaluateVarCoefficients(CgnsRoot.CgnsBase);
		
		time = time+dt;
		counter=counter+1;
		
		CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;
		
		FillNewsolToPrevsol(CgnsRoot.CgnsBase);
		
		if(CgnsRoot.CgnsBase->BaseParameters->flag2d==1)
		{
			//make_w_zero(CgnsRoot.CgnsBase);
			Solutions2d(CgnsRoot.CgnsBase,counter);
		}
		else
		{
			Solutions(CgnsRoot.CgnsBase,counter);
		}
		
		maxerror = CalculateError(CgnsRoot.CgnsBase);
		if(counter%100==0)
		{
			printf("maximum error in velocity time steps %12.8e\n",maxerror);
			printf("Time=%lf\n",CgnsRoot.CgnsBase->BaseParameters->Time);
		}
		
		if(CgnsRoot.CgnsBase->BaseParameters->TstepChoice==0)
		{
			CgnsRoot.CgnsBase->BaseParameters->TimeStep=CalculateTimeStep(CgnsRoot.CgnsBase);
			dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
			printf("deltaT=%lf\n",CgnsRoot.CgnsBase->BaseParameters->TimeStep);
		}
				
		for(vartag=1;vartag<5;vartag++)
			RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);

	}
	while(maxerror > *(CgnsRoot.CgnsBase->BaseParameters->convergenceCriterion));
		
	WriteCgnsFile(filename2,CgnsRoot);
	printsolver(CgnsRoot.CgnsBase);
	
	
	return 0;
}

