#include "APKsolverCoupled.h"
/**
The function solverCoupled() is the solver for three dimensional geometries which solves all the variables
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

int solverCoupled(Root CgnsRoot)
{
	int cntVar,timetag,vartag,CounterTstep=0;
	int total,index_var,iScal;
	double time=0.0,dt;
	FILE *save,*odt,*Q,*fp;
	char filename2[] = "SOLUTIONS/restart.cgns";
	int timeStepToken,insideTimeStepConverToken,variableIterTag,velocityIterTag,cntVariable,massVelocityIterTag,stop;
	double massVelocityResidual[3],velocityResidual[3],variableResidual[11];
	double ErrorTimeStepCoupled[27],ErrorTimeStep[27]; //Shyam, size changed to 16 to accommodate combustion variables.
	double residMassFlux;
	
	Q=fopen("SOLUTIONS/MassConservation.dat","w");
	dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
	printf("\nInside SolverCoupled.c\n");
	
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
			WriteCgnsFile(filename2,CgnsRoot);
			WriteCgnsFile_vertex("SOLUTIONS/vertex.cgns",CgnsRoot);
			printsolver(CgnsRoot.CgnsBase);
			dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
			odt=fopen("Settings/old_dt.dat","w");
			fprintf(odt,"%12.10lf %12.10lf",dt,CgnsRoot.CgnsBase->BaseParameters->Time);
			fclose(odt);                        
	                Write_Solution(CgnsRoot.CgnsBase,"SOLUTIONS/RESTART_DATA/solution.dat");
                        printf("Stopping AnuPravaha as per user request\n");
			//Write_gnuplot_style_solution(CgnsRoot.CgnsBase,"SOLUTIONS/gnuplot/gnu_solution.dat");
			exit(1);
		}
		CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;
		
		FillNewsolToPrevsol(CgnsRoot.CgnsBase);
		
		do/*iterates till everything converges inside the time loop*/
		{
			for(cntVar=0;cntVar<27;cntVar++)
				ErrorTimeStepCoupled[cntVar]=ErrorTimeStep[cntVar];
			if(CgnsRoot.CgnsBase->SolverScheme==1)
			{
				if(CounterTstep<1)
				{
					do/*mass flux loop*/
					{
						savevolumeflux(CgnsRoot.CgnsBase);
						
						CDcoeffCalculation(CgnsRoot.CgnsBase,1);
						
						do//Iterates mass velocities till convergence
						{
							massVelocityResidual[0] = massVelocityOneIter(CgnsRoot.CgnsBase, 1);//supply 1 for u, 2 for v and  3 for w and no other value
							
							massVelocityResidual[1] = massVelocityOneIter(CgnsRoot.CgnsBase, 2);//supply 1 for u, 2 for v and  3 for w and no other value
							
							massVelocityResidual[2] = massVelocityOneIter(CgnsRoot.CgnsBase, 3);//supply 1 for u, 2 for v and  3 for w and no other value
							
							massVelocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,massVelocityResidual);		
						}
						while(massVelocityIterTag==1);

						volumefluxstar(CgnsRoot.CgnsBase);//Calculates new Volumeflux from mass velocities

						PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux

						CalcNewVolumeFlux(CgnsRoot.CgnsBase);

						residMassFlux=FindVolumeFluxDifference(CgnsRoot.CgnsBase);

						//printf("mass flux residual =%e\n",residMassFlux);
					}
					while(residMassFlux > CgnsRoot.CgnsBase->BaseParameters->MassFluxConv);
				}
				else
				{
					CDcoeffCalculation(CgnsRoot.CgnsBase,1);
				
					do//Iterates mass velocities till convergence
					{
						massVelocityResidual[0] = massVelocityOneIter(CgnsRoot.CgnsBase, 1);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityResidual[1] = massVelocityOneIter(CgnsRoot.CgnsBase, 2);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityResidual[2] = massVelocityOneIter(CgnsRoot.CgnsBase, 3);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,massVelocityResidual);		
					}

					while(massVelocityIterTag==1);

					volumefluxstar(CgnsRoot.CgnsBase);//Calculates new Volumeflux from mass velocities

					PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux

					CalcNewVolumeFlux(CgnsRoot.CgnsBase);
				}
			}
			else if(CgnsRoot.CgnsBase->SolverScheme==0)
			{
				do/*mass flux loop*/
				{
					savevolumeflux(CgnsRoot.CgnsBase);
					
					CDcoeffCalculation(CgnsRoot.CgnsBase,1);

					do//Iterates mass velocities till convergence
					{
						massVelocityResidual[0] = massVelocityOneIter(CgnsRoot.CgnsBase, 1);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityResidual[1] = massVelocityOneIter(CgnsRoot.CgnsBase, 2);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityResidual[2] = massVelocityOneIter(CgnsRoot.CgnsBase, 3);//supply 1 for u, 2 for v and  3 for w and no other value

						massVelocityIterTag=FindVelocityIterToken(CgnsRoot.CgnsBase,massVelocityResidual);		
					}
					while(massVelocityIterTag==1);

					volumefluxstar(CgnsRoot.CgnsBase);//Calculates new Volumeflux from mass velocities

					PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux

					CalcNewVolumeFlux(CgnsRoot.CgnsBase);

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
				
				velocityResidual[2] = CalcVariableOneIter(CgnsRoot.CgnsBase,3);
				
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
					EvaluateEnergysource(CgnsRoot.CgnsBase);
					if(CgnsRoot.CgnsBase->flagConjugateProblem==0 && CgnsRoot.CgnsBase->flagSteadyState==0)/*Conjugate Problem and Steady State*/
					{
						//printf("i am in CalcVariableOneIter 6");
						variableResidual[0] = CalcVariableOneIterConjugate(CgnsRoot.CgnsBase,6);
						//printf("variableResidual[0]=%lf\n",variableResidual[0]);
						//ShearStress(CgnsRoot.CgnsBase,6);
					}
					else
					{
						//printf("i am in CalcVariableOneIterConjugate 6");
						variableResidual[0] = CalcVariableOneIter(CgnsRoot.CgnsBase,6);
						//ShearStress(CgnsRoot.CgnsBase,6);
						//printf("Temperature Residual = %e\n",variableResidual[0]);
					}
					//CalcNusseltNo(CgnsRoot.CgnsBase,6);/*This is not a generalised function*/
				}
				else
					variableResidual[0]=0.0;
				
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
				
				variableIterTag=FindVariableIterToken(variableResidual,CgnsRoot.CgnsBase);
			}
			while(variableIterTag==1);
			printf("\n");
			EvaluateVarCoefficients(CgnsRoot.CgnsBase);// Calculates density,viscosity etc. if they are functions.
			if(CgnsRoot.CgnsBase->flagVariableDensityProblem == 1)
			{
				//EvaluateVarCoefficients(CgnsRoot.CgnsBase);
				for(cntVariable=700;cntVariable<=713;cntVariable++)
					RefreshVariableCoefficients(cntVariable,CgnsRoot.CgnsBase);
			}
			
			FindTimeStepErrorArray(CgnsRoot.CgnsBase,ErrorTimeStep);


			insideTimeStepConverToken=FindInsideTimeStepConverToken(ErrorTimeStep,ErrorTimeStepCoupled);
		}//End of inside-time-step convergence loop
		while(insideTimeStepConverToken==1);

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
		CON=VarMassFluxConv(CgnsRoot.CgnsBase,6);
		fprintf(Q,"Energy fluxConv at all Boundary =%12.8e\n",CON);
	}
	for(iScal=0;iScal<10;iScal++)
	{
		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
		{
			CON=VarMassFluxConv(CgnsRoot.CgnsBase,600+iScal);
			fprintf(Q,"Scalar[%d] fluxConv at all Boundary =%12.8e\n",iScal,CON);
		}
	}*/
	fclose(Q);
	return 0;
}
double CalcVariableError(int variableTag,Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k,totalcells=0,SolveSolidTag=1;
	double *domain,*domainPrev;
	double rmsError,delta_t;
	double errorVariable=0.0;
	Zone *zone;
	
	if(variableTag==1||variableTag==2||variableTag==3 || (base->flagTurb==1 && (variableTag==600||variableTag==601)) ||variableTag==62||variableTag==63||variableTag==64||variableTag==65||variableTag==710||variableTag==68||variableTag==69) //variableTag==67|| as number-density equation is no longer solved for
	{
		SolveSolidTag = 0; 
		totalcells = base->TotalCellFluid;
	}
	else
	{
		SolveSolidTag = 1;
		totalcells = base->TotalCellDomain;
	}
		
	if(variableTag==600||variableTag==601)
		delta_t = base->BaseParameters->Turb->DualTimeStep;
	else
		delta_t = base->BaseParameters->TimeStep;
	
	numZones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			//totalcells = totalcells+(NI-1)*(NJ-1)*(NK-1);
		
			domain=MakeDomainArray(zone,variableTag);
			domainPrev=MakeDomainArray(zone,variableTag*10);
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				errorVariable = errorVariable + ((*(domain+index) - *(domainPrev+index))*(*(domain+index) - *(domainPrev+index)));///(delta_t*delta_t);
			}
		}
	}
	rmsError = sqrt(errorVariable/totalcells);
	rmsError = rmsError/delta_t;
	return(rmsError);
}

//##########################################################################################

int FindTimeStepErrorArray(Base *base,double rmsError[27])
{
	int indexVar,nop,nophase,npcount;
	
	nop = base->nphase;
	if(base->flagMPFV==1)
	{
		rmsError[0] = CalcVariableError(31,base);

		rmsError[1] = CalcVariableError(32,base);

		rmsError[2] = CalcVariableError(33,base); 

		npcount = 0;
		for(nophase=0;nophase<nop;nophase++)
		{	
			rmsError[3+3*nophase+npcount] = CalcVariableError(41+3*nophase,base);

			rmsError[4+3*nophase+npcount] = CalcVariableError(42+3*nophase,base);
			
			rmsError[5+3*nophase+npcount] = CalcVariableError(43+3*nophase,base); 
			rmsError[6+3*nophase+npcount] = CalcVariableError(34+nophase,base);
			npcount++;
		}
		//for(i=0;i<11;i++)
			//printf("\nrmsError[%d]is--%e",i,rmsError[i]);
		//Added by chandrasekhar to compute the residuals in case of turbulent MPFV
		for(indexVar=0 ;indexVar<10 ;indexVar++)
		{
			if(base->CgnsZone->Flowsol->Scal[indexVar]!=NULL)
			{	
				rmsError[indexVar+3+4*nop]=CalcVariableError(indexVar+600,base);
				//printf("\nrmsError[%d]is--%e",indexVar+3+4*nop,rmsError[indexVar+3+4*nop]);
			}
			else
				rmsError[indexVar+3+4*nop]=0.0;
		}
		

	}
	else
	{
	rmsError[0] = CalcVariableError(1,base);
	rmsError[1] = CalcVariableError(2,base);
	if(base->BaseParameters->flag2d==0)
	rmsError[2] = CalcVariableError(3,base); 
	else
	rmsError[2] = 0.0;
	
	

        if(base->flagTurbComb==0) //Shyam
        {
		if(base->CgnsZone->Flowsol->T!=NULL)
                {
			rmsError[3] = CalcVariableError(6,base);
			//printf(" Temp residuals = %e\n",rmsError[3]);
		}
        }
        else if(base->flagSolidification==1||base->flagDiffKinetCont==2) //Shyam
        {
                if(base->CgnsZone->Flowsol->enthalpy!=NULL)
                        rmsError[3] = CalcVariableError(7,base);
	}
        else
        rmsError[3]=0.0; //Shyam

	for(indexVar=0 ;indexVar<10 ;indexVar++)
	{
		if(base->CgnsZone->Flowsol->Scal[indexVar]!=NULL)
			rmsError[indexVar+4]=CalcVariableError(indexVar+600,base);
		else
			rmsError[indexVar+4]=0.0;
	}

        if(base->CgnsZone->Flowsol->MixFrac!=NULL) //Shyam
	rmsError[14] = CalcVariableError(14,base);
	else
	rmsError[14]=0.0;

        if(base->CgnsZone->Flowsol->FuelMassFrac!=NULL) //Shyam
	rmsError[15] = CalcVariableError(15,base);
	else
	rmsError[15]=0.0;
	}


	if(base->CgnsZone->Flowsol->U_Disperse!=NULL) //mrunalini
		rmsError[16] = CalcVariableError(62,base);
	else
		rmsError[16]=0.0;
	
	if(base->CgnsZone->Flowsol->V_Disperse!=NULL) //mrunalini
		rmsError[17] = CalcVariableError(63,base);
	else
		rmsError[17]=0.0;
	
	if(base->CgnsZone->Flowsol->W_Disperse!=NULL) //mrunalini
		rmsError[18] = CalcVariableError(64,base);
	else
		rmsError[18]=0.0;
	
	if(base->CgnsZone->Flowsol->VarthetaDisperse!=NULL) //adhi+shash
		rmsError[19] = CalcVariableError(65,base);
	else
		rmsError[19]=0.0;
	
/*	if(base->CgnsZone->Flowsol->NumDensity!=NULL) //adhi+shash
		rmsError[20] = CalcVariableError(67,base);
	else
		rmsError[20]=0.0;
*/	
	if((base->CgnsZone->Flowsol->Rho!=NULL)&&(base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1)) //shash
	{
		rmsError[21] = CalcVariableError(710,base);
	}
	else
		rmsError[21]=0.0;
	
	if(base->CgnsZone->Flowsol->DropFuelMassFrac!=NULL) //shash for fuel mass frac
 		rmsError[22] = CalcVariableError(68,base);
 	else
 		rmsError[22]=0.0;
	
	if(base->CgnsZone->Flowsol->Tdrop!=NULL) //shash for droplet temperature 
 		rmsError[23] = CalcVariableError(69,base);
 	else
 		rmsError[23]=0.0;
	if(base->flaginduced==1)
	{
		rmsError[24] = CalcVariableError(56,base);
		rmsError[25] = CalcVariableError(57,base);
		rmsError[26] = CalcVariableError(58,base); 
	}
	else
	{
		rmsError[24] = 0.0;
		rmsError[25] = 0.0;
		rmsError[26] = 0.0;
	}
	return 0;
}

//###################################################################################################

int FindInsideTimeStepConverToken(double ErrorNew[27],double ErrorBeforeIteration[27])
{
	int token=0,cntVar;
	//for(cntVar=0;cntVar<14;cntVar++)
	//	printf("\nVariable Count = %d Difference in ErrorNew and ErrorBeforeIteration is %lf",cntVar,(fabs(ErrorNew[cntVar]-ErrorBeforeIteration[cntVar])));

	for(cntVar=0;cntVar<27;cntVar++)
		if(fabs(ErrorNew[cntVar]-ErrorBeforeIteration[cntVar])>1.0e-2)
			token=1;
	return(token);
}

//###################################################################################################

int FindtimeStepToken(Base *base,double ErrorTimeStep[27],int CounterTstep) //Size changed to 16 to accomodate two combustion variables, Shyam
{
	int token=0,cntVar,nop,che;
	
	nop = base->nphase;
	if(base->flagMPFV==1)
	{
		che=3+4*nop;
		for(cntVar=0;cntVar<3+4*nop;cntVar++)
			if(ErrorTimeStep[cntVar]>(*(base->BaseParameters->convergenceCriterion+0)))
				token=1;
		
		//added by chandrasekhar to compute the residuals in case of turbulent MPFV
		if(base->flagTurb==1)
		{
			for(cntVar=che;cntVar<che+2;cntVar++)
				if(ErrorTimeStep[cntVar]>(*(base->BaseParameters->convergenceCriterion+cntVar-9)))
					token=1;
		}
	}
	else if(base->flagMPFV==0)
	{

		if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia ==1)//Mrunalini
		{
			for(cntVar=0;cntVar<3;cntVar++)
				if(ErrorTimeStep[cntVar]>1e-6)
					token=1;
			for(cntVar=4;cntVar<=4;cntVar++)
				if(ErrorTimeStep[cntVar]>1e-4)
					token=1;
		}
		else
		{
			for(cntVar=0;cntVar<3;cntVar++)
				if(ErrorTimeStep[cntVar]>(*(base->BaseParameters->convergenceCriterion+0)))
					token=1;  
			
			for(cntVar=3;cntVar<16;cntVar++)
			if(ErrorTimeStep[cntVar]>(*(base->BaseParameters->convergenceCriterion+cntVar-2)))
				token=1;  
		}
		
	}
	for(cntVar=16;cntVar<19;cntVar++)
		if(ErrorTimeStep[cntVar]>1e-6) //(*(base->BaseParameters->convergenceCriterion+0)) //disperse phase velocity components
			token=1;  
	for(cntVar=19;cntVar<20;cntVar++)
		if(ErrorTimeStep[cntVar]>1e-8)
			token=1;  
/*	for(cntVar=20;cntVar<=20;cntVar++)
		if(ErrorTimeStep[cntVar]>1e5)
			token=1;
*/
	for(cntVar=21;cntVar<=21;cntVar++)
		if(ErrorTimeStep[cntVar]>1e-5)
			token=1;  
	for(cntVar=22;cntVar<=22;cntVar++)
		if(ErrorTimeStep[cntVar]>1e-6)
			token=1;  
	for(cntVar=23;cntVar<=23;cntVar++)
		if(ErrorTimeStep[cntVar]>1e-4) //temperature of droplet
			token=1;  
	for(cntVar=24;cntVar<27;cntVar++)
		if(ErrorTimeStep[cntVar]>(*(base->BaseParameters->convergenceCriterion+0)))
			token=1;  
	
	if(CounterTstep%1 == 0)
		printTimeStepResid(base,ErrorTimeStep);
	//printf("cntVar_prime = %d \t ErrorTimeStep=%lf\n",cntVar_prime,ErrorTimeStep[cntVar_prime]);
	//printf("After printTimeStepResid\n");
	return(token);
}

//###################################################################################################

//########################################################################################################################################
int FindVelocityIterToken(Base *base,double velocityResidual[3])
{
	int token=0,cntVelocity;
	
	for(cntVelocity=0;cntVelocity<3;cntVelocity++)
	{
		if(velocityResidual[cntVelocity] > *(base->BaseParameters->convergenceCriterion))
			token=1;
	}
	return (token);
}

//#########################################################################################################################################
int FindVariableIterToken(double variableResidual[11],Base *base)
{
	int token=0,cntVariable;
	//FILE *convcri;
// 	convcri=fopen("Settings/convergenceCriteria.dat","r");
// 	for(cntVariable=0;cntVariable<12;cntVariable++)
// 		 fscanf(convcri,"%lf",(base->BaseParameters->convergenceCriterion+cntVariable));
// 	fclose(convcri);
	for(cntVariable=0;cntVariable<11;cntVariable++)
	{
		if(variableResidual[cntVariable] > *(base->BaseParameters->convergenceCriterion+1+cntVariable))
		{
			token=1;
		}
	}
	return(token);
}

//#########################################################################################################################################

int printTimeStepResid(Base *base,double ErrorTimeStep[27])
{
	int cntVariable,nop,nophase,npcount;
	FILE *errEvol;

	errEvol = fopen("SOLUTIONS/errEvolution.dat","a");
// 	fp = fopen("SOLUTIONS/VariableResidual.dat","w");

	printf("***************************\nCurrent Time = %lf\n",base->BaseParameters->Time);
// 	fprintf(fp,"Current Time = %lf\n",base->BaseParameters->Time);
	
	printf("U velocity residual= %e\nV velocity residual= %e\n",ErrorTimeStep[0],ErrorTimeStep[1]);
	fprintf(errEvol,"%12.8lf\t%e\t%e\t",base->BaseParameters->Time,ErrorTimeStep[0],ErrorTimeStep[1]);
	if(base->CgnsZone->Flowsol->w!=NULL)
	{
		printf("W velocity residual= %e\n",ErrorTimeStep[2]);
		fprintf(errEvol,"%e\t",ErrorTimeStep[2]);
	}
// 	fprintf(fp,"U velocity residual= %e\nV velocity residual= %e\nW velocity residual= %e\n",ErrorTimeStep[0],ErrorTimeStep[1],ErrorTimeStep[2]);
	nop =base->nphase;
	npcount =0;
	if(base->flagMPFV==1)
	{
		for(nophase=0;nophase<nop;nophase++)
		{
			printf("U velocity residual Phase[%d]= %e\nV velocity residualPhase[%d]= %e\nW velocity residualPhase[%d]= %e\nVartheta residualPhase[%d]=%e\n",nophase,ErrorTimeStep[3+3*nophase+npcount],nophase,ErrorTimeStep[4+3*nophase+npcount],nophase,ErrorTimeStep[5+3*nophase+npcount],nophase,ErrorTimeStep[6+3*nophase+npcount]);
			npcount++;
		}
	}
	if(base->flagTurbComb==0)//Shyam
	{
		if(base->CgnsZone->Flowsol->T!=NULL)
		{
			printf("Temperature residual= %e\n",ErrorTimeStep[3]);
			fprintf(errEvol,"%e\t",ErrorTimeStep[3]);
// 			fprintf(fp,"Temperature residual= %e\n",ErrorTimeStep[3]);
		}
	}
	if(base->flagTurbComb==0||base->flagDiffKinetCont==2)//Shyam
	{
		if(base->CgnsZone->Flowsol->enthalpy!=NULL)
		{
			printf("enthalpy residual= %e\n",ErrorTimeStep[3]);
			fprintf(errEvol,"%e\t",ErrorTimeStep[3]);
// 			fprintf(fp,"enthalpy residual= %e\n",ErrorTimeStep[3]);
		}
	}
	if(base->CgnsZone->Flowsol->MixFrac!=NULL)//Shyam
	{
		printf("Mixture Fraction residual= %e\n",ErrorTimeStep[14]);
		fprintf(errEvol,"%e\t",ErrorTimeStep[14]);
// 		fprintf(fp,"Mixture Fraction residual= %e\n",ErrorTimeStep[14]);
	}
	if(base->flagDiffKinetCont==2)//Shyam
	{
		if(base->CgnsZone->Flowsol->FuelMassFrac!=NULL)
		{
			printf("Fuel Mass Fraction residual= %e\n",ErrorTimeStep[15]);
			fprintf(errEvol,"%e\t",ErrorTimeStep[15]);
// 			fprintf(fp,"Mixture Fraction residual= %e\n",ErrorTimeStep[15]);
		}
	}
	for(cntVariable=0;cntVariable<10;cntVariable++)
	{
		if(base->CgnsZone->Flowsol->Scal[cntVariable]!=NULL)
		{
			if(base->flagMPFV==0)	
			{
				printf("Scalar[%d] residual= %e\n",cntVariable,ErrorTimeStep[cntVariable+4]);
				fprintf(errEvol,"%e\t",ErrorTimeStep[cntVariable+4]);
// 			fprintf(fp,"Scalar[%d] residual= %e\n",cntVariable,ErrorTimeStep[cntVariable+4]);
			}
			if(base->flagMPFV==1)
			{
				printf("Scalar[%d] residual= %e\n",cntVariable,ErrorTimeStep[cntVariable+3+4*nop]);
				fprintf(errEvol,"%e\t",ErrorTimeStep[cntVariable+3+4*nop]);
			}
		}
	}
	if((base->flagParticle==1)||(base->flagDroplet==1))
	{
		printf("U Disperse-Phase Velocity residual= %e\nV Disperse-Phase Velocity residual= %e\nW Disperse-Phase Velocity residual= %e\n",ErrorTimeStep[16],ErrorTimeStep[17],ErrorTimeStep[18]);
		fprintf(errEvol,"\t%e\t%e\t%e\t",ErrorTimeStep[16],ErrorTimeStep[17],ErrorTimeStep[18]);
// 	fprintf(fp,"U Disperse-Phase Velocity residual= %e\nV Disperse-Phase Velocity residual= %e\nW Disperse-Phase Velocity residual= %e\n",ErrorTimeStep[16],ErrorTimeStep[17],ErrorTimeStep[18]);
		printf("Volume Fraction Disperse residual= %e\n",ErrorTimeStep[19]);
		fprintf(errEvol,"\t%e\t",ErrorTimeStep[19]);
// 	fprintf(fp,"Volume Fraction Disperse residual= %e\n",ErrorTimeStep[19]);
	/*	if((base->flagDroplet==1))
		{
			printf("Number Density residual= %e\n",ErrorTimeStep[20]);
			fprintf(errEvol,"\t%e\t",ErrorTimeStep[20]);
// 			fprintf(fp,"Number Density residual= %e\n",ErrorTimeStep[20]);
		}
*/
		if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
		{
			printf("Carrier Density residual= %e\n",ErrorTimeStep[21]);
			fprintf(errEvol,"\t%e\t",ErrorTimeStep[21]);
// 			fprintf(fp,"Carrier Density residual= %e\n",ErrorTimeStep[21]);
		}
		if((base->BaseParameters->flagDropFuelMassFrac==1))
	    	{
	 		printf("Fuel Mass Frac residual= %e\n",ErrorTimeStep[22]);
			fprintf(errEvol,"\t%e\t",ErrorTimeStep[22]);
// 			fprintf(fp,"Carrier Density residual= %e\n",ErrorTimeStep[21]);
	      	}
		if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
		{
	 		printf("Droplet Temp residual= %e\n",ErrorTimeStep[23]);
			fprintf(errEvol,"\t%e\t",ErrorTimeStep[23]);
// 			fprintf(fp,"Carrier Density residual= %e\n",ErrorTimeStep[21]);
	      	}       
	}
	if(base->flaginduced==1)
	{
		printf("Bx residual= %e\nBy residual= %e\nBz residual= %e\n",ErrorTimeStep[24],ErrorTimeStep[25],ErrorTimeStep[26]);
		fprintf(errEvol,"%12.8lf\t%e\t%e\t%e\t",base->BaseParameters->Time,ErrorTimeStep[24],ErrorTimeStep[25],ErrorTimeStep[26]);
		//fprintf(fp,"Bx residual= %e\nBy residual= %e\nBz residual= %e\n",ErrorTimeStep[16],ErrorTimeStep[17],ErrorTimeStep[18]);
//		fprintf(fp,"Carrier Density residual= %e\n",ErrorTimeStep[21]);
	}
	fprintf(errEvol,"\n");
	printf("***************************\n");
	
	fclose(errEvol);
// 	fclose(fp);
	return 0;
}

int printTimeStepResidOneline(Base *base)
{
	int cntVariable,nop,nophase;
	FILE *errEvol;
	nop =base->nphase;

	errEvol = fopen("SOLUTIONS/errEvolution.dat","a");

//      fprintf(errEvol,"#variables = \n\"Time\"\t\t \"U rms\"\t \"V rms\"\t ");
        fprintf(errEvol,"#\"Time\"\t\t \"U rms\"\t \"V rms\"\t ");
	if(base->CgnsZone->Flowsol->w!=NULL) fprintf(errEvol,"\"W rms\"\t ");
	if(base->flagTurbComb==0)
	{
		if(base->CgnsZone->Flowsol->T!=NULL) { fprintf(errEvol,"\"Temperature rms\t"); }
	}
	if(base->flagTurbComb==0||base->flagDiffKinetCont==2)
	{
		if(base->CgnsZone->Flowsol->enthalpy!=NULL) { fprintf(errEvol,"\"Enthalpy rms\t"); }
	}
	if(base->CgnsZone->Flowsol->MixFrac!=NULL) { fprintf(errEvol,"\"Mix.Frac. rms\t"); }
	if(base->flagDiffKinetCont==2)
	{
		if(base->CgnsZone->Flowsol->FuelMassFrac!=NULL) { fprintf(errEvol,"\"Fuel Mass Frac. rms\t"); }
	}
	for(cntVariable=0;cntVariable<10;cntVariable++)
	{
		if(base->CgnsZone->Flowsol->Scal[cntVariable]!=NULL)
		{
			if(base->flagMPFV==0) { fprintf(errEvol,"\"Scalar[%d] rms\"\t",cntVariable+4); }
			if(base->flagMPFV==1) { fprintf(errEvol,"\"Scalar[%d] rms\"\t",cntVariable+3+4*nop); }
		}
	}
	if((base->flagParticle==1)||(base->flagDroplet==1))
	{
		fprintf(errEvol,"\"U Dis-Ph rms\" \t\"V Dis-Ph rms\" \t\"W Dis-Ph rms\"\t");
		fprintf(errEvol,"\"Vol.Frac.Dis rms\"\t");
		if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
		{
			fprintf(errEvol,"\"Carr.Den. rms\"\t");
		}
		if((base->BaseParameters->flagDropFuelMassFrac==1))
	    	{
			fprintf(errEvol,"\"FuelMass.Frac rms\"\t");
	      	}
		if((base->flagDroplet==1)&&(base->BaseParameters->flagVariableDia==1))
		{
			fprintf(errEvol,"\"DropletTemp rms\"\t");
	      	}       
	}
	if(base->flaginduced==1) { fprintf(errEvol,"\"Bx rms\"\t\"By rms\"\t\"Bz rms\"\t"); }

	fprintf(errEvol,"\n");
	fclose(errEvol);
}
//#########################################################################################################################################
