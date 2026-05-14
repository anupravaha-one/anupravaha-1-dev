#include "APKSolverElectricSemiCoupled.h"

/**
This function SolverElectricSemiCoupled() is the solver for three dimensional geometries which solves all electric variables in a semicoupled approach.
*/ 
int SolverElectricSemiCoupled(Root CgnsRoot)
{

	double ChargeError,PotentialError;
	double time=0.0,dt,ConvCriteria=1e-8;
	int stop,timeStepToken, count=0, val=1;
	int numZones, nZones;

	Zone *zone;
	FILE *save, *Tstep, *RMS;

	char InProcessFile[] = "SOLUTIONS/Output_Electric.cgns";
	char MidBreachFile[] = "SOLUTIONS/restart.cgns";
#ifdef FILE_OPEN
{

 	RMS = fopen("Settings/ESPOutput/ErrorEvolution_Electric.dat", "w");
 	fprintf(RMS,"VARIABLES = \"PotentialResidue\"\n\"ChargeResidue\"\n\"TimeSteps\"\n\n");
 	fclose(RMS);//This is done to remove any previous data
// 
 	Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","w");
 	fprintf(Tstep,"%d", 1);//This is done to remove any garbage value stored earlier
 	fclose(Tstep); 
}
#endif

	dt = CgnsRoot.CgnsBase->BaseParameters->EfieldTimeStep;
	ConvCriteria = CgnsRoot.CgnsBase->BaseParameters->EfieldConvCri;
	CgnsRoot.CgnsBase->BaseParameters->Time=0.0;
	
	CgnsRoot.CgnsBase->ChargeConv = 0;//This should always be equal to 1 for iterative boundary
	
	Coeffcalculation_ESP(CgnsRoot.CgnsBase);

	do
	{
		stop = 1;time+=dt;count++;
		
		//Setting ZoneWireTag to Zero (To be used while iterating charge density)
		nZones=(CgnsRoot.CgnsBase->Number_of_zones);
		for(numZones=0;numZones<nZones;numZones++)
		{
			zone = (CgnsRoot.CgnsBase->CgnsZone+numZones);
			zone->ZoneWireTag = 0;
			
		}
		
#ifdef FILE_OPEN
{		
		//print the number of timesteps of the problem each time
		Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","w");
		fprintf(Tstep,"%d", count);
		fclose(Tstep); 
		
		save=fopen("Settings/iwrite.dat","r");
		fscanf(save,"%d",&stop);
		fclose(save);
}
#endif
		if(stop==0)
		{
			WriteCgnsFile(MidBreachFile,CgnsRoot);
			printsolver(CgnsRoot.CgnsBase);
			exit(1);
		}
		else if(stop==2)
		{
			printf("\nI am Writing the solution into \"%s\" file\n", InProcessFile);
			WriteCgnsFile(InProcessFile,CgnsRoot);
			save=fopen("Settings/iwrite.dat","w");
			fprintf(save,"%d",1);
			fclose(save);
		}
		
		CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;
		
		FillNewsolToPrevsol(CgnsRoot.CgnsBase); 
		
		
		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Potential!=NULL)
		{
			CalcSource_opti(CgnsRoot.CgnsBase,8);
			CalcVariablePotential(CgnsRoot.CgnsBase,8);
			PotentialError = CalcErrorTimeStep(CgnsRoot.CgnsBase,8,dt);
		}
		else
		{
			PotentialError=0.0;
			printf("\nUnable to solve for Potential\nMemory Not Allocated, Exiting.......\n");
			exit(1);	
		}
		 
		//Calculate Electric Flux
		ElectricFluxXYZ(CgnsRoot.CgnsBase);
		
		CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,9);
		Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,9);
		CDcoeffCalculation(CgnsRoot.CgnsBase,9);
		CalcSource_opti(CgnsRoot.CgnsBase,9);
		
		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->RhoIon!=NULL)
		{
			CalcVariableCharge(CgnsRoot.CgnsBase,9);
			ChargeError = CalcErrorTimeStep(CgnsRoot.CgnsBase,9,dt);
		}
		else
		{
			ChargeError=0.0;
			printf("\nUnable to solve for Ion Charge Density\nMemory Not Allocated, Exiting.......\n");
			exit(1);
		}

		//if((count%1)==0)
		printf("\nmaximum error:\tPotential = %12.8e\tIon Charge Density = %12.8e\tTime Step=%d",PotentialError,ChargeError,count);
		
		//sav=fopen("Settings/ESPInput/WireChargeUpdateData.dat","r");
		//fscanf(sav,"%lf",&Omega);
		//fscanf(sav,"%lf",&Gamma);
		//fscanf(sav,"%d",&stp1);
		//fclose(sav);
		//if(count%1==0)
		//{
#ifdef FILE_OPEN
{
	
		RMS = fopen("Settings/ESPOutput/ErrorEvolution_Electric.dat", "a");
		fprintf(RMS,"%e\t%e\t%d\n",PotentialError,ChargeError,count);
		fclose(RMS);
}
#endif
		//}
		
		if((ChargeError>ConvCriteria) || (PotentialError>ConvCriteria))
			timeStepToken = 1;
		else
		{
			//conv=fopen("Settings/ESPInput/IonChargeConvergency.dat","r");
			//fscanf(conv,"%d",&val);
			//fclose(conv);
			val = CgnsRoot.CgnsBase->ChargeConv;
			if(val==0)
				timeStepToken = 0;
			else
				timeStepToken = 1;
			
		}
		
		SaveTimeStepData(CgnsRoot,count,CgnsRoot.CgnsBase->BaseParameters->Time);

	}
	while(timeStepToken==1);
	
	//Calculate Current density and total current on wire and plate
	WireCurrentDensity(CgnsRoot.CgnsBase,8);
	
	Free_Coefficients_ESP(CgnsRoot.CgnsBase);
	
	if(CgnsRoot.CgnsBase->flagEfield==1)//Only for fluid flow problems we need to calculate N-S source term.
		ESP_Source_Calculation(CgnsRoot.CgnsBase);
	
	return 0; 
}


