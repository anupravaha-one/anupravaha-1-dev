#include "APKinitializeRadiation.h"

/**
The Function initializeRadiation() allocated and creates memory  for radiation variable.
This function inturn calls
initializeRadiation_variable() and initializeRadiationSlab() for each variable to create the necessary memory space for that 
radiation variable in Flowsolution structure and also memory for radiation variable in slabs. The functions
initializeRadiationVarCoefficient() and initializeRadiationSlabVarCoeff() are used to initialize memory for 
refractive index in the zone and slab.*/

int initializeRadiation(Base *CgnsBase, int count)
	{
	 FILE *fp;
	 Zone *z;
	 int i, j, k, numZone, numIface, numBc, flag, flagvar, dir, varCoeffCount;
	 double conCri, refindex, intensity,abs,sca,constant,source;
	 
	 fp=fopen("Settings/Convergence.dat","w");
	 
	 z=CgnsBase->CgnsZone;
	 numZone=CgnsBase->Number_of_zones;
	 printf("Enter the Convergence Criterion for Radiative Intensity\n");
	 scanf("%lf",&conCri);
	 CgnsBase->BaseParameters->RadConCriterion=conCri;
	 fprintf(fp,"%e",conCri);
	 fclose(fp);
	 
	 printf("\nEnter Ratio of Diamond Scheme(1->Pure Diamond Scheme) to Step Scheme(0-> Pure Step Scheme):");
	 scanf("%lf",&constant);
	 CgnsBase->BaseParameters->RatioCenUpRadiation = constant;
	 
	 if(CgnsBase->BaseParameters->flagRadiation==1)
	 {
	 	for(i=0;i<numZone;i++)
	 	{
		(z+i)->Flowsol=(FlowSolution *)malloc(sizeof(FlowSolution));
		(z+i)->PrevFlowsol=(PrevFlowSolution *)malloc(sizeof(PrevFlowSolution));
		}
	}
		printf("\nInitialization for Refractive Index of the Medium\n");
		if( CgnsBase->Number_of_zones > 1)
		{
			printf("\nDo you want to specify Refractive Index for the WHOLE DOMAIN [0] or ZONE BY ZONE [1]:\n");
			scanf("%d",&flag);
		} else { flag = 0;}
		
		varCoeffCount = 16;  /*16->Refractive Index, not necessary currently*/ 
		if(flag==0)
			{
			 printf("If Refractive Index is constant [0] or variable [1]\n");
			 scanf("%d",&flagvar);
			 if(flagvar==0)
			 	{
				 printf("\nEnter the constant value of Refractive Index\n");
				 scanf("%lf",&refindex);
				 }
			else
				refindex=0.0;
				
		for(i=0;i<numZone;i++)
			{
		 	(z+i)->Flowsol->RefIndex=(double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+i), refindex);
		 	(z+i)->PrevFlowsol->RefIndex=(double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+i), refindex);
			(z+i)->Flowsol->flagRefIndex=flag; 
			numIface=(z+i)->Number_of_Interfaces;
			numBc=(z+i)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
					{
					for(j=0;j<numIface;j++)
					if(((z+i)->Iface+j)->Direction==dir)
					(((z+i)->Iface+j)->ZSlab)->RefIndex= (double *)initializeRadiationSlabVarCoeff(flagvar, CgnsBase, z+i, refindex, dir, ((z+i)->Iface+j)->Point_Range);
					
					for(k=0;k<numBc;k++)
					if(((z+i)->ZoneBc+k)->Direction==dir)
					(((z+i)->ZoneBc+k)->ZSlab)->RefIndex= (double *)initializeRadiationSlabVarCoeff(flagvar, CgnsBase, z+i, refindex,  dir, ((z+i)->ZoneBc+k)->Point_Range);
					}
			}
			}
		if(flag==1)
		{
		for(i=0;i<numZone;i++)
			{
			 printf("\nZone Nane is %s",(z+i)->zonename);
			 printf("\nNI=%d\tNJ=%d\tNK=%d",*((z+i)->Size+0),*((z+i)->Size+1),*((z+i)->Size+2));
			 printf("\nRefractive Index is constant for this Zone [0] or Variable [1]\n");
			 scanf("%d",&flagvar);
			 if(flagvar==0)
			 {
			 printf("\nEnter the Constant Value of Refractive Index for This Zone\n");
			 scanf("%lf",&refindex);
			 }
			 else
			 refindex=0.0;
			 (z+i)->Flowsol->RefIndex=(double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+i), refindex);
		 	(z+i)->PrevFlowsol->RefIndex=(double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+i), refindex);
			(z+i)->Flowsol->flagRefIndex=flag; 
			numIface=(z+i)->Number_of_Interfaces;
			numBc=(z+i)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
					{
					for(j=0;j<numIface;j++)
					if(((z+i)->Iface+j)->Direction==dir)
					(((z+i)->Iface+j)->ZSlab)->RefIndex= (double *)initializeRadiationSlabVarCoeff(flagvar, CgnsBase, z+i, refindex, dir, ((z+i)->Iface+j)->Point_Range);
					
					for(k=0;k<numBc;k++)
					if(((z+i)->ZoneBc+k)->Direction==dir)
					(((z+i)->ZoneBc+k)->ZSlab)->RefIndex= (double *)initializeRadiationSlabVarCoeff(flagvar, CgnsBase, z+i, refindex, dir, ((z+i)->ZoneBc+k)->Point_Range);
					} 
			}	
		}
	printf("\nDo you want to solve RTE for Absorption-Emission\n");
	printf("Enter 1->Yes,0->NO\n");
	scanf("%d",&(CgnsBase->BaseParameters->flagAbsorption));
	
	for(i=0;i<numZone;i++)
	(z+i)->Flowsol->kappa=(double *)malloc(sizeof(double));
	
	if((CgnsBase->BaseParameters->flagAbsorption)==1)
		{
	printf("\nEnter the Absorption Coefficient as we are assuming absorption coefficient is same for whole domain\n");
	scanf("%lf",&abs);
		 for(i=0;i<numZone;i++)
		 *((z+i)->Flowsol->kappa)=abs;
		}
	
	printf("\nDo you want to solve RTE for Scattering\n");
	printf("Enter 1->Yes,0->NO\n");
	scanf("%d",&(CgnsBase->BaseParameters->flagScattering));
	
	for(i=0;i<numZone;i++)
	(z+i)->Flowsol->sigma=(double *)malloc(sizeof(double));
	
	if((CgnsBase->BaseParameters->flagScattering)==1)
		{
	printf("\nEnter the Scattering Coefficient as we are assuming absorption coefficient is same for whole domain\n");
	scanf("%lf",&sca);
	for(i=0;i<numZone;i++)
		 *((z+i)->Flowsol->sigma)=sca;
		
	printf("\nEnter the Choice of Scattering Function\n");
	scanf("%d",&(CgnsBase->BaseParameters->Scatteringchoice));
		 }
	
	printf("\nInitialization for Radiative Intensity\n");
	i=1;//to restrict repetative memory location
	CgnsBase->VARTAG[count]=16;
	count=count+1;
	for(j=0;j<numZone;j++)
		{
		printf("\nEnter The Initialize value of Radiative Intensity\n");
		scanf("%lf",&intensity);
		(z+j)->Flowsol->RadIntensity=(double *)initializeRadiation_variable(i, CgnsBase, z+j, intensity);
		(z+j)->PrevFlowsol->RadIntensity=(double *)initializeRadiation_variable(i, CgnsBase, z+j, intensity);
		
		for(dir=1;dir<=6;dir++)
			{
			for(k=0;k<=numIface;k++)
				if(((z+j)->Iface+k)->Direction==dir)
				{
				((z+j)->Iface+k)->ZSlab->RadIntensity=(double *)initializeRadiationSlab(i, CgnsBase, dir,  ((z+j)->Iface+k)->Point_Range);
				//((z+j)->Iface+k)->ZSlab->RadIntensity=(double *)initializeRadiationSlab(i, CgnsBase, dir, ((z+j)->Iface+k)->Point_Range);
				}
			for(k=0;k<=numBc;k++)
				if(((z+j)->ZoneBc+k)->Direction==dir)
				{
				((z+j)->ZoneBc+k)->ZSlab->RadIntensity=(double *)initializeRadiationSlab(i,CgnsBase,  dir, ((z+j)->ZoneBc+k)->Point_Range);
				//((z+j)->ZoneBc+k)->ZSlab->RadIntensity=(double *)initializeRadiationSlab(i,CgnsBase, dir, ((z+j)->ZoneBc+k)->Point_Range);
				}
			}
		}
	if(CgnsBase->BaseParameters->flagRadiation==1)
		{
		printf("\nDo you want to Incorporate any Volumetric Source Term: Enter [0->No] or [1->Yes]\n");
		scanf("%d",&(CgnsBase->BaseParameters->flagRadSource));
		if((CgnsBase->BaseParameters->flagRadSource)==1)
		{
		 for(j=0;j<numZone;j++)
		 	{
			 printf("Zone Name=%s\n",(z+j)->zonename);
			 printf("Enter 0/1 for Volumetric Heat Generation in this Zone 1->Yes 0->No :\n");
			 scanf("%d",&((z+j)->flagRadHeatGeneration));
				
				if((z+j)->flagRadHeatGeneration==1)
				{
				printf("Enter whether source is constant[0] at each point or variable[1]: 0->Constant 1->Variable: \n");
				scanf("%d",&flagvar);
					
				(z+j)->flagRadHeatConVar=flagvar; 
				if (flagvar==0)
				{	
					printf("Enter the constant value: ");
					scanf("%lf",&source);
					
				}
				else
				source=1.0;
				((z+j)->Flowsol)->RadEnergysource= (double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+j), source);
				}
					
				else
				{
				(z+j)->flagRadHeatConVar=0;
				source=0.0;
				((z+j)->Flowsol)->RadEnergysource= (double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+j), source);
				}
			}
		}
		else
		{
			for(j=0;j<numZone;j++)
			{
				(z+j)->flagRadHeatGeneration = 0;
				(z+j)->flagRadHeatConVar   = 0;
				((z+j)->Flowsol)->RadEnergysource=(double *)malloc(sizeof(double));
			}
		}
		
		printf("\nDo you want to Enter the Fluid Temperature: Enter [0->No] or [1->Yes]\n");
		scanf("%d",&(CgnsBase->BaseParameters->flagTemp));
		if((CgnsBase->BaseParameters->flagTemp)==1)
		{
		 for(j=0;j<numZone;j++)
		 	{
			 printf("Zone Name=%s\n",(z+j)->zonename);
			 printf("Enter 0/1 for The Fluid Temperature in Kelvin in this Zone 1->Yes 0->No :\n");
			 scanf("%d",&((z+j)->flagZoneTemp));
				
				if((z+j)->flagZoneTemp==1)
				{
					printf("Enter whether Temperature is constant[0] at each point or variable[1]: 0->Constant 1->Variable: \n");
					scanf("%d",&flagvar);
					
					(z+j)->flagRadHeatConVar=flagvar; 
					if(flagvar==0)
					{
					printf("Enter the Temperature value in Kelvin: ");
					scanf("%lf",&source);
					}
					else
					source=0.0;
					//printf("Source=%e\n",source);
				((z+j)->Flowsol)->RadGasTemp= (double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+j), source);
				}
				else
				{
				(z+j)->flagZoneTemp=0;
				source=0.0;
				((z+j)->Flowsol)->RadGasTemp= (double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+j), source);
				}
			}
		}
		else
		{
			for(j=0;j<numZone;j++)
			{
				(z+j)->flagZoneTemp = 0;
				((z+j)->Flowsol)->RadGasTemp=(double *)malloc(sizeof(double));
			}
		}
	}
	if(CgnsBase->BaseParameters->flagflowRadiation==1)
		{
		for(j=0;j<numZone;j++)
		 	{
			 	flagvar=1;
				source=0.0;
				(z+j)->flagRadHeatConVar=1;
				CgnsBase->BaseParameters->flagTemp=1;
				((z+j)->Flowsol)->RadEnergysource= (double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+j), source);
				//((z+j)->Flowsol)->RadGasTemp= (double *)initializeRadiationVarCoefficient(flagvar, CgnsBase, (z+j), source);
			}
		}
	CgnsBase->total_var=count;
	printf("\nTotal number of variable=%d\n", count);
	return(count);
}
