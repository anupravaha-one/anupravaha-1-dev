#include "APKinitialize.h"

/**
The Function initialize() asks for which variables memory should be allocated and creates memory only for those
variables.

This function inturn calls
initialize_variable() and initializeslab() for each variable to create the necessary memory space for that 
variable in Flowsolution structure and also memory for that variable in slabs. The functions initializeFlux(),
initializeVarCoefficient() and initializeSlabVarCoeff() are used to initialize memory for volume flux and 
coefficients of variables like density viscosity thermal conductivity etc. This function also opens the files
variablelist.dat and coefficients.dat and writes all the user entered data like what variables to be evaluated and
what are there initialization values etc. The information stored in these files is read by the initializerestart
to allocate the same amount of memory and exactly the same initialization values during restart.
*/

int initialize(Base *CgnsBase)
{ 
	FILE *var,*var2,*var3,*convcri,*pressConvcri,*odt,*potentialMHDConvcri,*inducedMHDConvcri,*PsiMHDConvcri,*omega;//rajendra
	int i,j,k,p,flag,cntVariable,num_zone,count=0,vector,cntVar,scal,vortchoice,np;/* num_zone is the number of zones*/
	int numiface,numBc,iface,zBc,dir,varCoeffCount=0,Tchoice;
	Zone *z,*zNeighbour;
	Interface *IfacePtr;
	double constant,initViscosity,dt,time,initconductivity,initeffconductivity;
	double betastar=0.09,kappa=0.41,SigmaOmega1=0.5,Beta1=0.0750,SigmaOmega2=0.856,Beta2=0.0828;
	int conv,nop,cnt=0,nophase;
	int che;//checking purpose --- chandrasekhar.
	var=fopen("Settings/variablelist.dat","w");
	var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
	pressConvcri=fopen("Settings/pressureconvergCriteria.dat","w");
	if((CgnsBase->flagLMMHD==1)&&(CgnsBase->flaginduced==0)) //naveen
		potentialMHDConvcri=fopen("Settings/potentialMHDconvergCriteria.dat","w");
	if(CgnsBase->flaginduced==1) //rajendra
		{
			inducedMHDConvcri=fopen("Settings/inducedMHDconvergCriteria.dat","w");
			PsiMHDConvcri=fopen("Settings/PsiMHDconvergCriteria.dat","w");
			potentialMHDConvcri=fopen("Settings/potentialMHDconvergCriteria.dat","w");
		}
	convcri=fopen("Settings/convergenceCriteria.dat","r");
	if (convcri==NULL) {
	   printf("\nError opening file Settings/convergenceCriteria.dat to read....probably it does not exists....Trying to create....");
	   convcri=fopen("Settings/convergenceCriteria.dat","w");
	   if (convcri==NULL) { printf("Could not create....\n");
	   } else { printf("Created....\n"); }
	   for (i=0;i<=12;i++) { fprintf(convcri,"1e-8\n"); }
	   fprintf(convcri,"\n\n1. u,v,w,umass,vmass,wmass convergence criterion both in iterative loop and steady state convergence/All MPFV vairables\n");
	   fprintf(convcri,"2. Convergence criterion for Temperature/enthalpy,LMMHD variables,Flow_Potential,Droplet temperature for Droplet flows.\n");
	   fprintf(convcri,"3. From 3-12, the convergence criterion corresponds to 10 scalars.\n");
	   fclose(convcri);
	}
	omega=fopen("Settings/omega.dat","r");
	if (omega==NULL) {
	   printf("\nError opening file Settings/omega.dat to read....probably it does not exists....Trying to create....");
	   omega=fopen("Settings/omega.dat","w");
	   if (omega==NULL) { printf("Could not create....\n");
	   } else { printf("\nCreated....\n"); }
	   fprintf(omega,"0.6\n1.0\n1.0\n1.8\n\n0.7\n0.7\n0.7\n0.7\n0.7\n0.7\n0.7\n0.7\n0.7\n0.7\n\n\n");
	   fprintf(omega, "\n1. u,v,w,umass,vmass,wmass relaxation factors.  i.e., Momentum relaxation factors  \
			   \n2. Pressure relaxtion factor \
			   \n3. Energy relaxation factor for conv.diffusion equation. \
			   \n4. Relaxation factor for any diffusive equation without unsteady,convection terms. \
			   \n\
			   \nFirst four parameters are general and may not require modifications \
			   \n\
			   \nNext four parameters which corresponds to other variables apart from momentum and energy may need modification \
			   \naccording to the type of the problem and the variable solved. \
			   \n\
			   \n5. All other scalars/variables depending upon the module.  Observe the terminal while inputting the boundary conditions. \
  		           \nThe sequence of relaxation factors should be same as the sequence of boundary conditions.\n" );
	   fclose(omega);	
	}
	convcri=fopen("Settings/convergenceCriteria.dat","r");
	omega=fopen("Settings/omega.dat","r");
	if (convcri==NULL || pressConvcri==NULL || var==NULL || var2==NULL || omega==NULL)
	{
		printf("File opening if failed in APKinitialize.c\nExiting...........\n");
		exit(0);
	}
	
	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones;
	
	//Useful for conjugate heat transfer with turbulence
	//NOTE: This loop should be called before calling the function initialize_T().
		for(j=0;j<num_zone;j++)
	{
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		
		(z+j)->ConjugateZoneTag = 0;
		for(iface=0;iface<numiface;iface++)
		{
			IfacePtr = (z+j)->Iface+iface;
			zNeighbour = (Zone *) FindZonePtr(IfacePtr->Adjoining_Zone_Name, CgnsBase);
			IfacePtr->AdjoiningZonePtr = zNeighbour;//(Zone *) zNeighbour;
// 			printf("assigned=%u\tzNeighbour=%u\n", IfacePtr->AdjoiningZonePtr, zNeighbour);
			if((z+j)->flagSolidFluid - zNeighbour->flagSolidFluid !=0)//If two adjacent zones are of different phases
			{
				(z+j)->ConjugateZoneTag = 1;
					
				if((z+j)->flagSolidFluid  == 0)
					IfacePtr->ConjugateIfaceTag = 1;//Current Zone FLUID and Adjacent Zone SOLID
				else
					IfacePtr->ConjugateIfaceTag = 2;//Current Zone SOLID and Adjacent Zone FLUID
			}
			else
				IfacePtr->ConjugateIfaceTag = 0;//Current Zone does not have any Conjugate interfaces
		}
		}
	
/*#####################################################################################################################*/
	printf("\n########### ENTER THE FOLLOWING BASE PARAMETERS ############");
	
	//Bi-CgStab Initialization  //sandeep	
	//printf("\nBi-CgStab with Variable Density Case is under testing.\n");
	printf("\nEnter the Iterative Scheme ");
	printf("\nFor Gauss-Seidal(GS)-[Enter 0] ");
	printf("\nFor Bi-CGSTAB--------[Enter 1] ");
	printf("\nFor Optimised GS-----[Enter 2] : ");
	scanf("%d",&CgnsBase->BiCgflag);
	while(!(CgnsBase->BiCgflag==0 || CgnsBase->BiCgflag==1 || CgnsBase->BiCgflag==2))
	{
		printf("\nRefer the above options and enter the correct choice:");
		scanf("%d",&CgnsBase->BiCgflag);
	}
        fprintf(var2,"%d\n",CgnsBase->BiCgflag);        
        fprintf(var3,"%d \tEnter the Iterative Scheme 0->Gauss-Seidal, 1->Bi-CGSTAB, 2->Optimised GS: \n\n",CgnsBase->BiCgflag);
	
	if(CgnsBase->BiCgflag==2)//Impoved G-S method
	{
		for(j=0;j<num_zone;j++)
		{
			(z+j)->Residue= (double *)initialize_variable(1,(z+j),0.0);
		
			(z+j)->Correction= (double *)initialize_variable(1,(z+j),0.0);
			(z+j)->AP= (double *)initialize_variable(1,(z+j),0.0);
		
			(z+j)->Correction_Prev= (double *)initialize_variable(1,(z+j),0.0);
			(z+j)->AP_Prev= (double *)initialize_variable(1,(z+j),0.0);
			((z+j)->grid)->GRHS=(double *)initialize_variable(1,(z+j),0.0);
		}
		
	}
	else if(CgnsBase->BiCgflag==1)//Bi-CGStab
	{
		printf("Note: For BI-Cgstab Iterative scheme, enter stringent convergence criteria for Pressure, atleast as low as 1e-7.\n");
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Residue)= (double *)initialize_variable(1,(z+j),0.0);
			((z+j)->Residue_init)= (double *)initialize_variable(1,(z+j),0.0);	
			((z+j)->Pi)= (double *)initialize_variable(1,(z+j),0.0);
			((z+j)->Vi)= (double *)initialize_variable(1,(z+j),0.0);
			((z+j)->Si)= (double *)initialize_variable(1,(z+j),0.0);
			((z+j)->Ti)= (double *)initialize_variable(1,(z+j),0.0);
			((z+j)->grid)->GRHS=(double *)initialize_variable(1,(z+j),0.0);
		}
	}
	
	printf("\nEnter Convergence criteria  for PressureLoop : ");
	scanf("%lf",&constant);
	CgnsBase->BaseParameters->PressLoopConv = constant;
	fprintf(pressConvcri,"%12.8e\n",constant);
        fprintf(var2,"%12.8e\n",constant);        
        fprintf(var3,"%12.8e \tEnter Convergence criteria  for PressureLoop: \n\n",constant);
	
	if((CgnsBase->flagLMMHD==1)&&(CgnsBase->flaginduced==0)) //naveen
	{
		printf("\nEnter Convergence criteria for PotentialMHDLoop: ");
		scanf("%lf",&constant);
		CgnsBase->BaseParameters->PotentialMHDLoopConv = constant;
                fprintf(var2,"%12.8e\n",constant);        
                fprintf(var3,"%12.8e \tEnter Convergence criteria for PotentialMHDLoop: \n\n",constant);
                //fprintf(potentialMHDConvcri,"%12.8e\n",constant);
		fclose(potentialMHDConvcri);
	}
	if(CgnsBase->flaginduced==1) //rajendra
	{
		printf("\nEnter Convergence criteria for InducedMHDLoop: ");
		scanf("%lf",&constant);
                fprintf(var2,"%12.8e\n",constant);        
                fprintf(var3,"%12.8e \tEnter Convergence criteria for InducedMHDLoop: \n\n",constant);
		CgnsBase->BaseParameters->InducedMHDLoopConv = constant;
		fprintf(inducedMHDConvcri,"%12.8e\n",constant);
		fclose(inducedMHDConvcri);
		printf("\n Enter Convergence criteria for PsiMHDLoop: ");
		scanf("%lf",&constant);
                fprintf(var2,"%12.8e\n",constant);        
                fprintf(var3,"%12.8e \tEnter Convergence criteria for PsiMHDLoop: \n\n",constant);
		CgnsBase->BaseParameters->PsiMHDLoopConv = constant;
		fprintf(PsiMHDConvcri,"%12.8e\n",constant);
		fclose(PsiMHDConvcri);
		printf("\n Enter Convergence criteria for PotentialMHDLoop: ");
		scanf("%lf",&constant);                
                fprintf(var2,"%12.8e\n",constant);                
                fprintf(var3,"%12.8e \tEnter Convergence criteria for PotentialMHDLoop: \n\n",constant);
		CgnsBase->BaseParameters->PotentialMHDLoopConv = constant;
		fprintf(potentialMHDConvcri,"%12.8e\n",constant);
		fclose(potentialMHDConvcri);
	}
	printf("\nEnter Convergence criteria  for MassFlux Loop: ");
	scanf("%lf",&constant);	                
        fprintf(var2,"%12.8e\n",constant);        
        fprintf(var3,"%12.8e \tEnter Convergence criteria  for MassFlux Loop: \n\n",constant);
	CgnsBase->BaseParameters->MassFluxConv = constant;
	
	//printf("\nDo you want to enter the Time step[Enter 1] or let the Code fix it according to Stability Constraints[Enter 0]");
	//scanf("%d",&Tchoice);                
	Tchoice = 1;                
        //fprintf(var2,"%d\n",Tchoice);        
        //fprintf(var3,"%d \tDo you want to enter the Time step[Enter 1] or let the Code fix it according to Stability Constraints[Enter 0] \n\n",Tchoice);
	CgnsBase->BaseParameters->TstepChoice=Tchoice;
	
	if(Tchoice==1)
	{
		if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia==1)
		{
			printf("\nFor Dilute Gas-Droplet flows with Evaporation choose small time-steps\nFor inlet droplet diameter less than 100 microns choose time-step value as 1e-5\nFor inlet droplet diameter greater than 100 microns choose time step value as 1e-4");

		}

		printf("\nEnter TimeStep (dt):");
		scanf("%lf",&constant);                
                fprintf(var2,"%lf\n",constant);        
                fprintf(var3,"%lf \t Enter TimeStep (dt):\n\n",constant);
		CgnsBase->BaseParameters->TimeStep = constant;
	}
	else
	{
		CgnsBase->BaseParameters->TimeStep =0.001;
	}
	
	CgnsBase->BaseParameters->Turb->DualTimeStep = CgnsBase->BaseParameters->TimeStep / CgnsBase->BaseParameters->Turb->numDualTimeIter;

	//CgnsBase->BaseParameters=(Parameters *)malloc(sizeof(Parameters));
	CgnsBase->BaseParameters->flagGravity = 0;
	CgnsBase->BaseParameters->flagBuoyancy = 0;
	CgnsBase->BaseParameters->flagEnergySource = 0;
	CgnsBase->BaseParameters->flagDissipation = 0;
	CgnsBase->BaseParameters->TimeStepingSchemeTag=0;
/*#####################################################################################################################*/	
	printf("\nPlease Enter the time steping scheme ");
	printf("\nFor CrankNicolson scheme [Enter 1] : "); 
	printf("\nFor BTCS scheme          [Enter 0] : ");
	scanf("%d",&j);
        fprintf(var2,"%d\n",j);        
        fprintf(var3,"%d \tPlease Enter the time steping scheme [CrankNicolson->1/BTCS->0] \n\n",j);
	if(j==1)
		CgnsBase->BaseParameters->TimeStepingSchemeTag=1;

	printf("\nEnter the Convection Scheme(Deferred Corrected) QUICK [Enter=1] or Central Difference [Enter=0]: ");
	scanf("%d",&conv);
        fprintf(var2,"%d\n",conv);        
        fprintf(var3,"%d \tEnter the Convection Scheme(Deferred Corrected) QUICK [Enter=1] or Central Difference [Enter=0] \n\n",conv);
	CgnsBase->BaseParameters->ConvectionScheme=conv;
/*###########################################################################################*/	
	if(CgnsBase->BaseParameters->ConvectionScheme==0)
	{
		printf("\nEnter Ratio of CentralDiff(1->pure Centraldiff) to Upwind(0->pure Upwind):");
		scanf("%lf",&constant);
		fprintf(var2,"%lf\n",constant);        
                fprintf(var3,"%lf \tEnter Ratio of CentralDiff(1->pure Centraldiff) to Upwind(0->pure Upwind): \n\n",constant);
		CgnsBase->BaseParameters->RatioCenUP = constant;
	}
	else
	{
		printf("\nEnter Ratio of QUICK(1->pure QUICK) to Upwind(0->pure Upwind):");
		scanf("%lf",&constant);
		fprintf(var2,"%lf\n",constant);        
                fprintf(var3,"%lf \tEnter Ratio of QUICK(1->pure QUICK) to Upwind(0->pure Upwind): \n\n",constant);
		CgnsBase->BaseParameters->RatioCenUP = constant;
	}
	if(CgnsBase->flagTurb == 1)
	{
		if(CgnsBase->BaseParameters->ConvectionScheme==0)
		{
			printf("\nEnter Ratio of CentralDiff(1->pure CD2) to Upwind(0->pure Upwind) for TURBULENCE quantities:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);        
                        fprintf(var3,"%lf \tEnter Ratio of CentralDiff(1->pure CD2) to Upwind(0->pure Upwind) for TURBULENCE quantities: \n\n",constant);
			CgnsBase->BaseParameters->RatioCenUPTurb = constant;
		}
		else
		{
			printf("\nEnter Ratio of QUICK(1->pure QUICK) to Upwind(0->pure Upwind) for TURBULENCE quantities:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);        
                        fprintf(var3,"%lf \tEnter Ratio of QUICK(1->pure QUICK) to Upwind(0->pure Upwind) for TURBULENCE quantities: \n\n",constant);
			CgnsBase->BaseParameters->RatioCenUPTurb = constant;
		}
	}
	
	if(CgnsBase->flagConjugateProblem == 0)
	{
		printf("\nEnter Nature of Solution for conjugate problem (1->UnSteady State) (0->Steady State):");
		scanf("%lf",&constant);
		fprintf(var2,"%lf\n",constant);        
                fprintf(var3,"%lf \tEnter Nature of Solution for conjugate problem (1->UnSteady State) (0->Steady State): \n\n",constant);
		CgnsBase->flagSteadyState = constant;/*0->for steady state, 1->for Unsteady state*/
	}
	else
		CgnsBase->flagSteadyState = 0;/*0->for steady state, 1->for Unsteady state*/

	/*Allocation of memories for various Convergence Criteria*/
	CgnsBase->BaseParameters->convergenceCriterion=(double *)malloc(14*sizeof(double));
	for(cntVariable=0;cntVariable<14;cntVariable++)
		 fscanf(convcri,"%lf",(CgnsBase->BaseParameters->convergenceCriterion+cntVariable));

	CgnsBase->BaseParameters->OMEGA=(double *)malloc(14*sizeof(double));
	for(cntVariable=0;cntVariable<14;cntVariable++)
		 fscanf(omega,"%lf",(CgnsBase->BaseParameters->OMEGA+cntVariable));

	/*the first element is the criterion for velocities, second is for temperature and rest 10 are for the 10 scalars*/
	if(CgnsBase->BaseParameters->restart==1)
	{
		odt=fopen("Settings/old_dt.dat","r");
		fscanf(odt,"%lf %lf",&dt,&time);
		CgnsBase->BaseParameters->Time=time;
		fclose(odt);
		printf("The solution is restarted and time =%lf\n",time);
	}
	else if(CgnsBase->BaseParameters->restart==0)
	{
		CgnsBase->BaseParameters->Time=0.0;
	}
	else
	{
		printf("Wrong value in BaseParameters->restart ....... Exiting in initialize.c........\n");
		exit(0);
	}
	
/*#####################################################################################################################*/	
	for(j=0;j<num_zone;j++)
	{
		/* Allocation of memory for Flowsolution */
		(z+j)->Flowsol=(FlowSolution*)malloc(sizeof(FlowSolution));
		(z+j)->PrevFlowsol=(PrevFlowSolution*)malloc(sizeof(PrevFlowSolution));

                if(CgnsBase->flagTurbComb==1)
		(z+j)->Flowsol->flagTurbComb = 1; //Shyam
                else
		(z+j)->Flowsol->flagTurbComb = 0; //Shyam
	}
	
	/*Entering the value 1 for this question would enable the code 
	to allocate memory in the data structure for that particular variable  */
	if((CgnsBase->flagParticle==1)||(CgnsBase->flagDroplet==1))
		printf("\nNow you have to enter the property values for the carrier phase:\n");
	
	if(CgnsBase->flagDroplet==1)
		printf("\nIf you are considering Evaporation of droplet then choose density as variable:\n");
	
	if(CgnsBase->flagMPFV==0)
	{
	/*****************************MEMORY ALLOCATION FOR DENSITY*******************************************************/
		/*The memory allocation for Density for Flowsol node and Slabs at boundary and Interface*/
		printf("\nCreating memory for Density[kg/m3]:");
		if( CgnsBase->Number_of_zones > 1)
		{
			printf("\n Do you want to specify Density for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1):");
			scanf("%d",&flag);
			fprintf(var2,"%d\n",flag);        
                	fprintf(var3,"%d \tDo you want to specify Density for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1): \n\n",flag);
		} else { flag = 0;}
		if(flag==0) /*Initializing density for whole domain at a time*/
		{
			//printf("\n If Density is a CONSTANT, Enter 0");
			//printf("\n If Density is a VARIABLE, Enter 1:");
			//scanf("%d",&p);
			p = 0;
			//fprintf(var2,"%d\n",p);        
                        //fprintf(var3,"%d \tIf Density is a CONSTANT, Enter 0, VARIABLE Enter 1: \n\n",p);
			varCoeffCount = 1;  /*1->Density*/
			if(p==0)
			{
				printf("\nEnter the constant value of Density[kg/m3]:");
				scanf("%lf",&constant);
				fprintf(var2,"%12.8e\n",constant);        
                                fprintf(var3,"%12.8e \tEnter the constant value of Density[kg/m3]: \n\n",constant);
			}
			else
			{
				//adhi+shash
				if((CgnsBase->flagDroplet==1)&&(CgnsBase->BaseParameters->flagVariableDia==1))
				{
					printf("\nEnter the Initialization value of Density[kg/m3] for Carrier phase:");
					scanf("%lf",&constant);                
                                        fprintf(var2,"%12.8e\n",constant);        
                                        fprintf(var3,"%12.8e \tEnter the Initialization value of Density[kg/m3] for Carrier phase: \n\n",constant);
				}
				else
					constant=1.0;//adhi+shash
			}	
			for(j=0;j<num_zone;j++)
			{
				/*Memory allocation in the Flow Solution Node X Direction Velocity*/
				((z+j)->Flowsol)->Rho= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->PrevFlowsol)->Rho= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->Flowsol)->flagRho=p;
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							/*Memory allocation in the Interface Slab X Direction Velocity*/
							//printf("I am going to Iface\n");
							(((z+j)->Iface+iface)->ZSlab)->Rho=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}
	
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to zoneBc\n");
						(((z+j)->ZoneBc+zBc)->ZSlab)->Rho = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}
		}
		if(flag==1)
		{
			for(j=0;j<num_zone;j++)
			{ 
				printf("Zone Name is %s\n",(z+j)->zonename);
				printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				//printf("\n If Density is a CONSTANT, Enter 0");
				//printf("\n If Density is a VARIABLE, Enter 1:");
				//scanf("%d",&p);
				p = 0;
				//fprintf(var2,"%d\n",p);        
                                //fprintf(var3,"%d \tIf Density is a CONSTANT, Enter 0: If Density is a VARIABLE, Enter 1: \n\n",p);
				varCoeffCount = 1;  /*1->Density*/
				if(p==0)
				{
					printf("\nEnter the constant value of Density[kg/m3]:");
					scanf("%lf",&constant);
					fprintf(var2,"%12.8e\n",constant);        
                                        fprintf(var3,"%12.8e \tEnter the constant value of Density[kg/m3]: \n\n",constant);
				}
				else
					constant=1.0;
				
				((z+j)->Flowsol)->Rho= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->PrevFlowsol)->Rho= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->Flowsol)->flagRho=p;
				numiface=(z+j)->Number_of_Interfaces;
				//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							/*Memory allocation in the Interface Slab X Direction Velocity*/
							//printf("I am going to Iface\n");
							(((z+j)->Iface+iface)->ZSlab)->Rho=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}
	
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to zoneBc\n");
						(((z+j)->ZoneBc+zBc)->ZSlab)->Rho = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}
		}
	/***************************** END OF MEMORY ALLOCATION FOR DENSITY*******************************************************/
	
	/*****************************MEMORY ALLOCATION FOR VISCOSITY*******************************************************/	
		printf("\nCreating memory for VISCOSITY[kg/m-s]:");
		if( CgnsBase->Number_of_zones > 1)
		{
			printf("\n Do you want to specify Viscosity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
			scanf("%d",&flag);
			fprintf(var2,"%d\n",flag);        
                	fprintf(var3,"%d \tDo you want to specify Viscosity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",flag);
		} else { flag = 0;}
		if(flag==0) /*Initializing density for whole domain at a time*/
		{
			//printf("\n If Viscosity is a CONSTANT, Enter 0");
			//printf("\n If Viscosity is a VARIABLE, Enter 1:");
			//scanf("%d",&p);
			p = 0;
			//fprintf(var2,"%d\n",p);        
                        //fprintf(var3,"%d \tIf Viscosity is a CONSTANT, Enter 0: If Viscosity is a VARIABLE, Enter 1: \n\n",p);
			varCoeffCount = 2;  /*2->VISCOSITY*/
			if(p==0)
			{
				printf("\nEnter the constant value of VISCOSITY[kg/m-s]:");
				scanf("%lf",&constant);
				fprintf(var2,"%E\n",constant);        
                                fprintf(var3,"%E \tEnter the constant value of VISCOSITY[kg/m-s]: \n\n",constant);
			}
			else
				constant=1.0;
	
			for(j=0;j<num_zone;j++)
			{       
				/*Memory allocation in the Flow Solution Node X Direction Velocity*/
				((z+j)->Flowsol)->Viscoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->PrevFlowsol)->Viscoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->Flowsol)->flagViscoeff=p;
				numiface=(z+j)->Number_of_Interfaces;
				//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							/*Memory allocation in the Interface Slab X Direction Velocity*/
							//printf("I am going to Iface\n");
							(((z+j)->Iface+iface)->ZSlab)->Viscoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}
	
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to zoneBc\n");
						(((z+j)->ZoneBc+zBc)->ZSlab)->Viscoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}
		}
		if(flag==1)
		{
			for(j=0;j<num_zone;j++)
			{
				printf("Zone Name is %s\n",(z+j)->zonename);
				printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				//printf("\n If Viscosity is a CONSTANT, Enter 0");
				//printf("\n If Viscosity is a VARIABLE, Enter 1:");
				//scanf("%d",&p);
				p = 0;
				//fprintf(var2,"%d\n",p);        
                                //fprintf(var3,"%d \tIf Viscosity is a CONSTANT Enter 0: If Viscosity is a VARIABLE Enter 1 \n\n",p);
				varCoeffCount = 2;  /*2->VISCOSITY*/
				if(p==0)
				{
					printf("\nEnter the constant value of VISCOSITY[kg/m-s]:");
					scanf("%lf",&constant);
					fprintf(var2,"%12.8e\n",constant);        
                                        fprintf(var3,"%12.8e \tEnter the constant value of VISCOSITY[kg/m-s]: \n\n",constant);
				}
				else
					constant=1.0;
				
				((z+j)->Flowsol)->Viscoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->PrevFlowsol)->Viscoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->Flowsol)->flagViscoeff=p;
				numiface=(z+j)->Number_of_Interfaces;
				//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							/*Memory allocation in the Interface Slab X Direction Velocity*/
							//printf("I am going to Iface\n");
							(((z+j)->Iface+iface)->ZSlab)->Viscoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}
	
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to zoneBc\n");
						(((z+j)->ZoneBc+zBc)->ZSlab)->Viscoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}	
		}
		
	/*****************************END OF MEMORY ALLOCATION FOR VISCOSITY*******************************************************/
        //printf("\nEnter u_source:");
	//scanf("%lf",&constant);                
	constant = 0;                
        //fprintf(var2,"%lf\n",constant);        
        //fprintf(var3,"%lf \tEnter u_source: \n\n",constant);
	CgnsBase->BaseParameters->u_source = constant;

	//printf("\nEnter v_source:");
	//scanf("%lf",&constant);                
	constant = 0;                
        //fprintf(var2,"%lf\n",constant);        
        //fprintf(var3,"%lf \tEnter v_source: \n\n",constant);
	CgnsBase->BaseParameters->v_source = constant;


	for(j=0;j<num_zone;j++)
		{       
			//printf("\nInitialising ConjugateRms to zero\n");
			
			numiface=(z+j)->Number_of_Interfaces;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						((z+j)->Iface+iface)->ConjugateRms = 0.0;
					}
				}
			}
		}
	/*****************************MEMORY ALLOCATION FOR U,V,W VELOCITIES*******************************************************/
		/*The memory allocation for X Direction Velocity for Flowsol node and Slabs at boundary and Interface */
		
		
		
		CgnsBase->VARTAG[count]=1;//Variable Tag has to be stored bhans
		count=count+1;

		printf("\nCreating memory for U velocity and umass[m/s]");
		i=1;
		fprintf(var,"%d\n",i);

		if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {

		   if( CgnsBase->Number_of_zones > 1){
			printf("\nDo you want to initialize U velocity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);        
                        fprintf(var3,"%d \tDo you want to initialize U velocity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",p);
		   } 
		   if(p==0){
			printf("\nEnter the Initialization value of U velocity[m/s] for all Zones:");
			scanf("%lf",&constant);
			fprintf(var2,"%12.8e\n",constant);        
                	fprintf(var3,"%12.8e \tEnter the Initialization value of U velocity[m/s] for all Zones: \n\n",constant);
		   }
		}
		for(j=0;j<num_zone;j++)
		{
		   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
			if(p!=0){
				printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				printf("\nEnter the Initialization value of U velocity[m/s]:");
				scanf("%lf",&constant);
				fprintf(var2,"%12.8e\n",constant);        
                        	fprintf(var3,"%12.8e \tEnter the Initialization value of U velocity[m/s]: \n\n",constant);
			}
		   } else { constant = 0.0; }
			/*Memory allocation in the Flow Solution Node X Direction Velocity*/
			((z+j)->Flowsol)->u= (double *)initialize_variable(i,(z+j),constant);
                        ((z+j)->Flowsol)->utemp= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->u= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->umass= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->umass= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->u=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->uOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->umass = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->umassOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->u = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->umass = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->uOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->umassOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					}
				}
			}
		}

		
		/*The memory allocation for Y Direction Velocity for Flowsol node and Slabs at boundary and Interface */
		printf("\nCreating memory for V velocity and vmass[m/s]");
		i=1;
		fprintf(var,"%d\n",i);
		CgnsBase->VARTAG[count]=2;//Variable Tag has to be stored bhans
		count=count+1;
		
		if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
		   if( CgnsBase->Number_of_zones > 1){
			printf("\nDo you want to initialize V velocity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);        
                        fprintf(var3,"%d \tDo you want to initialize V velocity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",p);
		   } 
		   if(p==0){
			printf("\nEnter the Initialization value of V velocity[m/s] for all Zones:");
			scanf("%lf",&constant);
			fprintf(var2,"%12.8e\n",constant);                        
                	fprintf(var3,"%12.8e \tEnter the Initialization value of V velocity[m/s] for all Zones: \n\n",constant);
		   }
		}
		for(j=0;j<num_zone;j++)
		{
		   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
			if(p!=0){
				printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				printf("\nEnter the Initialization value of V velocity[m/s]:");
				scanf("%lf",&constant);
				fprintf(var2,"%12.8e\n",constant);                        
                        	fprintf(var3,"%12.8e \tEnter the Initialization value of V velocity[m/s]: \n\n",constant);
			}
		   } else { constant = 0.0;}
		
                        ((z+j)->Flowsol)->vtemp= (double *)initialize_variable(i,(z+j),constant);
	                ((z+j)->Flowsol)->v= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->v= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->vmass= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->vmass= (double *)initialize_variable(i,(z+j),constant); 
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->v = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->vOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->vmass = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->vmassOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
					}
				}
				
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->v = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->vmass = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->vOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->vmassOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					}
				}
			}
		}
	
		//The memory allocation for Z Direction Velocity for Flowsol node and Slabs at boundary and Interface 
		if(CgnsBase->BaseParameters->flag2d == 0) // means 3D case then print Note: In 2d case memory is NOT created for w Velocity
		{
			printf("\nCreating memory for W velocity and Wmass[m/s]");
			i=1;
			CgnsBase->VARTAG[count]=3;//Variable Tag has to be stored bhans
			count=count+1;
		}
		else
		{
			i=0;
		}
		fprintf(var,"%d\n",i);
		
		if(CgnsBase->BaseParameters->flag2d == 0)
		{
		   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
			if( CgnsBase->Number_of_zones > 1){
				printf("\nDo you want to initialize W velocity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
				scanf("%d",&p);
				fprintf(var2,"%d\n",p);        
                        	fprintf(var3,"%d \tDo you want to initialize W velocity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",p);
			} 
			if(p==0){
				printf("\nEnter the Initialization value of W velocity[m/s] for all Zones:");
				scanf("%lf",&constant);
				fprintf(var2,"%12.8e\n",constant);                        
                        	fprintf(var3,"%12.8e \tEnter the Initialization value of W velocity[m/s] for all Zones: \n\n",constant);
			}
		   }
		}
		else
			constant = 0.0;
			
		for(j=0;j<num_zone;j++)
		{
		   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
				if(p!=0){
					printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
					printf("\nEnter the Initialization value of W velocity[m/s]:");
					scanf("%lf",&constant);
					fprintf(var2,"%12.8e\n",constant);                        
                                	fprintf(var3,"%12.8e \tEnter the Initialization value of W velocity[m/s]: \n\n",constant);
				}
			}
			else
				constant = 0.0; 
		   } else {constant = 0.0;}
			
			//Memory allocation in the Flow Solution Node Z Direction Velocity
			((z+j)->Flowsol)->w= (double *)initialize_variable(i,(z+j),constant);
                        ((z+j)->Flowsol)->wtemp= (double *)initialize_variable(i,(z+j),constant);		
	                ((z+j)->PrevFlowsol)->w= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->wmass= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->wmass= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{	
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->w = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->wOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->wmass = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->wmassOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
					}
				}
	
				for(zBc=0;zBc<numBc;zBc++)
				{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->w = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->wmass = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->wOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->wmassOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					}
				}
			}
		}
	}
	else//Added by sushil 
	{	
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->u= NULL;
			((z+j)->Flowsol)->v= NULL;
			((z+j)->Flowsol)->w= NULL;
			((z+j)->PrevFlowsol)->u= NULL;
			((z+j)->PrevFlowsol)->v= NULL;
			((z+j)->PrevFlowsol)->w= NULL;
			((z+j)->Flowsol)->umass= NULL;
			((z+j)->Flowsol)->vmass= NULL;
			((z+j)->Flowsol)->wmass= NULL;
			((z+j)->PrevFlowsol)->umass= NULL;
			((z+j)->PrevFlowsol)->vmass= NULL;
			((z+j)->PrevFlowsol)->wmass= NULL;
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->u=NULL;
						(((z+j)->Iface+iface)->ZSlab)->v=NULL;
						(((z+j)->Iface+iface)->ZSlab)->w=NULL;
						(((z+j)->Iface+iface)->ZSlab)->uOld=NULL;
						(((z+j)->Iface+iface)->ZSlab)->vOld=NULL;
						(((z+j)->Iface+iface)->ZSlab)->wOld=NULL;
						(((z+j)->Iface+iface)->ZSlab)->umass = NULL;
						(((z+j)->Iface+iface)->ZSlab)->vmass = NULL;
						(((z+j)->Iface+iface)->ZSlab)->wmass = NULL;
						(((z+j)->Iface+iface)->ZSlab)->umassOld = NULL;
						(((z+j)->Iface+iface)->ZSlab)->vmassOld = NULL;
						(((z+j)->Iface+iface)->ZSlab)->wmassOld = NULL;
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->u = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->v = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->w = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->umass = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->vmass = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->wmass = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->uOld = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->vOld = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->wOld = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->umassOld = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->vmassOld = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->wmassOld = NULL;
					}
				}
			}
		}	
	}//end
/*****************************END OF MEMORY ALLOCATION FOR U,V,W VELOCITIES*******************************************************/
	
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->prevIterVariable = (double*)malloc((*((z+j)->Size+0)-1) * (*((z+j)->Size+1)-1) * (*((z+j)->Size+2)-1)*sizeof(double));
	}
	
/*****************************MEMORY ALLOCATION FOR MASSFLUXES*******************************************************/
/* Creatiing memory for massflux is necessary,hence no need to ask user to enter 0 or 1 for massflux*/
/* The function initializeFlux is defined in the file initialize_variable.c*/
	printf("\n Creating memory for MassFlux");
	i=1;
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->Massflux= (double *)initializeFlux(i,(z+j));
			((z+j)->PrevFlowsol)->Massflux= (double *)initializeFlux(i,(z+j)); 
		}
	//printf("\n Initializing default value of flagTurb to 0 (Laminar case)\n");
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->flagTurb=0;
	}
/*****************************END OF MEMORY ALLOCATION FOR MASSFLUXES*******************************************************/

/*****************************MEMORY ALLOCATION FOR PRESSURE*******************************************************/
	/*The memory allocation for Pressure for Flowsol node and Slabs at boundary and Interface */
	printf("\n Creating memory for Pressure");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[count]=4;//Variable Tag has to be stored bhans
	count=count+1;
		for(j=0;j<num_zone;j++)
		{
			//printf("\nEnter the Initialization value of Pressure:");
			//scanf("%lf",&constant);
			//fprintf(var2,"%12.8e\n",constant);        
                        //fprintf(var3,"%12.8e \tEnter the Initialization value of Pressure: \n\n",constant);
			constant = 0.0;
			/*Memory allocation in the Flow Solution Node For Pressure*/
			((z+j)->Flowsol)->P= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->P= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->PressureCoeff= (double *)malloc(sizeof(double));
			*(((z+j)->Flowsol)->PressureCoeff)=1.0;
			
			///To be deleted.  Added for Testing Purpose.
			((z+j)->Flowsol)->Pressure_Grad[0]= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->Pressure_Grad[1]= (double *)initialize_variable(i,(z+j),constant);
			
			if(CgnsBase->BaseParameters->flag2d==0)
				((z+j)->Flowsol)->Pressure_Grad[2]= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->Source_term= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->Diag_Source= (double *)initialize_variable(i,(z+j),constant);
			///Block ends here
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						//Memory allocation in the Interface Slab for Pressure*/
						(((z+j)->Iface+iface)->ZSlab)->P = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						//Memory allocation in the Boundary Slab for Pressure*/
						(((z+j)->ZoneBc+zBc)->ZSlab)->P = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					}
				}
			}
		}
/*****************************END OF MEMORY ALLOCATION FOR PRESSURE*******************************************************/
/****************************Initialiation for Potential flow solution****************************************************/
	printf("\n\nDo you want to solve for Potential flow before entering into solver [1->Yes|0-> No]:");
	scanf("%d",&j);                
        fprintf(var2,"%d\n",j);        
        fprintf(var3,"%d \tDo you want to solve for Potential flow before entering into solver [1->Yes|0-> No]: \n\n",j);

	if(j==1)
	{
		printf("\n*****NOTE: Potential flow solver WORKS ONLY incase of a grid with ORTHOGONAL BOUNDARIES at Inlet and Outlet******\n");
		
		printf("Please Note that the following boundary conditions have to be entered for potential flow problems\n\n");
		printf("Velocity specified boundaries..... Nuemann for Potential function \n");
		printf("for ex. a stationary wall has to be given Homogeneous Neumann \nand inlet boundary condition for");
		printf("potential should be given non-homogeneous neumann (with a value equal to the inlet velocity\n");
		printf("\nOutlet boundary condition for Flow_Potential should be homogeneous direchlet\nand all other boundaries can be specified as homogeneous Neumann boundary conditions");
	}

	CgnsBase->BaseParameters->flag_Potentialflow=j;
	i=j;
	if(j==1)
	{
		CgnsBase->VARTAG[count]=990;//Variable Tag has to be stored bhans
		count=count+1;
	}
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		/*Memory allocation in the Flow Solution Node For Potential*/
		((z+j)->Flowsol)->Flow_Potential= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->Flow_Potential= (double *)initialize_variable(i,(z+j),constant);

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
					//Memory allocation in the Interface Slab for Pressure*/
					(((z+j)->Iface+iface)->ZSlab)->Flow_Potential = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					//Memory allocation in the Boundary Slab for Pressure*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->Flow_Potential = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
				}
			}
		}
	}


	
	if(CgnsBase->flagDiffKinetCont==1) //Shyam, Boussinesq approximation case for diffusion controlled combustion problem
        {
                printf("\nDo you want to incorporate Buoyancy in Momentum Eqn.(Boussinesq Approximation)[1->Yes/0->No]:");
                scanf("%d",&j);
                fprintf(var,"%d\n",j);                
                fprintf(var2,"%d\n",j);        
                fprintf(var3,"%d \tDo you want to incorporate Buoyancy in Momentum Eqn.(Boussinesq Approximation)[1->Yes/0->No]: \n\n",j);
                if(j==1)
                {
                        CgnsBase->BaseParameters->flagBuoyancy = 1;
                        printf("\nEnter the value of Beta (Volumetric Expansion Coefficient[1/k])");
                        scanf("%lf",&constant);
                        fprintf(var,"%lf\n",constant);                        
                        fprintf(var2,"%lf\n",constant);        
                        fprintf(var3,"%lf \tEnter the value of Beta (Volumetric Expansion Coefficient[1/k]) \n\n",constant);
                        CgnsBase->BaseParameters->VolumetricExpCoeff = constant;

                        printf("\nEnter the reference Temperature (T reference):");
                        scanf("%lf",&constant);
                        fprintf(var,"%lf\n",constant);                        
                        fprintf(var2,"%lf\n",constant);        
                        fprintf(var3,"%lf \tEnter the reference Temperature (T reference): \n\n",constant);
                        CgnsBase->BaseParameters->TReference = constant;
                }
        }

/*****************************MEMORY ALLOCATION FOR TEMPERATURE*******************************************************/	
	/*The memory allocation for Temperature in the Flowsol node and Slabs at boundary and Interface */	
	//bhans START
        //printf("\nFor diffusion controlled combustion problem, please select '0'\n");
	printf("\nDO YOU WANT TO SOLVE ENERGY EQUATION (1 ->YES  0-> NO): ");
	scanf("%d",&i);
	fprintf(var,"%d\n",i);                        
        fprintf(var2,"%d\n",i);        
        fprintf(var3,"%d \tDO YOU WANT TO SOLVE ENERGY EQUATION (1 ->YES  0-> NO): \n\n",i);
	//bhans END        
	if(i==1)
	{
		count=count+1;
                if(CgnsBase->flagTurbComb==0) //Shyam
                {
			printf("\nCreating memory for Thermal Conductivity[w/m-k]:");
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify Thermal Conductivity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
				scanf("%d",&flag);
				fprintf(var2,"%d\n",flag);        
                        	fprintf(var3,"%d \tDo you want to specify Thermal Conductivity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",flag);
			} else { flag = 0;}
			if(flag==0) 
			{
				//printf("\n If Thermal Conductivity is a CONSTANT, Enter 0");
				//printf("\n If Thermal Conductivity is a VARIABLE, Enter 1:");
				//scanf("%d",&p);
				p = 0;
				//fprintf(var2,"%d\n",p);        
                                //fprintf(var3,"%d \tIf Thermal Conductivity is a CONSTANT Enter 0: If Thermal Conductivity is a VARIABLE Enter 1: \n\n",p);
				varCoeffCount = 3;  /*3->THERMAL CONDUCTIVITY*/
				if(p==0)
				{
					printf("\nEnter the constant value of THERMAL CONDUCTIVITY[w/m-k]:");
					scanf("%lf",&constant);
					fprintf(var2,"%12.8e\n",constant);        
                                        fprintf(var3,"%12.8e \tEnter the constant value of THERMAL CONDUCTIVITY[w/m-k]: \n\n",constant);
				}
				else
					constant=1.0;

				for(j=0;j<num_zone;j++)
				{       
					((z+j)->Flowsol)->Diffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->Diffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagDiffcoeff=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->Diffcoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
								//printf("I am leaving Iface\n");
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->Diffcoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
			if(flag==1)
			{
				for(j=0;j<num_zone;j++)
				{ 
					printf("Zone Name is %s\n",(z+j)->zonename);
					printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
					//printf("\n If Thermal Conductivity is a CONSTANT, Enter 0");
					//printf("\n If Thermal Conductivity is a VARIABLE, Enter 1:");
					//scanf("%d",&p);
					p = 0;
					//fprintf(var2,"%d\n",p);        
                                        //fprintf(var3,"%d \tIf Thermal Conductivity is a CONSTANT Enter 0: If Thermal Conductivity is a VARIABLE Enter 1:\n\n",p);
					varCoeffCount = 3;  /*3->THERMAL CONDUCTIVITY*/
					if(p==0)
					{
						printf("\nEnter the constant value of THERMAL CONDUCTIVITY[w/m-k]:");
						scanf("%lf",&constant);
						fprintf(var2,"%12.8e\n",constant);        
                                                fprintf(var3,"%12.8e \tEnter the constant value of THERMAL CONDUCTIVITY[w/m-k]: \n\n",constant);
					}
					else
						constant=1.0;

					((z+j)->Flowsol)->Diffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->Diffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagDiffcoeff=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->Diffcoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->Diffcoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
/**************************************END OF MEMORY FOR THERMAL CONDUCTIVITY*****************************************************************/

/**************************************CREATING MEMORY FOR SPECIFIC HEAT********************************************************************/
			printf("\nCreating memory for Specific Heat[J/kg-k]:");
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify Specific Heat for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1):");
				scanf("%d",&flag);
				fprintf(var2,"%d\n",flag);        
                        	fprintf(var3,"%d \tDo you want to specify Specific Heat for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1): \n\n",flag);
			} else { flag = 0;}
			if(flag==0) /*Initializing density for whole domain at a time*/
			{	//printf("\n If Specific Heat is a CONSTANT ,Enter 0");
				//printf("\n If Specific Heat is a VARIABLE,Enter 1:");
				//scanf("%d",&p);
				p = 0;
				//fprintf(var2,"%d\n",p);        
                                //fprintf(var3,"%d \tIf Specific Heat is a CONSTANT Enter 0: If Specific Heat is a VARIABLE Enter 1: \n\n",p);
				varCoeffCount = 4;  /*4->SPECIFIC HEAT*/
				if(p==0)
				{
					printf("\nEnter the constant value of SPECIFIC HEAT[J/kg-k]:");
					scanf("%lf",&constant);
					fprintf(var2,"%12.8e\n",constant);        
                                        fprintf(var3,"%12.8e \tEnter the constant value of SPECIFIC HEAT[J/kg-k]: \n\n",constant);
				}
				else
					constant=1.0;

				for(j=0;j<num_zone;j++)
				{       
					/*Memory allocation in the Flow Solution Node X Direction Velocity*/
					((z+j)->Flowsol)->SpecificHeat= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->SpecificHeat= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagSpecificHeat=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{
							if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->SpecificHeat=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
								//printf("I am leaving Iface\n");
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->SpecificHeat = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
			if(flag==1)
			{
				for(j=0;j<num_zone;j++)
				{ 
					printf("Zone Name is %s\n",(z+j)->zonename);
					printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
					//printf("\n If Specific Heat is a CONSTANT, Enter 0");
					//printf("\n If Specific Heat is a VARIABLE, Enter 1:");
					//scanf("%d",&p);
					p = 0;
					//fprintf(var2,"%d\n",p);        
                                        //fprintf(var3,"%d \tIf Specific Heat is a CONSTANT Enter 0: If Specific Heat is a VARIABLE Enter 1: \n\n",p);
					varCoeffCount = 4;  /*4->SPECIFIC HEAT*/
					if(p==0)
					{
						printf("\nEnter the constant value of SPECIFIC HEAT[J/kg-k]:");
						scanf("%lf",&constant);
						fprintf(var2,"%12.8e\n",constant);        
                                                fprintf(var3,"%12.8e \tEnter the constant value of SPECIFIC HEAT[J/kg-k]: \n\n",constant);
					}
					else
						constant=1.0;

					((z+j)->Flowsol)->SpecificHeat= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->SpecificHeat= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagSpecificHeat=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{
							if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->SpecificHeat=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
								//printf("I am leaving Iface\n");
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->SpecificHeat = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
		}//Shyam

/************************************Gravity********************************************************************************************/
	printf("\nDo you want to incorporate Gravity in Momentum.\n");
	printf("Use gravity for both Boussinesq and non-Boussinesq approx. for including Bouyancy Eqn.[1->Yes/0->No]:");
	scanf("%d",&i);
	fprintf(var,"%d\n",i);                        
        fprintf(var2,"%d\n",i);        
        fprintf(var3,"%d \tUse gravity for both Boussinesq and non-Boussinesq approx. for including Bouyancy Eqn.[1->Yes/0->No]: \n\n",i);
	if(i==1)
	{
		CgnsBase->BaseParameters->flagGravity = 1;
		printf("\nEnter the X-Component of Gravity:");
		scanf("%lf",&(CgnsBase->BaseParameters->gravityX));
		fprintf(var2,"%12.8e\n",(CgnsBase->BaseParameters->gravityX));        
                fprintf(var3,"%12.8e \tEnter the X-Component of Gravity: \n\n",CgnsBase->BaseParameters->gravityX);
		printf("\nEnter the Y-Component of Gravity:");
		scanf("%lf",&(CgnsBase->BaseParameters->gravityY));
		fprintf(var2,"%12.8e\n",(CgnsBase->BaseParameters->gravityY));                
                fprintf(var3,"%12.8e \tEnter the Y-Component of Gravity: \n\n",CgnsBase->BaseParameters->gravityY);
		printf("\nEnter the Z-Component of Gravity:");
		scanf("%lf",&(CgnsBase->BaseParameters->gravityZ));
		fprintf(var2,"%12.8e\n",(CgnsBase->BaseParameters->gravityZ));                
                fprintf(var3,"%12.8e \tEnter the Z-Component of Gravity: \n\n",CgnsBase->BaseParameters->gravityZ);
	}

/**************************************END OF MEMORY FOR SPECIFIC HEAT********************************************************************/
		printf("\nDo you want to incorporate Buoyancy in Momentum Eqn.(Boussinesq Approximation)[1->Yes/0->No]:");
		scanf("%d",&j);
		fprintf(var,"%d\n",j);                        
                fprintf(var2,"%d\n",j);        
                fprintf(var3,"%d \tDo you want to incorporate Buoyancy in Momentum Eqn.(Boussinesq Approximation)[1->Yes/0->No]: \n\n",j);
		if(j==1)
		{
			CgnsBase->BaseParameters->flagBuoyancy = 1;
			printf("\nEnter the value of Beta (Volumetric Expansion Coefficient[1/k])");
			scanf("%lf",&constant);
			fprintf(var,"%lf\n",constant);                        
                        fprintf(var2,"%lf\n",constant);        
                        fprintf(var3,"%lf \tEnter the value of Beta (Volumetric Expansion Coefficient[1/k]) \n\n",constant);
			CgnsBase->BaseParameters->VolumetricExpCoeff = constant;

			printf("\nEnter the reference Temperature (T reference):");
			scanf("%lf",&constant);
			fprintf(var,"%lf\n",constant);                        
                        fprintf(var2,"%lf\n",constant);                        
                        fprintf(var3,"%lf \tEnter the reference Temperature (T reference): \n\n",constant);
			CgnsBase->BaseParameters->TReference = constant;
		}
		printf("\nDo you want to incorporate Volumetric Heat Source in Energy Eqn.[1->Yes/0->No]:");
		scanf("%d",&j);
		fprintf(var,"%d\n",j);                        
                fprintf(var2,"%d\n",j);        
                fprintf(var3,"%d \tDo you want to incorporate Volumetric Heat Source in Energy Eqn.[1->Yes/0->No]: \n\n",j);
		if(j==1)
		{
			CgnsBase->BaseParameters->flagEnergySource = 1;
			for(j=0;j<num_zone;j++)
			{
				printf("Zone Name=%s\tSolid/Fluid Name=%s\n",(z+j)->zonename,(z+j)->SolidFluidName);
				printf("Enter 0/1 for Volumetric Heat Generation in this Zone 1->Yes 0->No :\n");
				scanf("%d",&p);
				fprintf(var,"%d\n",p);                        
                                fprintf(var2,"%d\n",p);        
                                fprintf(var3,"%d \tEnter 0/1 for Volumetric Heat Generation in this Zone 1->Yes 0->No : \n\n",p);
				if(p==0)
					(z+j)->flagHeatGeneration = 0;
				if(p==1)
					(z+j)->flagHeatGeneration = 1;
				if((z+j)->flagHeatGeneration==1)
				{
					printf("Enter whether source is constant[0] at each point or variable[1]: 0->Constant 1->Variable: \n");
					scanf("%d",&p);
					fprintf(var,"%d\n",p);                        
                                        fprintf(var2,"%d\n",p);        
                                        fprintf(var3,"%d \tEnter whether source is constant[0] at each point or variable[1]: 0->Constant 1->Variable: \n\n",p);
					(z+j)->flagHeat_con_var=p; 
					if (p==0)
					{	
						printf("Enter the constant value: ");
						scanf("%lf",&constant);
						fprintf(var,"%12.8e\n",constant);                        
                                                fprintf(var2,"%12.8e\n",constant);        
                                                fprintf(var3,"%12.8e \tEnter the constant value: \n\n",constant);
					}
					else
						constant=1.0;
					((z+j)->Flowsol)->Energysource= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				}
				else
				{
					p=0;
					(z+j)->flagHeat_con_var=0;
					constant=0.0;
					((z+j)->Flowsol)->Energysource= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				}
			}
		}
		else
		{
			CgnsBase->BaseParameters->flagEnergySource = 0;
			for(j=0;j<num_zone;j++)
			{
				(z+j)->flagHeatGeneration = 0;
				(z+j)->flagHeat_con_var   = 0;
			}
		}
		/*for(j=0;j<num_zone;j++)
		{
			printf("Zone Name=%s\tSolid/Fluid Name=%s\tflagHeatGeneration = %d\n",(z+j)->zonename,(z+j)->SolidFluidName,(z+j)->flagHeatGeneration);
		}*/

/* Heat Dissipation is not implemented yet  - So commenting it out 
		printf("\nDo you want to incorporate Dissipation in Energy Eqn.[1->Yes/0->No]:");
		scanf("%d",&j);
		fprintf(var,"%d\n",j);                        
                fprintf(var2,"%d\n",j);        
                fprintf(var3,"%d \tDo you want to incorporate Dissipation in Energy Eqn.[1->Yes/0->No]: \n\n",j);
*/

		j = 0;
		if(j==1)
			CgnsBase->BaseParameters->flagDissipation = 1;
		else
			CgnsBase->BaseParameters->flagDissipation = 0;
		if(CgnsBase->flagLMMHD==1){
		j = 0;
/* Dissipation was not implemented - so commenting out
                printf("\nDo you want to incorporate Joule Dissipation in Energy Eqn.[1->Yes/0->No]:");
		scanf("%d",&j);
		fprintf(var,"%d\n",j);                        
                fprintf(var2,"%d\n",j);        
                fprintf(var3,"%d \tDo you want to incorporate Joule Dissipation in Energy Eqn.[1->Yes/0->No]: \n\n",j);
*/
		if(j==1)
			CgnsBase->BaseParameters->flagJouledissipation = 1;
		else
			CgnsBase->BaseParameters->flagJouledissipation = 0;
		}	
		//bhans START
		if((CgnsBase->flagSolidification==1)||(CgnsBase->flagDiffKinetCont==2)||(CgnsBase->flagTwoPhase==1))
		{
			i=0;
			CgnsBase->VARTAG[count-1]=7;//Variable Tag has to be stored bhans
                        if(CgnsBase->flagDiffKinetCont==2)
                        initializeCombParameters(CgnsBase,var2);
			initialize_enthalpy(CgnsBase,var2);
		}
		else if(CgnsBase->flagTurbComb==0) //Shyam
		{
			i=1;
			CgnsBase->VARTAG[count-1]=6;//Variable Tag has to be stored bhans
        		fclose(var3);
			initialize_T(CgnsBase,var2);
                        printf("Temp. is initialized\n");
        		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		}
		//bhans END
	}
	else
	{
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->T= NULL;
			((z+j)->PrevFlowsol)->T= NULL;
			((z+j)->Flowsol)->enthalpy= NULL;
			((z+j)->PrevFlowsol)->enthalpy= NULL;
			((z+j)->Flowsol)->flagDiffcoeff=0;
			((z+j)->Flowsol)->flagSpecificHeat=0;
			((z+j)->Flowsol)->Diffcoeff=NULL;
			((z+j)->Flowsol)->SpecificHeat=NULL;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab for Temperature*/
					(((z+j)->Iface+iface)->ZSlab)->T = NULL;
					(((z+j)->Iface+iface)->ZSlab)->enthalpy=NULL;
					(((z+j)->Iface+iface)->ZSlab)->Diffcoeff=NULL;
					(((z+j)->Iface+iface)->ZSlab)->SpecificHeat=NULL;
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Boundary Slab for Temperature*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->T = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpy=NULL; 
					(((z+j)->ZoneBc+zBc)->ZSlab)->GradT = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->TOld = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpyOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->Diffcoeff=NULL; 
					(((z+j)->ZoneBc+zBc)->ZSlab)->SpecificHeat=NULL; 
					}
				}
			}
		}
	}
/*****************************EMD OF MEMORY ALLOCATION FOR TEMPERATURE******************************************************/
	//Initiation of combustion variables, Shyam
	if(CgnsBase->flagTurbComb==1)
        {
            if(CgnsBase->flagDiffKinetCont==1)
            initializeCombParameters(CgnsBase,var2);
            initializeCombVariables(CgnsBase,var2,&count);
        }

	fclose(var);fclose(var2);fclose(var3);
	
	//printf("\nCount before initializeElectric = %d\n", count);
// Memory Allocation for Electric Field Variables Begins Here
	if(CgnsBase->flagEfield==1 || CgnsBase->flagEfield==2)
	{
		cnt = count;
		count = initializeElectric(CgnsBase,cnt);
	}
	else
	{
		//Initialize electric variables to null
		for(j=0;j<num_zone;j++)
		{
			(z+j)->Flowsol->Potential=NULL;
			(z+j)->Flowsol->RhoIon=NULL;
			(z+j)->Flowsol->EfieldX=NULL;
			(z+j)->Flowsol->EfieldY=NULL;
			(z+j)->Flowsol->EfieldZ=NULL;
		}
	}
//	var=fopen("Settings/variablelist.dat","a");
//	var2=fopen("Settings/coefficients.dat","a");
//End of Memory Allocation for Electric Field Variables

/********************MEMORY ALLOCATION FOR LMMHD VARIABLES BEGINS HERE (naveen)***********************/
	if(CgnsBase->flagLMMHD==1)
	{
		cnt=count;
		if(CgnsBase->flaginduced==0)
		{
			for(j=0;j<num_zone;j++)
			{
				(z+j)->Flowsol->TotalMfieldX=NULL;//rajendra
				(z+j)->Flowsol->TotalMfieldY=NULL;//rajendra
				(z+j)->Flowsol->TotalMfieldZ=NULL;//rajendra
				(z+j)->Flowsol->Psi=NULL;
			}
		}
		else
		{
			for(j=0;j<num_zone;j++)
			{
				(z+j)->Flowsol->MfieldX=NULL;
				(z+j)->Flowsol->MfieldY=NULL;
				(z+j)->Flowsol->MfieldZ=NULL;
			}
		}
		count=initializeLMMHD(CgnsBase,cnt);
		
	}
	else
	{
		
		//Initialize LMMHD variables to null
		for(j=0;j<num_zone;j++)
		{
			(z+j)->Flowsol->PotentialMHD=NULL;
			(z+j)->Flowsol->Eleconductivity=NULL;
			(z+j)->Flowsol->Permeability=NULL;
			(z+j)->Flowsol->ElecurrentdensityX=NULL;
			(z+j)->Flowsol->ElecurrentdensityY=NULL;
			(z+j)->Flowsol->ElecurrentdensityZ=NULL;
			(z+j)->Flowsol->MfieldX=NULL;
			(z+j)->Flowsol->MfieldY=NULL;
			(z+j)->Flowsol->MfieldZ=NULL;
			(z+j)->Flowsol->TotalMfieldX=NULL;//rajendra
			(z+j)->Flowsol->TotalMfieldY=NULL;//rajendra
			(z+j)->Flowsol->TotalMfieldZ=NULL;//rajendra
			(z+j)->Flowsol->Psi=NULL;
		}
	}

/********************END OF MEMORY ALLOCATION FOR LMMHD VARIABLES***********************/

/*************************************Memory allocation for multiphase flow variables.*********************************/
	

	if(CgnsBase->flagMPFV==1)
	{	
		printf("\n Enter the no. of phase you want to solve\n");
		scanf("%d",&np);   
                
                var2=fopen("Settings/coefficients.dat","a");       	
                var3=fopen("Settings/ProgramReadInputValues.dat","a");              
                fprintf(var2,"%d\n",np);        
                fprintf(var3,"%d \tEnter the no. of phase you want to solve \n\n",np);
                fclose(var2);
                fclose(var3);
                
		CgnsBase->nphase=np;
		nop=CgnsBase->nphase;
		for(j=0;j<num_zone;j++)
		{	
			(z+j)->Nphase = np;
			nop = (z+j)->Nphase;
		}
		
		cnt=count;
		count=initializeMPFV(CgnsBase,cnt);
		
	}
	else
	{
		//Initialize multiphase flow variables to null
		CgnsBase->nphase=0;
		for(nophase=0;nophase<5;nophase++)
		{
			for(j=0;j<num_zone;j++)
			{
				
				
				(z+j)->Flowsol->RhoofPhase[nophase]=NULL;
				(z+j)->Flowsol->ViscosityofPhase[nophase]=NULL;
				(z+j)->Flowsol->U_phase[nophase]=NULL;
				(z+j)->Flowsol->V_phase[nophase]=NULL;
				(z+j)->Flowsol->W_phase[nophase]=NULL;
				(z+j)->Flowsol->Momentumofphasex[nophase]=NULL;
				(z+j)->Flowsol->Momentumofphasey[nophase]=NULL;
				(z+j)->Flowsol->Momentumofphasez[nophase]=NULL;
				(z+j)->Flowsol->Momentumx=NULL;
				(z+j)->Flowsol->Momentumxstar=NULL;
				(z+j)->Flowsol->Momentumy=NULL;
				(z+j)->Flowsol->Momentumystar=NULL;
				(z+j)->Flowsol->Momentumz=NULL;
				(z+j)->Flowsol->Momentumzstar=NULL;
				(z+j)->Flowsol->VarthetaPhase[nophase]=NULL;
				
				
				(z+j)->PrevFlowsol->RhoofPhase[nophase]=NULL;
				(z+j)->PrevFlowsol->ViscosityofPhase[nophase]=NULL;
				(z+j)->PrevFlowsol->U_phase[nophase]=NULL;
				(z+j)->PrevFlowsol->V_phase[nophase]=NULL;
				(z+j)->PrevFlowsol->W_phase[nophase]=NULL;
				(z+j)->PrevFlowsol->Momentumofphasex[nophase]=NULL;
				(z+j)->PrevFlowsol->Momentumofphasey[nophase]=NULL;
				(z+j)->PrevFlowsol->Momentumofphasez[nophase]=NULL;
				(z+j)->PrevFlowsol->Momentumx=NULL;
				(z+j)->PrevFlowsol->Momentumxstar=NULL;
				(z+j)->PrevFlowsol->Momentumy=NULL;
				(z+j)->PrevFlowsol->Momentumystar=NULL;
				(z+j)->PrevFlowsol->Momentumz=NULL;
				(z+j)->PrevFlowsol->Momentumzstar=NULL;
				(z+j)->PrevFlowsol->VarthetaPhase[nophase]=NULL;
			}
		}
	}
	
	if((CgnsBase->flagParticle==1)||(CgnsBase->flagDroplet==1))//adhi+shash
	{
		printf("\nNow enter the properties for disperse phase i.e. droplet or particle:\n");
		cnt=count;
		count=initializeDisperse(CgnsBase,cnt);
	}
	else
	{
		//Initialize Disperse flow variables to null
		for(j=0;j<num_zone;j++)
		{
				(z+j)->Flowsol->RhoDisperse=NULL;
				(z+j)->Flowsol->ViscosityDisperse=NULL;
				(z+j)->Flowsol->ViscoeffDummyDisperse=NULL;
				(z+j)->Flowsol->VarthetaDisperse=NULL;
//Commented by Mrunalini
/*
				(z+j)->Flowsol->MomentumDispersex=NULL;
				(z+j)->Flowsol->MomentumDispersey=NULL;
				(z+j)->Flowsol->MomentumDispersez=NULL;
				(z+j)->Flowsol->MomentumDispersexOld=NULL;
				(z+j)->Flowsol->MomentumDisperseyOld=NULL;
				(z+j)->Flowsol->MomentumDispersezOld=NULL;	
				(z+j)->Flowsol->NumDensity=NULL;
										*/
				(z+j)->Flowsol->U_Disperse=NULL;
				(z+j)->Flowsol->V_Disperse=NULL;
				(z+j)->Flowsol->W_Disperse=NULL;

				(z+j)->Flowsol->U_DisperseOld=NULL;
				(z+j)->Flowsol->V_DisperseOld=NULL;
				(z+j)->Flowsol->W_DisperseOld=NULL;

				((z+j)->Flowsol)->U_Dmass= NULL;
				((z+j)->Flowsol)->V_Dmass= NULL;
				((z+j)->Flowsol)->W_Dmass= NULL;

				(z+j)->Flowsol->ParticleDia=NULL;
			
				(z+j)->Flowsol->Evaporate= NULL;
				if(CgnsBase->flagMPFV==0)
				(z+j)->Flowsol->MassfluxDisperse=NULL;
				
				(z+j)->Flowsol->EvapoConst= NULL;
				
				(z+j)->Flowsol->Tdrop= NULL;
				(z+j)->Flowsol->CPdrop = NULL;
				(z+j)->Flowsol->flagCPdrop=0;
				(z+j)->Flowsol->DropFuelMassFrac= NULL;
				(z+j)->Flowsol->DropDiffcoeff= NULL;
				(z+j)->Flowsol->flagDropDiffcoeff=0;
				
				(z+j)->Flowsol->QDroptoGas= NULL;
				(z+j)->Flowsol->QLatent= NULL;
				
				(z+j)->PrevFlowsol->RhoDisperse=NULL;
				(z+j)->PrevFlowsol->ViscosityDisperse=NULL;
				(z+j)->PrevFlowsol->ViscoeffDummyDisperse=NULL;
				(z+j)->PrevFlowsol->VarthetaDisperse=NULL;
//Commented by Mrunalini
/*
				(z+j)->PrevFlowsol->MomentumDispersex=NULL;
				(z+j)->PrevFlowsol->MomentumDispersey=NULL;
				(z+j)->PrevFlowsol->MomentumDispersez=NULL;
				(z+j)->PrevFlowsol->MomentumDispersexOld=NULL;
				(z+j)->PrevFlowsol->MomentumDisperseyOld=NULL;
				(z+j)->PrevFlowsol->MomentumDispersezOld=NULL; 
				(z+j)->PrevFlowsol->NumDensity=NULL;
										*/
				(z+j)->PrevFlowsol->U_Disperse=NULL;
				(z+j)->PrevFlowsol->V_Disperse=NULL;
				(z+j)->PrevFlowsol->W_Disperse=NULL;

				(z+j)->PrevFlowsol->U_DisperseOld=NULL;
				(z+j)->PrevFlowsol->V_DisperseOld=NULL;
				(z+j)->PrevFlowsol->W_DisperseOld=NULL;

				((z+j)->PrevFlowsol)->U_Dmass= NULL;
				((z+j)->PrevFlowsol)->V_Dmass= NULL;
				((z+j)->PrevFlowsol)->W_Dmass= NULL;

				(z+j)->PrevFlowsol->ParticleDia=NULL;
				
				(z+j)->PrevFlowsol->Evaporate= NULL;
				if(CgnsBase->flagMPFV==0)
				(z+j)->PrevFlowsol->MassfluxDisperse=NULL;
				
				(z+j)->PrevFlowsol->EvapoConst= NULL;
				
				(z+j)->PrevFlowsol->Tdrop= NULL;
				(z+j)->PrevFlowsol->CPdrop = NULL;
				(z+j)->PrevFlowsol->DropFuelMassFrac= NULL;
				(z+j)->PrevFlowsol->DropDiffcoeff= NULL;
				
				(z+j)->PrevFlowsol->QDroptoGas= NULL;
				(z+j)->PrevFlowsol->QLatent= NULL;

				(z+j)->PrevFlowsol->RhoG_face= NULL;
				(z+j)->PrevFlowsol->VarthetaD_face= NULL;
				
				if(CgnsBase->flagMPFV==0)
				{	
					(z+j)->DragX=NULL;
			        	(z+j)->DragY=NULL;
					(z+j)->DragZ=NULL;
					(z+j)->PrevFlowsol->ParticleDia=NULL;
					(z+j)->PrevFlowsol->RhoDisperse=NULL;
					(z+j)->Flowsol->ParticleDia=NULL;
				}	
				(z+j)->t_p = NULL;
				
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab for Temperature*/
					(((z+j)->Iface+iface)->ZSlab)->RhoDisperse = NULL;
					(((z+j)->Iface+iface)->ZSlab)->ViscoeffDummyDisperse=NULL;
					(((z+j)->Iface+iface)->ZSlab)->ViscosityDisperse=NULL;
					(((z+j)->Iface+iface)->ZSlab)->VarthetaDisperse=NULL;

//Commented by Mrunalini
/*
					(((z+j)->Iface+iface)->ZSlab)->MomentumDispersex=NULL;
					(((z+j)->Iface+iface)->ZSlab)->MomentumDispersexOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->MomentumDispersey=NULL;
					(((z+j)->Iface+iface)->ZSlab)->MomentumDisperseyOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->MomentumDispersez=NULL;
					(((z+j)->Iface+iface)->ZSlab)->MomentumDispersezOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->NumDensity=NULL; 	*/
			
					(((z+j)->Iface+iface)->ZSlab)->U_Disperse=NULL;
					(((z+j)->Iface+iface)->ZSlab)->V_Disperse=NULL;
					(((z+j)->Iface+iface)->ZSlab)->W_Disperse=NULL;
				
					(((z+j)->Iface+iface)->ZSlab)->DropFuelMassFrac=NULL;
					(((z+j)->Iface+iface)->ZSlab)->DropDiffcoeff=NULL;
					(((z+j)->Iface+iface)->ZSlab)->Tdrop=NULL;
					(((z+j)->Iface+iface)->ZSlab)->CPdrop=NULL;

					(((z+j)->Iface+iface)->ZSlab)->U_DisperseOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->V_DisperseOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->W_DisperseOld=NULL;

					(((z+j)->Iface+iface)->ZSlab)->U_Dmass =NULL;
					(((z+j)->Iface+iface)->ZSlab)->V_Dmass =NULL;
					(((z+j)->Iface+iface)->ZSlab)->W_Dmass =NULL;

					(((z+j)->Iface+iface)->ZSlab)->U_DmassOld =NULL;
					(((z+j)->Iface+iface)->ZSlab)->V_DmassOld =NULL;
					(((z+j)->Iface+iface)->ZSlab)->W_DmassOld =NULL;

					(((z+j)->Iface+iface)->ZSlab)->VarthetaDisperseOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->DropFuelMassFracOld=NULL;
					(((z+j)->Iface+iface)->ZSlab)->TdropOld=NULL;

					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Boundary Slab for Temperature*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->RhoDisperse = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->ViscoeffDummyDisperse=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->ViscosityDisperse=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->VarthetaDisperse=NULL;
//Commented by Mrunalini
/*
					(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersex=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersexOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersey=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDisperseyOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersez=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersezOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->NumDensity=NULL; 			*/


					(((z+j)->ZoneBc+zBc)->ZSlab)->U_Disperse=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->V_Disperse=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->W_Disperse=NULL;
					
					(((z+j)->ZoneBc+zBc)->ZSlab)->DropFuelMassFrac=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->DropDiffcoeff=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->Tdrop=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->CPdrop=NULL; 

					(((z+j)->ZoneBc+zBc)->ZSlab)->U_DisperseOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->V_DisperseOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->W_DisperseOld=NULL;

					(((z+j)->ZoneBc+zBc)->ZSlab)->U_Dmass = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->V_Dmass = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->W_Dmass = NULL;

					(((z+j)->ZoneBc+zBc)->ZSlab)->U_DmassOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->V_DmassOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->W_DmassOld=NULL;

					(((z+j)->ZoneBc+zBc)->ZSlab)->VarthetaDisperseOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->DropFuelMassFracOld=NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->TdropOld=NULL;
					 
					 
					}
				}
			}
				
		}
	}	
	
	var=fopen("Settings/variablelist.dat","a");
	var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
/*******************************End of Memory allocation for multiphase flow variables.*********************************/




/*****************************MEMORY ALLOCATION FOR SCALARS*******************************************************/	
	/*Enter the number of the Scalars to be solved */
	if(CgnsBase->flagTurb == 0)
		printf("\nEnter number of Scalars to be solved [Max. of 10 Scalars is allowed]: ");
	else
	{
		printf("\nEnter Number of Scalars to be solved [Max. of 10 Scalars is allowed]:");

  		if(CgnsBase->flagTurbModel == 1 || CgnsBase->flagTurbModel == 3 || CgnsBase->flagTurbModel == 4 || CgnsBase->flagTurbModel == 5)
			printf("\nPlease note that Scalar1 and Scalar2 are considered as k and Epsilon");
		if(CgnsBase->flagTurbModel == 2 || CgnsBase->flagTurbModel == 6)
			printf("\nPlease note that Scalar1 and Scalar2 are considered as k and Omega");		  
		
		if(CgnsBase->flagTurbModel == 8)
		{
			printf("\nSpallart-Almaras is ONE-equation Model, Hence Enter number of scalars as 1 + desired number of scalars.");
		} else if (CgnsBase->flagTurbModel == 9) {
			printf("\nBaldwin-Lomax Model is ZERO-equation Model, Hence Enter number of scalars as 0 + desired number of scalars.");
		} else {
			printf("\nHence if you are solving for scalars, Enter number of scalars as 2 + desired number of scalars.");
		}

		printf("\n\nNow enter number of scalars: ");
	}
	printf("\n");
	scanf("%d",&scal);
	while(scal>10)
	{
		printf("\nMax allowed scalars is 10");
		printf("\nPlease enter the Number of Scalars(<= 10):");
		scanf("%d",&scal);
	}
	fprintf(var,"%d\n",scal);                        
        fprintf(var2,"%d\n",scal);        
        fprintf(var3,"%d \tPlease enter the Number of Scalars(<= 10): \n\n",scal);
	
	//CgnsBase->BaseParameters->Turb = (ParametersTurb *)calloc(1,sizeof(ParametersTurb));
	
	/*	if((CgnsBase->flagLogKEps == 0 || CgnsBase->flagTurbModel == 3) && CgnsBase->BaseParameters->restart == 0)//Is there any use of this..?-->Chandrasekhar
	{
		printf("***** Do you want to AUTO-INITIALIZE k and epsilon? *****\n");
		printf("NOTE: auto-initizalization is done using the assumption of equillibrium between production and dissipation rate of turbulent kinetic energy\n");
		printf("\nEnter 0->NO/1->YES: ");
		
		do
		{
			//scanf("%d", &flag);
			flag = 0;//Written for ad hoc purpose. <<< Narsimha 18th April, 2008
		}while(!(flag==0 || flag==1));
				
		CgnsBase->BaseParameters->Turb->flagAutoInit = flag;
		}*/
	
	//The looping over all the scalars starts here
	for(k=0;k<scal;k++)
	{
		i=1;
		CgnsBase->VARTAG[count+k]=600+k;//Variable Tag has to be stored bhans	
		
		if(CgnsBase->flagLMMHD==1)
			varCoeffCount = 7+k+1;  //4->SpecificHeat,5->Magnetic //rajendra Diffusivity,6->Scalcoeff1,.....15->Scalcoeff10.
		else
			varCoeffCount = 4+k+1;  //4->SpecificHeat,5->Scalcoeff1,6->Scalcoeff2,.....14->Scalcoeff10.
		

		printf("\nCreating memory for Scalar Coefficient[%d]:",k+1);

        	if( ( (k ==0) && (CgnsBase->flagTurb == 1) && (CgnsBase->flagTurbModel <= 8) ) ||
            	    ( (k ==1) && (CgnsBase->flagTurb == 1) && (CgnsBase->flagTurbModel <= 6) )	) {
			flag = 0;
		} else {
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify Scalar Coefficient[%d] for the WHOLE DOMAIN (Enter 0) or ZONE BY ZONE (Enter 1):",k+1);
				scanf("%d",&flag);
				fprintf(var2,"%d\n",flag);        
                		fprintf(var3,"%d \tDo you want to specify Scalar Coefficient[%d] for the WHOLE DOMAIN (Enter 0) or ZONE BY ZONE (Enter 1): \n\n",flag,k+1);
			} else { flag = 0;}	
		}


		if(flag==0) //Initializing for whole domain at a time
		{
        	if( ( (k ==0) && (CgnsBase->flagTurb == 1) && (CgnsBase->flagTurbModel <= 8) ) ||
              	    ( (k ==1) && (CgnsBase->flagTurb == 1) && (CgnsBase->flagTurbModel <= 6) )	) {
			p = 1;
			constant = 1.0;
		} else {
			printf("\n If Scalar Coefficient[%d] is CONSTANT ,Enter 0",k+1);
			printf("\n If Scalar Coefficient[%d] is a VARIABLE,Enter 1:",k+1);
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);        
                        fprintf(var3,"%d \tIf Scalar Coefficient[%d] is CONSTANT Enter 0: is a VARIABLE,Enter 1:\n\n",p,k+1);
			
			if(p==0)
			{
				printf("\nEnter the constant value of Scalar Coefficient[%d]:",k+1);
				scanf("%lf",&constant);
				fprintf(var2,"%12.8e\n",constant);        
                                fprintf(var3,"%12.8e \tEnter the constant value of Scalar Coefficient[%d]: \n\n",constant,k+1);
			}
			else
				constant=1.0;		
		}

        	if( (k ==0) && (CgnsBase->flagTurb == 1) ) {
			if(CgnsBase->flagTurbModel <= 6) {
				//printf("\nCreating memory for Scalar Coefficient[%d] - Turbulent Kinetic Energy:",k+1);
				printf(" - Turbulent Kinetic Energy:");
			}
			if(CgnsBase->flagTurbModel == 8) {
				printf(" - Spallart-Almaras is ONE-equation Model:");
			}
		} 
        	if( (k ==1) && (CgnsBase->flagTurb == 1) ) {
  			if(CgnsBase->flagTurbModel == 1 || CgnsBase->flagTurbModel == 3 || CgnsBase->flagTurbModel == 4 || CgnsBase->flagTurbModel == 5) {
				//printf("\nCreating memory for Scalar Coefficient[%d] - Dissipation of turbulent Kinetic Energy (Epsilon):",k+1);
				printf(" - Dissipation of turbulent Kinetic Energy (Epsilon):");
			}
			if(CgnsBase->flagTurbModel == 2 || CgnsBase->flagTurbModel == 6) {
				//printf("\nCreating memory for Scalar Coefficient[%d] - Omega:",k+1);
				printf(" - Omega (specific dissipation rate):");
			}
		} 

			for(j=0;j<num_zone;j++)
			{
				((z+j)->Flowsol)->Scalcoeff[k]= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->PrevFlowsol)->Scalcoeff[k]= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->Flowsol)->flagScalcoeff[k]=p;
				
				numiface=(z+j)->Number_of_Interfaces;
				//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{	if(((z+j)->Iface+iface)->Direction==dir)
						{
							(((z+j)->Iface+iface)->ZSlab)->Scalcoeff[k]=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Scalcoeff[k] = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}
		}
		
		if(flag==1)
		{
			for(j=0;j<num_zone;j++)
			{ 
				printf("Zone Name is %s\n",(z+j)->zonename);
				printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				printf("\n If Scalar Coefficient[%d] is CONSTANT,   Enter 0",k+1);
				printf("\n If Scalar Coefficient[%d] is a VARIABLE, Enter 1 :",k+1);
				scanf("%d",&p);
				fprintf(var2,"%d\n",p);        
                                fprintf(var3,"%d \tIf Scalar Coefficient[%d] is CONSTANT  Enter 0: is a VARIABLE, Enter 1 : \n\n",p,k+1);
				if(p==0)
				{
					printf("\nEnter the constant value of Scalar Coefficient[%d]:",k+1);
					scanf("%lf",&constant);
					fprintf(var2,"%12.8e\n",constant);        
                                        fprintf(var3,"%12.8e \tEnter the constant value of Scalar Coefficient[%d]: \n\n",constant,k+1);
				}
				else
					constant=1.0;

				((z+j)->Flowsol)->Scalcoeff[k]= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->PrevFlowsol)->Scalcoeff[k]= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				((z+j)->Flowsol)->flagScalcoeff[k]=p;
				numiface=(z+j)->Number_of_Interfaces;
				//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							(((z+j)->Iface+iface)->ZSlab)->Scalcoeff[k]=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Scalcoeff[k] = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}
		}
			
			
		if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
		   if( CgnsBase->Number_of_zones > 1){
			printf("\nDo you want to initialize Scal[%d] for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): ",k+1);
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);        
                        fprintf(var3,"%d \tDo you want to initialize Scal[%d] for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",p, k+1);
		   } 
		   if(p==0){
			printf("\nEnter the Initialization value of Scal[%d] for all Zones: ",k+1);
			scanf("%lf",&constant);
			fprintf(var2,"%12.8e\n",constant);        
                	fprintf(var3,"%12.8e \tEnter the Initialization value of Scal[%d] for all Zones: \n\n",constant,k+1);
		   }
		}
		for(j=0;j<num_zone;j++)
		{
// 			if(CgnsBase->BaseParameters->Turb->flagAutoInit == 1 && (k==0 || k==1))
// 			{
// 				constant = 0.0;
// 			}
// 			else
			//{
		   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
			if(p!=0){
				printf("\nEnter the Initialization value of Scal[%d] in zone %s:",k+1,(z+j)->zonename);
				scanf("%lf",&constant);
				fprintf(var2,"%12.8e\n",constant);        
                                fprintf(var3,"%12.8e \tEnter the Initialization value of Scal[%d] in zone %s \n\n",constant,k+1,(z+j)->zonename);
			}
		   } else { constant = 0.0;}
			//}

			//constant = 0.0;//Written for ad hoc purpose <<< Narsimha.
			
			if(CgnsBase->flagLogKEps == 1 && (k==0 || k==1))
			{
				if(constant != 0.0)
				{
					printf("This value [ %e ] is converted into its Logarithmic value", constant);
					constant = log(constant);
					printf(" [ %e ]\n", constant);
				}
				else
				{
					printf("Enter non-zero initial values for turbulence variables (Scalars 1 and 2) since Logarithmic form of k-epsilon model is choosen\n");
					exit(4); 
				}
			}

			/*Memory allocation in the Flow Solution Node For Scalars*/
			((z+j)->Flowsol)->Scal[k]= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->Scal[k]= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						//Memory allocation in the Interface Slab for Scalar
						(((z+j)->Iface+iface)->ZSlab)->Scal[k]= (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
						(((z+j)->Iface+iface)->ZSlab)->ScalOld[k]= (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						//Memory allocation in the Boundary Slab for Scalar
						(((z+j)->ZoneBc+zBc)->ZSlab)->Scal[k]= (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->ScalOld[k]= (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						(((z+j)->ZoneBc+zBc)->ZSlab)->GradScal[k]= (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					}
				}
			}
		}
	} /*END OF LOOPING OVER SOLVABLE SCALARS*/
	
	printf("\nDo you want to calculate VORTICITY?\nEnter:1(if YES)\nEnter:0(if NO)\n");
	scanf("%d",&vortchoice);
	fprintf(var,"%d\n",vortchoice);                        
        fprintf(var2,"%d\n",vortchoice);        
        fprintf(var3,"%d \tDo you want to calculate VORTICITY? Enter:1(if YES) Enter:0(if NO) \n\n",vortchoice);
	if(vortchoice==1)
		CgnsBase->BaseParameters->tagVorticity=1;
	else
		CgnsBase->BaseParameters->tagVorticity=0;
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->vorticity_X= (double *)initialize_variable(vortchoice,(z+j),0.5);
		((z+j)->Flowsol)->vorticity_Y= (double *)initialize_variable(vortchoice,(z+j),0.5);
		((z+j)->Flowsol)->vorticity_Z= (double *)initialize_variable(vortchoice,(z+j),0.5);
	}

	if(CgnsBase->BaseParameters->TimeStepingSchemeTag==1)
		for(j=0;j<num_zone;j++)
		{
			if((z+j)->Flowsol->u!=NULL)
				((z+j)->Flowsol)->PrevTStepFluxes_u = (double *)initialize_variable(1,(z+j),0.5);
			if((z+j)->Flowsol->v!=NULL)
				((z+j)->Flowsol)->PrevTStepFluxes_v = (double *)initialize_variable(1,(z+j),0.5);
			if((z+j)->Flowsol->w!=NULL)
				((z+j)->Flowsol)->PrevTStepFluxes_w = (double *)initialize_variable(1,(z+j),0.5);
			if((z+j)->Flowsol->T!=NULL)
			((z+j)->Flowsol)->PrevTStepFluxes_T = (double *)initialize_variable(1,(z+j),0.5);
			for(cntVar=0;cntVar<10;cntVar++)
				if((z+j)->Flowsol->Scal[cntVar]!=NULL)
					((z+j)->Flowsol)->PrevTStepFluxes_Scal[cntVar] = (double *)initialize_variable(1,(z+j),0.5);
		}

	
	/*The following loop would make the pointers for all the scalars which are not to be solved, point to NULL */
	i=0;
	for(k=scal;k<10;k++)
	{
		//fprintf(var,"%d\n",i);
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->Scal[k]= NULL;
			((z+j)->PrevFlowsol)->Scal[k]= NULL;
			((z+j)->Flowsol)->Scalcoeff[k]= NULL;
			((z+j)->PrevFlowsol)->Scalcoeff[k]= NULL;
			((z+j)->Flowsol)->flagScalcoeff[k]=0;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
					(((z+j)->Iface+iface)->ZSlab)->Scalcoeff[k]= NULL;
					(((z+j)->Iface+iface)->ZSlab)->Scal[k]= NULL;
					(((z+j)->Iface+iface)->ZSlab)->ScalOld[k]= NULL;
					}
				}
				for(zBc=0;zBc<numBc;zBc++)
				{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					(((z+j)->ZoneBc+zBc)->ZSlab)->Scalcoeff[k]= NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->Scal[k]= NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->ScalOld[k]= NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->GradScal[k]= NULL;
					}
				}
			}
		}
	}
	
	//Making combustion variable pointers to NULL, when combustion is not selected
	if(CgnsBase->flagTurbComb==0)//Shyam, START
	{
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->MixFrac= NULL;
			((z+j)->Flowsol)->FuelMassFrac= NULL;
			for(i=0;i<10;i++)
			{
				((z+j)->Flowsol)->InMassFrac[i]= NULL;
				((z+j)->PrevFlowsol)->InMassFrac[i]= NULL;
			}
			((z+j)->Flowsol)->OxiMassFrac= NULL;
			((z+j)->Flowsol)->ProdMassFrac= NULL;                        
			((z+j)->PrevFlowsol)->MixFrac= NULL;
			((z+j)->PrevFlowsol)->FuelMassFrac= NULL;
			((z+j)->PrevFlowsol)->OxiMassFrac= NULL;
			((z+j)->PrevFlowsol)->ProdMassFrac= NULL;
	
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
	
			for(iface=0;iface<numiface;iface++)
			{                                      
				(((z+j)->Iface+iface)->ZSlab)->OxiMassFrac= NULL;
				(((z+j)->Iface+iface)->ZSlab)->MixFrac= NULL;
				(((z+j)->Iface+iface)->ZSlab)->MixFracOld= NULL;
				(((z+j)->Iface+iface)->ZSlab)->FuelMassFrac= NULL;
				(((z+j)->Iface+iface)->ZSlab)->FuelMassFracOld= NULL;
				for(i=0;i<10;i++)
				{
					(((z+j)->Iface+iface)->ZSlab)->InMassFrac[i]= NULL;
					(((z+j)->Iface+iface)->ZSlab)->InMassFracOld[i]= NULL;
				}
				(((z+j)->Iface+iface)->ZSlab)->ProdMassFrac= NULL;
				(((z+j)->Iface+iface)->ZSlab)->ProdMassFracOld= NULL;
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				
				(((z+j)->ZoneBc+zBc)->ZSlab)->MixFrac= NULL;
				(((z+j)->ZoneBc+zBc)->ZSlab)->OxiMassFrac= NULL;
				(((z+j)->ZoneBc+zBc)->ZSlab)->MixFracOld= NULL;
				(((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFrac= NULL;
				(((z+j)->ZoneBc+zBc)->ZSlab)->FuelMassFracOld= NULL;
				for(i=0;i<10;i++)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->InMassFrac[i]= NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->InMassFracOld[i]= NULL;
				}
				(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFrac= NULL;
				(((z+j)->ZoneBc+zBc)->ZSlab)->ProdMassFracOld= NULL;
			}
		}
	}//Shyam, END.
	
	if(CgnsBase->flagTurb==1)
	{
		for(j=0;j<num_zone;j++)
		{
			if((z+j)->flagSolidFluid==0)
			{
// 				if(CgnsBase->BaseParameters->Turb->flagAutoInit == 1)
// 				{
// 					printf("Enter the initialization value of Turbulent Viscosity (No function)\n");
// 					scanf("%lf",&initViscosity);        
//                                      fprintf(var3,"%12.8e \tEnter the initialization value of Turbulent Viscosity (No function) \n\n",initViscosity);
// 				}
// 				else
				{
					initViscosity = 0.0;
				}
				initconductivity=0.0;
				initeffconductivity=0.0;
//				fprintf(var2,"%12.8e\n",initViscosity);
				((z+j)->Flowsol)->flagTurb=1;
				
				if(CgnsBase->flagTurbComb==0)//Shyam, If combustion is selected, 'Tcoeff' is initiated in initializeCombVariables(), Shyam
				{
					((z+j)->Flowsol)->Tcoeff=(TurbCoeff *)malloc(sizeof(TurbCoeff));
				}
				else if(CgnsBase->flagTurbComb==1)
				{
					((z+j)->Flowsol)->Tcoeff->EffScalcoeff=(double *)initialize_Effviscosity(1,(CgnsBase->CgnsZone+j),initViscosity); //Shyam
					((z+j)->Flowsol)->Tcoeff->sigmacombvar = 0.7;
				}
				
				((z+j)->Flowsol)->Tcoeff->TViscoeff=(double *)initialize_variable(1,(z+j),initViscosity);
				((z+j)->Flowsol)->Tcoeff->EffViscoeff=(double *)initialize_Effviscosity(1,(z+j),initViscosity);
				((z+j)->Flowsol)->Tcoeff->TProduction = (double *)initialize_variable(1,(z+j),0.0);
				((z+j)->Flowsol)->Tcoeff->f1 = (double *)initialize_variable(1,(z+j),1.0);
				((z+j)->Flowsol)->Tcoeff->Chi = 0.42;
				((z+j)->Flowsol)->Tcoeff->Cmu = 0.09;
				((z+j)->Flowsol)->Tcoeff->sigmaVolfrac =CgnsBase->number;//Turbulent prandtl number for turbulent dispersion--Turbulent MPFV
				//printf("\nTurbulent prandtl number is %lf",CgnsBase->number);
				((z+j)->Flowsol)->Tcoeff->sigma_f = 5.0;//urbulent prandtl number in case of Turbulent Gas-Particle flow
				((z+j)->Flowsol)->Tcoeff->sigma_r = 1.0e10;//Turbulent prandtl number for turbulent dispersion--Turbulent Gas-Particle flow
				//Added by chandrasekhar to Initialize the turbulent, effective , dummy effective viscosity of the phase in case of Turbulent MPFV
				if(CgnsBase->flagMPFV==1)
				{
					//Initiliazation of phase turbulent viscosity ,phase effective viscosity and phase dummy effective viscosity
					for(nophase=0;nophase<nop;nophase++)	
					{
						((z+j)->Flowsol)->Tcoeff->TViscoeffofphase[nophase]=(double *)initialize_variable(1,(z+j),initViscosity);	
						((z+j)->Flowsol)->Tcoeff->EffViscoeffofphase[nophase]=(double *)initialize_Effviscosityofphase(1,(z+j),initViscosity,nophase);
						((z+j)->Flowsol)->Tcoeff->DummyeffViscoeffofphase[nophase]=(double *)initialize_Dummyeffviscosityofphase(1,(z+j),initViscosity,nophase);
						((z+j)->Flowsol)->Tcoeff->DiffVolfrac[nophase]=(double *)initialize_variable(1,(z+j),initViscosity);	
					}
				}
				//End --- chandrasekhar
				
				//Added to initialize the turbulent, dummyviscosities and diffusion coeff. of volume fraction equation for disperse phase in case of Turbulent Gas-Particle flow
				if(CgnsBase->flagParticle==1)
				{
					((z+j)->Flowsol)->Tcoeff->TViscoeff_disperse=(double *)initialize_variable(1,(z+j),initViscosity);
					((z+j)->Flowsol)->Tcoeff->DummyViscoeff_disperse=(double *)initialize_variable(1,(z+j),initViscosity);	
					((z+j)->Flowsol)->Tcoeff->DiffVolfrac_disperse=(double *)initialize_variable(1,(z+j),initViscosity);
					((z+j)->Flowsol)->Tcoeff->T_e=(double *)initialize_variable(1,(z+j),initViscosity);	
				}
				
				if(((z+j)->Flowsol)->T!=NULL)//if turbulent energy equation is to be solved
				{
					((z+j)->Flowsol)->Tcoeff->Prt=1.0;
					((z+j)->Flowsol)->Tcoeff->TurbConductivity=(double *)initialize_variable(1,(z+j),initconductivity);
					((z+j)->Flowsol)->Tcoeff->EffTurbConductivity=(double *)initialize_variable(1,(z+j),initeffconductivity);
				}
				
				//Constants for k-Epsilon Model
				if(CgnsBase->flagTurbModel == 1)
				{
					((z+j)->Flowsol)->Tcoeff->C0epsi=(double *)initialize_variable(1,(z+j),0.0);
					((z+j)->Flowsol)->Tcoeff->C1epsi=(double *)initialize_variable(1,(z+j),1.44);
					((z+j)->Flowsol)->Tcoeff->C2epsi=(double *)initialize_variable(1,(z+j),1.92);
					((z+j)->Flowsol)->Tcoeff->TEprandtl=(double *)initialize_variable(1,(z+j),1.39);
					((z+j)->Flowsol)->Tcoeff->TKprandtl=(double *)initialize_variable(1,(z+j),1.0);
					((z+j)->Flowsol)->Tcoeff->sigmaK = 1.0;
					((z+j)->Flowsol)->Tcoeff->sigmaEpsilon = 1.39;
				}
				//Constants for k-Omega Model
				if(CgnsBase->flagTurbModel == 2 || CgnsBase->flagTurbModel == 6)
				{
					((z+j)->Flowsol)->Tcoeff->gamma		= 5.0/9.0;//(double *)initialize_variable(1,(z+j),5.0/9.0);
					((z+j)->Flowsol)->Tcoeff->gammastar 	= 1.0;//(double *)initialize_variable(1,(z+j),1.0);
					((z+j)->Flowsol)->Tcoeff->beta 		= 3.0/40.0;//(double *)initialize_variable(1,(z+j),3.0/40.0);
					((z+j)->Flowsol)->Tcoeff->betastar	= 0.09;//(double *)initialize_variable(1,(z+j),0.09);
					((z+j)->Flowsol)->Tcoeff->sigmaOmega 	= 0.5;//(double *)initialize_variable(1,(z+j),0.5);
					((z+j)->Flowsol)->Tcoeff->sigmaK 	= 0.5;//(double *)initialize_variable(1,(z+j),0.5);
				}
				//Constants for RNG K-epsilon Model
				if(CgnsBase->flagTurbModel == 3)
				{
					((CgnsBase->CgnsZone+j)->Flowsol)->flagTurb=1;
					((CgnsBase->CgnsZone+j)->Flowsol)->flagRNG_KEPS=1;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C2epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.68);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C1epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.42);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C0epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),0.0);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->TEprandtl=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.39);//Any use???
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->TKprandtl=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.0);//Any use???
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->sigmaK = 0.7179;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->sigmaEpsilon = 0.7179;
				}
				else
				{
					((CgnsBase->CgnsZone+j)->Flowsol)->flagRNG_KEPS=0;
				}

				if(CgnsBase->flagTurbModel == 4) // Launder and Sharma
				{
					((CgnsBase->CgnsZone+j)->Flowsol)->flagTurb=4;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C2epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.92);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C1epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.44);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C0epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),0.0);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->TEprandtl=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.39);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->TKprandtl=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.0);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->sigmaK = 1.0;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->sigmaEpsilon = 1.3;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->PhiK = (double *)initialize_variable(1,(CgnsBase->CgnsZone+j),0.0);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->PhiEps = (double *)initialize_variable(1,(CgnsBase->CgnsZone+j),0.0);
				}

				if(CgnsBase->flagTurbModel == 5)
				{
					((CgnsBase->CgnsZone+j)->Flowsol)->flagTurb=5;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C2epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.92);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C1epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.44);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->C0epsi=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),0.0);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->TEprandtl=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.39);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->TKprandtl=(double *)initialize_variable(1,(CgnsBase->CgnsZone+j),1.0);
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->sigmaK = 1.0;
					((CgnsBase->CgnsZone+j)->Flowsol)->Tcoeff->sigmaEpsilon = 1.3;
				}

				//Constants for SST Model
				if(CgnsBase->flagTurbModel == 6)
				{
					betastar	= 0.09;
					kappa 		= 0.41;
					//Set 1
					SigmaOmega1 	= 0.5;
					Beta1		= 0.0750;
					//Set 2
					SigmaOmega2 	= 0.856;
					Beta2		= 0.0828;

					((z+j)->Flowsol)->Tcoeff->betastar	= 0.09;
					((z+j)->Flowsol)->Tcoeff->a1		= 0.31;
					//Set 1 (Wilcox)
					((z+j)->Flowsol)->Tcoeff->SigmaK1	= 0.85;
					((z+j)->Flowsol)->Tcoeff->SigmaOmega1	= SigmaOmega1;
					((z+j)->Flowsol)->Tcoeff->Beta1		= Beta1;
					((z+j)->Flowsol)->Tcoeff->Gamma1	= Beta1/betastar - SigmaOmega1*kappa*kappa/sqrt(betastar);

					//Set 2 (Jones-Launder)
					((z+j)->Flowsol)->Tcoeff->SigmaK2	= 1.0;
					((z+j)->Flowsol)->Tcoeff->SigmaOmega2	= SigmaOmega2;
					((z+j)->Flowsol)->Tcoeff->Beta2		= Beta2;
					((z+j)->Flowsol)->Tcoeff->Gamma2	= Beta2/betastar - SigmaOmega2*kappa*kappa/sqrt(betastar);

					((z+j)->Flowsol)->Tcoeff->F1		= (double *)initialize_variable(1,(z+j),1.0);
				}
				else
				{
					((z+j)->Flowsol)->Tcoeff->F1		= NULL;
				}

				if(CgnsBase->flagTurbModel == 5 || CgnsBase->flagTurbModel == 6 || CgnsBase->flagTurbModel == 8 || CgnsBase->flagTurbModel == 9)
				{
					(z+j)->ZGeom->minWallDistance		= (double *)initialize_variable(1,(z+j),0.0);
				}
				else
				{
					(z+j)->ZGeom->minWallDistance		= NULL;
				}
				
				if(CgnsBase->flagTurbModel == 6 || CgnsBase->flagTurbModel == 8 || CgnsBase->flagTurbModel == 9)
				{
					(z+j)->Flowsol->Vorticity		= (double *)initialize_variable(1,(z+j),0.0);
				}
				else
				{
					(z+j)->Flowsol->Vorticity		= NULL;
				}
				
				//Parameters for Baldwin-Lomax Model
				if(CgnsBase->flagTurbModel == 9)
				{
					(z+j)->Flowsol->Tcoeff->BL		= (BaldwinLomax *)calloc(1,sizeof(BaldwinLomax));
					
					(z+j)->Flowsol->Tcoeff->BL->MixingLength= (double *)initialize_variable(1,(z+j),0.0);
					(z+j)->Flowsol->Tcoeff->BL->mu_ti	= (double *)initialize_variable(1,(z+j),0.0);
					(z+j)->Flowsol->Tcoeff->BL->mu_to	= (double *)initialize_variable(1,(z+j),0.0);
				}
				else
				{
					(z+j)->Flowsol->Tcoeff->BL		= NULL;
				}
				
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						IfacePtr = (z+j)->Iface+iface;
						if(IfacePtr->Direction==dir)
						{
							(IfacePtr->ZSlab)->EffViscoeff=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
							(IfacePtr->ZSlab)->TViscoeff=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
							
							//Added by chandrasekhar to Initialize the turbulent, effective , dummy effective viscosity of the phase in case of Turbulent MPFV at interfaces.
							if(CgnsBase->flagMPFV==1)
							{
								for(nophase=0;nophase<nop;nophase++)
								{
									(IfacePtr->ZSlab)->TViscoeffofphase[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);	
									(IfacePtr->ZSlab)->EffViscoeffofphase[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
									(IfacePtr->ZSlab)->DummyeffViscoeffofphase[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
									(IfacePtr->ZSlab)->DiffVolfrac[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
								}
							}
							//End --- chandrasekhar
							
							if(CgnsBase->flagParticle==1)
							{
								(IfacePtr->ZSlab)->TViscoeff_disperse=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);	
								(IfacePtr->ZSlab)->DummyViscoeff_disperse=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
								(IfacePtr->ZSlab)->DiffVolfrac_disperse=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
								(IfacePtr->ZSlab)->T_e=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
							}

							if(CgnsBase->flagTurbModel == 5 || CgnsBase->flagTurbModel == 9)
								(IfacePtr->ZSlab)->minWallDistance = (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
							else
								(IfacePtr->ZSlab)->minWallDistance = NULL;
							
							if((CgnsBase->flagTurbComb==0) && ((IfacePtr->ZSlab)->T!=NULL))//if turbulent energy equation is to be solved	
							{
								(IfacePtr->ZSlab)->TurbConductivity=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
								(IfacePtr->ZSlab)->EffTurbConductivity=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
							}
 							
							if(CgnsBase->flagTurbComb==1)//Shyam
								(IfacePtr->ZSlab)->EffScalcoeff=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,IfacePtr->Point_Range);
							
							if(CgnsBase->flagTurbModel == 9)
							{
								IfacePtr->ZSlab->BL			= (BaldwinLomax *)calloc(1,sizeof(BaldwinLomax));
								
								IfacePtr->ZSlab->BL->MixingLength 	= (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
								IfacePtr->ZSlab->BL->mu_ti		= (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
								IfacePtr->ZSlab->BL->mu_to		= (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
							}
							else
							{
								IfacePtr->ZSlab->BL = NULL;
							}
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
							(((z+j)->ZoneBc+zBc)->ZSlab)->EffViscoeff = (double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							(((z+j)->ZoneBc+zBc)->ZSlab)->TViscoeff = (double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							
							//Added by chandrasekhar to Initialize the turbulent, effective , dummy effective viscosity of the phase in case of Turbulent MPFV at boundaries
							if(CgnsBase->flagMPFV==1)
							{
								for(nophase=0;nophase<nop;nophase++)
								{
									(((z+j)->ZoneBc+zBc)->ZSlab)->TViscoeffofphase[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);	
									(((z+j)->ZoneBc+zBc)->ZSlab)->EffViscoeffofphase[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
									(((z+j)->ZoneBc+zBc)->ZSlab)->DummyeffViscoeffofphase[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
									(((z+j)->ZoneBc+zBc)->ZSlab)->DiffVolfrac[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);	
								}
							}
							//End --- chandrasekhar
							
							if(CgnsBase->flagParticle==1)
							{
								(((z+j)->ZoneBc+zBc)->ZSlab)->TViscoeff_disperse=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
								(((z+j)->ZoneBc+zBc)->ZSlab)->DummyViscoeff_disperse=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);	
								(((z+j)->ZoneBc+zBc)->ZSlab)->DiffVolfrac_disperse=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
								(((z+j)->ZoneBc+zBc)->ZSlab)->T_e=(double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);	
							}
							
							if(CgnsBase->flagTurbModel == 5 || CgnsBase->flagTurbModel == 9)
								(((z+j)->ZoneBc+zBc)->ZSlab)->minWallDistance = (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							else
								(((z+j)->ZoneBc+zBc)->ZSlab)->minWallDistance = NULL;
						
							if((((z+j)->ZoneBc+zBc)->ZSlab)->T!=NULL)//if turbulent energy equation is to be solved
							{
								(((z+j)->ZoneBc+zBc)->ZSlab)->TurbConductivity = (double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
								(((z+j)->ZoneBc+zBc)->ZSlab)->EffTurbConductivity = (double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							}
							
							if(CgnsBase->flagTurbComb==1)//Shyam
								(((z+j)->ZoneBc+zBc)->ZSlab)->EffScalcoeff = (double *)initializeSlabVarCoeff(1,(z+j),0.01,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						
							if(CgnsBase->flagTurbModel == 9)
							{
								((z+j)->ZoneBc+zBc)->ZSlab->BL			= (BaldwinLomax *)calloc(1,sizeof(BaldwinLomax));
								
								((z+j)->ZoneBc+zBc)->ZSlab->BL->MixingLength 	= (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
								((z+j)->ZoneBc+zBc)->ZSlab->BL->mu_ti		= (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
								((z+j)->ZoneBc+zBc)->ZSlab->BL->mu_to		= (double *)initializeSlabVarCoeff(1,(z+j),0.0,varCoeffCount,dir,IfacePtr->Point_Range);
							}
							else
							{
								((z+j)->ZoneBc+zBc)->ZSlab->BL = NULL;
							}
						}
					}
				}
			}
			else	//If the Zone is SOLID
			{
				((z+j)->Flowsol)->flagTurb=0;
				((z+j)->Flowsol)->Tcoeff = NULL;
				((z+j)->Flowsol)->Vorticity = NULL;
				(z+j)->ZGeom->minWallDistance = NULL;

				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							((z+j)->Iface+iface)->ZSlab->minWallDistance = NULL;
							((z+j)->Iface+iface)->ZSlab->BL = NULL;
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
							((z+j)->ZoneBc+zBc)->ZSlab->minWallDistance = NULL;
							((z+j)->ZoneBc+zBc)->ZSlab->BL = NULL;
						}
					}
				}
			}
		}
	}
	else 
	if(CgnsBase->flagTurb==0)
	{
		//This is must for Conjugate Heat Transfer Case
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->flagTurb=0;
			((z+j)->Flowsol)->Tcoeff = NULL;
			((z+j)->Flowsol)->Vorticity = NULL;
			(z+j)->ZGeom->minWallDistance = NULL;
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("\nNumber of Interfaces = %d\n", numiface);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->minWallDistance = NULL;
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->minWallDistance = NULL;
					}
				}
			}
		}
	}
	//initialize_solidification(CgnsBase,&count);
	initialize_TwoPhase(CgnsBase,&count);//Common function called for solidification-melting and flow boiling problems..created by anji
	
	//Nullyfying memory for Pressure head, water content, hydraulic conductivity-- added for groundwater, after upload
	i=0; 
	fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		(z+j)->Difftag=0;//This variable corresponds to ground water flows.  Since this is nowhere assigned a value for non-ground water flows it is assigned here.
		
		//Memory allocation in the Flow Solution Node For Pressure
		((z+j)->Flowsol)->H= NULL;
		((z+j)->PrevFlowsol)->H= NULL;
		 
		((z+j)->Flowsol)->WaterContent= NULL;
		((z+j)->PrevFlowsol)->WaterContent= NULL;
		
		((z+j)->Flowsol)->SatHydConductivity= NULL;
		((z+j)->Flowsol)->UnSatPermeability=  NULL;
		
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
					//(((z+j)->Iface+iface)->ZSlab)->H = NULL;
					(((z+j)->Iface+iface)->ZSlab)->WaterContent = NULL;
					(((z+j)->Iface+iface)->ZSlab)->SatHydConductivity = NULL;
					(((z+j)->Iface+iface)->ZSlab)->UnSatPermeability=  NULL;
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->H = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->WaterContent = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->SatHydConductivity = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->UnSatPermeability=  NULL;
				}
			}
		}
	}

	vector=scal+count;
	printf("\nTotal no of Solvable Variables:%d",vector);



	if(CgnsBase->flagTurb == 0)
		printf("\nTotal no of Scalars:%d",scal);
	else
	{
		if(CgnsBase->flagTurbModel == 1 || CgnsBase->flagTurbModel == 3 || CgnsBase->flagTurbModel == 4 || CgnsBase->flagTurbModel == 5)
			printf("\nTotal no of Scalars:%d, Scalar1 = k and Scalar2 = Epsilon",scal);
		if(CgnsBase->flagTurbModel == 2 || CgnsBase->flagTurbModel == 6)
			printf("\nTotal no of Scalars:%d, Scalar1 = k and Scalar2 = Omega",scal);
		
		if(CgnsBase->flagTurbModel == 8)
		{
			printf("\nTotal no of Scalars:%d, Scalar1 = Turbulent viscosity",scal);
		} else if (CgnsBase->flagTurbModel == 9) {
			printf("\nTotal no of Scalars:%d",scal);
		} 

	}
	//printf("\nTotal no of Scalars:%d",scal);


	CgnsBase->total_var    = vector;//Variable Tag has to be stored bhans
	for(i=0;i<vector;i++)
		printf("\nThe variable tag at %d\t is %d",i,CgnsBase->VARTAG[i]);
	
	fclose(var);
	fclose(var2);
        fclose(var3);
	fclose(convcri);
	fclose(pressConvcri);
	fclose(omega);	
        
        initialize_monitors(CgnsBase);

	return(vector);
}
