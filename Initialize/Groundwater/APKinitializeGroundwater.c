#include "APKinitializeGroundwater.h"

/**
The Function initialize() asks for which variables memory should be allocated and creates memory only for those
variables.

This function inturn calls initialize_variable() and initializeslab() for each variable to create the necessary memory space for that 
variable in Flowsolution structure and also memory for that variable in slabs. The functions initializeFlux(),
initializeVarCoefficient() and initializeSlabVarCoeff() are used to initialize memory for volume flux and 
coefficients of variables like density viscosity thermal conductivity etc. This function also opens the files
variablelist.dat and coefficients.dat and writes all the user entered data like what variables to be evaluated and
what are their initialization values etc. The information stored in these files is read by the initializerestart
to allocate the same amount of memory and exactly the same initialization values during restart.
*/
int initializeGroundwater(Base *CgnsBase)
{ 
	FILE *var,*var2,*convcri,*pressConvcri,*omega;;
	int i,j,k,p,flag,cntVariable,num_zone,count=0,vector,cntVar,scal;/* num_zone is the number of zones*/
	int numiface,numBc,iface,zBc,dir,varCoeffCount=0,Tchoice;//defo for groundwater
	Zone *z;
	double constant;
	int conv;
	int flag1;
	double *a,L,T;
	
	void *ptr;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *ptrZoneNeighbor;
	char tag,name[33];
	double con[9];
	
	int m,c,e;
	double theta_r,theta_s,alpha,n,a1,b1,d1,c1;
	
	var=fopen("Settings/variablelist.dat","a");
	var2=fopen("Settings/coefficients.dat","a");
	pressConvcri=fopen("Settings/pressureconvergCriteria.dat","a");
	convcri=fopen("Settings/convergenceCriteria.dat","r");
	omega=fopen("Settings/omega.dat","r");
	
	if (convcri==NULL || pressConvcri==NULL || var==NULL || var2==NULL || omega==NULL)
	{
		printf("File opening if failed in APKinitialize.c\nExiting...........\n");
		exit(0);
	}
	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones;
	
	//CgnsBase->BaseParameters=(Parameters *)malloc(sizeof(Parameters));
	CgnsBase->BaseParameters->flagGravity = 0;
	CgnsBase->BaseParameters->flagBuoyancy = 0;
	CgnsBase->BaseParameters->flagEnergySource = 0;
	CgnsBase->BaseParameters->flagDissipation = 0;
	CgnsBase->BaseParameters->TimeStepingSchemeTag=0;

/*###########################################################################################*/	
	printf("\n########### ENTER THE FOLLOWING BASE PARAMETERS ############");
	
	printf("\n If you want to solve the scalar with given velocity field Enter[1]");
	
	printf("\n If you do not want to solve the scalar \n or want to solve scalar from calculated velocity field Enter[0]  ");
	scanf("%d",&e);
	CgnsBase->flagVelocity = e;
	
	if(CgnsBase->flagVelocity==0)
	{
		printf("\nEnter Convergence criteria  for Pressure Head Loop:");
		scanf("%lf",&constant);
		CgnsBase->BaseParameters->HeadLoopConv = constant;
		fprintf(pressConvcri,"%12.8e\n",constant);//it store value of either pressure convergence or head convergence
	}
	
	printf("Do you want to enter the Time step[Enter 1] or let the Code fix it according to Stability Constraints[Enter 0]");
	scanf("%d",&Tchoice);
	CgnsBase->BaseParameters->TstepChoice=Tchoice;
	
	if(Tchoice==1)
	{
		printf("\nEnter TimeStep:");
		scanf("%lf",&constant);
		CgnsBase->BaseParameters->TimeStep = constant;
	}
	else
	{
		CgnsBase->BaseParameters->TimeStep =0.001;
	}
	
	CgnsBase->BaseParameters->Turb->DualTimeStep = CgnsBase->BaseParameters->TimeStep / CgnsBase->BaseParameters->Turb->numDualTimeIter;
	
	CgnsBase->flagSteadyState = 0;/*0->for steady state, 1->for Unsteady state*/
	/*Allocation of memories for various Convergence Criteria*/
	CgnsBase->BaseParameters->convergenceCriterion=(double *)malloc(12*sizeof(double));
	for(cntVariable=0;cntVariable<12;cntVariable++)
		 fscanf(convcri,"%lf",(CgnsBase->BaseParameters->convergenceCriterion+cntVariable));
	
	CgnsBase->BaseParameters->OMEGA=(double *)malloc(14*sizeof(double));
	for(cntVariable=0;cntVariable<14;cntVariable++)
		fscanf(omega,"%lf",(CgnsBase->BaseParameters->OMEGA+cntVariable));
	/*the first element is the criterion for velocities, second is for temperature and rest 10 are for the 10 scalars*/
	CgnsBase->BaseParameters->Time=0.0;
	
/*#####################################################################################################################*/	
	for(j=0;j<num_zone;j++)
	{
		/* Allocation of memory for Flowsolution */
		(z+j)->Flowsol=(FlowSolution*)malloc(sizeof(FlowSolution));
		(z+j)->PrevFlowsol=(PrevFlowSolution*)malloc(sizeof(PrevFlowSolution));
	}
	/*Entering the value 1 for this question would enable the code 
	to allocate memory in the data structure for that particular variable  */
	//**********************************Entering value of porosity*****************//////////
	
		
	for(j=0;j<num_zone;j++)
	{ 
		varCoeffCount =19;  /*15->Permeability*/
		constant=1;
		p=1;
		((z+j)->Flowsol)->UnSatPermeability= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
		((z+j)->Flowsol)->GradTheta= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
		
		numiface=(z+j)->Number_of_Interfaces;
		
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					(((z+j)->Iface+iface)->ZSlab)->UnSatPermeability=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					/*Memory allocation in the Interface Slab X Direction Velocity*/
					//printf("I am going to zoneBc\n");
					(((z+j)->ZoneBc+zBc)->ZSlab)->UnSatPermeability = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					//printf("I am leaving zoneBc\n");
				}
			}
		}
	}
		 
		
	//printf("creating memory for water content");
	i=1;
	fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		//  printf("enter initial value of mositure content");
		// scanf("%lf",&constant);
		constant=0.25;
		/*Memory allocation in the Flow Solution Node For Pressure*/
		((z+j)->Flowsol)->WaterContent= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->WaterContent= (double *)initialize_variable(i,(z+j),constant);

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
				(((z+j)->Iface+iface)->ZSlab)->WaterContent = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				(((z+j)->ZoneBc+zBc)->ZSlab)->WaterContent = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}
		
	//###################################################&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
	
	if((CgnsBase->flagSaturation==1)&&(CgnsBase->flagVelocity==0))
	{
		printf("\n MEDIUM IS HOMOGENOUS (Enter 0) or LAYERED SOIL (Enter 1):");//
		scanf("%d",&flag);
		if(flag==0) 
		{      
			printf("Enter the Type of Model for (K-Theta-PressureHead) porous medium (soil)");
			printf("\n ENTER 1 for VGM Model  &&  ENTER 2 for Haverkamp Model");
			scanf("%d",&flag1);

			for(j=0;j<num_zone;j++)
				((z+j)->flagtype_model)=flag1;

			if(flag1==1)// for VGM model
			{
				printf("Enter the value of various soil parameters for VGM model");

				printf("\n Enter the value of theta_residual(Residual water content)");
				scanf("%lf",&theta_r);

				printf("\nEnter the value of theta_saturation(Saturated water content which is equal to porosity)");
				scanf("%lf",&theta_s);

				printf("\nEnter the value of alpha(model parameter)");
				scanf("%lf",&alpha);	

				printf("\nEnter the value of  n(model parameter)");
				scanf("%lf",&n);

				for(j=0;j<num_zone;j++)
				{
					((z+j)->Theta_residual)=theta_r;
					((z+j)->Theta_saturation)=theta_s;

					((z+j)->Alpha)=alpha;
					((z+j)->n)=n;
				}

			}
			else if(flag1==2) // for Haverkamp model
			{

				printf("Enter the value of various soil parameters for Haverkamp model");

				printf("\n Enter the value of theta_residual(Residual water content)");
				scanf("%lf",&theta_r);

				printf("\nEnter the value of theta_saturation(Saturated water content which is equal to porosity)");
				scanf("%lf",&theta_s);

				printf("\nEnter the value of A, B, C, D(model parameter) insequence");
				scanf("%lf %lf %lf %lf",&a1,&b1,&c1, &d1);

				for(j=0;j<num_zone;j++)
				{
					((z+j)->Theta_residual)=theta_r;
					((z+j)->Theta_saturation)=theta_s;

					((z+j)->A)=a1;
					((z+j)->B)=b1;
					((z+j)->C)=c1;
					((z+j)->D)=d1;
				}
			}
		}
		if(flag==1)
		{
			for(j=0;j<num_zone;j++)
			{ 
				printf("Zone Name is %s\n",(z+j)->zonename);
				printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				printf("Enter the Type of Model for (K-Theta-PressureHead) porous medium (soil)");
				printf("\n ENTER 1 for VGM Model  &&  ENTER 2 for Haverkamp Model");
				scanf("%d",&flag1);

				printf("Enter the value of various soil parameters for VGM model");

				if(flag1==1)
				{
					printf("\n Enter the value of theta_residual(Residual water content)");
					scanf("%lf",&theta_r);

					printf("Enter the value of theta_saturation(Saturated water content which is equal to porosity)");
					scanf("%lf",&theta_s);

					printf("Enter the value of alpha(model parameters)");
					scanf("%lf",&alpha);

					printf("Enter the value of  n(model parameters)");
					scanf("%lf ",&n);

					((z+j)->Theta_residual)=theta_r;
					((z+j)->Theta_saturation)=theta_s;

					((z+j)->Alpha)=alpha;
					((z+j)->n)=n;
				}
				if(flag1==2)
				{
					printf("Enter the value of various soil parameters for Haverkamp model");

					printf("\n Enter the value of theta_residual(Residual water content)");
					scanf("%lf",&theta_r);

					printf("\nEnter the value of theta_saturation(Saturated water content which is equal to porosity)");
					scanf("%lf",&theta_s);

					printf("\nEnter the value of A, B, C, D(model parameter)");
					scanf("%lf %lf %lf %lf", &a1,&b1,&c1,&d1);


					((z+j)->Theta_residual)=theta_r;
					((z+j)->Theta_saturation)=theta_s;

					((z+j)->A)=a1;
					((z+j)->B)=b1;
					((z+j)->C)=c1;
					((z+j)->D)=d1;
				}
			}
		}
	}
	//@@}@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@////////////	
	
	if( CgnsBase->Number_of_zones > 1)
	{
		printf("\n Do you want to specify Porosity  for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1):");
		scanf("%d",&flag);
		fprintf(var2,"%d\n",flag);
	} else { flag = 0;}
	if(flag==1) 
	{
	        for(j=0;j<num_zone;j++)
		{ 
			printf("Zone Name is %s\n",(z+j)->zonename);
			printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
			
			printf("\nEnter the constant value of Porosity :");
			scanf("%lf",&constant);
			
			(z+j)->Porosity =constant;               
	        }
	}
	else 
	{
		printf("\nEnter the constant value of Porosity :");
		scanf("%lf",&constant);
	
		for(j=0;j<num_zone;j++)
		{ 
			(z+j)->Porosity =constant;  
		} 
	}
	
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@////////////	
	
	if(CgnsBase->flagVelocity==0)
	{
		if(CgnsBase->flagdeformable ==0)
		{
			for(j=0;j<num_zone;j++)
			(z+j)->Flowsol->Storativity =0.0;
		}
		else
		{
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify Storativity  for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1):");
				scanf("%d",&flag);
				fprintf(var2,"%d\n",flag);
			} else { flag = 0;}
			if(flag==1) 
			{
				for(j=0;j<num_zone;j++)
				{ 
					printf("Zone Name is %s\n",(z+j)->zonename);
					printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));

					printf("\nEnter the constant value of Storativity:");
					scanf("%lf",&constant);
					(z+j)->Flowsol->Storativity =constant;
				}
			}
			else if(flag==0) 
			{
				printf("\nEnter the constant value of Storativity:");
				scanf("%lf",&constant);

				for(j=0;j<num_zone;j++)
				(z+j)->Flowsol->Storativity =constant;
			}
		}
		
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@////////////	
		printf("\nCreating memory for Saturated Hydraulic conductivity:");
		if( CgnsBase->Number_of_zones > 1)
		{
			printf("\n Do you want to specify Hydraulic conductivity for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1):");
			scanf("%d",&flag);
			fprintf(var2,"%d\n",flag);
		} else { flag = 0;}
		if(flag==0) 
		{       printf("\n If Saturated Hydraulic conductivity is a CONSTANT, Enter 0");
			printf("\n If Saturated Hydraulic conductivity is a VARIABLE, Enter 1:");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);
			varCoeffCount = 17;  /*1->Permeability*/
			if(p==0)
			{
				printf("\nEnter the constant value of SatHydConductivity:");
				scanf("%lf",&constant);
				fprintf(var2,"%lf\n",constant);
			}
			else
				constant=0.0;
	
			for(j=0;j<num_zone;j++)
			{       
				if(CgnsBase->flagSaturation==0) 
				(z+j)->Flowsol->flagTestSaturation=0;// for saturated
				else
				(z+j)->Flowsol->flagTestSaturation=1;
				
				((z+j)->Flowsol)->SatHydConductivity= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				
				numiface=(z+j)->Number_of_Interfaces;
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							
							(((z+j)->Iface+iface)->ZSlab)->SatHydConductivity=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							
						}
					}
	
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						
						(((z+j)->ZoneBc+zBc)->ZSlab)->SatHydConductivity = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
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
				printf("\n If Saturated Hydraulic conductivity is a CONSTANT, Enter 0");
				printf("\n If Saturated Hydraulic conductivity is a VARIABLE, Enter 1:");
				scanf("%d",&p);
				fprintf(var2,"%d\n",p);
				varCoeffCount = 17;  
				
		
				if(CgnsBase->flagSaturation==0) 
					(z+j)->Flowsol->flagTestSaturation=0;
				else
					(z+j)->Flowsol->flagTestSaturation=1;
				
				if(p==0)
				{    
					(z+j)->Flowsol->flagDiffcoeff=p;
					printf("\nEnter the constant value of Saturated Hydraulic conductivity:");
					scanf("%lf",&constant);
					fprintf(var2,"%lf\n",constant);
				}
				else
				{
					(z+j)->Flowsol->flagDiffcoeff=p;
					constant=0.0;
				}
				
				((z+j)->Flowsol)->SatHydConductivity= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
				
				numiface=(z+j)->Number_of_Interfaces;
				//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
							
							(((z+j)->Iface+iface)->ZSlab)->SatHydConductivity=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							
						}
					}
	
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to zoneBc\n");
							(((z+j)->ZoneBc+zBc)->ZSlab)->SatHydConductivity = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						
						}
					}
				}
			}
		}
	}
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$///
	
	if(CgnsBase->flagVelocity==0)
	{
		printf("\n Creating memory for Pressure Head");
		CgnsBase->VARTAG[count]=19;//  for groundwater vartag to be decided Variable Tag has to be stored bhans 
		count=count+1;
		printf("\nEnter the Initialization value of Pressure Head:");
		scanf("%lf",&constant);	
	}
	else
		constant=1.0;
	
	i=1;
	fprintf(var,"%d\n",i);
		
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->H= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->H= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->HeadCoeff= (double *)malloc(sizeof(double));
		*(((z+j)->Flowsol)->HeadCoeff)=1.0;

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
				/*Memory allocation in the Interface Slab for Pressure*/
					(((z+j)->Iface+iface)->ZSlab)->H = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->H = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}
		
	
	// ######$$$$$$$parameter to be used in  contaminant transport for calculating diffusion ####////////
	if(CgnsBase->BaseParameters->flag2d==1) 
		c=4;	//Size of Delta inverse matrix for 2d
	else
		c=9;	//Size of Delta inverse matrix for 3d
		
	for(j=0;j<num_zone;j++)
	{
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					(((z+j)->Iface+iface)->ZSlab)->Delta = (double *)initializeslab_Delta(i,c,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->Delta = (double 	*)initializeslab_Delta(i,c,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}
	
	for(j=0;j<num_zone;j++)
	{
		for(iface=0;iface<6;iface++)
		{
			ptr = (void *)FindInterfaceBCPtrMod(z+j,iface+1,name,&tag);
			
			if(tag=='i')
			{
				facePtr=(Interface *)ptr;
				
				ptrZoneNeighbor = (Zone *)FindZonePtr(name,CgnsBase);
				
				FillerInterfaceDelta(c,facePtr,ptrZoneNeighbor,iface);
			}
			else if(tag=='b')
			{         
				bcPtr=(ZoneBoundary *)ptr;
				
				updateBCDeltaInverse(c,z+j,bcPtr,iface+1);
			}
		}
	}
		
// ###############################################////////
	for(j=0;j<num_zone;j++)
	{
		/*Memory allocation for Gradient of H in X,Y,Z directions */
		((z+j)->Flowsol)->GradU= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->GradV= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->GradW= (double *)initialize_variable(i,(z+j),constant);
	}
	
//##############################creating memory for u velocity#$$$$$$$$$$$$$$$$$$
	
	printf("\n Creating memory for U velocity(Darcy velocity)");
	i=1;
	fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		if(CgnsBase->flagVelocity==1)
		{
			printf("\nEnter the Initialization value of U velocity(Darcy velocity)");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		}//Memory allocation in the Flow Solution Node X Direction Velocity//
		else
		{
			constant=0.0;
			fprintf(var2,"%lf\n",constant);
		}
		//Memory allocation in the Flow Solution Node X Direction Velocity//
		((z+j)->Flowsol)->u= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->u= (double *)initialize_variable(i,(z+j),constant);
		//for groundwater flow no need for  umass but has been allocated memory as volumeflux cal need it in 
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
					(((z+j)->Iface+iface)->ZSlab)->u=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->umass = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->u = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->umass = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}
	
	/*The memory allocation for Y Direction Velocity for Flowsol node and Slabs at boundary and Interface */
	printf("\n Creating memory for V velocity (Darcy velocity)");
	i=1;
	fprintf(var,"%d\n",i);
	
	for(j=0;j<num_zone;j++)
	{
		if(CgnsBase->flagVelocity==1)
		{
			printf("\nEnter the Initialization value of V velocity(Darcy velocity)");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		}//Memory allocation in the Flow Solution Node X Direction Velocity//
		else
		{
			constant=0.0;
			fprintf(var2,"%lf\n",constant);
		}

		((z+j)->Flowsol)->v= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->v= (double *)initialize_variable(i,(z+j),constant);

		//for groundwater flow no need for  vmass but has been allocated memory as volumeflux cal need it in 
		((z+j)->Flowsol)->vmass= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->vmass= (double *)initialize_variable(i,(z+j),constant); 

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
					(((z+j)->Iface+iface)->ZSlab)->v = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					(((z+j)->Iface+iface)->ZSlab)->vmass = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->v = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->vmass = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}
	
	// ########################creating memory for w velocity########################////////
	
	if(CgnsBase->BaseParameters->flag2d == 0) 
	printf("\n Creating memory for W velocity (Darcy velocity) ");
	i=1;
        fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{       
		if((CgnsBase->BaseParameters->flag2d == 0)&&(CgnsBase->flagVelocity==1)) 
		{
			printf("\nEnter the Initialization value of W velocity(Darcy velocity)");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		}//Memory allocation in the Flow Solution Node X Direction Velocity//
		else
		{
			constant=0.0;
			fprintf(var2,"%lf\n",constant);
		}

		((z+j)->Flowsol)->w= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->w= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->wmass= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->wmass= (double *)initialize_variable(i,(z+j),constant);

		numiface=(z+j)->Number_of_Interfaces;
		//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
				(((z+j)->Iface+iface)->ZSlab)->w=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
				(((z+j)->Iface+iface)->ZSlab)->wmass=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					(((z+j)->ZoneBc+zBc)->ZSlab)->w = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					(((z+j)->ZoneBc+zBc)->ZSlab)->wmass = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
				}
			}
		}
	}
	
	 // end for loop for groundwater
	 

	printf("\nInitialising ConjugateRms to zero\n");
	for(j=0;j<num_zone;j++)
	{
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
	
	printf("\n Creating memory for MassFlux");
	i=1;
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->Massflux= (double *)initializeFlux(i,(z+j));
		((z+j)->PrevFlowsol)->Massflux= (double *)initializeFlux(i,(z+j)); 
	}
	printf("\n Initializing default value of flagTurb to 0 (Laminar case)\n");
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->flagTurb=0;
	}

	for(j=0;j<num_zone;j++)
		(z+j)->Difftag=0;
/*****************************MEMORY ALLOCATION FOR SCALARS*******************************************************/	
	/*Enter the number of the Scalars to be solved */
	
	printf("\n You  want to solve for contaminants transport -Enter number of Scalars(solute or contaminants) to be solved or if do not want to solve enter[0] [Max. of 10 Scalars is allowed]: ");
	printf("\n");
	scanf("%d",&scal);
	while(scal>10)
	{
		printf("\nMax allowed scalars is 10");
		printf("\nPlease enter the Number of Scalars(<= 10):");
		scanf("%d",&scal);
	}
	fprintf(var,"%d\n",scal);
	
	if(scal>=1)
	{
		//#####################################################################################################################//
		printf("Please Enter the time steping scheme\n[CrankNicolson->1/BTCS->0]");
		scanf("%d",&j);
		if(j==1)
			CgnsBase->BaseParameters->TimeStepingSchemeTag=1;
		//#####################################################################################################################//	
		printf("Enter the Convection Scheme(Deferred Corrected) QUICK[enter=1] or Central Difference[enter=0]\n");
		scanf("%d",&conv);
		CgnsBase->BaseParameters->ConvectionScheme=conv;

		//#####################################################################################################################//

		if(CgnsBase->BaseParameters->ConvectionScheme==0)
		{
			printf("\nEnter Ratio of CentralDiff(1->pure Centraldiff) to Upwind(0->pure Upwind):");
			scanf("%lf",&constant);
			//fprintf(var2,"%lf\n",constant);
			CgnsBase->BaseParameters->RatioCenUP = constant;
		}
		else
		{
			printf("\nEnter Ratio of QUICK(1->pure QUICK) to Upwind(0->pure Upwind):");
			scanf("%lf",&constant);
			//fprintf(var2,"%lf\n",constant);
			CgnsBase->BaseParameters->RatioCenUP = constant;
		}
	
	//#####################################################################################################################//
	}
	
	for(k=0;k<scal;k++)
	{
		for(j=0;j<num_zone;j++)
			(z+j)->Difftag=1;
		
		i=1;
		CgnsBase->VARTAG[count+k]=600+k;//Variable Tag has to be stored bhans	
		printf("\nCreating memory for Scalar Coefficient[%d]:",k+1);
		if( CgnsBase->Number_of_zones > 1)
		{
			printf("\n Do you want to specify Scalar Coefficient[%d] for the WHOLE DOMAIN (Enter 0) or ZONE BY ZONE (Enter 1):",k+1);
			scanf("%d",&flag);
			fprintf(var2,"%d\n",flag);
		} else { flag = 0;}
	//-------------------------------------------------------------------------------------------/////	
		if(flag==0) /*Initializing coeffcients for whole domain at a time*/
		{
			printf("\n  Enter '0' ,IF you want to give constant value of Dispersion Coefficient (This    value is   constant for given velocity field )");  
			printf("\nEnter '1' ,IF you want to give constant value of Dispersivity(longitudinal and transversal) (This value is dependent  for given velocity field )");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);
			a=(double *)malloc(9*sizeof(double)); 
			varCoeffCount = 4+k+1;   /*4->SpecificHeat,5->Scalcoeff1,6->Scalcoeff2,.....14->Scalcoeff10.*/
			if(p==0)
			{
				printf("\nEnter the constant value of Dispersion Coeffcient (given the component for Dispersion in sequence Dxx  Dxy Dxz, Dyx Dyy Dyz, Dzx Dzy Dzz:");
				for(m=0;m<9;m++)
				{
					scanf("%lf",&con[m]);
					*(a+m)=con[m];
				}
				
				for(j=0;j<num_zone;j++)
				{
					(z+j)->CalDispersion=p;   
				}
			}
			else
			{
				constant=0.0;
				printf("Enter the value of longitudinal dispersivity of medium");
				scanf("%lf",&L);

				printf("Enter the value of transversal dispersivity of medium ");
				scanf("%lf",&T);
				
				for(m=0;m<9;m++)
					*(a+m)=0.0;
				
				for(j=0;j<num_zone;j++)
				{
					((z+j)->alphaL)=L;
					((z+j)->alphaT)=T;
					(z+j)->CalDispersion=p;
				}
			}
			
			for(j=0;j<num_zone;j++)
			{
				((z+j)->Flowsol)->Scalcoeff[k]= (double *)initializeDispersionCoeff(p,(z+j),a,varCoeffCount);
				((z+j)->PrevFlowsol)->Scalcoeff[k]= (double *)initializeDispersionCoeff(p,(z+j),a,varCoeffCount);
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
						(((z+j)->Iface+iface)->ZSlab)->Scalcoeff[k]=(double *)initializeSlabDispersionCoeff(p,(z+j),a,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}
					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Scalcoeff[k] = (double *)initializeSlabDispersionCoeff(p,(z+j),a,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
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
				printf("Enter '0' ,IF you want to give constant value of Dispersion Coefficient (This value is constant for given velocity field )");  
				printf("Enter '1' ,IF you want to give constant value of Dispersivity(longitudinal and transversal) (This value is dependent  for given velocity field )");
				scanf("%d",&p);
				fprintf(var2,"%d\n",p);

				(z+j)->CalDispersion=p;

				varCoeffCount = 4+k+1;	
				if(p==0)
				{  
					printf("\nEnter the constant value of Dispersion Coeffcient (given the component for Dispersion in sequence Dxx Dxy Dxz, Dyx Dyy Dyz, Dzx Dzy Dzz:");
					for(m=0;m<9;m++)
					{ 
						scanf("%lf",&con[m]);
						*(a+m)=con[m];
					}
				}
				else
				{
					printf("Enter the value of longitudinal dispersivity of medium");
					scanf("%lf",&L);

					printf("Enter the value of transversal dispersivity of medium ");
					scanf("%lf",&T);

					(z+j)->alphaL=L;
					(z+j)->alphaT=T;

					for(m=0;m<9;m++)
						*(a+m)=0.0;
					
				} 
				((z+j)->Flowsol)->Scalcoeff[k]= (double *)initializeDispersionCoeff(p,(z+j),a,varCoeffCount);
				((z+j)->PrevFlowsol)->Scalcoeff[k]= (double *)initializeDispersionCoeff(p,(z+j),a,varCoeffCount);
				((z+j)->Flowsol)->flagScalcoeff[k]=p;
				numiface=(z+j)->Number_of_Interfaces;

				numBc=(z+j)->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{	if(((z+j)->Iface+iface)->Direction==dir)
						{
							(((z+j)->Iface+iface)->ZSlab)->Scalcoeff[k]=(double *)initializeSlabDispersionCoeff(p,(z+j),a,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							//printf("I am leaving Iface\n");
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Scalcoeff[k] = (double *)initializeSlabDispersionCoeff(p,(z+j),a,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//printf("I am leaving zoneBc\n");
						}
					}
				}
			}  
		}
		
	//---------------------- Gradient of scalar in X, Y, Z, direction------------------------------------------------/////
		for(j=0;j<num_zone;j++)
		{
			/*Memory allocation in the Flow Solution Node For Scalars*/
			((z+j)->Flowsol)->GradScalX[k]= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->GradScalY[k]= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->GradScalZ[k]= (double *)initialize_variable(i,(z+j),constant);
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						/*Memory allocation in the Interface Slab for Scalar*/
						(((z+j)->Iface+iface)->ZSlab)->GradScalX[k]= (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->GradScalY[k]= (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->GradScalZ[k]= (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{	
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Boundary Slab for Scalar*/
						(((z+j)->ZoneBc+zBc)->ZSlab)->GradScalX[k]= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->GradScalY[k]= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->GradScalZ[k]= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
		
		for(j=0;j<num_zone;j++)
		{
			printf("\nEnter the Initialization value of Scal[%d]:",k);
			scanf("%lf",&constant);                         
			fprintf(var2,"%lf\n",constant);
			
			/*Memory allocation in the Flow Solution Node For Scalars*/
			((z+j)->Flowsol)->Scal[k]= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->Scal[k]= (double *)initialize_variable(i,(z+j),constant);

			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						/*Memory allocation in the Interface Slab for Scalar*/
						(((z+j)->Iface+iface)->ZSlab)->Scal[k]= (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->ScalOld[k]= (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{	
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						/*Memory allocation in the Boundary Slab for Scalar*/
						(((z+j)->ZoneBc+zBc)->ZSlab)->Scal[k]= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->ScalOld[k]= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->GradScal[k]= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
	} /*END OF LOOPING OVER SOLVABLE SCALARS*/
	
	printf("Enter The required Time after which solution is Desired");
	scanf("%lf",&constant);
	CgnsBase->BaseParameters->RequiredTime= constant;	

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
				{	
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Scalcoeff[k]= NULL;
						(((z+j)->Iface+iface)->ZSlab)->Scal[k]= NULL;
						(((z+j)->Iface+iface)->ZSlab)->ScalOld[k]= NULL;
					}
				}
				for(zBc=0;zBc<numBc;zBc++)
				{	
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
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
	for(j=0;j<num_zone;j++)
		((CgnsBase->CgnsZone+j)->Flowsol)->flagTurb=0;
		
/**************************** START OF MAKING NULL for VARIABLES WHICH ARE NOT SOLVED *******************************************************/		
		/*The memory allocation for Pressure for Flowsol node and Slabs at boundary and Interface */
	//printf("\n Making Pressure NULL for Ground Water flow case, not PressureHead\n");
	i=0; // Not solving for normal Pressure.
	fprintf(var,"%d\n",i);
	//CgnsBase->VARTAG[count]=4;//Variable Tag has to be stored bhans
	//count=count+1;
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		/*Memory allocation in the Flow Solution Node For Pressure*/
		((z+j)->Flowsol)->P= NULL;
		((z+j)->PrevFlowsol)->P= NULL;
		((z+j)->Flowsol)->PressureCoeff= NULL;
		((z+j)->Flowsol)->T= NULL;
		((z+j)->PrevFlowsol)->T= NULL;
		((z+j)->Flowsol)->Potential= NULL;
		((z+j)->PrevFlowsol)->Potential = NULL;
		((z+j)->Flowsol)->enthalpy = NULL;
		((z+j)->Flowsol)->f = NULL;
		((z+j)->Flowsol)->Latent_Heat_Content = NULL;
	
		((z+j)->Flowsol)->Vorticity = NULL;
		((z+j)->Flowsol)->vorticity_X = NULL;
		((z+j)->Flowsol)->vorticity_Y = NULL;
		((z+j)->Flowsol)->vorticity_Z = NULL;
		
		((z+j)->Flowsol)->Gamma = NULL;
		((z+j)->Flowsol)->CrossDiff = NULL;
		((z+j)->Flowsol)->MixFrac = NULL;
		((z+j)->Flowsol)->FuelMassFrac = NULL;
		((z+j)->Flowsol)->OxiMassFrac = NULL;
		((z+j)->Flowsol)->InMassFrac[10] = NULL;
		((z+j)->Flowsol)->ProdMassFrac = NULL;
		((z+j)->Flowsol)->EfieldX = NULL;
		((z+j)->Flowsol)->EfieldY = NULL;
		((z+j)->Flowsol)->EfieldZ = NULL;
		((z+j)->Flowsol)->Electricflux = NULL;
		((z+j)->Flowsol)->RhoIon = NULL;
		((z+j)->Flowsol)->DiElecConst = NULL;
		((z+j)->Flowsol)->MobIon = NULL;
		((z+j)->Flowsol)->ChargeDiffcoeff = NULL;
		((z+j)->Flowsol)->PotentialMHD = NULL;
		((z+j)->Flowsol)->PotentialMHDCoeff = NULL;
		((z+j)->Flowsol)->MfieldX = NULL;
		((z+j)->Flowsol)->MfieldY = NULL;
		((z+j)->Flowsol)->MfieldZ = NULL;
		((z+j)->Flowsol)->PotentialMHDpoissonRHS = NULL;
		((z+j)->Flowsol)->Eleconductivity = NULL;
		((z+j)->Flowsol)->ElecurrentdensityX = NULL;
		((z+j)->Flowsol)->ElecurrentdensityY = NULL;
		((z+j)->Flowsol)->ElecurrentdensityZ = NULL;
		((z+j)->Flowsol)->UcrossBX = NULL;
		((z+j)->Flowsol)->UcrossBY = NULL;
		((z+j)->Flowsol)->UcrossBZ = NULL;

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	
				if(((z+j)->Iface+iface)->Direction==dir)
				{
				/*Memory allocation in the Interface Slab for Pressure*/
					(((z+j)->Iface+iface)->ZSlab)->P = NULL;
					(((z+j)->Iface+iface)->ZSlab)->T = NULL;
					(((z+j)->Iface+iface)->ZSlab)->Potential = NULL;
					(((z+j)->Iface+iface)->ZSlab)->enthalpy = NULL;
					(((z+j)->Iface+iface)->ZSlab)->f = NULL;
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{	
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				/*Memory allocation in the Boundary Slab for Pressure*/
					(((z+j)->ZoneBc+zBc)->ZSlab)->P = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->T = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->Potential = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpy = NULL;
					(((z+j)->ZoneBc+zBc)->ZSlab)->f = NULL;
				}
			}
		}
	}
/*****************************END OF MAKING NULL *******************************************************/
	
	vector=scal+count;
	printf("\nTotal no of Solvable Variables:%d",vector);
	printf("\nTotal no of Scalars:%d",scal);
	CgnsBase->total_var    = vector;//Variable Tag has to be stored bhans
	fclose(var);
	fclose(var2);
	fclose(convcri);
	fclose(pressConvcri);
	fclose(omega);
	return(vector);
}

