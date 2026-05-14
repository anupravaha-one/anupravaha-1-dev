#include "APKinitializeDisperse.h"
#define Cell_Dim 3
/**
The Function initializeDisperse() asks for which variables memory should be allocated and creates memory only for those variables. These variables are disperse phase flow variables that is either of the particle phase or droplet phase depending on the type of disperse phase choosen by the user. Same as Function initialize(), this function calls initialize_variable() and initializeslab() for disperse phase variable to create the necessary memory space for that
variable in Flowsolution structure and also memory for that variable in slabs. The functions initializeFlux(),
initializeVarCoefficient() and initializeSlabVarCoeff() are used to initialize memory for volume flux of the disperse phase and coefficients of variables like particle or droplet density, liquid droplet viscosity, mass diffusivity of the evaporated vapour, specific heat of liquid droplet etc. This function also opens the files variablelist.dat and coefficients.dat and writes all the user entered data like what variables to be evaluated and what are there initialization values etc. The information stored in these files is read by the initializerestart to allocate the same amount of memory and exactly the same initialization values during restart.
*/

int initializeDisperse(Base *CgnsBase,int cnt)
{
	FILE *var,*var2;
	int i,j,p,num_zone;/* num_zone is the number of zones*/
	int numiface,numBc,iface,zBc,dir,varCoeffCount=0;
	Zone *z;
	double  VolFrac,constantu,constantv,constantw,constantden,constantmu,cons=0.0,constantdia,constantn,constant;

	double *variable1,*variable2;//adhiraj
	int totalcells,k;//adhiraj

	int NI,NJ,NK,a,b,c,index;
	int *ptrSize;
	int flag;
	double CPvapour,epsilondrop;//shash

	var=fopen("Settings/variablelist.dat","a");
	var2=fopen("Settings/coefficients.dat","a");

	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones;
/*###############################################################################################################*///   Disperse PHASE FLOW VARIABLES  problem.
	printf("\nEnter the value of EPSILON (restriction in Volume Fraction) \n");
	printf("\nFor gas-particle flows take EPSILON <= 1e-5 and \nFor gas-droplet flows with evaporation take EPSILON >= 1e-4\n");
	scanf("%lf",&epsilondrop);
	CgnsBase->BaseParameters->EPSILONDrop = epsilondrop;
	/***************************************************************************************************************/

	/*The memory allocation for Density for Flowsol node and Slabs at boundary and Interface*/

		printf("\nCreating memory for Density[kg/m3] of Disperse phase :");
		printf("\n If Density is a CONSTANT, Enter 0");
		printf("\n If Density is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 1;  /*1->Density*/

		if(p==0)
		{
			printf("\nEnter the constant value of Density[kg/m3]for disperse phase :");
			scanf("%lf",&constantden);
			fprintf(var2,"%lf\n",constantden);
		}
		else
 			{
			cons=0.0;constantden=1.0;
			}
		for(j=0;j<num_zone;j++)
		{

		       /*Memory allocation in the Flow Solution Node */
			((z+j)->Flowsol)->RhoDisperse= (double *)initializeVarCoefficient(p,(z+j),constantden,varCoeffCount);
			((z+j)->PrevFlowsol)->RhoDisperse= (double *)initializeVarCoefficient(p,(z+j),constantden,varCoeffCount);

			CgnsBase->BaseParameters->flagRhoDisperse=p;

			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						/*Memory allocation in the Interface Slab */
						(((z+j)->Iface+iface)->ZSlab)->RhoDisperse=(double*)initializeSlabVarCoeff(p,(z+j),constantden,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						//}//printf("I am leaving Iface\n");
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{

					/*Memory allocation in the Interface Slab X Direction Velocity*/
					//printf("I am going to zoneBc\n");
					(((z+j)->ZoneBc+zBc)->ZSlab)->RhoDisperse = (double *)initializeSlabVarCoeff(p,(z+j),constantden,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					//printf("I am leaving zoneBc\n");
					}
				}
			}
		}
//---------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//exit(155);
	/*****************************MEMORY ALLOCATION FOR VISCOSITY*******************************************************/
	if((CgnsBase->flagDroplet==1))
	{
	printf("\nCreating memory for VISCOSITY[kg/m-s] of Disperse phase :");

		printf("\n If Viscosity is a CONSTANT, Enter 0");
		printf("\n If Viscosity is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 2;  /*2->VISCOSITY*/
		if(p==0)
		{     printf("\nEnter the constant value of VISCOSITY[kg/m-s]for disperse phase :");
			scanf("%lf",&constantmu);
			fprintf(var2,"%lf\n",constantmu);
		}
		else
			{
			cons=0.0;constantmu=1.0;
			}

		for(j=0;j<num_zone;j++)
		{
			//need to ask
			((z+j)->Flowsol)->ViscoeffDummyDisperse= (double *)initializeVarCoefficient(p,(z+j),cons,varCoeffCount);
			((z+j)->PrevFlowsol)->ViscoeffDummyDisperse= (double *)initializeVarCoefficient(p,(z+j),cons,varCoeffCount);


			/*Memory allocation in the Flow Solution Node */
			//for(nophase=0;nophase<nop;nophase++)
			//{
			((z+j)->Flowsol)->ViscosityDisperse= (double *)initializeVarCoefficient(p,(z+j),constantmu,varCoeffCount);
			((z+j)->PrevFlowsol)->ViscosityDisperse= (double *)initializeVarCoefficient(p,(z+j),constantmu,varCoeffCount);

			CgnsBase->BaseParameters->flagViscoeffDisperse=p;
			//}
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->ViscoeffDummyDisperse=(double *)initializeSlabVarCoeff(p,(z+j),cons,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to Iface\n");
						//for(nophase=0;nophase<nop;nophase++)
						//{
						(((z+j)->Iface+iface)->ZSlab)->ViscosityDisperse=(double *)initializeSlabVarCoeff(p,(z+j),constantmu,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						 //}
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{

					(((z+j)->ZoneBc+zBc)->ZSlab)->ViscoeffDummyDisperse = (double *)initializeSlabVarCoeff(p,(z+j),cons,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					/*Memory allocation in the Interface Slab */

					(((z+j)->ZoneBc+zBc)->ZSlab)->ViscosityDisperse = (double *)initializeSlabVarCoeff(p,(z+j),constantmu,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					//}
					}
				}
			}
		}

	}
/*****************************END OF MEMORY ALLOCATION FOR VISCOSITY*******************************************************/

/****************************************************************************************************************/

/*************Memory allocation of volume fraction for  disperse phase**********************************************/
	printf("\n Creating memory for volume fraction of Disperse phase");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=65;//Variable Tag has to be stored
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{
			printf("Zone Name is %s\n",(z+j)->zonename);
			printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
			printf("\n Enter the Initialization value of volume fraction  for disperse phase for this zone:");
			scanf("%lf",&VolFrac);
			fprintf(var2,"%lf\n",VolFrac);

			//Mrunalini
			printf("\nSince Normalized Volume-fraction is solved for in Dilute Gas-Droplet flows with evaporation solver, the flowsol and prevflowsol values are initialized to 1.0\n");
			printf("\n Change the value if necessary in Initialize/Disperse/APKInitializeDisperse.c");
			CgnsBase->BaseParameters->InitialVolumeFraction = VolFrac;

			/*Memory allocation in the Flow Solution Node for normalized volume fraction*/
			//Since, we are solving for normalized volume fraction, the flowsol and prevflowsol values are initialized to 1.
			((z+j)->Flowsol)->VarthetaDisperse= (double *)initialize_variable(i,(z+j),1.0);
			((z+j)->PrevFlowsol)->VarthetaDisperse= (double *)initialize_variable(i,(z+j),1.0);

			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->VarthetaDisperse=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,1.0);
						(((z+j)->Iface+iface)->ZSlab)->VarthetaDisperseOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,1.0);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->VarthetaDisperse = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,1.0);
						(((z+j)->ZoneBc+zBc)->ZSlab)->VarthetaDisperseOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,1.0);
					}
				}
			}



		}
	/***********End of Memory allocation of normalized-volume fraction of disperse phase**********************************************/

/**************************************************************************************/
	/*The memory allocation for Particle Diameter for Flowsol node */

		printf("\nCreating memory for Particle/Droplet Diameter (m) of Disperse phase :");

		printf("\nEnter the initialization value for Particle/Droplet Diameter (m):");
		scanf("%lf",&constantdia);

		p=CgnsBase->BaseParameters->flagVariableDia;

		//Mrunalini
		if(CgnsBase->BaseParameters->flagVariableDia==1)
			CgnsBase->BaseParameters->InitialDropletDiameter = constantdia;

		printf("\n CgnsBase->BaseParameters->flagVariableDia=%d",CgnsBase->BaseParameters->flagVariableDia);
		for(j=0;j<num_zone;j++)
		{

			if(p==0) /*If loop for checking whether memory is to be allocated or not*/
			{
				totalcells=1;

				variable1=(double*)malloc(totalcells*sizeof(double));
	                        variable2=(double*)malloc(totalcells*sizeof(double));

				for(k=0;k<totalcells;k++)
				{
					*(variable1+k)=constantdia;
					*(variable2+k)=constantdia;
				}
			}
			else /*p=1 case ie variable case*/
			{
				totalcells=1;

				for(k=0;k<Cell_Dim;k++)
				totalcells=totalcells*(*((z+j)->Size+k+Cell_Dim));

				variable1=(double*)malloc(totalcells*sizeof(double));
				variable2=(double*)malloc(totalcells*sizeof(double));

				if((variable1 == NULL)||(variable2 == NULL))
				{
				printf("Memory allocation not sufficient,Exiting from InitializeDisperse\n");
				exit(1);
				}
				for(k=0;k<totalcells;k++)
				    	{
				    	*(variable1+k)=constantdia;
				        *(variable2+k)=constantdia;
				        }
			}
				((z+j)->Flowsol)->ParticleDia=variable1;
				((z+j)->PrevFlowsol)->ParticleDia=variable2;
			}

/**********************Memory allocation for disperse phase velocity in x direction**********************************/

/********************************************************************************************************************************

NOTE: THE NEW ALGORITHM FOR GAS_PARTICLE AND GAS_DROPLET FLOWS SOLVES THE MOEMNTUM EQUATION FOR DIPSERSED PHASE FOR CALCULATING VELOCITY COMPONENTS AND NOT SPECIFIC MOMENTA!

********************************************************************************************************************************/

	printf("\n Creating memory for disperse Phase  velocity in x direction:");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=62;//Variable Tag has to be stored
	cnt=cnt+1;

		for(j=0;j<num_zone;j++)
		{

			printf("\nEnter the Initialization value of U[m/s] velocity for disperse phase :");
			scanf("%lf",&constantu);
			fprintf(var2,"%lf\n",constantu);
	//mrunalini
			/*Memory allocation in the Flow Solution Node X Direction momentum*/
		/*	((z+j)->Flowsol)->MomentumDispersex= (double *)initialize_variable(i,(z+j),constantu);
			((z+j)->PrevFlowsol)->MomentumDispersex= (double *)initialize_variable(i,(z+j),constantu);
			((z+j)->Flowsol)->MomentumDispersexOld= (double *)initialize_variable(i,(z+j),constantu);
			((z+j)->PrevFlowsol)->MomentumDispersexOld= (double *)initialize_variable(i,(z+j),constantu);*/

			//allocating memory for velocities of dispersed phase
			((z+j)->Flowsol)->U_Disperse= (double *)initialize_variable(i,(z+j),constantu);
			((z+j)->PrevFlowsol)->U_Disperse= (double *)initialize_variable(i,(z+j),constantu);

			((z+j)->Flowsol)->U_DisperseOld= (double *)initialize_variable(i,(z+j),constantu);
			((z+j)->PrevFlowsol)->U_DisperseOld= (double *)initialize_variable(i,(z+j),constantu);

			if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)//Mrunalini
			{
				//printf("\nINITIALIZING U-MASS VELOCITY IN FLOWSOL FOR EVAPORATING DROPLET FLOWS!!\n");
				((z+j)->Flowsol)->U_Dmass= (double *)initialize_variable(i,(z+j),constantu);
				//printf("\nINITIALIZING U-MASS VELOCITY IN PREVFLOWSOL FOR EVAPORATING DROPLET FLOWS!!\n");
				((z+j)->PrevFlowsol)->U_Dmass= (double *)initialize_variable(i,(z+j),constantu);
			}

			//allocating memory for drag
			((z+j)->DragX)= (double *)initialize_variable(i,(z+j),0.0);
			//allocating memory for particle relaxation time
			((z+j)->t_p)= (double *)initialize_variable(i,(z+j),0.0);

			//allocating memory diffusion of volume fraction in momentum equation coming form Turbulent dispersion force
			((z+j)->Flowsol->Diff_Disperse)= (double *)initialize_variable(i,(z+j),0.0);
			/********************************************************************/
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						//(((z+j)->Iface+iface)->ZSlab)->MomentumDispersex=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumDispersexOld=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//mrunalini
						(((z+j)->Iface+iface)->ZSlab)->U_Disperse=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantu);
						(((z+j)->Iface+iface)->ZSlab)->U_DisperseOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantu);

						if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
						{
							//printf("\nINITIALIZING U-MASS VELOCITY IN SLABS FOR EVAPORATING DROPLET FLOWS!!\n");
							(((z+j)->Iface+iface)->ZSlab)->U_Dmass = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantu);
							(((z+j)->Iface+iface)->ZSlab)->U_DmassOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantu);
						}
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersex = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersexOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//mrunalini
						(((z+j)->ZoneBc+zBc)->ZSlab)->U_Disperse=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantu);
						(((z+j)->ZoneBc+zBc)->ZSlab)->U_DisperseOld=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantu);

						if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
						{
							//printf("\nINITIALIZING U-MASS VELOCITY IN SLABS FOR EVAPORATING DROPLET FLOWS!!\n");
							(((z+j)->ZoneBc+zBc)->ZSlab)->U_Dmass = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantu);
							(((z+j)->ZoneBc+zBc)->ZSlab)->U_DmassOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantu);
						}
					}
				}
			}
		}

/**********************Memory allocation for disperse phase momentum in y direction**********************************/
	printf("\n Creating memory for disperse Phase  momentum  in y direction:");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=63;//Variable Tag has to be stored
	cnt=cnt+1;

		for(j=0;j<num_zone;j++)
		{

			printf("\nEnter the Initialization value of V[m/s] velocity for disperse phase :");
			scanf("%lf",&constantv);
			fprintf(var2,"%lf\n",constantv);

			/*Memory allocation in the Flow Solution Node Y Direction momentum*/
		/*	((z+j)->Flowsol)->MomentumDispersey= (double *)initialize_variable(i,(z+j),constantv);
			((z+j)->PrevFlowsol)->MomentumDispersey= (double *)initialize_variable(i,(z+j),constantv);
			((z+j)->Flowsol)->MomentumDisperseyOld= (double *)initialize_variable(i,(z+j),constantv);
			((z+j)->PrevFlowsol)->MomentumDisperseyOld= (double *)initialize_variable(i,(z+j),constantv);
	*/
			/************************************************************************/
			//allocating memory for velocities of dispersed phase
			((z+j)->Flowsol)->V_Disperse = (double *)initialize_variable(i,(z+j),constantv);
			((z+j)->PrevFlowsol)->V_Disperse = (double *)initialize_variable(i,(z+j),constantv);

			((z+j)->Flowsol)->V_DisperseOld = (double *)initialize_variable(i,(z+j),constantv);
			((z+j)->PrevFlowsol)->V_DisperseOld = (double *)initialize_variable(i,(z+j),constantv);

			if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
			{
				((z+j)->Flowsol)->V_Dmass= (double *)initialize_variable(i,(z+j),constantv);
				((z+j)->PrevFlowsol)->V_Dmass= (double *)initialize_variable(i,(z+j),constantv);
			}

			//allocating memory for drag
			((z+j)->DragY)= (double *)initialize_variable(i,(z+j),0.0);


			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						//(((z+j)->Iface+iface)->ZSlab)->MomentumDispersey=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumDisperseyOld=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//mrunalini
						(((z+j)->Iface+iface)->ZSlab)->V_Disperse=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantv);
						(((z+j)->Iface+iface)->ZSlab)->V_DisperseOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantv);

						if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
						{
							(((z+j)->Iface+iface)->ZSlab)->V_Dmass = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantv);
							(((z+j)->Iface+iface)->ZSlab)->V_DmassOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantv);
						}
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersey = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDisperseyOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//mrunalini
						(((z+j)->ZoneBc+zBc)->ZSlab)->V_Disperse=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantv);
						(((z+j)->ZoneBc+zBc)->ZSlab)->V_DisperseOld=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantv);

						if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
						{
							(((z+j)->ZoneBc+zBc)->ZSlab)->V_Dmass=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantv);
							(((z+j)->ZoneBc+zBc)->ZSlab)->V_DmassOld=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantv);
						}
					}
				}
			}
		}

/**************************Memory allocation of individual momentum in z direction*********************************/
	if(CgnsBase->BaseParameters->flag2d == 0)
        printf("\n Creating memory for disperse phase momentum  in z direction ");
	if(CgnsBase->BaseParameters->flag2d == 0)
	{
		i=1;
		CgnsBase->VARTAG[cnt]=64;//Variable Tag has to be stored bhans
		cnt=cnt+1;
	}
	else
	{
	        i=0;
	}
	  fprintf(var,"%d\n",i);
		for(j=0;j<num_zone;j++)
		{
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
				printf("\nEnter the Initialization value of w[m/s] velocity for disperse phase :");
				scanf("%lf",&constantw);
				fprintf(var2,"%lf\n",constantw);

			}
			else
			constantw=0.0;
			/*Memory allocation in the Flow Solution Node X Direction momentum*/
			//	((z+j)->Flowsol)->MomentumDispersez= (double *)initialize_variable(i,(z+j),constantw);
			//	((z+j)->PrevFlowsol)->MomentumDispersez= (double *)initialize_variable(i,(z+j),constantw);

			/************************************************************************/
			//allocating memory for velocities of dispersed phase
			((z+j)->Flowsol)->W_Disperse= (double *)initialize_variable(i,(z+j),constantw);
			((z+j)->PrevFlowsol)->W_Disperse= (double *)initialize_variable(i,(z+j),constantw);

			((z+j)->Flowsol)->W_DisperseOld = (double *)initialize_variable(i,(z+j),constantw);
			((z+j)->PrevFlowsol)->W_DisperseOld = (double *)initialize_variable(i,(z+j),constantw);

			if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
			{
				((z+j)->Flowsol)->W_Dmass= (double *)initialize_variable(i,(z+j),constantw);
				((z+j)->PrevFlowsol)->W_Dmass= (double *)initialize_variable(i,(z+j),constantw);
			}
			//allocating memory for drag
			((z+j)->DragZ)= (double *)initialize_variable(i,(z+j),0.0);


			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						//(((z+j)->Iface+iface)->ZSlab)->MomentumDispersez=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumDispersezOld=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//W_disperse For slabs by shashwat
						(((z+j)->Iface+iface)->ZSlab)->W_Disperse=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantw);
						(((z+j)->Iface+iface)->ZSlab)->W_DisperseOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantw);

						if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
						{
							(((z+j)->Iface+iface)->ZSlab)->W_Dmass = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantw);
							(((z+j)->Iface+iface)->ZSlab)->W_DmassOld = (double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantw);
						}
					}
				}
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersez = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumDispersezOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//W_disperse For slabs by shashwat
						(((z+j)->ZoneBc+zBc)->ZSlab)->W_Disperse=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantw);
						(((z+j)->ZoneBc+zBc)->ZSlab)->W_DisperseOld=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantw);
						if(CgnsBase->flagDroplet==1 && CgnsBase->BaseParameters->flagVariableDia == 1)
						{
							(((z+j)->ZoneBc+zBc)->ZSlab)->W_Dmass=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantw);
							(((z+j)->ZoneBc+zBc)->ZSlab)->W_DmassOld=(double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantw);
						}

					}
				}
			}
		}

/************End of memory allocation of velocity components in x,y,z direction*******************************/

/*****************************MEMORY ALLOCATION FOR MASSFLUXES*******************************************************/
/* Creatiing memory for massflux is necessary,hence no need to ask user to enter 0 or 1 for massflux*/
/* The function initializeFlux is defined in the file initialize_variable.c*/
	printf("\n Creating memory for MassFlux for the disperse phase");
	i=1;
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->MassfluxDisperse= (double *)initializeFlux(i,(z+j));
			((z+j)->PrevFlowsol)->MassfluxDisperse= (double *)initializeFlux(i,(z+j));
		}

	if(CgnsBase->flagDroplet == 1 && CgnsBase->BaseParameters->flagVariableDia ==1)
	{
		i=1;

		for(j=0;j<num_zone;j++)
		{
			((z+j)->PrevFlowsol->RhoG_face) = (double *)initializeFlux(i,(z+j));
			((z+j)->PrevFlowsol->VarthetaD_face) = (double *)initializeFlux(i,(z+j));
		}
	}
	else
	{
		for(j=0;j<num_zone;j++)
		{
			((z+j)->PrevFlowsol->RhoG_face) = NULL;
			((z+j)->PrevFlowsol->VarthetaD_face) = NULL;
		}
	}

/****
//NOTE: NUMBER DENSITY EQUATION IS NO LONGER SOLVED FOR IN THE GAS_DROPLET FLOWS WITH EVAPORATION MODULE
	///**************************Memory allocation of Number density of the droplet phase *********************************/
/*	if((CgnsBase->flagDroplet==1))
	{
	printf("\n Creating memory for Number Density of the disperse Phase :");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=67;//Variable Tag has to be stored
	cnt=cnt+1;

		for(j=0;j<num_zone;j++)
		{

			printf("\nEnter the Initialization value of Number density for Droplet phase :");
			scanf("%lf",&constantn);
			fprintf(var2,"%lf\n",constantn);

			//Memory allocation in the Flow Solution Node Number density

			((z+j)->Flowsol)->NumDensity= (double *)initialize_variable(i,(z+j),constantn);
			((z+j)->PrevFlowsol)->NumDensity= (double *)initialize_variable(i,(z+j),constantn);
		if(CgnsBase->BaseParameters->flagVariableDia == 1)
			{
			((z+j)->Flowsol)->Evaporate= (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->Evaporate= (double *)initialize_variable(i,(z+j),0.0);


			((z+j)->Flowsol)->EvapoConst= (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->EvapoConst= (double *)initialize_variable(i,(z+j),0.0);

			((z+j)->Flowsol)->QDroptoGas= (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->QDroptoGas= (double *)initialize_variable(i,(z+j),0.0);

			((z+j)->Flowsol)->QLatent = (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->QLatent = (double *)initialize_variable(i,(z+j),0.0);


			}
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->NumDensity=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);

					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->NumDensity = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
	}
	else
	{
	for(j=0;j<num_zone;j++)
		{

			((z+j)->Flowsol)->NumDensity= NULL;
			((z+j)->PrevFlowsol)->NumDensity= NULL;

			((z+j)->Flowsol)->Evaporate= NULL;
			((z+j)->PrevFlowsol)->Evaporate= NULL;


			((z+j)->Flowsol)->EvapoConst= NULL;
			((z+j)->PrevFlowsol)->EvapoConst= NULL;

			((z+j)->Flowsol)->QDroptoGas= NULL;
			((z+j)->PrevFlowsol)->QDroptoGas= NULL;

			((z+j)->Flowsol)->QLatent= NULL;
			((z+j)->PrevFlowsol)->QLatent= NULL;
			/************************************************************************/
	/***
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->NumDensity = NULL;

					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->NumDensity = NULL;
					}
				}
			}
		}
	}
*/

//Following lines written by mrunalini for evaporating gas-droplet flows

	if((CgnsBase->flagDroplet==1) && (CgnsBase->BaseParameters->flagVariableDia == 1))
	{
		i=1;
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->Evaporate= (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->Evaporate= (double *)initialize_variable(i,(z+j),0.0);


			((z+j)->Flowsol)->EvapoConst= (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->EvapoConst= (double *)initialize_variable(i,(z+j),0.0);

			((z+j)->Flowsol)->QDroptoGas= (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->QDroptoGas= (double *)initialize_variable(i,(z+j),0.0);

			((z+j)->Flowsol)->QLatent = (double *)initialize_variable(i,(z+j),0.0);
			((z+j)->PrevFlowsol)->QLatent = (double *)initialize_variable(i,(z+j),0.0);


		}
	}
	else
	{
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->Evaporate= NULL;
			((z+j)->PrevFlowsol)->Evaporate= NULL;


			((z+j)->Flowsol)->EvapoConst= NULL;
			((z+j)->PrevFlowsol)->EvapoConst= NULL;

			((z+j)->Flowsol)->QDroptoGas= NULL;
			((z+j)->PrevFlowsol)->QDroptoGas= NULL;

			((z+j)->Flowsol)->QLatent= NULL;
			((z+j)->PrevFlowsol)->QLatent= NULL;
		}
	}

	if((CgnsBase->flagDroplet==1)&&(CgnsBase->BaseParameters->flagVariableDia==1))
	{
	printf("\nMaking Density of the Carrier Phase as solvable variable:");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=710;//Variable Tag has to be stored
	cnt=cnt+1;
	}
        ///**************************Memory allocation of Fuel mass fraction of the droplet phase *********************************/
	if((CgnsBase->flagDroplet==1)&&(CgnsBase->BaseParameters->flagVariableDia==1)&&(CgnsBase->BaseParameters->flagDropFuelMassFrac==1))
	{
	printf("\n Creating memory for evaporated fuel mass fraction of the disperse Phase :");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=68;//Variable Tag has to be stored
	cnt=cnt+1;

		for(j=0;j<num_zone;j++)
		{

			printf("\nEnter the Initialization value of Evaporated Fuel mass fraction :");
			scanf("%lf",&constantn);
			fprintf(var2,"%lf\n",constantn);

			/*Memory allocation in the Flow Solution Node Number density*/

			((z+j)->Flowsol)->DropFuelMassFrac= (double *)initialize_variable(i,(z+j),constantn);
			((z+j)->PrevFlowsol)->DropFuelMassFrac= (double *)initialize_variable(i,(z+j),constantn);

			/************************************************************************/

			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->DropFuelMassFrac=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantn);
						(((z+j)->Iface+iface)->ZSlab)->DropFuelMassFracOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantn);

					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->DropFuelMassFrac = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantn);
						(((z+j)->ZoneBc+zBc)->ZSlab)->DropFuelMassFracOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantn);
					}
				}
			}
		}
			printf("\nCreating memory for Mass Diffusivity [m2/s]:");
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify mass diffusivity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
				scanf("%d",&flag);
				fprintf(var2,"%d\n",flag);
			} else { flag = 0;}
			if(flag==0)
			{
				printf("\n If Mass Diffusivity is a CONSTANT, Enter 0");
				printf("\n If Mass Diffusivity is a VARIABLE, Enter 1:");
				scanf("%d",&p);
				fprintf(var2,"%d\n",p);
				varCoeffCount = 155;
				if(p==0)
				{
					printf("\nEnter the constant value of mass diffusivity[m2/s]:");
					scanf("%lf",&constant);
					fprintf(var2,"%lf\n",constant);
				}
				else
					constant=0.0;

				for(j=0;j<num_zone;j++)
				{
					((z+j)->Flowsol)->DropDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->DropDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagDropDiffcoeff=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->DropDiffcoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
								//printf("I am leaving Iface\n");
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->DropDiffcoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
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
					printf("\n If Mass Diffusivity is a CONSTANT, Enter 0");
					printf("\n If Mass Diffusivity is a VARIABLE, Enter 1:");
					scanf("%d",&p);
					fprintf(var2,"%d\n",p);
					varCoeffCount = 155;
					if(p==0)
					{
						printf("\nEnter the constant value of Mass Diffusivity[m2/s]:");
						scanf("%lf",&constant);
						fprintf(var2,"%lf\n",constant);
					}
					else
						constant=0.0;

					((z+j)->Flowsol)->DropDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->DropDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagDropDiffcoeff=p;
					numiface=(z+j)->Number_of_Interfaces;

					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->DropDiffcoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->DropDiffcoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
/**************************************END OF MEMORY FOR THERMAL CONDUCTIVITY*****************************************************************/
	}
	else
	{
	for(j=0;j<num_zone;j++)
		{

			((z+j)->Flowsol)->DropFuelMassFrac= NULL;
			((z+j)->PrevFlowsol)->DropFuelMassFrac= NULL;
			((z+j)->Flowsol)->flagDropDiffcoeff=0;
			((z+j)->Flowsol)->DropDiffcoeff=NULL;
			((z+j)->PrevFlowsol)->DropDiffcoeff=NULL;
			numiface=(z+j)->Number_of_Interfaces;

			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->DropFuelMassFrac = NULL;
						(((z+j)->Iface+iface)->ZSlab)->DropFuelMassFracOld = NULL;
						(((z+j)->Iface+iface)->ZSlab)->DropDiffcoeff=NULL;
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->DropFuelMassFrac = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->DropFuelMassFracOld = NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->Diffcoeff=NULL;
					}
				}
			}
		}
	}
//######################################################################################//
///**************************Memory allocation of temperature of the droplet phase *********************************/
	if((CgnsBase->flagDroplet==1)&&(CgnsBase->BaseParameters->flagVariableDia==1))
	{
	printf("\n Creating memory for temperature of disperse Phase :");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=69;//Variable Tag has to be stored
	cnt=cnt+1;

		for(j=0;j<num_zone;j++)
		{

			printf("\nEnter the Initialization value of droplet temperature :");
			scanf("%lf",&constantn);
			fprintf(var2,"%lf\n",constantn);

			/*Memory allocation in the Flow Solution Node Number density*/

			((z+j)->Flowsol)->Tdrop = (double *)initialize_variable(i,(z+j),constantn);
			((z+j)->PrevFlowsol)->Tdrop = (double *)initialize_variable(i,(z+j),constantn);

			/************************************************************************/

			numiface=(z+j)->Number_of_Interfaces;

			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Tdrop=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantn);
						(((z+j)->Iface+iface)->ZSlab)->TdropOld=(double *)initSlab(i,dir,((z+j)->Iface+iface)->Point_Range,constantn);

					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Tdrop = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantn);
						(((z+j)->ZoneBc+zBc)->ZSlab)->TdropOld = (double *)initSlab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantn);
					}
				}
			}
		}
			printf("\nCreating memory for specific heat for liquid droplet [KJ/Kg K]:");
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify specific heat for liquid droplet for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
				scanf("%d",&flag);
				fprintf(var2,"%d\n",flag);
			} else { flag = 0;}
			if(flag==0)
			{
				printf("\n If Specific heat for liquid droplet is CONSTANT, Enter 0");
				printf("\n If Specific heat for liquid droplet is a VARIABLE, Enter 1:");
				scanf("%d",&p);
				fprintf(var2,"%d\n",p);
				varCoeffCount = 1555; //required or not, need to be checked   //shash
				if(p==0)
				{
					printf("\nEnter the constant value of specific heat for liquid dropletKJ/Kg K]:");
					scanf("%lf",&constant);
					fprintf(var2,"%lf\n",constant);
				}
				else
					constant=1.0;

				for(j=0;j<num_zone;j++)
				{
					((z+j)->Flowsol)->CPdrop = (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->CPdrop = (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagCPdrop = p;
					numiface=(z+j)->Number_of_Interfaces;

					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->CPdrop = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);

							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->CPdrop = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);

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

					printf("\n If Specific heat for liquid droplet is CONSTANT, Enter 0");
					printf("\n If Specific heat for liquid droplet is a VARIABLE, Enter 1:");
					scanf("%d",&p);
					fprintf(var2,"%d\n",p);
					varCoeffCount = 1555; //required or not, need to be checked  /*3->THERMAL CONDUCTIVITY*/ //shash
					if(p==0)
					{
					printf("\nEnter the constant value of specific heat for liquid dropletKJ/Kg K]:");
					scanf("%lf",&constant);
					fprintf(var2,"%lf\n",constant);
					}
					else
					constant=0.0;
					((z+j)->Flowsol)->CPdrop = (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->CPdrop = (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagCPdrop = p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->CPdrop=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->CPdrop = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
		printf("\nEnter the value of specific heat for evaporated vapour from droplets [KJ/Kg K]:");
		scanf("%lf",&CPvapour);
		CgnsBase->BaseParameters->CPvapour = CPvapour;
/**************************************END OF MEMORY FOR THERMAL CONDUCTIVITY*****************************************************************/
	}
	else
	{
	for(j=0;j<num_zone;j++)
		{

			((z+j)->Flowsol)->Tdrop = NULL;
			((z+j)->PrevFlowsol)->Tdrop= NULL;
			((z+j)->Flowsol)->flagCPdrop = 0;
			((z+j)->Flowsol)->CPdrop = NULL;
			((z+j)->PrevFlowsol)->CPdrop = NULL;
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Tdrop = NULL;
						(((z+j)->Iface+iface)->ZSlab)->TdropOld = NULL;
						(((z+j)->Iface+iface)->ZSlab)->CPdrop = NULL;
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Tdrop= NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->TdropOld= NULL;
						(((z+j)->ZoneBc+zBc)->ZSlab)->CPdrop=NULL;
					}
				}
			}
		}
	}



//#####################################################################################//
	for(j=0;j<num_zone;j++)
	{
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		ptrSize=(z+j)->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		//printf("\nNi %d NJ %d Nk %d\n",NI,NJ,NK);
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					for(c=0;c<NK-1;c++)
					for(b=0;b<NJ-1;b++)
					for(a=0;a<NI-1;a++)
					{
						if(dir==1||dir==2)
							index= c*(NJ-1)+b;
						if(dir==3||dir==4)
							index= c*(NI-1)+a;
						if(dir==5||dir==6)
							index= b*(NI-1)+a;
						//printf("\n Udisperse %lf",*((((zone)->Iface+iface)->ZSlab)->slabUDisperse+index));
						//*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);

					}

				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				//printf("\n----------------------------------------\n");
				//printf("\n----------------------------------------\n");
				//printf("\nNi %d NJ %d Nk %d\n",NI,NJ,NK);
					for(c=0;c<NK-1;c++)
					for(b=0;b<NJ-1;b++)
					for(a=0;a<NI-1;a++)
					{
						if(dir==1||dir==2)
							index= c*(NJ-1)+b;
						if(dir==3||dir==4)
							index= c*(NI-1)+a;
						if(dir==5||dir==6)
							index= b*(NI-1)+a;
						//printf("\n %d %lf\t %lf\t %lf\t %lf\t  \n",dir,*(((z+j)->ZoneBc+zBc)->ZSlab->MomentumDispersex+index),*(((z+j)->ZoneBc+zBc)->ZSlab->MomentumDispersey+index),*(((z+j)->ZoneBc+zBc)->ZSlab->RhoDisperse+index*(CgnsBase->BaseParameters->flagRhoDisperse)),*(((z+j)->ZoneBc+zBc)->ZSlab->VarthetaDisperse+index));
					//printf("\n###################################\n");
					}
				}
			}
		}

	}


	//initializeVarDisperse(CgnsBase); Commented as the momentum equations are solved for velocity components and not specific momenta

//printf("\n-------------------------------------------------------------------------\n");
for(j=0;j<num_zone;j++)
	{
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		ptrSize=(z+j)->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		//printf("\nNi %d NJ %d Nk %d\n",NI,NJ,NK);
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					for(c=0;c<NK-1;c++)
					for(b=0;b<NJ-1;b++)
					for(a=0;a<NI-1;a++)
					{
						if(dir==1||dir==2)
							index= c*(NJ-1)+b;
						if(dir==3||dir==4)
							index= c*(NI-1)+a;
						if(dir==5||dir==6)
							index= b*(NI-1)+a;
						//printf("\n Udisperse %lf",*((((zone)->Iface+iface)->ZSlab)->slabUDisperse+index));
						//*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);

					}

				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				//printf("\n----------------------------------------\n");
				//printf("\n----------------------------------------\n");
				//printf("\nNi %d NJ %d Nk %d\n",NI,NJ,NK);
					for(c=0;c<NK-1;c++)
					for(b=0;b<NJ-1;b++)
					for(a=0;a<NI-1;a++)
					{
						if(dir==1||dir==2)
							index= c*(NJ-1)+b;
						if(dir==3||dir==4)
							index= c*(NI-1)+a;
						if(dir==5||dir==6)
							index= b*(NI-1)+a;
						//printf("\n %d %lf\t %lf\t %lf\t %lf\t  \n",dir,*(((z+j)->ZoneBc+zBc)->ZSlab->MomentumDispersex+index),*(((z+j)->ZoneBc+zBc)->ZSlab->MomentumDispersey+index),*(((z+j)->ZoneBc+zBc)->ZSlab->RhoDisperse+index*(CgnsBase->BaseParameters->flagRhoDisperse)),*(((z+j)->ZoneBc+zBc)->ZSlab->VarthetaDisperse+index));
					//printf("\n###################################\n");
					}
				}
			}
		}

	}
 	//fclose(var);
 	//fclose(var2);
	//exit(155);
	return(cnt);
}
