#include "APKinitializeElectric.h"

/**
The Function initializeElectric() is used for creating memory for electric field variables

This function inturn calls
initialize_variable() and initializeslab() for each variable to create the necessary memory space for that 
variable in Flowsolution structure and also memory for that variable in slabs. For extra information see initialize().
*/
int initializeElectric(Base *CgnsBase, int cnt)
{ 
	FILE *var,*var2,*EfieldConv;
	int i,j,p,flag,num_zone;/* num_zone is the number of zones*/
	int numiface,numBc,iface,zBc,dir,varCoeffCount=0,Tchoice,conv;
	Zone *z;
	double constant;
	
	var=fopen("Settings/variablelist.dat","a");
	var2=fopen("Settings/coefficients.dat","a");
	EfieldConv=fopen("Settings/ElectricFieldconvergCriteria.dat","w");
	
	
	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones; 
	//base->BaseParameters
	if((CgnsBase->flagEfield)==1)
	{
		printf("\nEnter Scaling Factor for Electric Source Term in Momentum Equations:");
		scanf("%lf", &constant);
		CgnsBase->BaseParameters->ScaleFactor = constant;
	}
	printf("\nThe following information is specifically intended for ELECTRIC FIELD ONLY\n");
/*#####################################################################################################################*/
	printf("Please Enter the time steping scheme for ELECTRIC FIELD Equations\n[CrankNicolson->1/BTCS->0]");
	scanf("%d",&j);
	if(j==1)
		CgnsBase->BaseParameters->EfieldTSSTag=1; 
	else
		CgnsBase->BaseParameters->EfieldTSSTag=0;

/*#####################################################################################################################*/
	printf("Enter the Convection Scheme for ELECTRIC FIELD Equations(Deferred Corrected) QUICK[enter=1] or Central Difference[enter=0]\n");
	scanf("%d",&conv);
	CgnsBase->BaseParameters->EfieldConvScheme=conv;
/*###########################################################################################*/	
	printf("\n########### ENTER THE FOLLOWING BASE PARAMETERS ############");
	
	printf("\nEnter Convergence criteria for ELECTRIC FIELD Equations:");
	scanf("%lf",&constant);
	fprintf(EfieldConv,"%12.8e\n",constant);
	CgnsBase->BaseParameters->EfieldConvCri=constant;
	
	printf("Do you want to enter the Time step[Enter 1] or let the Code fix it according to Stability Constraints[Enter 0]\nFor ELECTRIC FIELD Equations:");
	scanf("%d",&Tchoice);
	CgnsBase->BaseParameters->EfieldTstepChoice=Tchoice;
	
	if(Tchoice==1)
	{
		printf("\nEnter TimeStep for ELECTRIC FIELD Equations:");
		scanf("%lf",&constant);
		CgnsBase->BaseParameters->EfieldTimeStep = constant;
	}
	else
	{
		CgnsBase->BaseParameters->EfieldTimeStep =1.0e-3;
	}
	
	if(CgnsBase->BaseParameters->EfieldConvScheme==0)
	{
		printf("\nEnter Ratio of CentralDiff(1->pure Centraldiff) to Upwind(0->pure Upwind):");
		scanf("%lf",&constant);
		CgnsBase->BaseParameters->EfieldRatioCenUP = constant;
	}
	else
	{
		printf("\nEnter Ratio of QUICK(1->pure QUICK) to Upwind(0->pure Upwind):");
		scanf("%lf",&constant);
		CgnsBase->BaseParameters->EfieldRatioCenUP = constant;
	}
	
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MEMORY ALLOCATION FOR DiElectric Constant &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	printf("\nCreating memory for DiElectric Coefficient[F/m]:");
	printf("\n\n********************************************************************************\n");
	printf("\nYou have to enter the value of DiElectric constant for a DIMENSIONAL problem and\n");
	printf("INVERSE of the charge density coefficient in the POISSON EQUATION for a NON-DIMENSIONAL problem\n");
	printf("\n***********************************************************************************************\n");
	if( CgnsBase->Number_of_zones > 1)
	{
		printf("\n Do you want to specify Permittivity for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
		scanf("%d",&flag);
		fprintf(var2,"%d\n",flag);
	} else { flag = 0;}
	if(flag==0) /*Initializing Dielectric constant for whole domain at a time*/
	{
		printf("\n If Permittivity is a CONSTANT, Enter 0");
		printf("\n If Permittivity is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 5;  /*5->DiElectric Constant*/
		if(p==0)
		{
			printf("\nEnter the constant value of DiElectric Coefficient[F/m]:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		}
		else
			constant=0.0;

		for(j=0;j<num_zone;j++)
		{       
			((z+j)->Flowsol)->DiElecConst= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->DiElecConst= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flagDiElec = p;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->DiElecConst=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->DiElecConst = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
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
			printf("\n If Permittivity is a CONSTANT, Enter 0");
			printf("\n If Permittivity is a VARIABLE, Enter 1:");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);
			varCoeffCount = 5;  /*5->DiElectric Constant*/
			if(p==0)
			{
				printf("\nEnter the constant value of DiElectric Constant[F/m]:");
				scanf("%lf",&constant);
				fprintf(var2,"%lf\n",constant);
			}
			else
				constant=0.0;
			
			((z+j)->Flowsol)->DiElecConst= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->DiElecConst= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flagDiElec = p;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->DiElecConst=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					(((z+j)->ZoneBc+zBc)->ZSlab)->DiElecConst = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}	
	}
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& END OF MEMORY ALLOCATION FOR DiElectric Constant &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     MEMORY ALLOCATION FOR Ion Mobility Rate      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	printf("\nCreating memory for Ion Mobility Rate[Sq.m/v-s]:");
	printf("\n\n*************************************************************************\n");
	printf("\nYou have to enter the value of Ion Mobility for a DIMENSIONAL problem and\n");
	printf("The convective coefficient in the CHARGE DENSITY EQUATION for a NON-DIMENSIONAL problem\n");
	printf("\n***************************************************************************************\n");
	if( CgnsBase->Number_of_zones > 1)
	{
		printf("\n Do you want to specify Mobility for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
		scanf("%d",&flag);
		fprintf(var2,"%d\n",flag);
	} else { flag = 0;}
	if(flag==0) /*Initializing Dielectric constant for whole domain at a time*/
	{
		printf("\n If Mobility is a CONSTANT, Enter 0");
		printf("\n If Mobility is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 6;  /*6->Ion Mobility Rate*/
		if(p==0)
		{
			printf("\nEnter the constant value of Ion Mobility Rate[Sq.m/v-s]:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		}
		else
			constant=0.0;

		for(j=0;j<num_zone;j++) 
		{       
			((z+j)->Flowsol)->MobIon= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->MobIon= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flagMobIon = p;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->MobIon=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->MobIon = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
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
			printf("\n If Mobility is a CONSTANT, Enter 0");
			printf("\n If Mobility is a VARIABLE, Enter 1:");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);
			varCoeffCount = 6;  /*6->Ion Mobility Rate*/
			if(p==0)
			{
				printf("\nEnter the constant value of Ion Mobility Rate[Sq.m/v-s]:");
				scanf("%lf",&constant);
				fprintf(var2,"%lf\n",constant);
			}
			else
				constant=0.0;
			
			((z+j)->Flowsol)->MobIon= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->MobIon= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flagMobIon = p;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->MobIon=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					(((z+j)->ZoneBc+zBc)->ZSlab)->MobIon = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}	
	}
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   END OF MEMORY ALLOCATION FOR Ion Mobility Rate   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MEMORY ALLOCATION FOR Charge Diffusion Coefficient &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	printf("\nCreating memory for Ion Charge Diffusion Coefficient[Sq.m/s]:");
	printf("\n\n*************************************************************************\n");
	printf("\nYou have to enter the value of Charge Diffusion for a DIMENSIONAL problem and\n");
	printf("The diffusion coefficient in the CHARGE DENSITY EQUATION for a NON-DIMENSIONAL problem\n");
	printf("\n**************************************************************************************\n");
	if( CgnsBase->Number_of_zones > 1)
	{
		printf("\n Do you want to specify Ion Charge Diffusion for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
		scanf("%d",&flag);
		fprintf(var2,"%d\n",flag);
	} else { flag = 0;}
	if(flag==0) /*Initializing Dielectric constant for whole domain at a time*/
	{
		printf("\n If Ion Charge Diffusion is a CONSTANT, Enter 0");
		printf("\n If Ion Charge Diffusion is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 7;  /*7->Charge Diffusion Coefficient*/
		if(p==0)
		{
			printf("\nEnter the constant value of Charge Diffusion Coefficient[Sq.m/s]:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
		}
		else
			constant=0.0;

		for(j=0;j<num_zone;j++) 
		{       
			((z+j)->Flowsol)->ChargeDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->ChargeDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flagChargeDiffCoeff = p;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->ChargeDiffcoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->ChargeDiffcoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
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
			printf("\n If Ion Charge Diffusion is a CONSTANT, Enter 0");
			printf("\n If Ion Charge Diffusion is a VARIABLE, Enter 1:");
			scanf("%d",&p);
			fprintf(var2,"%d\n",p);
			varCoeffCount = 7;  /*7->Charge Diffusion Coefficient*/
			if(p==0)
			{
				printf("\nEnter the constant value of Charge Diffusion Coefficient[Sq.m/s]:");
				scanf("%lf",&constant);
				fprintf(var2,"%lf\n",constant);
			}
			else
				constant=0.0;
			
			((z+j)->Flowsol)->ChargeDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->ChargeDiffcoeff= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flagChargeDiffCoeff = p;
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->ChargeDiffcoeff=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					(((z+j)->ZoneBc+zBc)->ZSlab)->ChargeDiffcoeff = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}	
	}
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& END OF MEMORY ALLOCATION FOR Charge Diffusion Coefficient &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      END OF MEMORY ALLOCATION FOR Electric Potential      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
/*The memory allocation for Electrical Parameters for Flowsol node and Slabs at boundary and Interface */
	printf("\n Creating memory for Electric Potential[Volts]");
	i=1;
	fprintf(var,"%d\n",i);
	//printf("\nCount before potential = %d\n", cnt);
	CgnsBase->VARTAG[cnt]=8;//Variable Tag has to be stored
	cnt=cnt+1;
	//CgnsBase->VARTAG[cnt]=11;//EfieldX
	//cnt=cnt+1;
	//CgnsBase->VARTAG[cnt]=12;//EfieldY
	//cnt=cnt+1;
	//CgnsBase->VARTAG[cnt]=13;//EfieldZ
	//cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{
			printf("\nEnter the Initialization value of Electric Potential[Volts]:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
			((z+j)->Flowsol)->Potential= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->Potential= (double *)initialize_variable(i,(z+j),constant);
			if(CgnsBase->flagEfield==1)
			{
				/** These two variables are added by bhanesh to improve the efficiency of souce calculation for ESP problems.  The source were calculated unnecessilarily every time.  This has been avoided by storing the sources in these variables rather than calculating them agin*/
				
				((z+j)->Flowsol)->ESP_Source[0]=(double *)initialize_variable(i,(z+j),0.0);
				((z+j)->Flowsol)->ESP_Source[1]=(double *)initialize_variable(i,(z+j),0.0);
				
				if(CgnsBase->BaseParameters->flag2d == 0)
					((z+j)->Flowsol)->ESP_Source[2]=(double *)initialize_variable(i,(z+j),0.0);
			}	
			/*
			//Allocating Memory For Electric Field in different directions
			((z+j)->Flowsol)->EfieldX= (double *)initialize_variable(i,(z+j),0);
			((z+j)->PrevFlowsol)->EfieldX= (double *)initialize_variable(i,(z+j),0);
			fprintf(var,"%d\n",i);
			

			((z+j)->Flowsol)->EfieldY= (double *)initialize_variable(i,(z+j),0);
			((z+j)->PrevFlowsol)->EfieldY= (double *)initialize_variable(i,(z+j),0);
			fprintf(var,"%d\n",i);
			

			((z+j)->Flowsol)->EfieldZ= (double *)initialize_variable(i,(z+j),0);
			((z+j)->PrevFlowsol)->EfieldZ= (double *)initialize_variable(i,(z+j),0);
			fprintf(var,"%d\n",i);
			*/
			
			//Allocating Memory For Electric Flux
			((z+j)->Flowsol)->Electricflux= (double *)initializeFlux(i,(z+j));
			((z+j)->PrevFlowsol)->Electricflux= (double *)initializeFlux(i,(z+j));
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Potential = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->PotentialOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->EfieldX = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->EfieldY = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->EfieldZ = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Potential =    (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->PotentialOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->EfieldX =    (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->EfieldY =    (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->EfieldZ =    (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}//End of Electric Potential
	
	printf("\n Creating memory for Ion Charge Density[Coulombs per metercube]");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=9;//Variable Tag has to be stored
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{
			printf("\nEnter the Initialization value of Ion Space Charge Density [Coulombs per metercube]:");
			scanf("%lf",&constant);
			fprintf(var2,"%lf\n",constant);
			((z+j)->Flowsol)->RhoIon= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->RhoIon= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->RhoIon    = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->RhoIonOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->RhoIon    = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->RhoIonOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}//End of Ion Charge Density
		
	printf("\nDo you want to calculate Electric Field in Three Directions?\nEnter:1(if YES)\nEnter:0(if NO)");
	scanf("%d",&Tchoice);
	fprintf(var,"%d\n",Tchoice);
	if(Tchoice==1)
		CgnsBase->BaseParameters->TagEfield=1;
	else
		CgnsBase->BaseParameters->TagEfield=0;
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->EfieldX= (double *)initialize_variable(Tchoice,(z+j),0);
		((z+j)->Flowsol)->EfieldY= (double *)initialize_variable(Tchoice,(z+j),0);
		((z+j)->Flowsol)->EfieldZ= (double *)initialize_variable(Tchoice,(z+j),0);
	}
		
	fclose(EfieldConv);

	return(cnt);
}
