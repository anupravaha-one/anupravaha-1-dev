#include "APKinitializeMPFV.h"
/**
 This subroutine is developed for multiphase flow(bubbly ), where all the flow variables which is being used in case of mutiphase flow has been initialized and their memory allocation has been done here itself.For individual phase momentum directly velocity is stored in the memory which has been alterd in while deriving global momentum in the InitializeVarMPFV().
 */

int initializeMPFV(Base *CgnsBase,int cnt)
{ 
	FILE *var,*var2;
	int i,j,p,num_zone;/* num_zone is the number of zones*/
	int numiface,numBc,iface,zBc,dir,varCoeffCount=0,nophase,nop;
	Zone *z;
	double  *VolFrac,constantu[5],constantv[5],constantw[5],constantden[5],constantmu[5],cons=0.0,sumvol,constant;
	int che;	
	var=fopen("Settings/variablelist.dat","a");
	var2=fopen("Settings/coefficients.dat","a");
	
	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones;
	nop=CgnsBase->nphase;
	VolFrac = (double *)malloc(nop*num_zone*sizeof(double));
	/*#####################################################################################################################*/	//   MULTI PHASE FLOW VARIABLES  problem.	
	
/***************************************************************************************************************/	/*The memory allocation for Density for Flowsol node and Slabs at boundary and Interface*/
	
	printf("\n\nAlways give first phase as a continuous phase and 2nd phase as disperse phase:\n ");
	printf("\n i.e, Fisrt phase should be heavier phase and 2nd phase must be lighter one\n");
	for(nophase=0;nophase<nop;nophase++)	
	{	
		
		printf("\nCreating memory for Density[kg/m3] of  phase %d:",nophase);
		printf("\n If Density is a CONSTANT, Enter 0");
		printf("\n If Density is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 1;  /*1->Density*/
		
		if(p==0)
		{
			
	
			printf("\nEnter the constant value of Density[kg/m3]for phase %d:",nophase);
			scanf("%lf",&constantden[nophase]);
			fprintf(var2,"%lf\n",constantden[nophase]);
			
		
		}
		
		else
 			cons=0.0;

		for(j=0;j<num_zone;j++)
		{
		 
		
			((z+j)->Flowsol)->Rho= (double *)initializeVarCoefficient(1,(z+j),cons,varCoeffCount);
			((z+j)->PrevFlowsol)->Rho= (double *)initializeVarCoefficient(1,(z+j),cons,varCoeffCount);
			((z+j)->Flowsol)->flagRho=1;//default value
		 
		 
		      //  for(nophase=0;nophase<nop;nophase++)
			//{
			/*Memory allocation in the Flow Solution Node */
			((z+j)->Flowsol)->RhoofPhase[nophase]= (double *)initializeVarCoefficient(p,(z+j),constantden[nophase],varCoeffCount);
			((z+j)->PrevFlowsol)->RhoofPhase[nophase]= (double *)initializeVarCoefficient(p,(z+j),constantden[nophase],varCoeffCount);
			
			CgnsBase->BaseParameters->flagRhoPhase[nophase]=p;
			//}
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Rho=(double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						
						
											
						/*Memory allocation in the Interface Slab */
						
  						//for(nophase=0;nophase<nop;nophase++)
						//{
						(((z+j)->Iface+iface)->ZSlab)->RhoofPhase[nophase]=(double*)initializeSlabVarCoeff(p,(z+j),constantden[nophase],varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						//}//printf("I am leaving Iface\n");
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					(((z+j)->ZoneBc+zBc)->ZSlab)->Rho = (double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					
					/*Memory allocation in the Interface Slab X Direction Velocity*/
					//printf("I am going to zoneBc\n");
					//for(nophase=0;nophase<nop;nophase++)
					//{
					(((z+j)->ZoneBc+zBc)->ZSlab)->RhoofPhase[nophase] = (double *)initializeSlabVarCoeff(p,(z+j),constantden[nophase],varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					
					//}//printf("I am leaving zoneBc\n");
					}
				}
			}
		}
		
		
	/*****************************MEMORY ALLOCATION FOR VISCOSITY*******************************************************/	
		printf("\nCreating memory for VISCOSITY[kg/m-s] of  phase %d:",nophase);
	
		printf("\n If Viscosity is a CONSTANT, Enter 0");
		printf("\n If Viscosity is a VARIABLE, Enter 1:");
		scanf("%d",&p);
		fprintf(var2,"%d\n",p);
		varCoeffCount = 2;  /*2->VISCOSITY*/
		if(p==0)
		{     //  for(nophase=0;nophase<nop;nophase++)
			//{
			printf("\nEnter the constant value of VISCOSITY[kg/m-s]for phase %d:",nophase);
			scanf("%lf",&constantmu[nophase]);
			fprintf(var2,"%lf\n",constantmu[nophase]);
			//}
		}
		else
			cons=0.0;

		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->Viscoeff= (double *)initializeVarCoefficient(1,(z+j),cons,varCoeffCount);
			((z+j)->PrevFlowsol)->Viscoeff= (double *)initializeVarCoefficient(1,(z+j),cons,varCoeffCount);
			((z+j)->Flowsol)->flagViscoeff=1;//default value
			
			((z+j)->Flowsol)->ViscoeffDummy[nophase]= (double *)initializeVarCoefficient(1,(z+j),cons,varCoeffCount);
			((z+j)->PrevFlowsol)->ViscoeffDummy[nophase]= (double *)initializeVarCoefficient(1,(z+j),cons,varCoeffCount);
			/*Memory allocation in the Flow Solution Node */
			//for(nophase=0;nophase<nop;nophase++)
			//{
			((z+j)->Flowsol)->ViscosityofPhase[nophase]= (double *)initializeVarCoefficient(p,(z+j),constantmu[nophase],varCoeffCount);
			((z+j)->PrevFlowsol)->ViscosityofPhase[nophase]= (double *)initializeVarCoefficient(p,(z+j),constantmu[nophase],varCoeffCount);
			
			CgnsBase->BaseParameters->flagViscoeffPhase[nophase]=p;
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
						(((z+j)->Iface+iface)->ZSlab)->Viscoeff=(double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);	
						
						/*Memory allocation in the Interface Slab X Direction Velocity*/
						//printf("I am going to Iface\n");
						//for(nophase=0;nophase<nop;nophase++)
						//{
						
						(((z+j)->Iface+iface)->ZSlab)->ViscoeffDummy[nophase]=(double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->ViscosityofPhase[nophase]=(double *)initializeSlabVarCoeff(p,(z+j),constantmu[nophase],varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						 //}
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					
					(((z+j)->ZoneBc+zBc)->ZSlab)->Viscoeff = (double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					
					(((z+j)->ZoneBc+zBc)->ZSlab)->ViscoeffDummy[nophase] = (double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					/*Memory allocation in the Interface Slab */
					
					//for(nophase=0;nophase<nop;nophase++)
					//{
					(((z+j)->ZoneBc+zBc)->ZSlab)->ViscosityofPhase[nophase] = (double *)initializeSlabVarCoeff(p,(z+j),constantmu[nophase],varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					//}
					}
				}
			}
		}
		
	}
	
	printf("\nCreating memory for bubble diameter[m] of disperse phase :");
	
	printf("\n If bubble diameter is a CONSTANT, Enter 0");
	printf("\n If bubble diameter is a VARIABLE, Enter 1:");
	scanf("%d",&p);
	fprintf(var2,"%d\n",p);
	varCoeffCount = 3;  
	if(p==0)
	{     //  for(nophase=0;nophase<nop;nophase++)
			//{
		printf("\nEnter the constant value of bubble diameter[m]for disperse phase :");
		scanf("%lf",&constant);
		fprintf(var2,"%lf\n",constant);
			//}
	}
	else
		cons=0.0;
	
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->ParticleDia= (double *)initialize_variable(1,(z+j),constant);
		((z+j)->PrevFlowsol)->ParticleDia= (double *)initialize_variable(1,(z+j),constant);
		
		printf("\n\ndia=%e",*(((z+j)->Flowsol)->ParticleDia+0));
		CgnsBase->BaseParameters->flagVariableDia = p;
		
			
		
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
					//(((z+j)->Iface+iface)->ZSlab)->ParticleDia=(double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);	
						
					
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					
					//(((z+j)->ZoneBc+zBc)->ZSlab)->ParticleDia = (double *)initializeSlabVarCoeff(1,(z+j),cons,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					
					
				}
			}
		}
	}
		

	
	
/****************************************************************************************************************/


/*****************************END OF MEMORY ALLOCATION FOR VISCOSITY*******************************************************/

/**************************************************************************************************************/		/*Memory allocation for individual and overall momentum*/

/******************************************************************************************************************/
	printf("\n Creating memory for overall momentum  in x direction");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=31;//Variable Tag has to be stored 
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{       
			
			/*Memory allocation in the Flow Solution Node X Direction momentum*/
			((z+j)->Flowsol)->Momentumx= (double *)initialize_variable(i,(z+j),cons);
			((z+j)->PrevFlowsol)->Momentumx= (double *)initialize_variable(i,(z+j),cons);
			((z+j)->Flowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),cons);
			((z+j)->PrevFlowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),cons);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Momentumx=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumxOld=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->Momentumxstar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumxstarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumx = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumxstar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumxOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumxstarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
					
/**************************Overall momentum in y direction*********************************************************/
		
	printf("\n Creating memory for overall momentum  in y direction");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=32;//Variable Tag has to be stored 
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{       
			/*Memory allocation in the Flow Solution Node Y Direction momentum*/
			((z+j)->Flowsol)->Momentumy= (double *)initialize_variable(i,(z+j),cons);
			((z+j)->PrevFlowsol)->Momentumy= (double *)initialize_variable(i,(z+j),cons);
			((z+j)->Flowsol)->Momentumystar= (double *)initialize_variable(i,(z+j),cons);
			((z+j)->PrevFlowsol)->Momentumystar= (double *)initialize_variable(i,(z+j),cons);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Momentumy=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumyOld=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->Momentumystar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumystarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumy= (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumystar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumyOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumystarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
		
/*********************Overall momentum in Z direction*************************************************************/
	if(CgnsBase->BaseParameters->flag2d == 0)
	
	printf("\n Creating memory for overall momentum  in Z direction");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=33;//Variable Tag has to be stored 
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{       
			if(CgnsBase->BaseParameters->flag2d == 0)	
			{	cons=0.0;
				
// 				((z+j)->Flowsol)->Momentumz= (double *)initialize_variable(i,(z+j),cons);
// 				((z+j)->PrevFlowsol)->Momentumz= (double *)initialize_variable(i,(z+j),cons);
// 				((z+j)->Flowsol)->Momentumzstar= (double *)initialize_variable(i,(z+j),cons);
// 				((z+j)->PrevFlowsol)->Momentumzstar= (double *)initialize_variable(i,(z+j),cons);
			 
			}
				
 			else
 			cons=0.0;	
				/*Memory allocation in the Flow Solution Node X Direction Velocity*/
 			((z+j)->Flowsol)->Momentumz= (double *)initialize_variable(i,(z+j),cons);
 			((z+j)->PrevFlowsol)->Momentumz= (double *)initialize_variable(i,(z+j),cons);
 			((z+j)->Flowsol)->Momentumzstar= (double *)initialize_variable(i,(z+j),cons);
 			((z+j)->PrevFlowsol)->Momentumzstar= (double *)initialize_variable(i,(z+j),cons);
 			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Momentumz=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumzOld=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->Momentumzstar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumzstarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumz = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumzstar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumzOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumzstarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
		
/*****************************END OF MEMORY ALLOCATION FOR Overall MOMENTUM IN X,Y,Z**************************** *******************************************************/


/*************Memory allocation of volume fraction for  each phase**********************************************/
	sumvol=0.0;
	
	
	for(nophase=0;nophase<nop;nophase++)
	{
	printf("\n Creating memory for volume fraction of phases");
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=34+nophase;//Variable Tag has to be stored 
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{
			printf("\n Enter the Initialization value of volume fraction  for phase %d in zone %s:",nophase,(z+j)->zonename);
			scanf("%lf",&VolFrac[j*nop+nophase]);
			fprintf(var2,"%lf\n",VolFrac[j*nop+nophase]);	
			
			/*Memory allocation in the Flow Solution Node y Direction momentum*/
			((z+j)->Flowsol)->VarthetaPhase[nophase]= (double *)initialize_variable(i,(z+j),VolFrac[j*nop+nophase]);
			((z+j)->PrevFlowsol)->VarthetaPhase[nophase]= (double *)initialize_variable(i,(z+j),VolFrac[j*nop+nophase]);
			//((z+j)->Flowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			//((z+j)->PrevFlowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->VarthetaPhase[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->VarthetaPhaseOld[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->Momentumystar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumystarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->VarthetaPhase[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumystar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->VarthetaPhaseOld[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumystarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		
				
			
		}
	}
		
	for(j=0;j<num_zone;j++)
	{
		sumvol=0.0;
		
		for(nophase=0;nophase<nop;nophase++)
		{	
			sumvol = sumvol+VolFrac[j*nop+nophase];
		}
		if(sumvol!=1)
		{
			printf("\n\n Give proper value of volume fraction so that the sum of volume fraction should be 1");
			exit(0);
		}
	}
	

	
	
/***********End of Memory allocation of volume fraction of each phase**********************************************/

/**********************Memory allocation for individual momentum in x direction**********************************/
	for(nophase=0;nophase<nop;nophase++)
	{
	printf("\n Creating memory for individual momentum  in x direction for phase %d:",nophase);
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=40+3*nophase+1;//Variable Tag has to be stored 
	cnt=cnt+1;
	
		for(j=0;j<num_zone;j++)
		{       
			
			printf("\nEnter the Initialization value of U[m/s] velocity for phase %d in zone %s:",nophase,(z+j)->zonename);
			scanf("%lf",&constantu[nophase]);
			fprintf(var2,"%lf\n",constantu[nophase]);
			
			/*Memory allocation in the Flow Solution Node X Direction momentum*/
			((z+j)->Flowsol)->Momentumofphasex[nophase]= (double *)initialize_variable(i,(z+j),constantu[nophase]);
			((z+j)->PrevFlowsol)->Momentumofphasex[nophase]= (double *)initialize_variable(i,(z+j),constantu[nophase]);
			//((z+j)->Flowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			//((z+j)->PrevFlowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Momentumofphasex[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumofphasexOld[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->Momentumxstar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumxstarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumofphasex[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumxstar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumofphasexOld[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumxstarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
					
/**********************Memory allocation for individual momentum in y direction**********************************/
	
	printf("\n Creating memory for individual momentum  in y direction for phase %d:",nophase);
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=40+3*nophase+2;//Variable Tag has to be stored 
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{       
			printf("\nEnter the Initialization value of V[m/s] velocity for phase %d in zone %s:",nophase,(z+j)->zonename);
			scanf("%lf",&constantv[nophase]);
			fprintf(var2,"%lf\n",constantv[nophase]);
			
			
			
			/*Memory allocation in the Flow Solution Node y Direction momentum*/
			((z+j)->Flowsol)->Momentumofphasey[nophase]= (double *)initialize_variable(i,(z+j),constantv[nophase]);
			((z+j)->PrevFlowsol)->Momentumofphasey[nophase]= (double *)initialize_variable(i,(z+j),constantv[nophase]);
			//((z+j)->Flowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			//((z+j)->PrevFlowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Momentumofphasey[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumofphaseyOld[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->Momentumystar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumystarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumofphasey[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumystar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumofphaseyOld[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumystarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}		
	
/**************************Memory allocation of individual momentum in z direction*********************************/ 	if(CgnsBase->BaseParameters->flag2d == 0)

	printf("\n Creating memory for individual momentum  in z direction for phase %d",nophase);
	i=1;
	fprintf(var,"%d\n",i);
	CgnsBase->VARTAG[cnt]=40+3*nophase+3;//Variable Tag has to be stored 
	cnt=cnt+1;
		for(j=0;j<num_zone;j++)
		{       
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
				printf("\nEnter the Initialization value of w[m/s] velocity for phase %d in zone %s:",nophase,(z+j)->zonename);
				scanf("%lf",&constantw[nophase]);
				fprintf(var2,"%lf\n",constantw[nophase]);
				
			}	
			else
			constantw[nophase]=0.0;
			/*Memory allocation in the Flow Solution Node X Direction momentum*/
				((z+j)->Flowsol)->Momentumofphasez[nophase]= (double *)initialize_variable(i,(z+j),constantw[nophase]);
				((z+j)->PrevFlowsol)->Momentumofphasez[nophase]= (double *)initialize_variable(i,(z+j),constantw[nophase]);
				//((z+j)->Flowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
				//((z+j)->PrevFlowsol)->Momentumxstar= (double *)initialize_variable(i,(z+j),constant);
			
			numiface=(z+j)->Number_of_Interfaces;
			//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	if(((z+j)->Iface+iface)->Direction==dir)
					{
						(((z+j)->Iface+iface)->ZSlab)->Momentumofphasez[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						(((z+j)->Iface+iface)->ZSlab)->MomentumofphasezOld[nophase]=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->Momentumzstar = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
						//(((z+j)->Iface+iface)->ZSlab)->MomentumzstarOld = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumofphasez[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->Momentumzstar = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumofphasezOld[nophase] = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						//(((z+j)->ZoneBc+zBc)->ZSlab)->MomentumzstarOld = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}	
	}
	
// 	printf("\n Creating memory for diffFlux\n\n");
// 	i=1;
// 	for(j=0;j<num_zone;j++)
// 	{	for(nophase=0;nophase<nop;nophase++)
// 		{
// 			((z+j)->Flowsol)->diffu[nophase] = (double *)initializeDiffFlux(i,(z+j));
// 			((z+j)->Flowsol)->diffv[nophase] = (double *)initializeDiffFlux(i,(z+j));
// 			((z+j)->Flowsol)->diffw[nophase] = (double *)initializeDiffFlux(i,(z+j));
// 		}	
// 	}
	
	
	//checking
	printf("\n Creating memory for diffFlux\n\n");
	i=1;
	for(j=0;j<num_zone;j++)
	{	
		
		((z+j)->Flowsol)->diffu = (double *)initializeDiffFlux(i,(z+j));
		((z+j)->Flowsol)->diffv = (double *)initializeDiffFlux(i,(z+j));
		((z+j)->Flowsol)->diffw = (double *)initializeDiffFlux(i,(z+j));
		
		(z+j)->DragX = (double *)initialize_variable(i,(z+j),0.0);
		(z+j)->DragY = (double *)initialize_variable(i,(z+j),0.0);
		(z+j)->DragZ = (double *)initialize_variable(i,(z+j),0.0);
		(z+j)->Dragcoeff = (double *)initialize_variable(i,(z+j),0.0);
		(z+j)->Bx = (double *)initialize_variable(i,(z+j),0.0);
			
		(z+j)->LiftX = (double *)initialize_variable(i,(z+j),0.0);
		(z+j)->LiftY = (double *)initialize_variable(i,(z+j),0.0);
		(z+j)->LiftZ = (double *)initialize_variable(i,(z+j),0.0);
		//(z+j)->Flowsol->ParticleDia =(double *)initialize_variable(i,(z+j),0.0);
		((z+j)->Flowsol)->MassfluxDisperse= (double *)initializeFlux(i,(z+j));
		((z+j)->PrevFlowsol)->MassfluxDisperse= (double *)initializeFlux(i,(z+j)); 
	}
	
	//CgnsBase->BaseParameters->flagVariableDia = 0;//This is default value to keep bubble dia constant
	//end 
/************End of memory allocation of individual momentum in x,y,z direction*******************************/	



	
	
	initializeVarMPFV(CgnsBase);
	
// 	fclose(var);
// 	fclose(var2);
	return(cnt);
}


