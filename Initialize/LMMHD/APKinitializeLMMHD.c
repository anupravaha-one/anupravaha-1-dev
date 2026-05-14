#include "APKinitializeLMMHD.h"


double *non_uniform(int Findex,Zone *z,int Fnum);
double *non_uniform_slab(int j,int dir,Interface *Iface,Base *CgnsBase,int Findex);
double *non_uniform_slabBC(int j,int dir,ZoneBoundary *ZoneBc,Base *CgnsBase,int Findex);
void *UFnc(int Findex,double *Field,int B,int i,int j,int k,Zone *z,double centre[3]);

int initializeLMMHD(Base *CgnsBase, int cnt)
{
	FILE *var,*var2,*var3;
	int i,j,num_zone,p,flag;
	int MHDconst,Findex; 					//Added by Kuldhir
	int numiface,numBc,dir,iface,zBc,varCoeffCount=0;
	Zone *z;
	double constant,constantX,constantY,constantZ,sigma,sigmawall,t_wall,L,constantB;//rajendra
	
	var=fopen("Settings/variablelist.dat","a");
	
	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones; 
	printf("\nCreating memory for Electric potentialMHD[Volts]:");
	CgnsBase->VARTAG[cnt]=17;//Variable Tag has to be stored 
	cnt=cnt+1;
	i=1;
	fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		/*Memory allocation in the Flow Solution Node For PotentialMHD*/
		((z+j)->Flowsol)->PotentialMHD= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->PrevFlowsol)->PotentialMHD= (double *)initialize_variable(i,(z+j),constant);
		((z+j)->Flowsol)->PotentialMHDCoeff= (double *)malloc(sizeof(double));
		*(((z+j)->Flowsol)->PotentialMHDCoeff)=1.0;

		
				/** These two variables are added by rajendra to improve the efficiency of souce calculation for LMMHD problems.  The source were calculated unnecessilarily every time.  This has been avoided by storing the sources in these variables rather than calculating them agin*/
				
				((z+j)->Flowsol)->LMMHD_Source[0]=(double *)initialize_variable(i,(z+j),0.0);
				((z+j)->Flowsol)->LMMHDmass_Source[0]=(double *)initialize_variable(i,(z+j),0.0);
				((z+j)->Flowsol)->LMMHD_Source[1]=(double *)initialize_variable(i,(z+j),0.0);
				((z+j)->Flowsol)->LMMHDmass_Source[1]=(double *)initialize_variable(i,(z+j),0.0);
				
				if(CgnsBase->BaseParameters->flag2d == 0)
				{	
					((z+j)->Flowsol)->LMMHD_Source[2]=(double *)initialize_variable(i,(z+j),0.0);
					((z+j)->Flowsol)->LMMHDmass_Source[2]=(double *)initialize_variable(i,(z+j),0.0);
				}	
		
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
				/*Memory allocation in the Interface Slab for Potential*/
				  (((z+j)->Iface+iface)->ZSlab)->PotentialMHD = (double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
				    (((z+j)->Iface+iface)->ZSlab)->PotentialMHDold = (double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
				}
			}
		
			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				/*Memory allocation in the Boundary Slab for Potential*/
				  (((z+j)->ZoneBc+zBc)->ZSlab)->PotentialMHD = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
				    (((z+j)->ZoneBc+zBc)->ZSlab)->PotentialMHDold = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
				}
			}
		}
	}
	



/*This LMMHD case is modified by Kuldhir on 11 Jan 2010 to incorporate variable magnetic field case*/
if(CgnsBase->flaginduced==0)
{
	printf("\nCreating memory for X-component of Imposed Magnetic field:");
	CgnsBase->VARTAG[cnt]=25;//Variable Tag has to be stored 
	//cnt=cnt+1;
	//i=1;
	for(j=0;j<num_zone;j++)
	{       
		if (num_zone > 1) printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		printf("\n Do you want to enter constant(0) Magnetic field or function (1)? Enter your option here:");
		scanf("%d",&MHDconst);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",MHDconst);
		fprintf(var3,"%d \tDo you want to enter constant(0) Magnetic field or function (1)? Enter your option here:\n\n",MHDconst);
		fclose(var2);
		fclose(var3);

		if(MHDconst==0)
		{
			printf("\n Enter X-component of Imposed Magnetic field[Tesla]:");
			scanf("%lf",&constantX);
			((z+j)->Flowsol)->MfieldX= (double *)initialize_variable(i,(z+j),constantX);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constantX);
			fprintf(var3,"%lf \tEnter X-component of Imposed Magnetic field[Tesla]: \n\n",constantX);
			fclose(var2);
			fclose(var3);
		}
		else if(MHDconst==1)
		{
			printf("\n Enter function index for X-component of Imposed Magnetic field[Tesla]:");
			scanf("%d",&Findex);
			((z+j)->Flowsol)->MfieldX= (double *)non_uniform(Findex,z+j,1);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",Findex);
			fprintf(var3,"%d \tEnter function index for X-component of Imposed Magnetic field[Tesla]: \n\n",Findex);
			fclose(var2);
			fclose(var3);
		} 
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					if(MHDconst==0)
					  (((z+j)->Iface+iface)->ZSlab)->MfieldX=(double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constantX);
					else if(MHDconst==1)
					{
						(((z+j)->Iface+iface)->ZSlab)->MfieldX=(double *)non_uniform_slab(j,dir,(z+j)->Iface+iface,CgnsBase,Findex);
					}
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					if(MHDconst==0)
					  (((z+j)->ZoneBc+zBc)->ZSlab)->MfieldX = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantX);
					else if(MHDconst==1)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->MfieldX= (double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
					}
				}
			}
		}
	}
	
	printf("\nCreating memory for Y-component of Imposed Magnetic field:");
	CgnsBase->VARTAG[cnt]=26;//Variable Tag has to be stored 
	//cnt=cnt+1;
	//i=1;
	for(j=0;j<num_zone;j++)
	{       
		if (num_zone > 1) printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		printf("\n Do you want to enter constant(0) Magnetic field or function (1)? Enter your option here:");
		scanf("%d",&MHDconst);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",MHDconst);
		fprintf(var3,"%d \tDo you want to enter constant(0) Magnetic field or function (1)? Enter your option here:\n\n",MHDconst);
		fclose(var2);
		fclose(var3);

		if(MHDconst==0)
		{
			printf("\n Enter Y-component of Imposed Magnetic field[Tesla]:");
			scanf("%lf",&constantY);
			((z+j)->Flowsol)->MfieldY= (double *)initialize_variable(i,(z+j),constantY);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constantY);
			fprintf(var3,"%lf \tEnter Y-component of Imposed Magnetic field[Tesla]: \n\n",constantY);
			fclose(var2);
			fclose(var3);
		}
		else if(MHDconst==1)
		{
			printf("\n Enter function index for Y-component of Imposed Magnetic field[Tesla]:");
			scanf("%d",&Findex);
			((z+j)->Flowsol)->MfieldY= (double *)non_uniform(Findex,z+j,2);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",Findex);
			fprintf(var3,"%d \tEnter function index for Y-component of Imposed Magnetic field[Tesla]: \n\n",Findex);
			fclose(var2);
			fclose(var3);
		} 
		
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					if(MHDconst==0)
					  (((z+j)->Iface+iface)->ZSlab)->MfieldY=(double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constantY);
					else if(MHDconst==1)
					{
						(((z+j)->Iface+iface)->ZSlab)->MfieldY=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
					}
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					if(MHDconst==0)
					  (((z+j)->ZoneBc+zBc)->ZSlab)->MfieldY = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantY);
					else if(MHDconst==1)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->MfieldY = (double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
					}
				}
			}
		}
	}
			
	printf("\nCreating memory for Z-component of Imposed Magnetic field:");
	CgnsBase->VARTAG[cnt]=27;//Variable Tag has to be stored 
	//cnt=cnt+1;
	//i=1;
	for(j=0;j<num_zone;j++)
	{       
		if (num_zone > 1) printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		printf("\n Do you want to enter constant(0) Magnetic field or function (1)? Enter your option here:");
		scanf("%d",&MHDconst);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",MHDconst);
		fprintf(var3,"%d \tDo you want to enter constant(0) Magnetic field or function (1)? Enter your option here:\n\n",MHDconst);
		fclose(var2);
		fclose(var3);

		if(MHDconst==0)
		{
			printf("\n Enter Z-component of Imposed Magnetic field[Tesla]:");
			scanf("%lf",&constantZ);
			((z+j)->Flowsol)->MfieldZ= (double *)initialize_variable(i,(z+j),constantZ);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constantZ);
			fprintf(var3,"%lf \tEnter Z-component of Imposed Magnetic field[Tesla]: \n\n",constantZ);
			fclose(var2);
			fclose(var3);
		}
		else if(MHDconst==1)
		{
			printf("\n Enter function index for Z-component of Imposed Magnetic field[Tesla]:");
			scanf("%d",&Findex);
			((z+j)->Flowsol)->MfieldZ= (double *)non_uniform(Findex,z+j,3);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",Findex);
			fprintf(var3,"%d \tEnter function index for Z-component of Imposed Magnetic field[Tesla]: \n\n",Findex);
			fclose(var2);
			fclose(var3);
		} 
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					if(MHDconst==0)
					  (((z+j)->Iface+iface)->ZSlab)->MfieldZ=(double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constantZ);
					else if(MHDconst==1)
					{
						(((z+j)->Iface+iface)->ZSlab)->MfieldZ=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
					}
				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					if(MHDconst==0)
					  (((z+j)->ZoneBc+zBc)->ZSlab)->MfieldZ = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constantZ);
					else if(MHDconst==1)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->MfieldZ = (double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
					}
				}
			}
		}
	}
}
/***************************MEMORY ALLOCATION FOR ELECTRICAL CONDUCTIVITY*******************************************************/
	/*The memory allocation for ELECTRICAL CONDUCTIVITY for Flowsol node and Slabs at boundary and Interface*/
	printf("\nCreating memory for Electrical conductivity of liquid metal:");
	if( CgnsBase->Number_of_zones > 1)
	{
		printf("\n Do you want to specify electrical conductivity for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1):");
		scanf("%d",&flag);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",flag);
		fprintf(var3,"%d \tDo you want to specify electrical conductivity for the WHOLE DOMAIN(Enter 0) or ZONE BY ZONE(Enter 1): \n\n",flag);
		fclose(var2);
		fclose(var3);
	} else { flag = 0;}

	if(flag==0) /*Initializing electrical conductivity for whole domain at a time*/
	{
		printf("\n If electrical conductivity is a CONSTANT, Enter 0");
		printf("\n If electrical conductivity is a VARIABLE, Enter 1:");
		scanf("%d",&p);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",p);
		fprintf(var3,"%d \tIf electrical conductivity is a CONSTANT, Enter 0 / is a VARIABLE, Enter 1: \n\n",p);
		fclose(var2);
		fclose(var3);

		varCoeffCount = 5;  /*5->electrical conductivity*/
		if(p==0)
		{
			printf("\nEnter the constant value of electrical conductivity(mho/m):");
			scanf("%lf",&constant);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constant);
			fprintf(var3,"%lf \tEnter the constant value of electrical conductivity(mho/m): \n\n",constant);
			fclose(var2);
			fclose(var3);
		}
		else
			constant=0.0;

		for(j=0;j<num_zone;j++)
		{       
			
			((z+j)->Flowsol)->Eleconductivity= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->Eleconductivity= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flageleccond=p;
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
						(((z+j)->Iface+iface)->ZSlab)->Eleconductivity=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						//printf("I am leaving Iface\n");
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab X Direction Velocity*/
					//printf("I am going to zoneBc\n");
					(((z+j)->ZoneBc+zBc)->ZSlab)->Eleconductivity = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					
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
			printf("\n If Electrical conductivity is a CONSTANT, Enter 0");
			printf("\n If Electrical conductivity is a VARIABLE, Enter 1:");
			scanf("%d",&p);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",p);
			fprintf(var3,"%d \tIf Electrical conductivity is a CONSTANT, Enter 0 / is a VARIABLE, Enter 1: \n\n",p);
			fclose(var2);
			fclose(var3);

			varCoeffCount = 5;  /*5->Density*/
			if(p==0)
			{
				printf("\nEnter the constant value of electrical conductivity(mho/m):");
				scanf("%lf",&constant);

				var2=fopen("Settings/coefficients.dat","a");
				var3=fopen("Settings/ProgramReadInputValues.dat","a");
				fprintf(var2,"%lf\n",constant);
				fprintf(var3,"%lf \tEnter the constant value of electrical conductivity(mho/m): \n\n",constant);
				fclose(var2);
				fclose(var3);
			}
			else
				constant=0.0;
			
			((z+j)->Flowsol)->Eleconductivity= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->Eleconductivity= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
			((z+j)->Flowsol)->flageleccond=p;
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
						(((z+j)->Iface+iface)->ZSlab)->Eleconductivity=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						//printf("I am leaving Iface\n");
					}
				}

				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
					/*Memory allocation in the Interface Slab X Direction Velocity*/
					//printf("I am going to zoneBc\n");
					(((z+j)->ZoneBc+zBc)->ZSlab)->Eleconductivity = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					//printf("I am leaving zoneBc\n");
					}
				}
			}
		}
	}
/***************************** END OF MEMORY ALLOCATION FOR ELECTRICAL CONDUCTIVITY*******************************************************/



/*This LMMHD case is modified by Kuldhir on 11 Jan 2010 to incorporate variable magnetic field case*/
	
if(CgnsBase->flaginduced==1)
{
	printf("\nCreating memory for X-component of Total Magnetic field:");
	CgnsBase->VARTAG[cnt]=56;//Variable Tag has to be stored //rajendra
	cnt=cnt+1;
	//i=1;
	for(j=0;j<num_zone;j++)
	{       
		if (num_zone > 1) printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		printf("\n Do you want to enter constant(0) Magnetic field or function (1)? Enter your option here:");
		scanf("%d",&MHDconst);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",MHDconst);
		fprintf(var3,"%d \tDo you want to enter constant(0) Magnetic field or function (1)? Enter your option here:\n\n",MHDconst);
		fclose(var2);
		fclose(var3);

		if(MHDconst==0)
		{
			printf("\nEnter the Initialization value of TotalMfieldX[Tesla]:");//rajendra
			scanf("%lf",&constantB);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constantB);
			fprintf(var3,"%lf \tEnter the Initialization value of TotalMfieldX[Tesla]: \n\n",constantB);
			fclose(var2);
			fclose(var3);

			((z+j)->Flowsol)->TotalMfieldX= (double *)initialize_variable(i,(z+j),constantB);
		 	((z+j)->PrevFlowsol)->TotalMfieldX= (double *)initialize_variable(i,(z+j),constantB);//rajendra
		}
		else if(MHDconst==1)
		{
			printf("\n Enter function index for X-component of Magnetic field[Tesla]:");
			scanf("%d",&Findex);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",Findex);
			fprintf(var3,"%d \tEnter function index for X-component of Imposed Magnetic field[Tesla]: \n\n",Findex);
			fclose(var2);
			fclose(var3);

			((z+j)->Flowsol)->TotalMfieldX= (double *)non_uniform(Findex,z+j,1);
			((z+j)->PrevFlowsol)->TotalMfieldX= (double *)non_uniform(Findex,z+j,1);
		}

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					if(MHDconst==0)
					{
					  (((z+j)->Iface+iface)->ZSlab)->TotalMfieldX=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					    (((z+j)->Iface+iface)->ZSlab)->TotalMfieldXold=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);//rajendra
					}
					else if(MHDconst==1)
					{
						(((z+j)->Iface+iface)->ZSlab)->TotalMfieldX=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
						(((z+j)->Iface+iface)->ZSlab)->TotalMfieldXold=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
					}
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					if(MHDconst==0)
					{
					  (((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldX = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					    (((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldXold = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
					else if(MHDconst==1)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldX=(double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
						(((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldXold=(double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
					}			
				}
			}
		}
	}
			
	printf("\nCreating memory for Y-component of Total Magnetic field:");
	CgnsBase->VARTAG[cnt]=57;//Variable Tag has to be stored //rajendra
	cnt=cnt+1;
	//i=1;
	for(j=0;j<num_zone;j++)
	{       
		if (num_zone > 1) printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		printf("\n Do you want to enter constant(0) Magnetic field or function (1)? Enter your option here:");
		scanf("%d",&MHDconst);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",MHDconst);
		fprintf(var3,"%d \tDo you want to enter constant(0) Magnetic field or function (1)? Enter your option here:\n\n",MHDconst);
		fclose(var2);
		fclose(var3);

		if(MHDconst==0)
		{
			printf("\nEnter the Initialization value of TotalMfieldY[Tesla]:");//rajendra	
			scanf("%lf",&constantB);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constantB);
			fprintf(var3,"%lf \tEnter the Initialization value of TotalMfieldY[Tesla]: \n\n",constantB);
			fclose(var2);
			fclose(var3);
		
			((z+j)->Flowsol)->TotalMfieldY= (double *)initialize_variable(i,(z+j),constantB);
		 	((z+j)->PrevFlowsol)->TotalMfieldY= (double *)initialize_variable(i,(z+j),constantB);//rajendra
		}
		else if(MHDconst==1)
		{
			printf("\n Enter function index for Y-component of Magnetic field[Tesla]:");
			scanf("%d",&Findex);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",Findex);
			fprintf(var3,"%d \tEnter function index for Y-component of Imposed Magnetic field[Tesla]: \n\n",Findex);
			fclose(var2);
			fclose(var3);

			((z+j)->Flowsol)->TotalMfieldY= (double *)non_uniform(Findex,z+j,1);
			((z+j)->PrevFlowsol)->TotalMfieldY= (double *)non_uniform(Findex,z+j,1);
		}
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					if(MHDconst==0)
					{
					  (((z+j)->Iface+iface)->ZSlab)->TotalMfieldY=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					    (((z+j)->Iface+iface)->ZSlab)->TotalMfieldYold=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
					else if(MHDconst==1)
					{
						(((z+j)->Iface+iface)->ZSlab)->TotalMfieldY=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
						(((z+j)->Iface+iface)->ZSlab)->TotalMfieldYold=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
					}
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					if(MHDconst==0)
					{
					  (((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldY = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					    (((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldYold = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);//rajendra
					}
					else if(MHDconst==1)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldY=(double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
						(((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldYold=(double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
					}
				}
			}
		}
	}
			
	printf("\nCreating memory for Z-component of Total Magnetic field:");
	CgnsBase->VARTAG[cnt]=58;//Variable Tag has to be stored //rajendra
	cnt=cnt+1;
	//i=1;
	for(j=0;j<num_zone;j++)
	{       
		if (num_zone > 1) printf("Zonename :%s \t NI=%d NJ=%d NK=%d\n",(z+j)->zonename,*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		printf("\n Do you want to enter constant(0) Magnetic field or function (1)? Enter your option here:");
		scanf("%d",&MHDconst);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%d\n",MHDconst);
		fprintf(var3,"%d \tDo you want to enter constant(0) Magnetic field or function (1)? Enter your option here:\n\n",MHDconst);
		fclose(var2);
		fclose(var3);

		if(MHDconst==0)
		{       
			printf("\nEnter the Initialization value of TotalMfieldZ[Tesla]:");//rajendra
			scanf("%lf",&constantB);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%lf\n",constantB);
			fprintf(var3,"%lf \tEnter the Initialization value of TotalMfieldZ[Tesla]: \n\n",constantB);
			fclose(var2);
			fclose(var3);
	
			((z+j)->Flowsol)->TotalMfieldZ= (double *)initialize_variable(i,(z+j),constantB);
	 		((z+j)->PrevFlowsol)->TotalMfieldZ= (double *)initialize_variable(i,(z+j),constantB);//rajendra
		}
		else if(MHDconst==1)
		{
			printf("\n Enter function index for Z-component of Magnetic field[Tesla]:");
			scanf("%d",&Findex);

			var2=fopen("Settings/coefficients.dat","a");
			var3=fopen("Settings/ProgramReadInputValues.dat","a");
			fprintf(var2,"%d\n",Findex);
			fprintf(var3,"%d \tEnter function index for Z-component of Imposed Magnetic field[Tesla]: \n\n",Findex);
			fclose(var2);
			fclose(var3);

			((z+j)->Flowsol)->TotalMfieldZ= (double *)non_uniform(Findex,z+j,1);
			((z+j)->PrevFlowsol)->TotalMfieldZ= (double *)non_uniform(Findex,z+j,1);
		}
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	
				if(((z+j)->Iface+iface)->Direction==dir)
				{
					if(MHDconst==0)
					{
					  (((z+j)->Iface+iface)->ZSlab)->TotalMfieldZ=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					    (((z+j)->Iface+iface)->ZSlab)->TotalMfieldZold=(double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);//rajendra
					}
					else if(MHDconst==1)
					{
						(((z+j)->Iface+iface)->ZSlab)->TotalMfieldZ=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
						(((z+j)->Iface+iface)->ZSlab)->TotalMfieldZold=(double *)non_uniform_slab(j,dir,((z+j)->Iface+iface),CgnsBase,Findex);
					}
				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					if(MHDconst==0)
					{
					  (((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldZ = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					    (((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldZold = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);//rajendra
					}
					else if(MHDconst==1)
					{
						(((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldZ=(double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
						(((z+j)->ZoneBc+zBc)->ZSlab)->TotalMfieldZold=(double *)non_uniform_slabBC(j,dir,((z+j)->ZoneBc+zBc),CgnsBase,Findex);
					}
				}
			}
		}
	}			
	
			printf("\nCreating memory for magnetic permeability of the liquid metal:");//rajendra
			if( CgnsBase->Number_of_zones > 1)
			{
				printf("\n Do you want to specify magnetic permeability for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
				scanf("%d",&flag);

				var2=fopen("Settings/coefficients.dat","a");
				var3=fopen("Settings/ProgramReadInputValues.dat","a");
				fprintf(var2,"%d\n",flag);
				fprintf(var3,"%d \tDo you want to specify magnetic permeability for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",flag);
				fclose(var2);
				fclose(var3);
			} else { flag = 0;}

			if(flag==0) 
			{
				printf("\n If magnetic permeability is a CONSTANT, Enter 0");
				printf("\n If magnetic permeability is a VARIABLE, Enter 1:");
				scanf("%d",&p);

				var2=fopen("Settings/coefficients.dat","a");
				var3=fopen("Settings/ProgramReadInputValues.dat","a");
				fprintf(var2,"%d\n",p);
				fprintf(var3,"%d \tIf magnetic permeability is a CONSTANT, Enter 0 / is a VARIABLE, Enter 1: \n\n",p);
				fclose(var2);
				fclose(var3);

				varCoeffCount = 6;  /*6->permeability*/
				if(p==0)
				{
					printf("\nEnter the constant value of magnetic permeability[VS/Am]:");
					scanf("%lf",&constant);
					
					var2=fopen("Settings/coefficients.dat","a");
					var3=fopen("Settings/ProgramReadInputValues.dat","a");
					fprintf(var2,"%lf\n",constant);
					fprintf(var3,"%lf \tEnter the constant value of magnetic permeability[VS/Am]: \n\n",constant);
					fclose(var2);
					fclose(var3);
				}
				else
					constant=0.0;

				for(j=0;j<num_zone;j++)
				{       
					((z+j)->Flowsol)->Permeability= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->Permeability= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagpermeability=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->Permeability=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
								//printf("I am leaving Iface\n");
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->Permeability = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
			//printf("\n Permeability=%e",*Permeability);
			if(flag==1)
			{
				for(j=0;j<num_zone;j++)
				{ 
					printf("Zone Name is %s\n",(z+j)->zonename);
					printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
					printf("\n If magnetic permeability is a CONSTANT, Enter 0");
					printf("\n If magnetic permeability is a VARIABLE, Enter 1:");
					scanf("%d",&p);

					var2=fopen("Settings/coefficients.dat","a");
					var3=fopen("Settings/ProgramReadInputValues.dat","a");
					fprintf(var2,"%d\n",p);
					fprintf(var3,"%d \tIf magnetic permeability is a CONSTANT, Enter 0 / is a VARIABLE, Enter 1: \n\n",p);
					fclose(var2);
					fclose(var3);

					varCoeffCount = 6;  /*6->permeability*/
					if(p==0)
					{
						printf("\nEnter the constant value of magnetic permeability[VS/Am]:");
						scanf("%lf",&constant);

						var2=fopen("Settings/coefficients.dat","a");
						var3=fopen("Settings/ProgramReadInputValues.dat","a");
						fprintf(var2,"%lf\n",constant);
						fprintf(var3,"%lf \tEnter the constant value of magnetic permeability[VS/Am]: \n\n",constant);
						fclose(var2);
						fclose(var3);
					}
					else
						constant=0.0;

					//printf("sigma in this zone =%lf\n",constant);
					((z+j)->Flowsol)->Permeability= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->PrevFlowsol)->Permeability= (double *)initializeVarCoefficient(p,(z+j),constant,varCoeffCount);
					((z+j)->Flowsol)->flagpermeability=p;
					numiface=(z+j)->Number_of_Interfaces;
					//printf("NO. of interface=%d No. of Boundary=%d\n",numiface,numBc);
					numBc=(z+j)->Number_of_BC_Patches;
					for(dir=1;dir<=6;dir++)
					{
						for(iface=0;iface<numiface;iface++)
						{	if(((z+j)->Iface+iface)->Direction==dir)
							{
								(((z+j)->Iface+iface)->ZSlab)->Permeability=(double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
							}
						}

						for(zBc=0;zBc<numBc;zBc++)
						{
							if(((z+j)->ZoneBc+zBc)->Direction==dir)
							{
							(((z+j)->ZoneBc+zBc)->ZSlab)->Permeability = (double *)initializeSlabVarCoeff(p,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
							//printf("I am leaving zoneBc\n");
							}
						}
					}
				}
			}
	printf("\n Creating memory for MagneticFlux");//rajendra
		i=1;
		for(j=0;j<num_zone;j++)
		{
			((z+j)->Flowsol)->Magneticflux= (double *)initializeFlux(i,(z+j));
			((z+j)->PrevFlowsol)->Magneticflux= (double *)initializeFlux(i,(z+j)); 
		}
		printf("\nCreating memory for Scalar Psi");//rajendra
		CgnsBase->VARTAG[cnt]=59;//Variable Tag has to be stored 
		cnt=cnt+1;
		i=1;
		for(j=0;j<num_zone;j++)
		{
			constant = 0.0;
			//Memory allocation in the Flow Solution Node For Psi
			((z+j)->Flowsol)->Psi= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->PrevFlowsol)->Psi= (double *)initialize_variable(i,(z+j),constant);
			((z+j)->Flowsol)->PsiCoeff= (double *)malloc(sizeof(double));
			*(((z+j)->Flowsol)->PsiCoeff)=1.0;
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{	
					if(((z+j)->Iface+iface)->Direction==dir)
					{
						//Memory allocation in the Interface Slab for Potential
					  (((z+j)->Iface+iface)->ZSlab)->Psi = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					    (((z+j)->Iface+iface)->ZSlab)->Psiold = (double *)initializeslab(i,dir,((z+j)->Iface+iface)->Point_Range);
					}
				}
		
				for(zBc=0;zBc<numBc;zBc++)
				{	
					if(((z+j)->ZoneBc+zBc)->Direction==dir)
					{
						//Memory allocation in the Boundary Slab for Potential
					  (((z+j)->ZoneBc+zBc)->ZSlab)->Psi = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					    (((z+j)->ZoneBc+zBc)->ZSlab)->Psiold = (double *)initializeslab(i,dir,((z+j)->ZoneBc+zBc)->Point_Range);
					}
				}
			}
		}
// 	
	printf("\nCreating memory for magnetic diffusivity of the liquid metal:");//rajendra
	
		{
         	for(j=0;j<num_zone;j++)
	        {
		  if(((z+j)->Flowsol)->flagpermeability == 0 && 	((z+j)->Flowsol)->flageleccond == 0 )
		    {
                      ((z+j)->Flowsol)->flagMagdiff=0;
		    }
		  else
                      ((z+j)->Flowsol)->flagMagdiff=1;
		
        	      varCoeffCount = 7; 
		
		      constant=0.0;
				
			((z+j)->Flowsol)->Magdiff= (double *)initializeVarCoefficient(0,(z+j),constant,varCoeffCount);
			((z+j)->PrevFlowsol)->Magdiff= (double *)initializeVarCoefficient(0,(z+j),constant,varCoeffCount);

			*((z+j)->Flowsol)->Magdiff= 1.0/(*(((z+j)->Flowsol)->Permeability)* (*(((z+j)->Flowsol)->Eleconductivity)));
			*((z+j)->PrevFlowsol)->Magdiff= 1.0/(*(((z+j)->PrevFlowsol)->Permeability)* (*(((z+j)->PrevFlowsol)->Eleconductivity)));
			
			numiface=(z+j)->Number_of_Interfaces;
			numBc=(z+j)->Number_of_BC_Patches;
			
				for(dir=1;dir<=6;dir++)
				{
					for(iface=0;iface<numiface;iface++)
					{
						if(((z+j)->Iface+iface)->Direction==dir)
						{
						(((z+j)->Iface+iface)->ZSlab)->Magdiff=(double *)initializeSlabVarCoeff(0,(z+j),constant,varCoeffCount,dir,((z+j)->Iface+iface)->Point_Range);
						  *(((z+j)->Iface+iface)->ZSlab)->Magdiff = 1.0/(*((((z+j)->Iface+iface)->ZSlab)->Permeability)* (*((((z+j)->Iface+iface)->ZSlab)->Eleconductivity)));
						}
					}

					for(zBc=0;zBc<numBc;zBc++)
					{
						if(((z+j)->ZoneBc+zBc)->Direction==dir)
						{
						(((z+j)->ZoneBc+zBc)->ZSlab)->Magdiff = (double *)initializeSlabVarCoeff(0,(z+j),constant,varCoeffCount,dir,((z+j)->ZoneBc+zBc)->Point_Range);
						  *(((z+j)->ZoneBc+zBc)->ZSlab)->Magdiff = 1.0/(*((((z+j)->ZoneBc+zBc)->ZSlab)->Permeability)* (*((((z+j)->ZoneBc+zBc)->ZSlab)->Eleconductivity)));
						}
					}
				}
		  }
		}
			
	
		
				
}	
	printf("\nCreating memory for UcrossBX:");
	CgnsBase->VARTAG[cnt]=22;//Variable Tag has to be stored 
	//cnt=cnt+1;
	i=1;
	//fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		((z+j)->Flowsol)->UcrossBX= (double *)initialize_variable(i,(z+j),constant);
			
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
				  (((z+j)->Iface+iface)->ZSlab)->UcrossBX= (double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				  (((z+j)->ZoneBc+zBc)->ZSlab)->UcrossBX = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
				}
			}
		}
	}
	printf("\nCreating memory for UcrossBY:");
	CgnsBase->VARTAG[cnt]=23;//Variable Tag has to be stored 
	//cnt=cnt+1;
	i=1;
	//fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		((z+j)->Flowsol)->UcrossBY= (double *)initialize_variable(i,(z+j),constant);
			
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
				  (((z+j)->Iface+iface)->ZSlab)->UcrossBY= (double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
				}
			}
		
			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				  (((z+j)->ZoneBc+zBc)->ZSlab)->UcrossBY = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
				}
			}
		}
	}
	printf("\nCreating memory for UcrossBZ:");
	CgnsBase->VARTAG[cnt]=24;//Variable Tag has to be stored 
	//cnt=cnt+1;
	i=1;
	//fprintf(var,"%d\n",i);
	for(j=0;j<num_zone;j++)
	{
		constant = 0.0;
		((z+j)->Flowsol)->UcrossBZ= (double *)initialize_variable(i,(z+j),constant);
		
		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{	if(((z+j)->Iface+iface)->Direction==dir)
				{
				  (((z+j)->Iface+iface)->ZSlab)->UcrossBZ= (double *)initslabmag(i,dir,((z+j)->Iface+iface)->Point_Range,constant);
				}
			}
	
			for(zBc=0;zBc<numBc;zBc++)
			{	if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
				  (((z+j)->ZoneBc+zBc)->ZSlab)->UcrossBZ = (double *)initslabmag(i,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
				}
			}
		}
	}
	//	printf("\nCreating memory for X component of electrical current density:");
	//	printf("\n Enter the initial value of X component of  electrical current density [A/m^2]:");
	//	scanf("%lf",&constant);
	for(j=0;j<num_zone;j++)
	 ((z+j)->Flowsol)->ElecurrentdensityX= (double *)initialize_variable(i,(z+j),0.0);
        
//printf("\nCreating memory for Y component of electrical current density:");
//	printf("\n Enter the initial value of Y component of  electrical current density [A/m^2]:");
//	scanf("%lf",&constant);
	for(j=0;j<num_zone;j++)
	 	((z+j)->Flowsol)->ElecurrentdensityY= (double *)initialize_variable(i,(z+j),0.0);	
      
//	printf("\nCreating memory for Z component of electrical current density:");
//	printf("\n Enter the initial value of Z component of  electrical current density [A/m^2]:");
//	scanf("%lf",&constant);
	for(j=0;j<num_zone;j++)
	 	((z+j)->Flowsol)->ElecurrentdensityZ= (double *)initialize_variable(i,(z+j),0.0);	
         
	printf("\n Creating memory for PotentialMHDpoissonRHS");
	for(j=0;j<num_zone;j++)
	{
		((z+j)->Flowsol)->PotentialMHDpoissonRHS= (double *)initializeFlux(i,(z+j));
		((z+j)->PrevFlowsol)->PotentialMHDpoissonRHS= (double *)initializeFlux(i,(z+j)); 
	}
	/*	printf("\n Do you want to enter conductivity and thickness of wall and characteristic length (YES-1,NO-0):");
	scanf("%d",&i);
	if(i==1)
	{
		printf("\n Enter the condctivity of wall:");
		scanf("%lf",&sigmawall);
		printf("\n Enter the thickness of wall:");
		scanf("%lf",&t_wall);
		printf("\n Enter the characteristic length:");
		scanf("%lf",&L);

		var2=fopen("Settings/coefficients.dat","a");
		var3=fopen("Settings/ProgramReadInputValues.dat","a");
		fprintf(var2,"%lf\n",sigmawall);
		fprintf(var3,"%lf \tEnter the condctivity of wall: \n\n",sigmawall);
		fprintf(var2,"%lf\n",t_wall);
		fprintf(var3,"%lf \tEnter the thickness of wall: \n\n",t_wall);
		fprintf(var2,"%lf\n",L);
		fprintf(var3,"%lf \tEnter the characteristic length:  \n\n",L);
		fclose(var2);
		fclose(var3);
		
		(CgnsBase->BaseParameters)->wallconductanceratio = (sigmawall * t_wall)/(sigma * L);
		}*/
			/*for(j=0;j<num_zone;j++)
			{
				if(constantX==0&&constantY==0&&constantZ==0)
					CgnsBase->flagLMMHD=0;
			}*/
	//cnt = cnt-6;//As we are not applying BC's for B and UXB variables we have to do this.
printf("\n count is %d",cnt);//rajendra
	return(cnt);
}





#include "APKFindCellCentre.h"
#include "APKCall_user_func.h"
double *non_uniform(int Findex,Zone *z,int Fnum)
{
	int c1,c2,c3,c4,c5;
	int ni,nj,nk,i,j,k,*ptrsize;
	double cen_x,cen_y,cen_z,value;
	double *crdX,*crdY,*crdZ,*Field;
	char name1[]="CoordinateX";
	char name2[]="CoordinateY";
      	char name3[]="CoordinateZ";
	ZoneBoundary *zbc;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3];

	ni=*(z->Size);
	nj=*(z->Size+1);
	nk=*(z->Size+2); 
	ptrsize=z->Size;
	zbc=z->ZoneBc;	

	Field=(double*)malloc((ni-1)*(nj-1)*(nk-1)*sizeof(double));
        if(Field==NULL)
	{
		printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
		exit(1);
	}
	
	c1=strcmp(z->CgnsGrid->Array_Name1,name1);
	c2=strcmp(z->CgnsGrid->Array_Name1,name2);
	c3=strcmp(z->CgnsGrid->Array_Name1,name3);
	c4=strcmp(z->CgnsGrid->Array_Name2,name2);
	c5=strcmp(z->CgnsGrid->Array_Name2,name3);
	if(c1==0)
	{
		crdX=z->CgnsGrid->X;    
		if(c4==0)
		{
			crdY=z->CgnsGrid->Y;
			crdZ=z->CgnsGrid->Z;    
		}
		else
		{
			crdY=z->CgnsGrid->Z;
			crdZ=z->CgnsGrid->Y;    
		}
	}
	if(c2==0)
	{
		crdX=z->CgnsGrid->Y;
		if(c5==0)
		{
			crdY=z->CgnsGrid->Z;
			crdZ=z->CgnsGrid->X;
		}
		else
		{
			crdY=z->CgnsGrid->X;
			crdZ=z->CgnsGrid->Z;
		}
	}
	if(c3==0)
	{
		crdX=z->CgnsGrid->Z;
		if(c4==0)
		{
			crdY=z->CgnsGrid->Y;
			crdZ=z->CgnsGrid->X;
		}
		else
		{
			crdY=z->CgnsGrid->X;
			crdZ=z->CgnsGrid->Y;
		}
	}

//printf("Findex=%d\t Fnum=%d\t c1to5=[%d,%d,%d,%d,%d]\t ni,nj,nk=[%d,%d,%d]\t",Findex,Fnum,c1,c2,c3,c4,c5,ni,nj,nk);
	
		
	for(k=0;k<nk-1;k++)
		for(j=0;j<nj-1;j++)
			for(i=0;i<ni-1;i++)
			{
				ptrPoint1X=crdX+ni*nj*k+ni*j+i;
				ptrPoint1Y=crdY+ni*nj*k+ni*j+i;
				ptrPoint1Z=crdZ+ni*nj*k+ni*j+i;
				FindCellCentre(ptrsize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
				cen_x=centroid[0];
				cen_y=centroid[1];
				cen_z=centroid[2];
				switch(Findex)
				{
					case 1:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_1(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 2:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_2(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 3:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_3(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 4:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_4(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 5:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_5(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 6:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_6(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 7:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_7(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 8:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_8(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 9:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_9(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 10:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_10(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 11:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_11(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 12:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_12(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 13:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_13(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 14:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_14(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 15:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_15(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
					case 16:
						*(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i)=user_func_16(cen_x,cen_y,cen_z,0.0,z,i,j,k,zbc);
						break;
				}
				//printf("\n\n i,j,k=[%d,%d,%d]\t centroid=[%.2lf,%.2lf,%.2lf]\t",i,j,k,cen_x,cen_y,cen_z);
				//printf("\n crdX,Y,Z=[%lf,%.2lf,%.2lf]\t Px,Py,Pz=[%.2lf,%.2lf,%.2lf]\t Field=%.2lf\n\n",*crdX,*crdY,*crdZ,*ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z, *(Field+(ni-1)*(nj-1)*k+(ni-1)*j+i));
			}
	return(Field);
}

#include "APKFindZonePtr.h"
#include "APKFindFaceCenterBC.h"
#include "APKGetSurfacePoints.h"
#include "APKTransformCellIndices1.h"
#include "APKAdjoiningZonePointerXYZ.h"
double *non_uniform_slab(int numZones,int dir,Interface *Iface,Base *CgnsBase,int Findex)
{
	Zone *ptrZone,*CurrentZone;
	ZoneBoundary *zbc;
	double surfPoints[4][3],*coordX,*coordY,*coordZ,*Field,centNeighbour[3],*ptrPoint[3];
	int i,j,k,Tag,temp[2][3],NI,NJ,NK,indices1[3],indices2[3];

	CurrentZone=(CgnsBase)->CgnsZone+numZones;

	NI=*(CurrentZone->Size+0);
	NJ=*(CurrentZone->Size+1);
	NK=*(CurrentZone->Size+2);

	coordX=(CurrentZone->CgnsGrid)->X;
	coordY=(CurrentZone->CgnsGrid)->Y; 
	coordZ=(CurrentZone->CgnsGrid)->Z;

	ptrZone=(Zone *)FindZonePtr(Iface->Adjoining_Zone_Name,CgnsBase);
	Tag=strcmp((CurrentZone->SolidFluidName),(ptrZone->SolidFluidName));

//printf("\n\n[NI,NJ,NK]=[%d,%d,%d]\tcorX,corY,corZ=[%.2lf,%.2lf,%.2lf]\t AdZoneName=%s\n\n" ,NI,NJ,NK, *coordX, *coordY, *coordZ,Iface->Adjoining_Zone_Name);

	switch(dir)
	{
		case 1://East 
			Field=(double*)malloc((NK-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][2]=NI-1;
			indices1[0]=NI;

			for(k=0;k<(NK-1);k++)
				for(j=0;j<(NJ-1);j++)
				{
					temp[1][0]=k,temp[1][1]=j;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					if(Tag!=0)
						FindFaceCenterBC(surfPoints,centNeighbour);
					else
					{
                				indices1[1]=j+1,indices1[2]=k+1;
                				TransformCellIndices1(indices1,Iface,indices2);
	        				AdjoiningZonePointerXYZ(ptrZone, indices2, ptrPoint);
	        				FindCellCentre(ptrZone->Size,ptrPoint[0], ptrPoint[1], ptrPoint[2], centNeighbour);
        				}
					UFnc(Findex,Field,NJ,NI-1,j,k, CurrentZone, centNeighbour);
//printf("\teast: [j,k]=[%d,%d]\t",j,k);
				}
			break;
		case 2://West
			Field=(double*)malloc((NK-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][2]=0;
			indices1[0]=0;

			for(k=0;k<(NK-1);k++)
				for(j=0;j<(NJ-1);j++)
				{
					temp[1][0]=k,temp[1][1]=j;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					if(Tag!=0)
						FindFaceCenterBC(surfPoints,centNeighbour);
					else
					{
                				indices1[1]=j+1,indices1[2]=k+1;
                				TransformCellIndices1(indices1,Iface,indices2);
	        				AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
	        				FindCellCentre(ptrZone->Size,ptrPoint[0],ptrPoint[1],ptrPoint[2], centNeighbour);
        				}
					UFnc(Findex,Field,NJ,0,j,k, CurrentZone, centNeighbour);
//printf("\twest: [j,k]=[%d,%d]\t",j,k);
				}
			break;
		case 3://North
			Field=(double*)malloc((NK-1)*(NI-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][2]=NJ-1;
			indices1[1]=NJ;

			for(k=0;k<(NK-1);k++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=k,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					if(Tag!=0)
						FindFaceCenterBC(surfPoints,centNeighbour);
					else
					{
                				indices1[0]=i+1,indices1[2]=k+1;
                				TransformCellIndices1(indices1,Iface,indices2);
	        				AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
	        				FindCellCentre(ptrZone->Size,ptrPoint[0],ptrPoint[1],ptrPoint[2], centNeighbour);
        				}
					UFnc(Findex,Field,NI,NJ-1,i,k, CurrentZone, centNeighbour);
//printf("\tnorth: [i,k]=[%d,%d]\t",i,k);
				}
			break;
		case 4://South
			Field=(double*)malloc((NK-1)*(NI-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][2]=0;
			indices1[1]=0;

			for(k=0;k<(NK-1);k++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=k,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					if(Tag!=0)
						FindFaceCenterBC(surfPoints,centNeighbour);
					else
					{
                				indices1[0]=i+1,indices1[2]=k+1;
                				TransformCellIndices1(indices1,Iface,indices2);
	        				AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
	        				FindCellCentre(ptrZone->Size,ptrPoint[0],ptrPoint[1],ptrPoint[2], centNeighbour);
        				}
					UFnc(Findex,Field,NI,0,i,k, CurrentZone, centNeighbour);
//printf("\tsouth: [i,k]=[%d,%d]\t",i,k);
				}
			break;
		case 5://Top
			Field=(double*)malloc((NI-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][2]=NK-1;
			indices1[2]=NK;

			for(j=0;j<(NJ-1);j++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=j,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					if(Tag!=0)
						FindFaceCenterBC(surfPoints,centNeighbour);
					else
					{
                				indices1[0]=i+1,indices1[1]=j+1;
                				TransformCellIndices1(indices1,Iface,indices2);
	        				AdjoiningZonePointerXYZ(ptrZone, indices2, ptrPoint);
	        				FindCellCentre(ptrZone->Size,ptrPoint[0], ptrPoint[1], ptrPoint[2], centNeighbour);
        				}
					UFnc(Findex,Field,NI,NK-1,i,j, CurrentZone, centNeighbour);
//printf("\ttop: [i,j]=[%d,%d]\t",i,j);
				}
			break;
		case 6://Bottom
			Field=(double*)malloc((NI-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][2]=0;
			indices1[2]=0;

			for(j=0;j<(NJ-1);j++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=j,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					if(Tag!=0)
						FindFaceCenterBC(surfPoints,centNeighbour);
					else
					{
                				indices1[0]=i+1,indices1[1]=j+1;
                				TransformCellIndices1(indices1,Iface,indices2);
	        				AdjoiningZonePointerXYZ(ptrZone, indices2, ptrPoint);
	        				FindCellCentre(ptrZone->Size,ptrPoint[0], ptrPoint[1], ptrPoint[2], centNeighbour);
        				}
					UFnc(Findex,Field,NI,0,i,j, CurrentZone, centNeighbour);
//printf("\tbottom: [i,j]=[%d,%d]\t",i,j);
				}
			break;
	}
	return(Field);

}

double *non_uniform_slabBC(int numZones,int dir,ZoneBoundary *ZoneBc,Base *CgnsBase,int Findex)
{
	Zone *ptrZone,*CurrentZone;
	ZoneBoundary *zbc;

	double surfPoints[4][3],*coordX,*coordY,*coordZ,*Field,centNeighbour[3],*ptrPoint[3];
	int i,j,k,Tag,temp[2][3],NI,NJ,NK,indices1[3],indices2[3];

	CurrentZone=(CgnsBase)->CgnsZone+numZones;

	NI=*(CurrentZone->Size+0);
	NJ=*(CurrentZone->Size+1);
	NK=*(CurrentZone->Size+2);

	coordX=(CurrentZone->CgnsGrid)->X;
	coordY=(CurrentZone->CgnsGrid)->Y; 
	coordZ=(CurrentZone->CgnsGrid)->Z;

//printf("\n\n[NI,NJ,NK]=[%d,%d,%d]\tcorX,corY,corZ=[%.2lf,%.2lf,%.2lf]\t \"Boundary\"\n\n" ,NI,NJ,NK, *coordX, *coordY, *coordZ);

	switch(dir)
	{
		case 1://East 
			Field=(double*)malloc((NK-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][2]=NI-1;
			indices1[0]=NI;

			for(k=0;k<(NK-1);k++)
				for(j=0;j<(NJ-1);j++)
				{
					temp[1][0]=k,temp[1][1]=j;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					
						FindFaceCenterBC(surfPoints,centNeighbour);
					
					UFnc(Findex,Field,NJ,NI-1,j,k, CurrentZone, centNeighbour);
//printf("\teast: [j,k]=[%d,%d]\t",j,k);
				}
			break;
		case 2://West
			Field=(double*)malloc((NK-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][2]=0;
			indices1[0]=0;

			for(k=0;k<(NK-1);k++)
				for(j=0;j<(NJ-1);j++)
				{
					temp[1][0]=k,temp[1][1]=j;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					
						FindFaceCenterBC(surfPoints,centNeighbour);
					
					UFnc(Findex,Field,NJ,0,j,k, CurrentZone, centNeighbour);
//printf("\twest: [j,k]=[%d,%d]\t",j,k);
				}
			break;
		case 3://North
			Field=(double*)malloc((NK-1)*(NI-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][2]=NJ-1;
			indices1[1]=NJ;

			for(k=0;k<(NK-1);k++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=k,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					
						FindFaceCenterBC(surfPoints,centNeighbour);
					
					UFnc(Findex,Field,NI,NJ-1,i,k, CurrentZone, centNeighbour);
//printf("\tnorth: [i,k]=[%d,%d]\t",i,k);
				}
			break;
		case 4://South
			Field=(double*)malloc((NK-1)*(NI-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][2]=0;
			indices1[1]=0;

			for(k=0;k<(NK-1);k++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=k,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					
						FindFaceCenterBC(surfPoints,centNeighbour);
					
					UFnc(Findex,Field,NI,0,i,k, CurrentZone, centNeighbour);
//printf("\tsouth: [i,k]=[%d,%d]\t",i,k);
				}
			break;
		case 5://Top
			Field=(double*)malloc((NI-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][2]=NK-1;
			indices1[2]=NK;

			for(j=0;j<(NJ-1);j++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=j,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					
						FindFaceCenterBC(surfPoints,centNeighbour);
					
					UFnc(Findex,Field,NI,NK-1,i,j, CurrentZone, centNeighbour);
//printf("\ttop: [i,j]=[%d,%d]\t",i,j);
				}
			break;
		case 6://Bottom
			Field=(double*)malloc((NI-1)*(NJ-1)*sizeof(double));
			if(Field==NULL)
			{
				printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
				exit(1);
			}

			temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][2]=0;
			indices1[2]=0;

			for(j=0;j<(NJ-1);j++)
				for(i=0;i<(NI-1);i++)
				{
					temp[1][0]=j,temp[1][1]=i;
					GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
					
						FindFaceCenterBC(surfPoints,centNeighbour);
					
					UFnc(Findex,Field,NI,0,i,j, CurrentZone, centNeighbour);
//printf("\tbottom: [i,j]=[%d,%d]\t",i,j);
				}
			break;
	}
	return(Field);

}

void *UFnc(int Findex,double *Field,int B,int a,int b,int c,Zone *z,double *centre)
{
	ZoneBoundary *zbc;
	double cen_x,cen_y,cen_z;

	cen_x=*(centre+0);
	cen_y=*(centre+1);
	cen_z=*(centre+2);

	switch(Findex)
	{
		case 1:
			*(Field+(B-1)*c+b)=user_func_1(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 2:
			*(Field+(B-1)*c+b)=user_func_2(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 3:
			*(Field+(B-1)*c+b)=user_func_3(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 4:
			*(Field+(B-1)*c+b)=user_func_4(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 5:
			*(Field+(B-1)*c+b)=user_func_5(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 6:
			*(Field+(B-1)*c+b)=user_func_6(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 7:
			*(Field+(B-1)*c+b)=user_func_7(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 8:
			*(Field+(B-1)*c+b)=user_func_8(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 9:
			*(Field+(B-1)*c+b)=user_func_9(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 10:
			*(Field+(B-1)*c+b)=user_func_10(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 11:
			*(Field+(B-1)*c+b)=user_func_11(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 12:
			*(Field+(B-1)*c+b)=user_func_12(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 13:
			*(Field+(B-1)*c+b)=user_func_13(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 14:
			*(Field+(B-1)*c+b)=user_func_14(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 15:
			*(Field+(B-1)*c+b)=user_func_15(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
		case 16:
			*(Field+(B-1)*c+b)=user_func_16(cen_x,cen_y,cen_z,0.0,z,a,b,c,zbc);
		break;
	}
	//printf("\n\n Findex=%d\t B=%d\t a,b,c=[%d,%d,%d]\t centroid=[%.2lf,%.2lf,%.2lf]\t Field=%.2lf\n\n",Findex,B,a,b,c,cen_x,cen_y,cen_z,*(Field+(B-1)*c+b));
}
