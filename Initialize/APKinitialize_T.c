#include "APKinitialize_T.h"
int initialize_T(Base *CgnsBase,FILE *var2)
{
	FILE *var3;
	int j,num_zone,dir,iface,numiface,numBc,zBc,flag,flag_wh1=0;
	double constant;
	Zone *z;
	Interface *IfacePtr;
	var3=fopen("Settings/ProgramReadInputValues.dat","a");

		
	z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
	num_zone=(CgnsBase)->Number_of_zones;

	printf("\nCreating memory for Temperature ");
	if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
	   if( CgnsBase->Number_of_zones > 1)
	   {
		printf("\n Do you want to Initialize Temperature for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1):");
		scanf("%d",&flag);
		fprintf(var2,"%d\n",flag);
		fprintf(var3,"%d \tDo you want to Initialize Temperature for the WHOLE DOMAIN(Enter 0) OR ZONE BY ZONE(Enter 1): \n\n",flag);
	   } else { flag = 0;} 
	} else { flag = 0;}
	
	for(j=0;j<num_zone;j++)
	{
	  	(z+j)->ConjugateZoneTag = 0;

	   if ( (CgnsBase->BaseParameters->restart == 0) && (CgnsBase->BaseParameters->read_initialdata_cgns == 0) ) {
		if (flag != 0) {
			printf("\nZone Name is %s\n",(z+j)->zonename);
			printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
			printf("Enter the Initialization value of Temperature:\n");
			scanf("%lf",&constant);
			fprintf(var2,"%lf \n",constant);
			fprintf(var3,"%lf \tEnter the Initialization value of Temperature: \n\n",constant);
			//fprintf(var2,"%lf \t Enter the Initialization value of Temperature\n",constant);
		} else {
			if( CgnsBase->Number_of_zones == 1) {
				printf("\nThere is only one Zone in the Mesh. Zone Name is %s\n",(z+j)->zonename);
				printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
				printf("Enter the Initialization value of Temperature:\n");
				scanf("%lf",&constant);
				fprintf(var2,"%lf \n",constant);
				fprintf(var3,"%lf \tEnter the Initialization value of Temperature: \n\n",constant);
			} else {
				if (flag_wh1 == 0) {
					printf("Enter the Initialization value of Temperature for the Entire Domain:\n");
					scanf("%lf",&constant);
					fprintf(var2,"%lf \n",constant);
					fprintf(var3,"%lf \tEnter the Initialization value of Temperature for the Entire Domain: \n\n",constant);
					flag_wh1 = 1;
				}
			}
		} 
	   } else { constant = 0.0;}

		/*Memory allocation in the Flow Solution Node For Temperature*/
		((z+j)->Flowsol)->T= (double *)initialize_variable(1,(z+j),constant);
		((z+j)->PrevFlowsol)->T= (double *)initialize_variable(1,(z+j),constant);
		((z+j)->Flowsol)->enthalpy= (double *)initialize_variable(0,(z+j),constant);
		((z+j)->PrevFlowsol)->enthalpy= (double *)initialize_variable(0,(z+j),constant);

		numiface=(z+j)->Number_of_Interfaces;
		numBc=(z+j)->Number_of_BC_Patches;
			
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				IfacePtr = (z+j)->Iface+iface;
				if(IfacePtr->Direction==dir)
				{
					//Memory allocation in the Interface Slab for Temperature
					(IfacePtr->ZSlab)->T = (double *)initSlab(1,dir,IfacePtr->Point_Range,constant);
					(IfacePtr->ZSlab)->TOld = (double *)initSlab(1,dir,IfacePtr->Point_Range,constant);
					(IfacePtr->ZSlab)->enthalpy = (double *)initSlab(0,dir,IfacePtr->Point_Range,constant);
						
						if(IfacePtr->ConjugateIfaceTag == 1 && (CgnsBase->flagTurbModel == 1 || CgnsBase->flagTurbModel == 3))
						(z+j)->dummySlabTemp[dir-1] = (double *)initSlab(1,dir,IfacePtr->Point_Range,constant);
					else
					(z+j)->dummySlabTemp[dir-1] = (IfacePtr->ZSlab)->T;
				}
			}

			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((z+j)->ZoneBc+zBc)->Direction==dir)
				{
					//Memory allocation in the Boundary Slab for Temperature
					(((z+j)->ZoneBc+zBc)->ZSlab)->T = (double *)initSlab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpy = (double *)initSlab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					(((z+j)->ZoneBc+zBc)->ZSlab)->GradT = (double *)initSlab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					(((z+j)->ZoneBc+zBc)->ZSlab)->TOld = (double *)initSlab(1,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
					(((z+j)->ZoneBc+zBc)->ZSlab)->enthalpyOld = (double *)initSlab(0,dir,((z+j)->ZoneBc+zBc)->Point_Range,constant);
						
						(z+j)->dummySlabTemp[dir-1] = (((z+j)->ZoneBc+zBc)->ZSlab)->T;
				}
			}
		}
	}

        fclose(var3);
	return 0;
}



