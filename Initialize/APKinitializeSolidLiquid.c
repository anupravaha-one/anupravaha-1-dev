#include "APKinitializeSolidLiquid.h"

/**
This function asks the user to identify which zones are solids and which zones are fluid.
 
If any zone is named as solid then this function initializes the flag for ConjugateProblem in the base to 0,indicating that
this is a conjugate problem requring special interface conditions.
*/
int initializeSolidLiquid(Base *CgnsBase)
{
        FILE *var2,*var3;
	char name[33];
	int num_zone,j,flag,i,totalcellFluid,*ptrSize,totalcellSolid;
	Zone *z;
	
        var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
	do
	{
		printf("\nIF ALL THE ZONES ARE FLUID,THEN [ENTER 0] OTHERWISE [ENTER 1]:");
		scanf("%d",&flag);
	}while(!(flag==1 || flag==0));
	fprintf(var2,"%d\n",flag);        
        fprintf(var3,"%d \tIF ALL THE ZONES ARE FLUID,THEN ENTER 0 OTHERWISE ENTER 1: \n\n",flag);

	if(flag==0)
	{
		z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
		num_zone=(CgnsBase)->Number_of_zones;
		
		for(j=0;j<num_zone;j++)
		{
			(z+j)->flagSolidFluid = 0;
			//printf("Zone name=%s\n",(z+j)->zonename);
			strcpy((z+j)->SolidFluidName,"FLUID") ;
			printf("Zonename: %s\tFluid/Solid type: %s\n",(z+j)->zonename,(z+j)->SolidFluidName);
			printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
		}
	}
	else
	{
		z=CgnsBase->CgnsZone; /* z is the pointer to first Zone*/
		num_zone=(CgnsBase)->Number_of_zones;
		
		for(j=0;j<num_zone;j++)
		{
			printf("Zone Name is %s\n",(z+j)->zonename);
			printf("NI=%d NJ=%d NK=%d\n",*((z+j)->Size+0),*((z+j)->Size+1),*((z+j)->Size+2));
			printf("IS THIS ZONE A SOLID OR FLUID(Enter 0->Fluid,1->Solid)\n");
			do
			{
				scanf("%d",&i);
			}while(!(i==1 || i==0));
			fprintf(var2,"%d\n",i);        
                        fprintf(var3,"%d \tIS THIS ZONE A SOLID OR FLUID(Enter 0->Fluid,1->Solid) \n\n",i);
			(z+j)->flagSolidFluid = i;
			
			if(i==0)
			{
				printf("THIS ZONE IS FLUID. ENTER THE TYPE OF IT AS F1,F2 etc\n");
				scanf("%s",name);			
                                fprintf(var2,"%s\n",name);           
                                fprintf(var3,"%s \tTHIS ZONE IS FLUID. ENTER THE TYPE OF IT AS F1,F2 etc \n\n",name);     
				printf("Fluid name=%s\n",name);
				strcpy((z+j)->SolidFluidName,name);
				printf("Zonename: %s\tFluid/Solid type: %s\n",(z+j)->zonename,(z+j)->SolidFluidName);
			}
			if(i==1)
			{
				printf("THIS ZONE IS SOLID. ENTER THE TYPE OF IT AS S1,S2 etc\n");
				scanf("%s",name);                                
                                fprintf(var2,"%s\n",name);           
                                fprintf(var3,"%s \tTHIS ZONE IS SOLID. ENTER THE TYPE OF IT AS S1,S2 etc \n\n",name);    
				printf("Fluid name=%s\n",name);
				strcpy((z+j)->SolidFluidName,name);
				printf("Zonename: %s\tFluid/Solid type: %s\n",(z+j)->zonename,(z+j)->SolidFluidName);
			}
		}
	}
	totalcellFluid=0;
	num_zone=(CgnsBase)->Number_of_zones;
	//printf("num_zone=%d\n",num_zone);
	z = (CgnsBase->CgnsZone);
	for(j=0;j<num_zone;j++)
	{
		ptrSize= (z+j)->Size;
		if((z+j)->flagSolidFluid==0)
		{
			totalcellFluid = totalcellFluid+(*(ptrSize+3)) * (*(ptrSize+4)) * (*(ptrSize+5));/*NcellI*NcellJ*NcellK*/
			//printf("Zname=%s\tNI=%d\tNJ=%d\tNK=%d\n",(z+j)->zonename,(*(ptrSize+3)),(*(ptrSize+4)),(*(ptrSize+5)));
			//printf("TotalCell=%d\n",totalcells);
		}
	}
	(CgnsBase->TotalCellFluid)=totalcellFluid;
	printf("TotalCellFluid=%d\n",(CgnsBase->TotalCellFluid));
	totalcellSolid=0;
	num_zone=(CgnsBase)->Number_of_zones;
	z = (CgnsBase->CgnsZone);
	for(j=0;j<num_zone;j++)
	{
		ptrSize= (z+j)->Size;
		if((z+j)->flagSolidFluid==1)
		{
			totalcellSolid = totalcellSolid+(*(ptrSize+3)) * (*(ptrSize+4)) * (*(ptrSize+5));/*NcellI*NcellJ*NcellK*/
			//printf("Zname=%s\tNI=%d\tNJ=%d\tNK=%d\n",(z+j)->zonename,(*(ptrSize+3)),(*(ptrSize+4)),(*(ptrSize+5)));
			//printf("TotalCell=%d\n",totalcells);
		}
	}
	(CgnsBase->TotalCellSolid)=totalcellSolid;
	printf("TotalcellSolid=%d\n",(CgnsBase->TotalCellSolid));
	(CgnsBase->TotalCellDomain)=totalcellSolid+totalcellFluid;
	printf("TotalCellDomain=%d\n",(CgnsBase->TotalCellDomain));
	//printf("TotalCell=%d\tafter filling=%d\n",totalcells,(CgnsBase->TotalCellDomain));
	if(((CgnsBase->TotalCellDomain)-(CgnsBase->TotalCellFluid))==0)
		(CgnsBase->flagConjugateProblem)=1;/*0->ConjugateProblem 1->Normal Flow Problem*/
	else
		(CgnsBase->flagConjugateProblem)=0;/*0->ConjugateProblem 1->Normal Flow Problem*/
	//printf("TotalCell=%d\tafter filling=%d\n",totalcells,(CgnsBase->TotalCellDomain));
	
        fclose(var2);
        fclose(var3);
	return (0);
}
