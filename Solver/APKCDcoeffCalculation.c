#include "APKCDcoeffCalculation.h"

/**
It calls the functions which calculate convection-diffusion equation coefficients and deffered correction coefficients
depending time stepping scheme viz., BTCS, Crank Nicholson.
*/

int CDcoeffCalculation(Base *base,int variableTag)
{
//	printf("Entered CDcoeffCalculation\n");
	Zone *zone;
	int nZones,iZone,i;

	int DEFCOR_TAG=0;

	nZones=(base->Number_of_zones);
	

	if(variableTag==1 || variableTag==2 || variableTag==3)
		DEFCOR_TAG=1;
	else if(variableTag==31 || variableTag==32 || variableTag==33)
		DEFCOR_TAG=1;
	
	if(variableTag==9)
		DEFCOR_TAG=1;
	
	for(iZone=0;iZone<nZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		i = 0;
		if(base->BaseParameters->TimeStepingSchemeTag == 0)
		{
			if(zone->flagSolidFluid==0)
			{
				//FineGridCoefficients1(base,zone,i,variableTag);
				Conv_Diff_Assembly(base,zone,i,variableTag);

				if(DEFCOR_TAG==1)//This ensures that deferred correction is not called for variables that do not need it
				{
				
					if(base->BaseParameters->ConvectionScheme==1)
					{
						FineGridDefQVKCoefficients1(base,zone,variableTag);
					}
					else if(base->BaseParameters->ConvectionScheme==0)
					{
						FineGridDefCDCoefficients1(base,zone,variableTag);
					}
				}
			}
			else
				FineGridCoefficients1_solid(base,zone,i,variableTag);
			
		}
		else if(base->BaseParameters->TimeStepingSchemeTag == 1)
		{
			if(zone->flagSolidFluid==0)
			{
				//FineGridCoefficients2(base,zone,i,variableTag);
				Conv_Diff_Assembly(base,zone,i,variableTag);
		
				if(DEFCOR_TAG==1)//This ensures that deferred correction is not called for variables that do not need it
				{
				
					if(base->BaseParameters->ConvectionScheme==1)
					{
						FineGridDefQVKCoefficients1(base,zone,variableTag);
					}
					else if(base->BaseParameters->ConvectionScheme==0)
					{
						FineGridDefCDCoefficients1(base,zone,variableTag);
					}
				}
			}
			else
				FineGridCoefficients1_solid(base,zone,i,variableTag);
		}
		else
		{
			printf("correct tag is not stored in BaseParameters->TimeStepingSchemeTag\n");
			printf("exiting in UpdateVar\n");
			exit(1);
		}
	}
	return 0;
}


