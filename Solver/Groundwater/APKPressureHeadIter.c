#include "APKPressureHeadIter.h"

/**
This function solves for pressurehead in the current iteration . It calculate the residual in all zone and 
check whther connvergence is reached or not ..It call UpdateWaterContent and RefreshVariableCoefficients, to update value of water content and unsaturated permeabilty. Similar to this for 3D function PressureHeadIter3D() is called
*/ // these have written for pressure head
int PressureHeadIter(Base *base)
{
	FILE *pressConvcri;
	Zone *zone;
	int pressureItercase;
	
	double residBase,residZone;
	
	int totalBaseCells=0,iZone,numZones;	
	
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	pressureItercase = 0;  
	do
	{	
		residBase=0.0;

		pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
		
		fscanf(pressConvcri,"%lf",&(base->BaseParameters-> HeadLoopConv));// to make changes here
		fclose(pressConvcri);
		
		for(iZone=0;iZone<numZones;iZone++)
		{      
			zone=base->CgnsZone+iZone;
			ExtractVariables(19,base,zone);
			
			residZone=UpdatePressureHead2d(base,zone);

			residBase=residBase+residZone;
		}
	    
	        if((base->flagSaturation==1||base->flagSaturation==2))
		{
			UpdateWaterContent(base);
			RefreshVariableCoefficients(719,base);
		}
	
				
		residBase=(residBase/numZones);
		
		printf("Residual Pressure Head =%12.8e\n",residBase);
		
	}while(residBase > base->BaseParameters->HeadLoopConv);
	
	RefreshSurfaceCells(19,base);
		
	return 0;
  }
  
  
int PressureHeadIter3D(Base *base)
{
	FILE *pressConvcri;
	Zone *zone;
	
	double residBase,residZone;
	
	int totalBaseCells=0,iZone,numZones;	
		
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
		
	
	do
	{
		residBase=0.0;

		pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
		fscanf(pressConvcri,"%lf",&(base->BaseParameters-> HeadLoopConv));// to make changes here
		fclose(pressConvcri);
		
		for(iZone=0;iZone<numZones;iZone++)
		{       
			zone=base->CgnsZone+iZone;

			ExtractVariables(19,base,zone);
			
			residZone=UpdatePressureHead3D(base,zone);

			residBase=residBase+residZone;
		}
	       if(base->flagSaturation==1)
		{
			UpdateWaterContent(base);
			RefreshVariableCoefficients(719,base);
		}
						
		residBase=(residBase/numZones);
		printf("\n Residual Pressure Head = %12.8e",residBase);
		
	}while(residBase > base->BaseParameters->HeadLoopConv);
		
	RefreshSurfaceCells(19,base);
		
	return 0;
}


    
  
