#include "APKPotentialMHDIter.h"

/**
This function solves for potentialMHD in the current iteration after the velocities have been converged. This function loops over all the zones and solves for potentialMHD until the total residue is less than the specified value.This function calls other functions for solving potential for different cases of geometry (2D or 3D) in which the actual G-S iterations are executed. */

double PotentialMHDIter(Base *base)
{
	//FILE *potentialMHDConvcri;
	Zone *zone;
	int potentialMHDIterCount=0;
	int potentialMHDItercase;
	double residBase,residZone,PotentialMHDLoopConv;
	int totalBaseCells=0,iZone,numZones;
	numZones=base->Number_of_zones;

        totalBaseCells = base->TotalCellDomain;

	PotentialMHDLoopConv=base->BaseParameters->PotentialMHDLoopConv;

	if(base->BaseParameters->flag2d == 1)
		potentialMHDItercase = 0; /*case for calling UpdatePotentialMHD2d*/
	
	else if(base->BaseParameters->flag2d == 0)
		potentialMHDItercase = 1; /*case for calling UpdatePotentialMHD, this is for normal 3d*/
	
	do
	{
		residBase=0.0;
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			
			{
				ExtractVariables(17,base,zone);

				switch(potentialMHDItercase)
				{
					case 0:
                                            	residZone=UpdateVarCD2d(base,zone,17);
                                              
					break;
					case 1: 
                                              
						residZone=UpdateVarCD(base,zone,17);
                                             
					break;
				}
				residBase+= residZone;
				
			}
		}
		
		//residBase=(residBase/numZones);
		residBase = sqrt(residBase/totalBaseCells);
		potentialMHDIterCount+= 1;
		if(potentialMHDIterCount%10000 == 0)
			printf("residual potentialMHD=%12.8e\n",residBase);
		
	}
	while(residBase > PotentialMHDLoopConv);
//base->BaseParameters->PotentialMHDLoopConv 
	//Older fixed criterion 1e-6 removed by Kuldhir
	
	RefreshSurfaceCells(17,base);
		
	return (residBase);
}
