#include "APKmassVelocityOneIter.h"

/**
This function loops over all zones and calculates the mass velocities for only one iteration
(because of the coupling between the equations).This function calls the findResidual function 
to find residual for the G-S iterations. 
Totalresidual is used to find RMS residual for entire geometry rather than adding RMS 
residual of individual zones (as in APKmassVelocity.c).
*/


double massVelocityOneIter(Base *base,int variableTag)//supply 1 for u, 2 for v and  3 for w and no other value
{
	Zone *zone;
	int totalcells,numZones,nZones;
	double totalresidual,iterresidual;
	nZones=(base->Number_of_zones);
	totalcells=(base->TotalCellFluid);
	totalresidual = 0.0;
	
	if(base->BaseParameters->ConvectionScheme==1)//QUICK
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);
				if(zone->flagSolidFluid==0)
				{
					ExtractMassVelocities(variableTag*100,base,zone);
					
					iterresidual = findresidualQVK2d(base,zone,variableTag);
					
					totalresidual = totalresidual+iterresidual;
				}
			}
			totalresidual = sqrt(totalresidual/totalcells);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);
				if(zone->flagSolidFluid==0)
				{
					ExtractMassVelocities(variableTag*100,base,zone);
				
					iterresidual = findresidualQVK(base,zone,variableTag);
					
					totalresidual = totalresidual+iterresidual;
				}
			}
			totalresidual = sqrt(totalresidual/totalcells);
		}
	}
	else if(base->BaseParameters->ConvectionScheme==0)//CD
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);
				if(zone->flagSolidFluid==0)
				{
					ExtractMassVelocities(variableTag*100,base,zone);
					
					iterresidual = findresidualCD2d(base,zone,variableTag);
					
					totalresidual = totalresidual+iterresidual;
				}
			}
			totalresidual = sqrt(totalresidual/totalcells);
			
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);
				if(zone->flagSolidFluid==0)
				{
					ExtractMassVelocities(variableTag*100,base,zone);
					
					iterresidual = findresidualCD(base,zone,variableTag);
					
					totalresidual = totalresidual+iterresidual;
				}
			}
			totalresidual = sqrt(totalresidual/totalcells);
		}
	}
	
	RefreshMassVelocity(variableTag*100,base);
	
	return (totalresidual);
}

