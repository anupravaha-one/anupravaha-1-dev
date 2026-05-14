#include "APKmassVelocity.h"
/**
This function loops over all zones and calculates the mass velocities for
the current iteration.This function calls the findResidual function to find residual for
the G-S iterations. 
This function is called until the massflux is converged in every time step.
*/

int massVelocity(Base *base,int variableTag)/*supply 1 for u, 2 for v and  3 for w and no other value*/
{
	Zone *zone;
	int totalcells,numZones,nZones;
	double totalresidual,iterresidual;
	double convcri;
	
	nZones=(base->Number_of_zones);
	totalcells = base->TotalCellFluid;
	
	//NOTE: This will avoid problems even if one forgets inclusion of this statement before solving for the varialbe through massVelocity() <<< Narsimha 19th April, 2008
	//CDcoeffCalculation(CgnsRoot.CgnsBase,variableTag);
	
	base->masstag=1;//Added by chandrasekhar --- Turbulent MPFV
	convcri =find_convcri(variableTag,base);
	if(base->BaseParameters->ConvectionScheme==1)//QUICK
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				totalresidual = 0.0;
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
			while(totalresidual > convcri);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				totalresidual = 0.0;
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
			while(totalresidual > convcri);
		}
	}
	else if(base->BaseParameters->ConvectionScheme==0)//CD
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				totalresidual = 0.0;
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
			while(totalresidual > convcri);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				totalresidual = 0.0;
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
			while(totalresidual > convcri);
		}
	}
	
	RefreshMassVelocity(variableTag*100,base);
	
	return 0;
}
