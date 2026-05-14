#include "APKRefreshMassVelocities.h"
/**
This function is called only once during restart for implementing the boundary conditions for mass velocities.
*/

int RefreshMassVelocities(Base *base)
{
	int cntZone,numZones;
	Zone *zone;
	
	numZones = base->Number_of_zones;
	for(cntZone=0;cntZone < numZones;cntZone++)
	{
		zone = base->CgnsZone+cntZone;
		//if( (zone->flagSolidFluid==0) ) <-- Condition not necessary, By default all slabs especially massflux, are initialized to 0.5.Hence excluding solids zones will make solid zone slabs to have a default value of 0.5
		if(base->flagMPFV==1)
		{
			ExtractMassVelocities(3100,base,zone);
			ExtractMassVelocities(3200,base,zone);
			if(base->BaseParameters->flag2d==0)
				ExtractMassVelocities(3300,base,zone);
		}
		else if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
		{
			ExtractMassVelocities(100,base,zone);
			ExtractMassVelocities(200,base,zone);
			ExtractMassVelocities(6200,base,zone);
			ExtractMassVelocities(6300,base,zone);
			if(base->BaseParameters->flag2d==0)
			{
				ExtractMassVelocities(300,base,zone);
				ExtractMassVelocities(6400,base,zone);
			}
		}
		else
		{
			ExtractMassVelocities(100,base,zone);
			ExtractMassVelocities(200,base,zone);
			if(base->BaseParameters->flag2d==0)
				ExtractMassVelocities(300,base,zone);
		}	
	}
	
	return 0;
}

/**
This function is called in every Gauss-Siedel iteration for implementing the boundary conditions for mass velocities.
NOTE: variableTag supplied to this function is 100 for Umass, 200 for Vmass and 300 for Wmass.
*/
//variableTag: 100->Umass, 200->Vmass and 300->Wmass 
int RefreshMassVelocity(int vartag,Base *base)
{
	int cntZone,numZones;
	
	numZones=base->Number_of_zones;
	for(cntZone=0;cntZone < numZones;cntZone++)
	{
		ExtractMassVelocities(vartag,base,base->CgnsZone+cntZone);
	}
	
	return 0;
}




//----------------------------------------------------------
//------------------ Earlier Version -----------------------
//----------------------------------------------------------
// #include "APKRefreshMassVelocities.h"
// 
// /**
// 
// This function is called only once during restart for implementing the boundary conditions
// for mass velocities.
// 
// */
// int RefreshMassVelocities(Base *base)/*supply 1 for u, 2 for v and  3 for w and no other value*/
// {
// 	int BC_index,nZones,numBC,vartag;
// 	Zone *zone;
// 	int numZones,num;
// 	
// 	ZoneBoundary *zonebc;
// 	nZones=(base->Number_of_zones);
// 	
// 	if(base->BaseParameters->flag2d==0)
// 		num=3;
// 	else
// 		num=2;
// 	for(vartag=1;vartag<=num;vartag++)
// 	{
// 		for(numZones=0;numZones<nZones;numZones++)
// 		{
// 			zone = (base->CgnsZone+numZones);
// 			//if( (zone->flagSolidFluid==0) ) <-- Condition not necessary, By default all slabs especially massflux, are initialized to 0.5.Hence excluding solids zones will make solid zone slabs to have a default value of 0.5
// 			{
// 				numBC = zone->Number_of_BC_Patches;
// 				for(BC_index=0;BC_index<numBC;BC_index++)
// 				{
// 					zonebc = (zone->ZoneBc+BC_index);
// 					velocitymassBc(base,vartag,zone,zonebc);
// 				}
// 			}
// 		}
// 	}
// 	return 0;
// }



