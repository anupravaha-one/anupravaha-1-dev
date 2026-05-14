#include "APKRefreshSurfaceCells.h"

/**
This function is used for applying boundary conditions and/or updating the interface values of
adjoining zones. This function inturn calls the "ExtractVariables" that will actually implement 
BC for a boundary or update interface values for an interface.
*/

int RefreshSurfaceCells(int variableTag,Base *base)
{
	int cntZone,numZones;
	
	//Added by chandrasekhar --- To update the boundary conditions from velocity to momentum in case of Turbulent MPFV
	int BC_index;
	ZoneBoundary *zbc;
	Zone *z;
	//End--chandrasekhar
	/*printf("base=%u\n",base);*/
	
	numZones=base->Number_of_zones;
	for(cntZone=0;cntZone < numZones;cntZone++)
	{
	  //	printf("\n i m in zone number %d\n",cntZone+1);
		ExtractVariables(variableTag,base,base->CgnsZone+cntZone);
	}
	
	//Added by chandrasekhar --- To update the boundary conditions from velocity to momentum in case of Turbulent MPFV
	for(cntZone=0;cntZone < numZones;cntZone++)
	{
		z=base->CgnsZone+cntZone;
		zbc=(base->CgnsZone+cntZone)->ZoneBc;
		//The updation from velocity to momentum is done only after updating the velocities and volume fractions,that's why this is called only after updating the vartag 46
		if((base->flagMPFV==1)&&(base->flagTurb==1)&&(variableTag==46)&&(base->masstag==0))
		{
			
			//printf("I am entering the updation of momentum at the boundaries for turbulence\n");
			for(BC_index=0;BC_index<(z->Number_of_BC_Patches);BC_index++)
				update_boundaryDMPFV(base,z,zbc+BC_index);
		}
	}
	//end--chandrasekhar
	return 0;
}


//Note:We can use the same approch that used to update the momentum from velocities with in the time loop. But here the 
//values of phase momentum has not been updated before updating the global momentum(this is important especially at inlet)
//And in the time loop as already we have updated the phase momentums at boundaries in APKMain.c itself we can update 
//the momentums immediately in boundary_d_t_sTurb using APKupdate_boundaryTurbDMPFV.c
