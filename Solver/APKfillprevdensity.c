#include "APKfillprevdensity.h"

/**
This function updates the old (previous time step) values to new (current time step) values
of all variables to be solved. This function loops over all the cell centers and updates the
old values with the recent values of the corresponding variable.
*/

int fillprevdensity(Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k;
	Zone *zone;
	
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);		
		/**
		For Variable density formulation in solidification problems.
		As Solversemicoupled is used for solidification problems,
		The density of the previous time step is same as the current step.
		So density has been stored for 'n-1' time step also.**/
		if((zone->Flowsol->flagRho)==1 && (base->flagSolidification==1))
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*((zone->PrevFlowsol)->Rho+index) =*((zone->Flowsol)->Rho+index);
					*((zone->PrevFlowsol)->Rho_prev+index) =*((zone->Flowsol)->Rho+index);
			}
		}
	}
	return 0;
}
