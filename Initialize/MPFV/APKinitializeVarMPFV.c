#include "APKinitializeVarMPFV.h"
/**
 This subroutine is used specially for the global momentum i.e, for the overall momentum of phases by taking sum over all phases of the multiphase flow individual momentum and  the value  for the global momentum is filled in the memory  alloted for the same. Individual momentum value is also calculated from density ,velocity and volume fraction of each phase .Initially in individual momentum velocity  is filled which is altered in this function  in order to get phase momentum and global momentum
 */
double initializeVarMPFV(Base *CgnsBase)

{	
	int i,j,k,NI,NJ,NK,numZones,cntZ,*ptrSize,nop,che;
	Zone *zone;
	//double *variable;
	int totalcells,index,nophase,flagrhophase=0;
	numZones= CgnsBase->Number_of_zones;
	
	nop=CgnsBase->nphase;
	
	EvaluateVarCoefficients(CgnsBase);
	RefreshVariableCoefficients(710,CgnsBase);
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (CgnsBase->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);

	
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
			
        	/*************Inividual momentum calculation in x,y,z direction*************************************/
//Here u,v,w is stored in individual momentum memory itself so no need to multiply all this equation with velocities of each phase.			
			for(nophase=0;nophase<nop;nophase++)
			{
			flagrhophase=(CgnsBase->BaseParameters->flagRhoPhase[nophase]);
			
			*((zone->Flowsol)->Momentumofphasex[nophase]+index)=(*((zone->Flowsol)->VarthetaPhase[nophase]+index)*(*((zone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((zone->Flowsol)->Momentumofphasex[nophase]+index)));
			
			*((zone->Flowsol)->Momentumofphasey[nophase]+index)=(*((zone->Flowsol)->VarthetaPhase[nophase]+index)*(*((zone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((zone->Flowsol)->Momentumofphasey[nophase]+index)));
			
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
					*((zone->Flowsol)->Momentumofphasez[nophase]+index)=(*((zone->Flowsol)->VarthetaPhase[nophase]+index)*(*((zone->Flowsol)->RhoofPhase[nophase]+index*flagrhophase))*(*((zone->Flowsol)->Momentumofphasez[nophase]+index)));
			}
			
			
			}
			
			*((zone->Flowsol)->Momentumx+index)=0.0;
			*((zone->Flowsol)->Momentumy+index)=0.0;
			
			
/************************overall momentum calculation from individual momentum***************************************/ 
			
			for(nophase=0;nophase<nop;nophase++)
			{
			*((zone->Flowsol)->Momentumx+index)=*((zone->Flowsol)->Momentumx+index)+(*((zone->Flowsol)->Momentumofphasex[nophase]+index));
			
			
			*((zone->Flowsol)->Momentumy+index)=*((zone->Flowsol)->Momentumy+index)+(*((zone->Flowsol)->Momentumofphasey[nophase]+index));
			}
			
			if(CgnsBase->BaseParameters->flag2d == 0)
			{	*((zone->Flowsol)->Momentumz+index)=0.0;
				
				for(nophase=0;nophase<nop;nophase++)
				
				*((zone->Flowsol)->Momentumz+index)=*((zone->Flowsol)->Momentumz+index)+(*((zone->Flowsol)->Momentumofphasez[nophase]+index));
			
			}
			
		}
				
	}
		printf("\n coming out from initializevar mPFV FUN\n\n");
return(0);
}	


