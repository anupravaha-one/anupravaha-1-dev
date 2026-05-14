#include "APKinitializeVarDisperse.h"
/**
The Function initializeVarDisperse() derive the boundary conditions for solvable variable that is specific momenta from those of the volume fraction and density of the disperse phase. 
*/

double initializeVarDisperse(Base *CgnsBase)

{	
	int i,j,k,NI,NJ,NK,numZones,cntZ,*ptrSize;
	Zone *zone;
	//double *variable;
	int totalcells,index,flagrhoDisperse=0;
	numZones= CgnsBase->Number_of_zones;
	
	
	//EvaluateVarCoefficients(CgnsBase);
	//RefreshVariableCoefficients(710,CgnsBase);
	
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
			flagrhoDisperse=(CgnsBase->BaseParameters->flagRhoDisperse);
			
			*((zone->Flowsol)->MomentumDispersex+index)=(*((zone->Flowsol)->VarthetaDisperse+index)*(*((zone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->Flowsol)->MomentumDispersex+index)));
			
			*((zone->Flowsol)->MomentumDispersey+index)=(*((zone->Flowsol)->VarthetaDisperse+index)*(*((zone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->Flowsol)->MomentumDispersey+index)));
			
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
					*((zone->Flowsol)->MomentumDispersez+index)=(*((zone->Flowsol)->VarthetaDisperse+index)*(*((zone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->Flowsol)->MomentumDispersez+index)));
			}
			
			//flagrhoDisperse=(CgnsBase->BaseParameters->flagRhoDisperse);
			
			*((zone->PrevFlowsol)->MomentumDispersex+index)=(*((zone->PrevFlowsol)->VarthetaDisperse+index)*(*((zone->PrevFlowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->PrevFlowsol)->MomentumDispersex+index)));
			
			*((zone->PrevFlowsol)->MomentumDispersey+index)=(*((zone->PrevFlowsol)->VarthetaDisperse+index)*(*((zone->PrevFlowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->PrevFlowsol)->MomentumDispersey+index)));
			
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
					*((zone->PrevFlowsol)->MomentumDispersez+index)=(*((zone->PrevFlowsol)->VarthetaDisperse+index)*(*((zone->PrevFlowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->PrevFlowsol)->MomentumDispersez+index)));
			}
		//for momentumDisperseOLd values
		*((zone->Flowsol)->MomentumDispersexOld+index)=(*((zone->Flowsol)->VarthetaDisperse+index)*(*((zone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->Flowsol)->MomentumDispersexOld+index)));
			
			*((zone->Flowsol)->MomentumDisperseyOld+index)=(*((zone->Flowsol)->VarthetaDisperse+index)*(*((zone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->Flowsol)->MomentumDisperseyOld+index)));
			
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
					*((zone->Flowsol)->MomentumDispersezOld+index)=(*((zone->Flowsol)->VarthetaDisperse+index)*(*((zone->Flowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->Flowsol)->MomentumDispersezOld+index)));
			}
			
			//flagrhoDisperse=(CgnsBase->BaseParameters->flagRhoDisperse);
			
			*((zone->PrevFlowsol)->MomentumDispersexOld+index)=(*((zone->PrevFlowsol)->VarthetaDisperse+index)*(*((zone->PrevFlowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->PrevFlowsol)->MomentumDispersexOld+index)));
			
			*((zone->PrevFlowsol)->MomentumDisperseyOld+index)=(*((zone->PrevFlowsol)->VarthetaDisperse+index)*(*((zone->PrevFlowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->PrevFlowsol)->MomentumDisperseyOld+index)));
			
			if(CgnsBase->BaseParameters->flag2d == 0)
			{
					*((zone->PrevFlowsol)->MomentumDispersezOld+index)=(*((zone->PrevFlowsol)->VarthetaDisperse+index)*(*((zone->PrevFlowsol)->RhoDisperse+index*flagrhoDisperse))*(*((zone->PrevFlowsol)->MomentumDispersezOld+index)));
			}	
		//printf("\n %lf\t %lf \t %lf \t %lf \t %lf   \n",*((zone->Flowsol)->MomentumDispersex+index),*((zone->Flowsol)->MomentumDispersey+index),*((zone->PrevFlowsol)->MomentumDispersex+index),*((zone->PrevFlowsol)->MomentumDispersey+index),*((zone->Flowsol)->VarthetaDisperse+index)); 	
		}
	}
			
			//*((zone->Flowsol)->Momentumx+index)=0.0;
			//*((zone->Flowsol)->Momentumy+index)=0.0;
			
			
/************************overall momentum calculation from individual momentum***************************************/ 
/*			
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
			
	}*/
		printf("\n coming out from initializeVarDisperse\n\n");
return(0);
}	


