#include "APKLMMHD_Source_Calculation.h"

/**
This function is used for calculation of pressure gradient of each cell at once and storing the value in an array.
This value of Pressure gradient will be used for subsequent calculation of Navier-Stokes equations.
*/


int LMMHD_Source_Calculation(Base *base)
{
	int NI,NJ,NK,index[3],*ptrSize,indexP;
	int k,j,i,numZones,nZones,flagrho,flagelectcond;
	double *variable[3],rho,sigma,Dphi[3];
	double B_x, B_y, B_z, velPu, velPv, velPw,*volume;
	nZones = base->Number_of_zones;
	Zone *zone;
	
	
	//printf("Calculating LMMHD Momentum Source values in N-S Equations.....\n\n");
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		volume= zone->ZGeom->Volume;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		variable[0] = zone->Flowsol->LMMHD_Source[0];
		variable[1] = zone->Flowsol->LMMHD_Source[1];
		
		if(base->BaseParameters->flag2d==0)
			variable[2] = zone->Flowsol->LMMHD_Source[2];
		
		indexP=0;
		flagrho = zone->Flowsol->flagRho;
		
		
		flagelectcond=zone->Flowsol->flageleccond;//rajendra

		
				
		
		if(base->BaseParameters->flag2d==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				sigma = *(zone->PrevFlowsol->Eleconductivity+indexP*flagelectcond);
				index[0]=k;index[1]=j;index[2]=i;
				
				FindgradientPhiMHD(base,zone,index,Dphi);
				
				velPu = *(zone->Flowsol->u+indexP);
				velPv = *(zone->Flowsol->v+indexP);
				//velPw = *(zone->Flowsol->w+indexP);
				velPw=0.0;
				
				
					
				if(base->flagLMMHD==1)
				{
				
					{
						B_x = *(zone->Flowsol->MfieldX+indexP);
						B_y = *(zone->Flowsol->MfieldY+indexP);
						B_z = *(zone->Flowsol->MfieldZ+indexP);
					}
				}
					


				*(variable[0]+indexP)=(*(volume+indexP))*(sigma/rho)*((B_y*Dphi[2]/(*(volume+indexP)) - B_z*Dphi[1]/(*(volume+indexP)) + B_x*(velPv*B_y+velPw*B_z)));

				*(variable[1]+indexP)=(*(volume+indexP))*(sigma/rho)*((B_z*Dphi[0]/(*(volume+indexP)) - B_x*Dphi[2]/(*(volume+indexP)) + B_y*(velPu*B_x+velPw*B_z)));
				
				
 				
				indexP++;
				
			}
		}
		else if(base->BaseParameters->flag2d==0)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				sigma = *(zone->PrevFlowsol->Eleconductivity+indexP*flagelectcond);
				index[0]=k;index[1]=j;index[2]=i;
 				

				FindgradientPhiMHD(base,zone,index,Dphi);

				velPu = *(zone->Flowsol->u+indexP);
				velPv = *(zone->Flowsol->v+indexP);
				velPw = *(zone->Flowsol->w+indexP);

				
				if(base->flagLMMHD==1)
				{
				
				
					{
						B_x = *(zone->Flowsol->MfieldX+indexP);
						B_y = *(zone->Flowsol->MfieldY+indexP);
						B_z = *(zone->Flowsol->MfieldZ+indexP);
					}
				}

 				

				*(variable[0]+indexP)=(*(volume+indexP))*(sigma/rho)*((B_y*Dphi[2]/(*(volume+indexP)) - B_z*Dphi[1]/(*(volume+indexP)) + B_x*(velPv*B_y+velPw*B_z)));

				*(variable[1]+indexP)=(*(volume+indexP))*(sigma/rho)*((B_z*Dphi[0]/(*(volume+indexP)) - B_x*Dphi[2]/(*(volume+indexP)) + B_y*(velPu*B_x+velPw*B_z)));


				*(variable[2]+indexP)=(*(volume+indexP))*(sigma/rho)*((B_x*Dphi[1]/(*(volume+indexP)) - B_y*Dphi[0]/(*(volume+indexP)) + B_z*(velPu*B_x+velPv*B_y)));

				indexP++;
			}
		}
	}
	return(0);
}


