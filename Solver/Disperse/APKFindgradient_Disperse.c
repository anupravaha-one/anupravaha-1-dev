#include "APKFindgradient_Disperse.h"
/**
The function Findgradient_Disperse() calls the functions  CalcGradVar()and CalcGradVar_Volfrac and to compute the gradients
of volume fraction of disperse phase and that of velocity components multiplied by density of disperse phase and Diffusion coefficient of disperse phase respectively in case of Turbulent gas particle flow, which arise as source terms in the momentum equations of particles due to turbulent dispersion force term.
 */
double Findgradient_Disperse(Base *base,Zone *zone,int index[3],int variableTag)
{
	double source;
	double *gradVariable1Slab[6],*gradVariable2Slab[6],*gradVariable3Slab[6];
	double *variable1,*variable2,*variable3;
	double DV[3],DW[3],D_Chi[3];
	
	if(base->flagMPFV==0)
	{
		if(base->flagTurb == 1)
		{
			if((variableTag==62)||(variableTag==63)||(variableTag==64))//disperse phase momentum
			{
				variable3 = zone->PrevFlowsol->VarthetaDisperse;//
				MakeSlabPtrArray(zone,gradVariable3Slab,650);
			}
		}
	}
	
	CalcGradVar(base,zone,variable3,gradVariable3Slab,index,D_Chi);
	
	if(base->flagMPFV==0)
	{
		if(base->BaseParameters->flag2d==0)
		{
			switch(variableTag)
			{
				case 62:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersex;
					MakeSlabPtrArray(zone,gradVariable1Slab,620);
					CalcGradVar_Volfrac(base,zone,variable1,gradVariable1Slab,index,DV,variableTag);
					source=D_Chi[0]*DV[0]+D_Chi[1]*DV[1]+D_Chi[2]*DV[2];
				}
				case 63:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersey;
					MakeSlabPtrArray(zone,gradVariable1Slab,630);
					CalcGradVar_Volfrac(base,zone,variable1,gradVariable1Slab,index,DV,variableTag);
					source=D_Chi[0]*DV[0]+D_Chi[1]*DV[1]+D_Chi[2]*DV[2];
				}
				case 64:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersex;
					MakeSlabPtrArray(zone,gradVariable1Slab,640);
					CalcGradVar_Volfrac(base,zone,variable1,gradVariable1Slab,index,DV,variableTag);
					source=D_Chi[0]*DV[0]+D_Chi[1]*DV[1]+D_Chi[2]*DV[2];
				}
			}
		}
		else
		{
			switch(variableTag)
			{
				
				case 62:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersex;
					
					MakeSlabPtrArray(zone,gradVariable1Slab,620);
					
					CalcGradVar_Volfrac(base,zone,variable1,gradVariable1Slab,index,DV,variableTag);
	
					source=D_Chi[0]*DV[0]+D_Chi[1]*DV[1];
				}
				break;
				
				case 63:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersey;
					
					MakeSlabPtrArray(zone,gradVariable1Slab,630);
					
					CalcGradVar_Volfrac(base,zone,variable1,gradVariable1Slab,index,DV,variableTag);
	
					source=D_Chi[0]*DV[0]+D_Chi[1]*DV[1];
				}
				break;
			}
		}
	}	
	
	
	return (source);
}


