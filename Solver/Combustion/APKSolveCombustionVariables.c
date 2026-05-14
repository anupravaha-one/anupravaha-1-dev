#include "APKSolveCombustionVariables.h"
/** 
This function Calculates the combustion variables.  After flow is solved for that particular timestep, the combustion variables are computed.  Different codes for Diffusion controlled and kinematicall controlled cases.
Remember that premixed case is not tested yet.
**/



int SolveCombustionVariables(Base *CgnsBase)
{

	int GStoken;
	double variableResidual[4];//Timeerror[4];
	//CalcMuEff(CgnsBase,14);
	
	if(CgnsBase->flagDiffKinetCont==1) //Diffusion controlled reaction
	{
		CallCalcCoefficientsForVariables(CgnsBase,0,14);
		CDcoeffCalculation(CgnsBase,14);
		CalcSource_opti(CgnsBase,14);
		CalcVariable(CgnsBase,14);
	}
        else if(CgnsBase->flagDiffKinetCont==2) //Kinetically Controlled
        {

                if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                {
			
			CallCalcCoefficientsForVariables(CgnsBase,0,7);
			CDcoeffCalculation(CgnsBase,7);
			CalcSource_opti(CgnsBase,7);
			CalcVariable(CgnsBase,7);
			
			CallCalcCoefficientsForVariables(CgnsBase,0,14);
			CDcoeffCalculation(CgnsBase,14);
			CalcSource_opti(CgnsBase,14);
			CalcVariable(CgnsBase,14);
			
			CallCalcCoefficientsForVariables(CgnsBase,0,15);
			CDcoeffCalculation(CgnsBase,15);
			CalcSource_opti(CgnsBase,15);
			CalcVariable(CgnsBase,15);
                }
                else //Premixed case
                {
			do
                	{
				variableResidual[0] = CalcVariableOneIter(CgnsBase,7);
                		variableResidual[1] = CalcVariableOneIter(CgnsBase,14);
                                variableResidual[2] = CalcVariableOneIter(CgnsBase,15);
				variableResidual[3] = CalcVariableOneIter(CgnsBase,16);

                                GStoken = FindGSResidualCombVariables(variableResidual);

                	}while(GStoken == 1);

                }
        }
return(0);
}

int FindGSResidualCombVariables(double *variableResidual)
{
	int token=0,i;
	for(i=0;i<4;i++)
        {
        	if(*(variableResidual+i) > 1.0e-08)
                {
                	token = 1;
                }

        }
	return(token);
}


