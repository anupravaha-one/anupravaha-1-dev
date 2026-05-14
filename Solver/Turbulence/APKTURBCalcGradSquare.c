#include "APKTURBCalcGradSquare.h"
/**This Function is used to Calculate the sum of squares of the gradients of a Variable in X, Y, Z directions.
Application:  Used in kturb() function [LINE: 123, 142] in the case of Turbulent Kinetic Energy (K) and its Dissipation rate (Epsilon) in logarithmic formulation.
*/
double CalcGradSquare(Base *base,Zone *zone,int index[3],int variableTag)
{
	double source = 0;
	double *VariableSlab[6];
	double *variable;
	double DK[3],DEps[3];
	
	switch(variableTag)
	{
		case 600:
		{
			/*variable = zone->Flowsol->Scal[0];
			MakeSlabPtrArray(zone,VariableSlab,600);*/
			
			//variable = MakeDomainArray(zone, 6000);
			variable = zone->PrevFlowsol->Scal[0];
			MakeSlabPtrArray(zone,VariableSlab,6000);
			
			CalcGradVar(base,zone,variable,VariableSlab,index,DK);
			source = DK[0]*DK[0] + DK[1]*DK[1] + DK[2]*DK[2];
		}
		break;
		
		case 601:
		{
			//variable = zone->Flowsol->Scal[1];
			//MakeSlabPtrArray(zone,VariableSlab,601);
			
			//variable = MakeDomainArray(zone, 6010);
			variable = zone->PrevFlowsol->Scal[1];
			MakeSlabPtrArray(zone,VariableSlab,6010);
			
			CalcGradVar(base,zone,variable,VariableSlab,index,DEps);
			source = DEps[0]*DEps[0] + DEps[1]*DEps[1] + DEps[2]*DEps[2];
		}
		break;
		
		default:
			printf("\nImproper Variable Tag in CalcGradSquare.c\n");
		break;
	}
	
	return (source);
}

