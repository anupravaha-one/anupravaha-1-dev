#include "APKTURBCalcVariableOneIterOmega.h"
/**
The function CalcVariableOneIterOmega() is supposed to loop over all the zones and calls some other function which does the Gauss Seidel
iterations of the Omega equation.  This function will be written when k-Omega model is implemented.
*/
double CalcVariableOneIterOmega(Base *base,int variableTag)
{
	double baseResidual = 0.0;
	int totalBaseCells;
	double zoneResidual;
	int *ptrSize,nZones;
	Zone *zone;
	int numZones;
	double *tempVar,*prevTstepVar,*variableDomain;
	
	nZones=(base->Number_of_zones);
	
	totalBaseCells = 0;

	// baseResidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones*/
	RefreshSurfaceCells(variableTag,base);/*supply tag for Variable*/
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		
		variableDomain=MakeDomainArray(zone,variableTag); //variableDomain stores variable pointer
		prevTstepVar=MakeDomainArray(zone,variableTag*10);   //prevTstepVar stores prevFlowsol->variable 
		
		//printf("variable6=%lf\tprevtstep6=%lf\n",*(variableDomain+1234),*(prevTstepVar+1232));
		
		tempVar=(double*)malloc((*(ptrSize+3)) * (*(ptrSize+4)) * (*(ptrSize+5))*sizeof(double));
		
		totalBaseCells = totalBaseCells+(*(ptrSize+3)) * (*(ptrSize+4)) * (*(ptrSize+5));/*NcellI*NcellJ*NcellK*/

		ExtractVariables(variableTag,base,zone);
		
		fillnewtoprev(variableDomain,tempVar,ptrSize);/*tempVar stores variable of previous iteration (not time step)*/
		//printf("I am going to Omega");
		//zoneResidual = UpdateVar(base,zone,variableDomain,prevTstepVar,variableTag,tempVar);
		zoneResidual = UpdateVarCD2d(base,zone,variableTag);
		//ExtractVariables(variableTag,base,zone);
		baseResidual = baseResidual+zoneResidual;
		free(tempVar);
	}

	baseResidual = sqrt(baseResidual/totalBaseCells);
	//printf("Residue variable Epsilon = %12.8e tag=%d\n",baseResidual,variableTag);
	//exit(0);
	return (baseResidual);
}
