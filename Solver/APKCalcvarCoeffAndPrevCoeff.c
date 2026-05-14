#include "APKCalcvarCoeffAndPrevCoeff.h"
/**
The function CalcvarCoefficient() returns the current value of "coefficient" corresponding to a particular variable that
arises in the convective terms. For example, in the case of temperature we have Cp, specific heat arising in the convective
terms.
*/
double CalcvarCoefficient(Base *base,Zone *zone,int variableTag,int INDEX)
{
	double varCoefficient;
	int flagrho,flagheat;
	
	if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==7)||(variableTag==56)||(variableTag==57)||(variableTag==58)||(base->flagMPFV==1)||variableTag==62||variableTag==63||variableTag==64||variableTag==65||variableTag==710||(variableTag==68))//rajendra//added disperse phase....adhi+shash)//rajendra //||variableTag==67
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= 1.0;
	}
	else if(variableTag==69)
	{
		flagrho=(zone->Flowsol)->flagRho;
		flagheat=(zone->Flowsol)->flagCPdrop;
		varCoefficient= *(zone->Flowsol->CPdrop+(INDEX*flagheat));/* 1.0*Cp */
	}
	else if(variableTag==6)
	{
		flagrho=(zone->Flowsol)->flagRho;
		flagheat=(zone->Flowsol)->flagSpecificHeat;
		varCoefficient= *(zone->Flowsol->SpecificHeat+(INDEX*flagheat));/* 1.0*Cp */
	}
	else if((variableTag>= 600) && (variableTag<=609))
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= 1.0;/*This is correct orelse PLEASE CHANGE THIS ACCORDINGLY*/
	}
        else if((variableTag==7)||(variableTag==14)||(variableTag==15)) //Shyam
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= 1.0;
	}
	else if(variableTag==8)
	{
		varCoefficient= 1.0;
	}
	else if(variableTag==9)
	{
		varCoefficient= 1.0;
	}
	return(varCoefficient);
}

/**
The function CalcvarCoefficientPrev() returns the previous value of "coefficient" corresponding to a particular variable that
arises in the convective terms. For example, in the case of temperature we have Cp, specific heat arising in the convective
terms.
*/
double CalcvarCoefficientPrev(Base *base,Zone *zone,int variableTag,int INDEX)
{
	double varCoefficient;
	int flagrho,flagheat;
	
	if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==7)||(variableTag==56)||(variableTag==57)||(variableTag==58)||(base->flagMPFV==1)||variableTag==62||variableTag==63||variableTag==64||variableTag==65||variableTag==710||(variableTag==68))//rajendra//added disperse phase...adhi+shash)//rajendra //||variableTag==67
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= 1.0;//*(zone->PrevFlowsol->Rho+(INDEX*flagrho));
	}
	else if(variableTag==69)
	{
		flagrho=(zone->Flowsol)->flagRho;
		flagheat=(zone->Flowsol)->flagCPdrop;
		varCoefficient= *(zone->PrevFlowsol->CPdrop+(INDEX*flagheat));/* 1.0*Cp */
	}
	else if(variableTag==6)
	{
		flagrho=(zone->Flowsol)->flagRho;
		flagheat=(zone->Flowsol)->flagSpecificHeat;
		varCoefficient= *(zone->PrevFlowsol->SpecificHeat+(INDEX*flagheat));/* 1.0*Cp */
		//varCoefficient= *(zone->PrevFlowsol->Rho+(INDEX*flagrho))*(*(zone->PrevFlowsol->SpecificHeat+(INDEX*flagheat)));/* Rho*Cp */
	}
	else if((variableTag>= 600) && (variableTag<=609))
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= 1.0;/*PLEASE CHANGE THIS ACCORDINGLY*/
		//varCoefficient= *(zone->PrevFlowsol->Rho+(INDEX*flagrho));/*PLEASE CHANGE THIS ACCORDINGLY*/
	}
        else if((variableTag==7)||(variableTag==14)||(variableTag==15)) //Shyam
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= 1.0;
	}
	else if(variableTag==8)
	{
		varCoefficient= 1.0;
	}
	else if(variableTag==9)
	{
		varCoefficient= 1.0;
	}

	return(varCoefficient);
}

	

