#include "APKCalcvarCoeffAndPrevCoeffConjugate.h"
/**
The function CalcvarCoefficientConjugate() is similar to CalcvarCoefficient() and is used only incase of conjugate heat transfer.
This function returns the current value of "coefficient" corresponding to a particular variable that
arises in the convective terms. For example, in the case of temperature we have Cp, specific heat arising in the convective
terms.
*/
double CalcvarCoefficientConjugate(Base *base,Zone *zone,int variableTag,int INDEX)
{
	double varCoefficient;
	int flagrho,flagheat;
	
	if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==7))
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= *(zone->Flowsol->Rho+(INDEX*flagrho));
	}
	else if(variableTag==6)
	{
		flagrho=(zone->Flowsol)->flagRho;
		flagheat=(zone->Flowsol)->flagSpecificHeat;
		varCoefficient= *(zone->Flowsol->SpecificHeat+(INDEX*flagheat)) * *(zone->Flowsol->Rho+(INDEX*flagrho));/* rho*Cp */
	}
	else if((variableTag>= 600) && (variableTag<=609))
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= *(zone->Flowsol->Rho+(INDEX*flagrho));/*This is correct orelse PLEASE CHANGE THIS ACCORDINGLY*/
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
The function CalcvarCoefficientConjugatePrev() is similar to CalcvarCoefficientPrev() and is used only incase of conjugate heat transfer.
This function returns the previous value of "coefficient" corresponding to a particular variable that
arises in the convective terms. For example, in the case of temperature we have Cp, specific heat arising in the convective
terms.
*/
double CalcvarCoefficientPrevConjugate(Base *base,Zone *zone,int variableTag,int INDEX)
{
	double varCoefficient;
	int flagrho,flagheat;
	
	if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==7))
	{
		flagrho=(zone->Flowsol)->flagRho;
		varCoefficient= *(zone->PrevFlowsol->Rho+(INDEX*flagrho));
	}
	else if(variableTag==6)
	{
		flagrho=(zone->Flowsol)->flagRho;
		flagheat=(zone->Flowsol)->flagSpecificHeat;
		//varCoefficient= *(zone->PrevFlowsol->SpecificHeat+(INDEX*flagheat));/* 1.0*Cp */
		varCoefficient= *(zone->PrevFlowsol->Rho+(INDEX*flagrho))*(*(zone->PrevFlowsol->SpecificHeat+(INDEX*flagheat)));/* Rho*Cp */
	}
	else if((variableTag>= 600) && (variableTag<=609))
	{
		flagrho=(zone->Flowsol)->flagRho;
		//varCoefficient= 1.0;/*PLEASE CHANGE THIS ACCORDINGLY*/
		varCoefficient= *(zone->PrevFlowsol->Rho+(INDEX*flagrho));/*PLEASE CHANGE THIS ACCORDINGLY*/
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

	

