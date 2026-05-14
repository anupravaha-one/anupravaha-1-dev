#include "APKModifiedDiffCoeff.h"

double ModifiedDiffCoeff(Zone *z,Interface *facePtr,double diffCoef,int indArray[3],int domainCellIndex)
{
	double deltaPlus,deltaP,uTau,Kmod;
	double Pr,Prt,a1,a2,beta,Tplus,y1plus,y2plus;
	int NI,NJ,NK,i,j,k,flagrho,flagviscoeff,flagSpecificHeat,flagDiffcoeff,slabCellIndex;
	double rho,mu,Cp,K,Kt,Chi;

	Prt = 1.0;
	Kt = diffCoef;
	Chi = z->Flowsol->Tcoeff->Chi;
	
	NI = *(z->Size+0);
	NJ = *(z->Size+1);
	NK = *(z->Size+2);
	
	i = indArray[2];
	j = indArray[1];
	k = indArray[0];
	
	switch(facePtr->Direction)
	{
		case 1:
		case 2:
			slabCellIndex = (NJ-1)*k+j;
			break;
		case 3:
		case 4:
			slabCellIndex = (NI-1)*k+i;
			break;
		case 5:
		case 6:
			slabCellIndex = (NI-1)*j+i;
			break;
	}
	
	flagrho = z->Flowsol->flagRho;
	flagviscoeff = z->Flowsol->flagViscoeff;
	flagSpecificHeat = z->Flowsol->flagSpecificHeat;
	flagDiffcoeff = z->Flowsol->flagDiffcoeff;

	uTau = *(facePtr->uTau+slabCellIndex);
	deltaP = *(facePtr->deltaP+slabCellIndex);
	
	rho = *(z->Flowsol->Rho+flagrho*domainCellIndex);
	mu = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
	Cp = *(z->Flowsol->SpecificHeat+flagSpecificHeat*domainCellIndex);
	K = *(z->Flowsol->Diffcoeff+flagDiffcoeff*domainCellIndex);
	
	deltaPlus = rho * deltaP * uTau / mu;
	if(deltaPlus < 11.067)	//<---- Scallable Wall Function Treatment; Ref: Heat Transfer Predictions using advanced turbulence models, [CFX], Mentor
	{
		deltaPlus =  11.067;
		deltaP = deltaPlus * mu/(rho * uTau);
	}
	
	Pr = mu*Cp/K;
	y1plus = 10.0/pow(Pr,1.0/3.0);
	a1 = 0.001 * Prt;
	a2 = 15.0 * pow(Pr,(2.0/3.0));
	y2plus = sqrt(Chi/a1);
	beta = a2 - Prt/(2.0*Chi) * (1.0 + log(Chi/a1));
	
	if(deltaPlus < y1plus)
	{
		Tplus = Pr * deltaPlus;
	}
	else 
	{
		if(deltaPlus < y2plus)
			Tplus = a2 - Prt/(2.0*a1*deltaPlus*deltaPlus);
		else
			Tplus = (Prt/Chi) * log(deltaPlus) + beta;
	}
	
	Kmod = Kt * Prt/(Chi*Tplus);
	
	return(Kmod);
}
