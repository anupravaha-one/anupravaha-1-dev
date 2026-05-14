#include"APKSpalartAllmarasSource.h"

double SpalartAllmarasSource(Base *base,Zone *zone,double *source,int index[3])
{
	int flagrho,flagViscoeff,i,j,k,NI,NJ,NK,indexP;
	double diacoeff,Production,Destruction,CrossDiffusion;
	double Cb1=0.1355,Cb2=0.622,Cw1,Cw2=0.3,kappa=0.41,sigma=2.0/3.0;
	//double Cv1=7.1,Cw3=2.0;
	double fv1,fv2,fw,X,d,r,g,temp;
	double vol,rho,nu_t,nu_tnew,nu,mu,S,Stilda;
	double *NuOld,*slabNuOld[6],Dnu[3];
	
	Cw1 = Cb1/(kappa*kappa) + (1.0+Cb2)/sigma;
	
	//Put zone loop here
	flagrho = zone->Flowsol->flagRho;
	flagViscoeff = zone->Flowsol->flagViscoeff;
	
	NuOld	= zone->PrevFlowsol->Scal[0];
	MakeSlabPtrArray(zone,slabNuOld,6000);
	
	NI=*(zone->Size+0);
	NJ=*(zone->Size+1);
	NK=*(zone->Size+2);
	k=index[0],j=index[1],i=index[2];
	indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
	if(flagrho == 1)
	{
		rho = *(zone->PrevFlowsol->Rho+indexP);
	}
	else
	{
		rho = *(zone->PrevFlowsol->Rho);
	}
	
	if(flagViscoeff == 1)
	{
		mu = *(zone->Flowsol->Viscoeff+indexP);
	}
	else
	{
		mu = *(zone->Flowsol->Viscoeff);
	}
	
	nu  = mu/rho;
	nu_t = *(NuOld+indexP);
	nu_tnew = *(zone->Flowsol->Scal[0]+indexP);
	X = nu_t/nu;
	
	temp = X*X*X;
	fv1 = temp/(temp+357.911);	//fv1 = temp/(temp+Cv1*Cv1*Cv1);
	fv2 = 1.0 - X/(1.0+X*fv1);
	
	S = *(zone->Flowsol->Vorticity+indexP);
	d  = *(zone->ZGeom->minWallDistance+indexP);
	
	temp = kappa*d;
	temp = temp*temp;
	Stilda = S + nu_t/temp*fv2;	//Stilda = S + nu_t/(kappa*kappa*d*d)*fv2;
	
	//Calculating the function `fw' used in Destruction term
	r = nu_t/(Stilda*temp);		//r = nu_t/(Stilda*kappa*kappa*d*d);
	temp = r*r*r;
	temp = temp*temp;
	g = r + Cw2*(temp-r);		//g = r + Cw2*(r*r*r*r*r*r-r);
	temp = g*g*g;
	temp = temp*temp;
	fw = g * pow(65.0/(temp+64.0), 1.0/6.0);	//fw = g * pow((1.0+pow(Cw3,6))/(pow(g,6)+pow(Cw3,6)),1.0/6.0);
	
	//Calculating the Production term
	Production = Cb1*Stilda*nu_t;
	
	//Calculating the Cross Diffusion term
	CalcGradVar(base,zone,NuOld,slabNuOld,index,Dnu);
	temp = Dnu[0]*Dnu[0] + Dnu[1]*Dnu[1] + Dnu[2]*Dnu[2];
	CrossDiffusion = Cb2/sigma * temp;
	
	vol = *(zone->ZGeom->Volume+indexP);
	
	//Calculating the Destruction term
	temp = Cw1*fw/(d*d)*nu_t;
	Destruction = temp * (-nu_t);
	//Destruction = temp * (2.0*nu_tnew - nu_t);
	
	diacoeff = vol*temp*2.0;
	*source = vol * (Production - Destruction + CrossDiffusion);
	//*source = *source + vol * (Production - Destruction + CrossDiffusion);
	
	return(diacoeff);
}
