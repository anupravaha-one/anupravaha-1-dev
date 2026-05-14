#include "APKTURBCalcEpsilonWallTerm.h"

/**
CalcTurbulentProduction function calculates the production term which comes in the K and Epsilon equations. (See Report for further details
regarding Production term). This function calls TurbProduction() function which loops over the cells and according to
2D or 3D problem, calls grad2D() or grad3D() which actually calculates the production term.
*/
int CalcEpsWallTerm(Base *base)
{
	int nZones;
	Zone *zone;
	int numZones;
	
	nZones=(base->Number_of_zones);
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		
		if(zone->flagSolidFluid==0)
		{
			EpsilonWall(base,zone,zone->Flowsol->Tcoeff->PhiEps);
		}
	}
	
	return 0;
}

int EpsilonWall(Base *base,Zone *zone,double *GEps)
{
	int i,j,k,NI,NJ,NK,index[3],INDEX,face;
	int flagviscosity,flagrho,*ptrsize;
	double rho,*vols,*ptrVolume,*beta,*variable;
	double *slabVariable[6],*slabVolume[6],Fdiff[6];
	double Fsum,Viscosity,mue_t;
	
	ptrsize = zone->Size;
	ptrVolume= zone->ZGeom->Volume;
	flagrho = zone->Flowsol->flagRho; 
	flagviscosity = zone->Flowsol->flagViscoeff;
	NI=*(ptrsize+0);
	NJ=*(ptrsize+1);
	NK=*(ptrsize+2);
	variable=MakeDomainArray(zone,1);
	MakeSlabPtrArray(zone,slabVariable,1);/*Called for making slabs of u*/
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	beta= zone->ZGeom->Beta;
	vols = zone->ZGeom->Volume;
	
	if(base->BaseParameters->flag2d==1)
	{
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*INDEX);
			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);
			
			DoubleDerivative2d(index,variable,vols,ptrsize,slabVariable,slabVolume,beta,Fdiff);
			Fsum=0.0;
			
			for(face=0;face<4;face++)
				Fsum=Fsum+Fdiff[face];
			*(GEps+INDEX)=2.0*Viscosity/rho*mue_t/rho*Fsum*Fsum/(*(vols+INDEX));
		}
		variable=MakeDomainArray(zone,2);
		MakeSlabPtrArray(zone,slabVariable,2);/*Called for making slabs of u*/
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*INDEX);
			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);
			
			DoubleDerivative2d(index,variable,vols,ptrsize,slabVariable,slabVolume,beta,Fdiff);
			Fsum=0.0;
			
			for(face=0;face<4;face++)
				Fsum = Fsum+Fdiff[face];
			*(GEps+INDEX) = *(GEps+INDEX) + 2.0*Viscosity/rho*mue_t/rho*Fsum*Fsum/(*(vols+INDEX));
		}
	}
	else
	{
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*INDEX);
			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);
			
			DoubleDerivative3d(index,variable,vols,ptrsize,slabVariable,slabVolume,beta,Fdiff);
			Fsum=0.0;
			
			for(face=0;face<6;face++)
				Fsum=Fsum+Fdiff[face];
			*(GEps+INDEX)=2.0*Viscosity/rho*mue_t/rho*Fsum*Fsum/(*(vols+INDEX));
		}
		variable=MakeDomainArray(zone,2);
		MakeSlabPtrArray(zone,slabVariable,2);/*Called for making slabs of u*/
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*INDEX);
			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);
			
			DoubleDerivative3d(index,variable,vols,ptrsize,slabVariable,slabVolume,beta,Fdiff);
			Fsum=0.0;
			
			for(face=0;face<6;face++)
				Fsum = Fsum+Fdiff[face];
			*(GEps+INDEX) = *(GEps+INDEX) + 2.0*Viscosity/rho*mue_t/rho*Fsum*Fsum/(*(vols+INDEX));
		}
		variable=MakeDomainArray(zone,3);
		MakeSlabPtrArray(zone,slabVariable,3);/*Called for making slabs of u*/
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*INDEX);
			rho = *(zone->PrevFlowsol->Rho+INDEX*flagrho);
			mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+INDEX);
			
			DoubleDerivative3d(index,variable,vols,ptrsize,slabVariable,slabVolume,beta,Fdiff);
			Fsum=0.0;
			
			for(face=0;face<6;face++)
				Fsum = Fsum+Fdiff[face];
			*(GEps+INDEX) = *(GEps+INDEX) + 2.0*Viscosity/rho*mue_t/rho*Fsum*Fsum/(*(vols+INDEX));
		}
	}
	
return 0;
}
