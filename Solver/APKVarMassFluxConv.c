#include "APKVarMassFluxConv.h"

/**
This function is exactly similar to "MassFluxConv.c" but is used for a variable massflux case.
This function checks for the massflux conservation in the domain, i.e. total amount coming in and
total amount leaving the domain. This function sums up the massflux at all the faces in all zones
and returns the total difference between the incoming and outgoing fluxes.
*/

double VarMassFluxConv(Base *base,int vartag)
{
	int iface,NI,NJ,NK,*ptrSize,nZone;
	int i,j,k,count,dir,scheme;
	char tag,name[33];
	void *ptr;
	double Temflux,*F,Fdiff[6];
	double *slabVariable[6],*slabVolume[6];
	double *deferCorrFlux,*upwindFlux,*diagConvect,*tempVar;
	double *vols,*beta,*variable,diffdia,*Fc,*Fc_n,alpha=0.5,Kcond;
	int index[3],INDEX;
	Zone *zone;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	
	Temflux=0.0;
	count=base->Number_of_zones;
	for(nZone=0;nZone<base->Number_of_zones;nZone++)
	{
		zone=(base->CgnsZone+nZone);
		ptrSize=zone->Size;
		F = zone->Flowsol->Massflux; 
		scheme = base->SolverScheme;
		
		NI = *(ptrSize+0);
		NJ = *(ptrSize+1);
		NK = *(ptrSize+2);
		
		deferCorrFlux = (double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
		diagConvect = (double *)malloc((NK-1)*(NJ-1)*(NI-1)*sizeof(double));
		upwindFlux = (double *)malloc(6*(NK-1)*(NJ-1)*(NI-1)*sizeof(double));
		
		MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of u */
		MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
		variable=MakeDomainArray(zone,vartag); //variableDomain stores variable pointer
		tempVar=MakeDomainArray(zone,vartag*10); //variableDomain stores variable pointer
		
		
		for(iface=0;iface<6;iface++)
		{
			ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
			if(tag=='i')
			{
				facePtr=(Interface *)ptr;
			}
			else if(tag=='b')
			{
				bcPtr=(ZoneBoundary *)ptr;
				dir=bcPtr->Direction;
				vols= zone->ZGeom->Volume;
				beta= zone->ZGeom->Beta;
				switch(dir)
				{
					case 1:
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							index[0]=k,index[1]=j,index[2]=NI-2;
							convfluxC2D(tempVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							i=NI-2;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							Fc= (upwindFlux+INDEX*6);
							Fc_n= (deferCorrFlux+INDEX*6);
							Temflux=Temflux +Kcond*Fdiff[1]+*(Fc+1)+*(Fc_n+1) ;
						}
					break;
					case 2:
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							index[0]=k,index[1]=j,index[2]=0;
							convfluxC2D(tempVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							i=0;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							Fc= (upwindFlux+INDEX*6);
							Fc_n= (deferCorrFlux+INDEX*6);
							Temflux=Temflux +Kcond*Fdiff[0]+*(Fc+0)+*(Fc_n+0)  ;
						}
					break;
					case 3:
					for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=k,index[1]=NJ-2,index[2]=i;
							convfluxC2D(tempVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							j=NJ-2;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							Fc= (upwindFlux+INDEX*6);
							Fc_n= (deferCorrFlux+INDEX*6);
							Temflux=Temflux +Kcond*Fdiff[3]+*(Fc+3)+*(Fc_n+3) ;
						}
					break;
					case 4:
					for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=k,index[1]=0,index[2]=i;
							convfluxC2D(tempVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							j=0;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							Fc= (upwindFlux+INDEX*6);
							Fc_n= (deferCorrFlux+INDEX*6);
							Temflux=Temflux +Kcond*Fdiff[2]+*(Fc+2)+*(Fc_n+2)  ;
						}
					break;
					case 5:
					for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=NK-2,index[1]=j,index[2]=i;
							convfluxC2D(tempVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							k=NK-2;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							Fc= (upwindFlux+INDEX*6);
							Fc_n= (deferCorrFlux+INDEX*6);
							Temflux=Temflux +Kcond*Fdiff[5]+*(Fc+5)+*(Fc_n+5) ;
						}
					break;
					case 6:
					for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=0,index[1]=j,index[2]=i;
							convfluxC2D(tempVar,variable,alpha,zone,deferCorrFlux,upwindFlux,diagConvect,vartag,scheme,index);
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							k=0;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							Fc= (upwindFlux+INDEX*6);
							Fc_n= (deferCorrFlux+INDEX*6);
							Temflux=Temflux +Kcond*Fdiff[4]+*(Fc+4)+*(Fc_n+4);
						}
					break;
				}
			}
		}
		free(deferCorrFlux);
		free(upwindFlux);
		free(diagConvect);
	}
	return (Temflux);
}
