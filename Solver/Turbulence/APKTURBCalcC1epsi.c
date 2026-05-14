#include "APKTURBCalcC1epsi.h"

/**
This function calculates coefficient C1 for epsilon equation for both standard k-epsilon and RNG k-epsilon. It loops over all the 
cells and calculates coefficient C1 using cell values of k, epsilon and turbulent production.
 */

int CalcC1epsi(Base *base)
{
	Zone *zone;
	int i,j,k,NI,NJ,NK,nZones,numZones,cellindex,*ptrSize,flagrho,flagViscoeff;
	double *Tproduction,*domainptr,*domainC0,*TKE,*EPS,*volume;
	double etta,a,b,viscosity,*rho,f2,f1,fmue,Rt,Ry,y;
	
	nZones = base->Number_of_zones;
	//printf("i m in new function called...\n");
	for(numZones=0; numZones<nZones; numZones++)
	{
		zone = base->CgnsZone+numZones;
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			
			Tproduction = zone->Flowsol->Tcoeff->TProduction;
			domainC0=zone->Flowsol->Tcoeff->C0epsi;
			
			TKE = zone->Flowsol->Scal[0];
			EPS = zone->Flowsol->Scal[1];
			//printf("zone->Flowsol->flagRNG_KEPS == %d & zone->Flowsol->flagTurb ==%d\n",zone->Flowsol->flagRNG_KEPS,zone->Flowsol->flagTurb);
			if((zone->Flowsol->flagRNG_KEPS == 1) && (zone->Flowsol->flagTurb == 1)) // check for RNG k-epsilon
			{
				if(base->flagParticle==1)
					domainptr = zone->Flowsol->Tcoeff->C2epsi;
				else	
					domainptr = zone->Flowsol->Tcoeff->C1epsi;
				
				rho = zone->Flowsol->Rho;
				flagrho = zone->Flowsol->flagRho;
				volume = zone->ZGeom->Volume;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					etta = (*(TKE+cellindex) / *(EPS+cellindex))*sqrt(*(Tproduction+cellindex)* (*(rho+cellindex*flagrho)/ *(volume+cellindex)));
					if(base->flagParticle==1)
						a = (1.0 - etta/4.38)*pow(etta,3)*0.09;
					else	
						a = (1.0 - etta/4.38)*etta;
					
					b = 1.0 + 0.015*pow(etta,3);
					if(base->flagParticle==0)
					{
					*(domainptr + cellindex) = 1.42 - a/b;
					*(domainC0 + cellindex) = 0.0;
					if(*(domainptr + cellindex) > 1.42)
						*(domainptr + cellindex) = 1.42;
					if(*(domainptr + cellindex) < 0.0)
						*(domainptr + cellindex) = 0.0;
					}
					else if(base->flagParticle==1)
					{
						*(domainptr + cellindex)=1.68+a/b;
						if(*(domainptr + cellindex) > 1.68)
							*(domainptr + cellindex) = 1.68;
						if(*(domainptr + cellindex) < 0.0)
							*(domainptr + cellindex) = 0.0;
					}
				/*	a = (1.0 - etta/4.38)*pow(etta,3);
					b = 1.0 + 0.015*pow(etta,3);
				//	*(domainptr + cellindex) = 1.42;
					*(domainC0 + cellindex) = a/b;
					if(*(domainC0 + cellindex) > 1.42) // upper limit for coefficient
					*(domainC0 + cellindex) = 1.42;
					if(*(domainC0 + cellindex) < 0.0) // lower limit for coefficient
				*(domainC0 + cellindex) = 0.0; */
					//printf("C1=%lf\ti=%d\tj=%d\tk=%d\n",*(domainptr + cellindex),i,j,k);
				}
			}
			else if((zone->Flowsol->flagRNG_KEPS == 0) && (zone->Flowsol->flagTurb==1)) // check for standard k-epsilon
			{
				domainptr = zone->Flowsol->Tcoeff->C1epsi;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*(domainptr + cellindex) = 1.44;
					*(domainC0 + cellindex) = 0.0;
				}
			}	
			else if((zone->Flowsol->flagRNG_KEPS == 0) && (zone->Flowsol->flagTurb==4)) // Launder and Sharma
			{
				domainptr = zone->Flowsol->Tcoeff->C2epsi;
				rho = zone->Flowsol->Rho;
				flagrho = zone->Flowsol->flagRho;
				flagViscoeff = zone->Flowsol->flagViscoeff;
			
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*(domainC0 + cellindex) = 0.0;
					viscosity = *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff);
				
					if(*(EPS+cellindex)!=0.0)
					{
						f2 = 1.0 - 0.3 * exp(-(pow((*(rho+flagrho*cellindex) * pow(*(TKE+cellindex),2.0)/(*(EPS+cellindex) * viscosity)),2.0)));
						*(domainptr+cellindex) = 1.92 * f2;
					}
					else if((*(TKE+cellindex)==0.0) &&(*(EPS+cellindex)==0.0))
					{
						f2=0.7;
						*(domainptr+cellindex) = 0.7 * 1.92;
					}
					else
						printf("f2 is not defined\n");
				
				}
			}	
			else if((zone->Flowsol->flagRNG_KEPS == 0) && (zone->Flowsol->flagTurb==5)) // check for Lam Bremhorst
			{
				domainptr = zone->Flowsol->Tcoeff->C2epsi;
				rho = zone->Flowsol->Rho;
				flagrho = zone->Flowsol->flagRho;
				flagViscoeff = zone->Flowsol->flagViscoeff;
			
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*(domainC0 + cellindex) = 0.0;
					viscosity = *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff);
					Rt = *(rho+flagrho*cellindex) * pow(*(TKE+cellindex),2.0)/(*(EPS+cellindex) * viscosity);
					y = *(zone->ZGeom->minWallDistance+cellindex);
					Ry= *(rho+flagrho*cellindex) * y * sqrt(*(TKE+cellindex))/viscosity;
					if(Ry>1.0)
						fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
					else 
						fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
					if((Rt<1.0)&&(fmue>1.0))                                                                                             			fmue=0.01116225;
					fmue=0.01116225;
					
					f2 = 1.0 - exp(-(pow(Rt,2.0)));
					if(f2 < 0.01)  // limiter for f2 in lam bremhorst model
						f2 = 0.01;
					*(domainptr+cellindex) = 1.92 * f2;
					
						
					f1 = 1.0 + pow(0.05/fmue,3.0);
					*(zone->Flowsol->Tcoeff->C1epsi+cellindex)=1.44*f1;
					
				}
			}	
		}
	}
	return 0;
}
