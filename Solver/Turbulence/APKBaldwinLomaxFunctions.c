double CalcMixingLength(Base *base)
{
	int NI,NJ,NK,numZones,cellindex;
	int flagLogKEps,flagViscoeff=0;
	int i,j,k,*ptrSize,flagrho,numIface,numZbc,dir;
	int iface,zBc,*ptrRange,index,UpperLimitSlabCount=0,UpperLimitCount=0,LowerLimitSlabCount=0,LowerLimitCount=0;
	Zone *zone;
	double *domainPtr,*slabIptr,*slabBcptr,gammastar;
	double rho,*slabrho,*KE,*Epsilon,*Omega,*slabKE,*slabEpsilon,*slabOmega,fmue,viscosity;
	double mu,nu,y;
	Interface *IfacePtr;
	ZoneBoundary *bcPtr;
	BaldwinLomax *BL;

	for(numZones=0;numZones < base->Number_of_zones;numZones++)
	{
		zone = base->CgnsZone+numZones;

		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			BL = zone->Flowsol->Tcoeff->BL;
			flagViscoeff = zone->Flowsol->flagViscoeff;
			flagrho = zone->Flowsol->flagRho;
			
			mu  = *(zone->Flowsol->Viscoeff);
			rho = *(zone->Flowsol->Rho);
			
			cellindex = 0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				if(flagViscoeff == 1)
				{
					mu = *(zone->Flowsol->Viscoeff+cellindex);
				}
				if(flagrho == 1)
				{
					rho = *(zone->Flowsol->Rho+cellindex);
				}
				Utau = ;
				vorticity = *(zone->Flowsol->Vorticity+cellindex);
				y  = *(zone->ZGeom->minWallDistance+cellindex);
				lmix  = MixingLength(y,Utau,mu,rho);
				mu_ti = lmix * lmix * vorticity;
				
				*(BL->MixingLength+cellindex) 	= lmix;
				*(BL->mu_t_i+cellindex) 	= mu_ti;
				cellindex++;
			}

			numZbc=(zone)->Number_of_BC_Patches;

			for(zBc=0;zBc<numZbc;zBc++)
			{
				bcPtr = zone->ZoneBc+zBc;
				dir = bcPtr->Direction;

				mu = *(bcPtr->ZSlab->Viscoeff);
				rho = *(bcPtr->ZSlab->Rho);

				index = 0;
				if(dir==1||dir==2)
				{
					for(k=0; k<NK-1; k++) 
					for(j=0; j<NJ-1; j++)
					{
						if(flagViscoeff==1)
							mu = *(bcPtr->ZSlab->Viscoeff+index);
						if(flagrho==1)
							rho = *(bcPtr->ZSlab->Rho+index);
						Utau = ;
						vorticity = ;
						y  = *(bcPtr->ZSlab->minWallDistance+index);
						lmix  = MixingLength(y,Utau,mu,rho);

						index++;
					}
				}
				if(dir==3||dir==4)
				{
					for(k=0; k<NK-1; k++) 
					for(i=0; i<NI-1; i++)
					{
						*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
						index++;
					}
				}
				if(dir==5||dir==6)
				{
					for(j=0; j<NJ-1; j++)
					for(i=0; i<NI-1; i++)
					{
						*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
						index++;
					}
				}
			}
		}
	}
	break;
}


double MixingLength(double y, double Utau, double mu, double rho)
{
	double lmix,yPlus,nu;
	
	nu = mu/rho;
	yPlus = y*Utau/nu;
	
	lmix = 0.40 * y * (1.0 - exp(-yPlus/26.0));
	
	return(lmix);
}

double F_Kleb(double y, double delta)
{
	double F_Kleb;
	
	F_Kleb = y/delta;
	F_Kleb = F_Kleb * F_Kleb;
	F_Kleb = F_Kleb * F_Kleb * F_Kleb;//This is done (instead of using pow() function) to reduce the number of product operations
	
	F_Kleb = 1.0/(1.0 + 5.5 * F_Kleb);//F_Kleb is equal to (y/delta)^6.0
	return(F_Kleb);
}

double F_wake(double y, double delta)
{
}
double F_max(double y, double delta)
{
}
