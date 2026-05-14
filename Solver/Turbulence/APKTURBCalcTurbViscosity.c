#include "APKTURBCalcTurbViscosity.h"

/**
This function calculates Turbulent viscosity using the values of k and epsilon. It loops over all the cells and calculates
the turbulent viscosity from cell valus of k, epsilon and rho. This function also loops over all the interface and zone boundaries
and calculates the slab values of Turbulent viscosity.
*/
int CalcTurbViscosity(Base *base)
{
	int NI,NJ,NK,nZones,numZones,cellindex,totalCells;
	int flagLogKEps,flagViscoeff=0,flagLimTimeScale=0,flagLimTVisc=0;
	int i,j,k,*ptrSize,flagrho,numIface,numZbc,dir;
	int iface,zBc,*ptrRange,index,UpperLimitSlabCount=0,UpperLimitCount=0,LowerLimitSlabCount=0,LowerLimitCount=0;
	Zone *zone;
	double *domainPtr,*slabBcptr,gammastar;
	double *rho,*slabrho,*KE,*Epsilon,*Omega,*slabKE,*slabEpsilon,*slabOmega,fmue,viscosity;
	double *slabViscoeffPtr,Vorticity,arg2,F2=0.0,a1,nu,y,Rho,fv1,nu_t,X;
	double Rt,*f1,Ry,upperlimit=1.0e4,mu,Tau,mu_tmax;
	int nophase;//Added by chandrasekhar --- Turbulent MPFV
	//FILE *fp;
	//double *slabIptr;
	Interface *IfacePtr;
	ZoneBoundary *zBcPtr;
	Zone *ptrZoneNeighbor;

	nZones=((base)->Number_of_zones);
	flagLogKEps = base->flagLogKEps;
	flagLimTimeScale = base->BaseParameters->Turb->flagLimTimeScale;
	flagLimTVisc = base->BaseParameters->Turb->flagLimTVisc;
	
	switch(base->flagTurbModel)
	{
		case 1:	// For Case 1 and Case3 the conditionds are same, hence it goes into Case3 //standard k-Epsilon Model
		
		case 3:		//RNG k-epsilon model
			if(flagLogKEps == 0)
			{
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = base->CgnsZone+numZones;

					if(zone->flagSolidFluid==0)
					{
						flagViscoeff = zone->Flowsol->flagViscoeff;
						ptrSize = zone->Size;
						NI=*(ptrSize+0);
						NJ=*(ptrSize+1);
						NK=*(ptrSize+2);


						domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
						flagrho = zone->Flowsol->flagRho;
						rho = zone->Flowsol->Rho;
						KE = zone->Flowsol->Scal[0];
						Epsilon = zone->Flowsol->Scal[1];
						mu = *(zone->Flowsol->Viscoeff);
						
						cellindex = 0;
						
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							if(flagLimTimeScale == 1 && base->BaseParameters->Time != 0)
							{
								Tau = 1.0/(*(zone->Flowsol->Gamma+cellindex));
								*(domainPtr+cellindex) = *(rho+flagrho*cellindex) * (zone->Flowsol->Tcoeff->Cmu) * *(KE+cellindex) * Tau;
							}
							else
								*(domainPtr+cellindex) = *(rho+flagrho*cellindex) * (zone->Flowsol->Tcoeff->Cmu) * (*(KE+cellindex) * *(KE+cellindex))/(*(Epsilon+cellindex));
							//printf("\nzone-%s $$ *(domainPtr+cellindex)=%e",zone->zonename,*(domainPtr+cellindex));
							if(flagViscoeff == 1)
							{
								mu = *(zone->Flowsol->Viscoeff+cellindex);
							}
							
							if((*(domainPtr+cellindex)) < 0.0)
							{
								printf("k=%e eps=%e\n", *(KE+cellindex),*(Epsilon+cellindex));
								printf("[%s k-%d j-%d i-%d] cellindex=%d Eddy Viscosity=%e, limiting it to %e\n", zone->zonename,k,j,i,cellindex,*(domainPtr+cellindex),1.0e-4*mu);
								*(domainPtr+cellindex) = (1.0e-4* mu);
// 								*(Epsilon+cellindex) = 1.0e-12;// * Uref * Uref * Uref * Uref * *(rho+flagrho*cellindex)/ mu;
								LowerLimitCount++;
							}
							
							if(flagLimTVisc==1)
							{
								//mu_tmax = *(KE+cellindex) * sqrt(1.0/3.0 * *(rho+flagrho*cellindex) * *(zone->ZGeom->Volume+cellindex)/ *(zone->Flowsol->Tcoeff->TProduction+cellindex));
								mu_tmax = upperlimit* mu;
								if(*(domainPtr+cellindex) > mu_tmax)
								{
									//printf("%d [%s %d %d %d] rho=%e k=%e eps=%e Tau=%e Eddy Viscosity=%e, limiting it to 10^5*mu[%e]\n",cellindex,zone->zonename,k,j,i,*(rho+flagrho*cellindex),*(KE+cellindex),*(Epsilon+cellindex),Tau,*(domainPtr+cellindex),upperlimit* mu);
									*(domainPtr+cellindex) = mu_tmax;
									*(KE+cellindex) = sqrt(mu_tmax * *(Epsilon+cellindex)/(0.09 * *(rho+flagrho*cellindex)));
									UpperLimitCount++;
								}
							}
							if(base->flagMPFV==1)
							{
								for(nophase=0;nophase<(base->nphase);nophase++)
								{
									*((zone->Flowsol)->Tcoeff->TViscoeffofphase[nophase]+cellindex)=*(domainPtr+cellindex);//*(*((zone->Flowsol)->RhoofPhase[nophase]+cellindex*0)/(*(rho+flagrho*cellindex)) );
									//printf("\n*((zone->Flowsol)->Tcoeff->TViscoeffofphase[%d]+cellindex)=%e",nophase,*((zone->Flowsol)->Tcoeff->TViscoeffofphase[nophase]+cellindex));
								}
							}
							cellindex++;
						}

						numZbc=(zone)->Number_of_BC_Patches;

						for(zBc=0;zBc<numZbc;zBc++)
						{
							index = 0;
							dir = (zone->ZoneBc+zBc)->Direction;
							ptrRange=(zone->ZoneBc+zBc)->Point_Range;

							slabViscoeffPtr = (zone->ZoneBc+zBc)->ZSlab->Viscoeff;
							slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
							slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
							slabKE =  (zone->ZoneBc+zBc)->ZSlab->Scal[0];
							slabEpsilon =  (zone->ZoneBc+zBc)->ZSlab->Scal[1];
							domainPtr = zone->Flowsol->Tcoeff->TViscoeff;

							if(dir==1||dir==2)
							{
								if(dir==1)
									i=NI-1;
								else
									i=1;
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								{
									if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
									{
										cellindex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
										*(slabBcptr+index) = *(domainPtr+cellindex);
										if(base->flagParticle==1)
										*((zone->ZoneBc+zBc)->ZSlab->T_e+index)=0.125*(*(KE+cellindex)/ *(Epsilon+cellindex) );
										//printf("%e\t{%d %d %d}[%d %d %s %d %d %d]Boundary Slab Turb Visc = %e\n", *(domainPtr+cellindex),NI,NJ,NK,dir,cellindex,zone->zonename,k,j,i,*(slabBcptr+index));
									}
									else
									{
// 										if(flagLimTimeScale == 1 && base->BaseParameters->Time != 0)
// 										{
// 											cellindex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
// 											Tau = 1.0/(*(zone->Flowsol->Gamma+cellindex));
// 											*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * *(slabKE+index) * Tau;
// 										}
// 										else
										*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * (*(slabKE+index) * *(slabKE+index))/(*(slabEpsilon+index));
										//printf("%e\t{%d %d %d}[%d %d %s %d %d %d]Boundary Slab Turb Visc = %e\n", *(domainPtr+cellindex),NI,NJ,NK,dir,cellindex,zone->zonename,k,j,i,*(slabBcptr+index));
									}
									
									if(*(slabBcptr+index) < 0.0 )
									{
										printf("k=%e eps=%e\n", *(slabKE+index),*(slabEpsilon+index));
										printf("dir=%d [%s k-%d j-%d i-%d] index=%d Eddy Viscosity=%e, limiting it to %e\n",dir,zone->zonename,k,j,i,index,*(slabBcptr+index),1.0e-4* *(slabViscoeffPtr+index*flagViscoeff));
										*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
										//*(slabEpsilon+index) = 1.0e-12;
										LowerLimitSlabCount++;
									}
									if(flagLimTVisc==1)
									{
										mu_tmax = upperlimit * *(slabViscoeffPtr+index*flagViscoeff);
										if( *(slabBcptr+index) > mu_tmax)
										{
											//printf("dir=%d [%s k-%d j-%d i-%d] Eddy Viscosity=%e, limiting it to %e\n",dir,zone->zonename,k,j,i,*(slabBcptr+index),upperlimit* *(slabViscoeffPtr+index*flagViscoeff));
											*(slabBcptr+index) = mu_tmax;
											*(slabKE+index) = sqrt(mu_tmax * *(slabEpsilon+index)/(0.09 * *(slabrho+flagrho*index)));
											UpperLimitSlabCount++;
										}
									}
									if(base->flagMPFV==1)
									{
										for(nophase=0;nophase<(base->nphase);nophase++)
										{
										*((zone->ZoneBc+zBc)->ZSlab->TViscoeffofphase[nophase]+index)=*(slabBcptr+index);//*(*((zone->ZoneBc+zBc)->ZSlab->RhoofPhase[nophase]+index*0)*(*(slabrho+flagrho*index)));
										}
									}
									index++;
								}
							}
							if(dir==3||dir==4)
							{
								if(dir==3)
									j=NJ-1;
								else
									j=1;
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
									{
										cellindex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
										*(slabBcptr+index) = *(domainPtr+cellindex);
										if(base->flagParticle==1)
										*((zone->ZoneBc+zBc)->ZSlab->T_e+index)=0.125*(*(KE+cellindex)/ *(Epsilon+cellindex) );
										//printf("%e\t{%d %d %d}[%d %d %s %d %d %d]Boundary Slab Turb Visc = %e\n", *(domainPtr+cellindex),NI,NJ,NK,dir,cellindex,zone->zonename,k,j,i,*(slabBcptr+index));
									}
									else
									{
// 										if(flagLimTimeScale == 1 && base->BaseParameters->Time != 0)
// 										{
// 											cellindex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
// 											Tau = 1.0/(*(zone->Flowsol->Gamma+cellindex));
// 											*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * *(slabKE+index) * Tau;
// 										}
// 										else
										*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * (*(slabKE+index) * *(slabKE+index))/(*(slabEpsilon+index));
										//printf("%e\t{%d %d %d}[%d %d %s %d %d %d]Boundary Slab Turb Visc = %e\n", *(domainPtr+cellindex),NI,NJ,NK,dir,cellindex,zone->zonename,k,j,i,*(slabBcptr+index));	
									}
									
									if( *(slabBcptr+index) < 0.0 )
									{
										printf("k=%e eps=%e\n", *(slabKE+index),*(slabEpsilon+index));
										printf("dir=%d [%s k-%d j-%d i-%d] index=%d Eddy Viscosity=%e, limiting it to %e\n",dir,zone->zonename,k,j,i,index,*(slabBcptr+index),1.0e-4* *(slabViscoeffPtr+index*flagViscoeff));
										*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
										//*(slabEpsilon+index) = 1.0e-12;
										LowerLimitSlabCount++;
									}
									if(flagLimTVisc==1)
									{
										mu_tmax = upperlimit * *(slabViscoeffPtr+index*flagViscoeff);
										if( *(slabBcptr+index) > mu_tmax)
										{
											//printf("dir=%d [%s k-%d j-%d i-%d] Eddy Viscosity=%e, limiting it to %e\n",dir,zone->zonename,k,j,i,*(slabBcptr+index),mu_tmax);
											*(slabBcptr+index) = mu_tmax;
											*(slabKE+index) = sqrt(mu_tmax * *(slabEpsilon+index)/(0.09 * *(slabrho+flagrho*index)));
											UpperLimitSlabCount++;
										}
									}
									if(base->flagMPFV==1)
									{
										for(nophase=0;nophase<(base->nphase);nophase++)
										{
										*((zone->ZoneBc+zBc)->ZSlab->TViscoeffofphase[nophase]+index)=*(slabBcptr+index);//*(*((zone->ZoneBc+zBc)->ZSlab->RhoofPhase[nophase]+index*0)*(*(slabrho+flagrho*index)));
										}
									}
									index++;
								}
							}
							if(dir==5||dir==6)
							{
								if(dir==5)
									k=NK-1;
								else
									k=1;
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
									{
										cellindex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
										*(slabBcptr+index) = *(domainPtr+cellindex);
										if(base->flagParticle==1)
										*((zone->ZoneBc+zBc)->ZSlab->T_e+index)=0.125*(*(KE+cellindex)/ *(Epsilon+cellindex) );
										//printf("%e\t{%d %d %d}[%d %d %s %d %d %d]Boundary Slab Turb Visc = %e\n", *(domainPtr+cellindex),NI,NJ,NK,dir,cellindex,zone->zonename,k,j,i,*(slabBcptr+index));
									}
									else
									{
// 										if(flagLimTimeScale == 1 && base->BaseParameters->Time != 0)
// 										{
// 											cellindex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
// 											Tau = 1.0/(*(zone->Flowsol->Gamma+cellindex));
// 											*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * *(slabKE+index) * Tau;
// 										}
// 										else
										*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * (*(slabKE+index) * *(slabKE+index))/(*(slabEpsilon+index));
										//printf("%e\t{%d %d %d}[%d %d %s %d %d %d]Boundary Slab Turb Visc = %e\n", *(domainPtr+cellindex),NI,NJ,NK,dir,cellindex,zone->zonename,k,j,i,*(slabBcptr+index));	
									}
									
									if( *(slabBcptr+index) < 0.0 )
									{
										printf("k=%e eps=%e\n", *(slabKE+index),*(slabEpsilon+index));
										printf("dir=%d [%s k-%d j-%d i-%d] index=%d Eddy Viscosity=%e, limiting it to %e\n",dir,zone->zonename,k,j,i,index,*(slabBcptr+index),1.0e-4* *(slabViscoeffPtr+index*flagViscoeff));
										*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
										//*(slabEpsilon+index) = 1.0e-12;
										LowerLimitSlabCount++;
									}
									if(flagLimTVisc==1)
									{
										mu_tmax = upperlimit * *(slabViscoeffPtr+index*flagViscoeff);
										if( *(slabBcptr+index) > mu_tmax)
										{
											//printf("dir=%d [%s k-%d j-%d i-%d] Eddy Viscosity=%e, limiting it to %e\n",dir,zone->zonename,k,j,i,*(slabBcptr+index),mu_tmax);
											*(slabBcptr+index) = mu_tmax;
											*(slabKE+index) = sqrt(mu_tmax * *(slabEpsilon+index)/(0.09 * *(slabrho+flagrho*index)));
											UpperLimitSlabCount++;
										}
									}
									if(base->flagMPFV==1)
									{
										for(nophase=0;nophase<(base->nphase);nophase++)
										{
										*((zone->ZoneBc+zBc)->ZSlab->TViscoeffofphase[nophase]+index)=*(slabBcptr+index);//*(*((zone->ZoneBc+zBc)->ZSlab->RhoofPhase[nophase]+index*0)*(*(slabrho+flagrho*index)));
										}
									}
									index++;
								}
							}
						}
					}
				}
				if(UpperLimitCount!=0 || UpperLimitSlabCount!=0)
					printf("*** Max. Turb Viscosity Ratio is limited to %e in %d Cell(s) and %d Slab(s) ***\n", upperlimit,UpperLimitCount,UpperLimitSlabCount);
				if(LowerLimitCount!=0 || LowerLimitSlabCount!=0)
					printf("*** Min. Turb Viscosity Ratio is limited to %e in %d Cell(s) and %d Slab(s) ***\n", 1.0e-4,LowerLimitCount,LowerLimitSlabCount);
			}
			else if(flagLogKEps == 1)
			{	//Calculating Turbulent Viscosity using Logarithmic values of Turbulent Kinetic Energy and its Dissipation
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base)->CgnsZone+numZones;

					if(zone->flagSolidFluid==0)
					{
						ptrSize = zone->Size;
						NI=*(ptrSize+0);
						NJ=*(ptrSize+1);
						NK=*(ptrSize+2);

						domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
						flagrho = zone->Flowsol->flagRho;
						rho = zone->Flowsol->Rho;
						KE = zone->Flowsol->Scal[0];
						Epsilon = zone->Flowsol->Scal[1];

						cellindex = 0;
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							*(domainPtr+cellindex) = *(rho+flagrho*cellindex) * (zone->Flowsol->Tcoeff->Cmu) * exp(*(KE+cellindex)*2.0 - *(Epsilon+cellindex));
							cellindex++;
						}

						numZbc=(zone)->Number_of_BC_Patches;

						for(zBc=0;zBc<numZbc;zBc++)
						{
							index = 0;
							dir = (zone->ZoneBc+zBc)->Direction;
							ptrRange=(zone->ZoneBc+zBc)->Point_Range;

							slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
							slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
							slabKE =  (zone->ZoneBc+zBc)->ZSlab->Scal[0];
							slabEpsilon =  (zone->ZoneBc+zBc)->ZSlab->Scal[1];

							if(dir==1||dir==2)
							{
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								{
									*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * exp(*(slabKE+index)*2.0 - *(slabEpsilon+index));
									index++;
								}
							}
							if(dir==3||dir==4)
							{
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * exp(*(slabKE+index)*2.0 - *(slabEpsilon+index));
									index++;
								}
							}
							if(dir==5||dir==6)
							{
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
									{
										*(slabBcptr+index) = *(slabrho+flagrho*index) * (zone->Flowsol->Tcoeff->Cmu) * exp(*(slabKE+index)*2.0 - *(slabEpsilon+index));
										index++;
									} 
							}
						}
					}
				}
			}
		break;

		case 2:		//k-Omega Model
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = base->CgnsZone+numZones;

				if(zone->flagSolidFluid==0)
				{
					flagViscoeff = zone->Flowsol->flagViscoeff;
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);

					domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
					flagrho = zone->Flowsol->flagRho;
					rho = zone->Flowsol->Rho;
					KE = zone->Flowsol->Scal[0];
					Omega = zone->Flowsol->Scal[1];
					gammastar = zone->Flowsol->Tcoeff->gammastar;

					cellindex = 0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						*(domainPtr+cellindex) = gammastar * *(rho+flagrho*cellindex) * *(KE+cellindex)/(*(Omega+cellindex));

						if((*(domainPtr+cellindex)) < 0.0)
						{
							*(domainPtr+cellindex) = (1.0e-4* *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff));
						}
						cellindex++;
					}

					numZbc=(zone)->Number_of_BC_Patches;

					for(zBc=0;zBc<numZbc;zBc++)
					{
						index = 0;
						dir = (zone->ZoneBc+zBc)->Direction;
						ptrRange=(zone->ZoneBc+zBc)->Point_Range;

						slabViscoeffPtr = (zone->ZoneBc+zBc)->ZSlab->Viscoeff;
						slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
						slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
						slabKE =  (zone->ZoneBc+zBc)->ZSlab->Scal[0];
						slabOmega =  (zone->ZoneBc+zBc)->ZSlab->Scal[1];

						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							{
								*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
								if( *(slabBcptr+index) < 0.0 )
								{
									*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
								}
								index++;
							}
						}
						if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
								if( *(slabBcptr+index) < 0.0 )
								{
									*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
								}
								index++;
							}
						}
						if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
								if( *(slabBcptr+index) < 0.0 )
								{
									*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
								}
								index++;
							}
						}
					}
				}
			}
		break;

		case 4:		// launder and sharma
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = base->CgnsZone+numZones;

				if(zone->flagSolidFluid==0)
				{
					flagViscoeff = zone->Flowsol->flagViscoeff;
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);

					domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
					flagrho = zone->Flowsol->flagRho;
					rho = zone->Flowsol->Rho;
					KE = zone->Flowsol->Scal[0];
					Epsilon = zone->Flowsol->Scal[1];		
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
						viscosity = *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff);
						if(*(Epsilon+cellindex)!=0.0)
						{
							fmue=exp(-3.4/(pow((1.0+0.02 * *(rho+flagrho*cellindex) * pow(*(KE+cellindex),2.0)/(*(Epsilon+cellindex) * viscosity)),2.0)));
							*(domainPtr+cellindex) = *(rho+flagrho*cellindex) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(KE+cellindex),2.0)/(*(Epsilon+cellindex));
						}
						else if((*(KE+cellindex)==0.0) &&(*(Epsilon+cellindex)==0.0))
						{
							*(domainPtr+cellindex) = 0.0;
					//		printf("i=%d\tj=%d\tk=%d\tmue_t=%lf\n",i,j,k,*(domainPtr+cellindex));
						}
						else
							printf("mu_t is not defined in domain\n");
						if((*(domainPtr+cellindex)) < 0.0)//(1.0e-4* mu))
						{
							*(domainPtr+cellindex) = (1.0e-4* viscosity);
						}
						//	printf("[%d %d %d] fmue=%e\tmue_t=%e\n",i,j,k,fmue,*(domainPtr+cellindex));
					}	

//					numIface=(zone)->Number_of_Interfaces;
					numZbc=(zone)->Number_of_BC_Patches;

					for(dir=1;dir<=6;dir++)
					{
// 						for(iface=0;iface<numIface;iface++)
// 						{
// 							if(((zone->Iface+iface)->Direction)==dir)
// 							{
// 					//			printf("printing the interface values for Zone=%d\n",numZones);
// 								ptrRange=(zone->Iface+iface)->Point_Range;
// 								slabIptr=(zone->Iface+iface)->ZSlab->TViscoeff;
// 								slabrho = (zone->Iface+iface)->ZSlab->Rho;
// 								slabKE =  (zone->Iface+iface)->ZSlab->Scal[0];
// 								slabEpsilon =  (zone->Iface+iface)->ZSlab->Scal[1];
// 								//printf("i m in interface loop for MU_T\n");
// 								if(dir==1||dir==2)
// 								{
// 									for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
// 									for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
// 									{
// 										index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
// 										if(*(slabEpsilon+index)!=0.0)
// 										{
// 											viscosity = *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 											fmue=exp(-3.4/(pow((1.0+0.02* *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity)),2.0)));
// 											*(slabIptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
// 										}
// 										else if((*(slabEpsilon+index)==0.0)&&(*(slabKE+index)==0.0))
// 											*(slabIptr+index)=0.0;
// 										else
// 											printf("mu_t is not defined at interface\n");
// 										if((*(slabIptr+index)) < 0.0)//(1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff)))
// 										{
// 											*(slabIptr+index) = 1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 										}
// 					//					printf("inteface dir=%d\tindex=%d\tfmue=%lf\tmue_t=%lf\n",dir,index,fmue,*(slabIptr+index));
// 									}
// 								}
// 								if(dir==3||dir==4)
// 								{
// 									for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
// 									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
// 									{
// 										index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
// 										if(*(slabEpsilon+index)!=0.0)
// 										{
// 											viscosity = *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 											fmue=exp(-3.4/(pow((1.0+0.02* *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity)),2.0)));
// 											*(slabIptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
// 										}
// 										else if((*(slabEpsilon+index)==0.0)&&(*(slabKE+index)==0.0))
// 											*(slabIptr+index)=0.0;
// 										else
// 											printf("mu_t is not defined at interface\n");
// 										if((*(slabIptr+index)) < 0.0)// < (1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff)))
// 										{
// 											*(slabIptr+index) = 1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 										}
// 					//					printf("inteface dir=%d\tindex=%d\tfmue=%lf\tmue_t=%lf\n",dir,index,fmue,*(slabIptr+index));
// 									//	printf("dir=%d\tindex=%d\tmue_t=%lf\n",dir,index,*(slabIptr+index));
// 									}
// 								}
// 								if(dir==5||dir==6)
// 								{
// 									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
// 									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
// 									{
// 										index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
// 										if(*(slabEpsilon+index)!=0.0)
// 										{
// 											viscosity = *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 											fmue=exp(-3.4/(pow((1.0+0.02* *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity)),2.0)));
// 											*(slabIptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
// 										}
// 										else if((*(slabEpsilon+index)==0.0)&&(*(slabKE+index)==0.0))
// 											*(slabIptr+index)=0.0;
// 										else
// 											printf("mu_t is not defined at interface\n");
// 										if((*(slabIptr+index)) < 0.0)// < (1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff)))
// 										{
// 											*(slabIptr+index) = 1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 										}
// 					//					printf("inteface dir=%d\tindex=%d\tfmue=%lf\tmue_t=%lf\n",dir,index,fmue,*(slabIptr+index));
// 									//	printf("dir=%d\tindex=%d\tmue_t=%lf\n",dir,index,*(slabIptr+index));
// 									}
// 								}							
// 							}
// 						}
						for(zBc=0;zBc<numZbc;zBc++)
						{
							if((zone->ZoneBc+zBc)->Direction==dir)
							{
								//printf("printing the Bc values for Zone=%d\n",numZones);
								ptrRange=(zone->ZoneBc+zBc)->Point_Range;
								slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
								//slabBcptr = (zone->ZoneBc+zBc)->ZSlab->EffViscoeff;
								slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
								slabKE =  (zone->ZoneBc+zBc)->ZSlab->Scal[0];
								slabEpsilon =  (zone->ZoneBc+zBc)->ZSlab->Scal[1];

								if(dir==1||dir==2)
								{
									for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									{
										index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
										if(*(slabEpsilon+index)!=0.0)
										{
											viscosity = *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff);
											fmue=exp(-3.4/(pow((1.0+0.02* *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity)),2.0)));
											*(slabBcptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
										}
										else if((*(slabEpsilon+index)==0.0)&&(*(slabKE+index)==0.0))
										{
											*(slabBcptr+index) = 0.0;
											//printf("dir=%d\tj=%d\tk=%d\tmue_t=%lf\n",dir,j,k,*(slabBcptr+index));
										}
										else 
											printf("mu_t is not defined at zone BC\n");
										if((*(slabBcptr+index))  < 0.0)//< (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff)))
										{
											*(slabBcptr+index) = (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff));
										}
										//printf("Zone BC dir=%d\tindex=%d\tfmue=%lf\tmue_t=%lf\n",dir,index,fmue,*(slabBcptr+index));
									}
								}
								if(dir==3||dir==4)
								{
									for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
									{
										index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
										if(*(slabEpsilon+index)!=0.0)
										{
											viscosity = *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff);
											fmue=exp(-3.4/(pow((1.0+0.02* *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity)),2.0)));
											*(slabBcptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
										}
										else if((*(slabEpsilon+index)==0.0)&&(*(slabKE+index)==0.0))
										{
											*(slabBcptr+index) = 0.0;
											//printf("dir=%d\ti=%d\tk=%d\tmue_t=%lf\n",dir,i,k,*(slabBcptr+index));
										}
										else 
											printf("mu_t is not defined at zone BC\n");
										if((*(slabBcptr+index)) < 0.0)// < (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff)))
										{
											*(slabBcptr+index) = (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff));
										}
										//printf("Zone BC dir=%d\tindex=%d\tfmue=%lf\tmue_t=%lf\n",dir,index,fmue,*(slabBcptr+index));
									}
								}
								if(dir==5||dir==6)
								{
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
									{
										index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
										if(*(slabEpsilon+index)!=0.0)
										{
											viscosity = *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff);
											fmue=exp(-3.4/(pow((1.0+0.02* *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity)),2.0)));
											*(slabBcptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
										}
										else if((*(slabEpsilon+index)==0.0)&&(*(slabKE+index)==0.0))
										{
											*(slabBcptr+index) = 0.0;
											//printf("dir=%d\ti=%d\tj=%d\tmue_t=%lf\n",dir,i,j,*(slabBcptr+index));
										}
										else 
											printf("mu_t is not defined at Zone BC\n");
										if((*(slabBcptr+index))  < 0.0)//< (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff)))
										{
											*(slabBcptr+index) = (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff));
										}
										//printf("Zone BC dir=%d\tindex=%d\tfmue=%lf\tmue_t=%lf\n",dir,index,fmue,*(slabBcptr+index));
									}
								}
							}
						}
					}
				}
			}
		break;

		case 5://Lam-Bremhorst(modified by santo)
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = base->CgnsZone+numZones;
				f1 = zone->Flowsol->Tcoeff->f1;
				if(zone->flagSolidFluid==0)	
				{
					flagViscoeff = zone->Flowsol->flagViscoeff;
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);
					domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
					flagrho = zone->Flowsol->flagRho;
					rho = zone->Flowsol->Rho;
					KE = zone->Flowsol->Scal[0];
					Epsilon = zone->Flowsol->Scal[1];
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
						viscosity = *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff);
						if(*(Epsilon+cellindex)!=0.0)
						{
							Rt = *(rho+flagrho*cellindex) * pow(*(KE+cellindex),2.0)/(*(Epsilon+cellindex) * viscosity);
							y = *(zone->ZGeom->minWallDistance+cellindex);
							Ry= *(rho+flagrho*cellindex) * y * sqrt(fabs(*(KE+cellindex)))/viscosity;
							if(Ry>1.0)
							fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
							else 
							fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
							if((Rt<1.0)&&(fmue>1.0))
								fmue=0.01116225;
							*(domainPtr+cellindex) = *(rho+flagrho*cellindex) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(KE+cellindex),2.0)/(*(Epsilon+cellindex));
						}
						else 
							*(domainPtr+cellindex) =0.0;
						if((*(domainPtr+cellindex)) < 0.0)
							*(domainPtr+cellindex) = (1.0e-4* viscosity);
						
					
					}

					numIface=(zone)->Number_of_Interfaces;
					numZbc=(zone)->Number_of_BC_Patches;

					for(dir=1;dir<=6;dir++)
					{
// 						for(iface=0;iface<numIface;iface++)
// 						{
// 							if(((zone->Iface+iface)->Direction)==dir)
// 							{
// 								ptrRange=(zone->Iface+iface)->Point_Range;
// 
// 								slabIptr=(zone->Iface+iface)->ZSlab->TViscoeff;
// 								slabrho = (zone->Iface+iface)->ZSlab->Rho;
// 								slabKE =  (zone->Iface+iface)->ZSlab->Scal[0];
// 								slabEpsilon =  (zone->Iface+iface)->ZSlab->Scal[1];
// 								
// 								if(dir==1||dir==2)
// 								{
// 									for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
// 									for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
// 									{
// 										index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
// 										if(*(slabEpsilon+index)!=0.0)
// 										{
// 											viscosity = *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 											Rt = *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity);
// 											y = *((zone->Iface+iface)->ZSlab->minWallDistance+index);
// 											Ry= *(slabrho+flagrho*index) * y * sqrt(*(slabKE+index))/viscosity;
// 											if(Ry>1.0)
// 											fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
// 											else 
// 											fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
// 											if((Rt<1.0)&&(fmue>1.0))
// 												fmue=0.01116225;
// 											*(slabIptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
// 										}
// 										else 
// 											*(slabIptr+index) = 0.0;
// 										if((*(slabIptr+index)) < 0.0)
// 											*(slabIptr+index) = 1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 										
// 										
// 				
// 									}
// 								}
// 								if(dir==3||dir==4)
// 								{
// 									for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
// 									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
// 									{
// 										index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
// 										if(*(slabEpsilon+index)!=0.0)
// 										{
// 											viscosity = *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 											Rt = *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity);
// 											y = *((zone->Iface+iface)->ZSlab->minWallDistance+index);
// 											Ry= *(slabrho+flagrho*index) * y * sqrt(*(slabKE+index))/viscosity;
// 											if(Ry>1.0)
// 											fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
// 											else 
// 											fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
// 											if((Rt<1.0)&&(fmue>1.0))
// 												fmue=0.01116225;
// 											*(slabIptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
// 										}
// 										else 
// 											*(slabIptr+index) = 0.0;
// 											
// 										if((*(slabIptr+index)) < 0.0)
// 											*(slabIptr+index) = 1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 										
// 									}
// 								}
// 								if(dir==5||dir==6)
// 								{
// 									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
// 									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
// 									{
// 										index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
// 										if(*(slabEpsilon+index)!=0.0)
// 										{
// 											viscosity = *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 											Rt = *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity);
// 											y = *((zone->Iface+iface)->ZSlab->minWallDistance+index);
// 											Ry= *(slabrho+flagrho*index) * y * sqrt(*(slabKE+index))/viscosity;
// 											if(Ry>1.0)
// 											fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
// 											else 
// 											fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
// 											if((Rt<1.0)&&(fmue>1.0))
// 												fmue=0.01116225;
// 											*(slabIptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
// 										}
// 										else 
// 											*(slabIptr+index) = 0.0;
// 										if((*(slabIptr+index)) < 0.0)
// 											*(slabIptr+index) = 1.0e-4* *((zone->Iface+iface)->ZSlab->Viscoeff+index*flagViscoeff);
// 										
// 									}
// 								}
// 							}
// 						}
						for(zBc=0;zBc<numZbc;zBc++)//calculating turbulent viscsity at the boundaries
						{
							if((zone->ZoneBc+zBc)->Direction==dir)
							{
								//printf("printing the Bc values for Zone=%d\n",numZones);
								ptrRange=(zone->ZoneBc+zBc)->Point_Range;
								slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
								//slabBcptr = (zone->ZoneBc+zBc)->ZSlab->EffViscoeff;
								slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
								slabKE =  (zone->ZoneBc+zBc)->ZSlab->Scal[0];
								slabEpsilon =  (zone->ZoneBc+zBc)->ZSlab->Scal[1];


								if(dir==1||dir==2)
								{
									for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									{
										index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
										viscosity = *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff);
										if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
										{
											*(slabBcptr+index) =0.0;
										}
										else
										{
											if(*(Epsilon+cellindex)!=0.0)
											{
												Rt = *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity);
												y = *((zone->ZoneBc+zBc)->ZSlab->minWallDistance+index);
												Ry= *(slabrho+flagrho*index) * y * sqrt(*(slabKE+index))/viscosity;
												if(Ry>1.0)
												fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
												else 
												fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
												if((Rt<1.0)&&(fmue>1.0))
													fmue=0.01116225;
												*(slabBcptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
											}
											else 
												*(slabBcptr+index) = 0.0;
											
											if((*(slabBcptr+index))  < 0.0)
												*(slabBcptr+index) = (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff));
											
											
											
										}
								
									}
								}
								if(dir==3||dir==4)
								{
									for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
									{
										index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
										viscosity = *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff);
										if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
										{
											*(slabBcptr+index) =0.0;
										}
										else
										{
											if(*(Epsilon+cellindex)!=0.0)
											{
												Rt = *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity);
												y = *((zone->ZoneBc+zBc)->ZSlab->minWallDistance+index);
												Ry= *(slabrho+flagrho*index) * y * sqrt(*(slabKE+index))/viscosity;
												if(Ry>1.0)
												fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
												else 
												fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
												if((Rt<1.0)&&(fmue>1.0))
													fmue=0.01116225;
												*(slabBcptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
											}
											else 
												*(slabBcptr+index) =0.0;
																						
											if((*(slabBcptr+index))  < 0.0)
												*(slabBcptr+index) = (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff));
											
										}
								
									}
								}
								if(dir==5||dir==6)
								{
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
									{
										index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
										viscosity = *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff);
										if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
										{
											*(slabBcptr+index) =0.0;
										}
										else
										{
											if(*(Epsilon+cellindex)!=0.0)
											{
												Rt = *(slabrho+flagrho*index) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index) * viscosity);
												y = *((zone->ZoneBc+zBc)->ZSlab->minWallDistance+index);
												Ry= *(rho+flagrho*cellindex) * y * sqrt(*(KE+cellindex))/viscosity;
												if(Ry>1.0)
												fmue = pow((1.0 - exp(-0.0165*Ry)),2.0)*(1.0 + 20.5/Rt);
												else 
												fmue=(0.0165*Ry)*(0.0165*Ry)*(1.0 + 20.5/Rt);
												if((Rt<1.0)&&(fmue>1.0))
													fmue=0.01116225;
												*(slabBcptr+index) = *(slabrho+flagrho*index) * fmue * (zone->Flowsol->Tcoeff->Cmu) * pow(*(slabKE+index),2.0)/(*(slabEpsilon+index));
											}
											else 
												*(slabBcptr+index) = 0.0;
											if((*(slabBcptr+index))  < 0.0)
												*(slabBcptr+index) = (1.0e-4* *((zone->ZoneBc+zBc)->ZSlab->Viscoeff+index*flagViscoeff));
												
											
										}
								
									} 
								}
							}
						}
					}
				}
			}
			
		break;
		case 6:		//SST Turbulence Model
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = base->CgnsZone+numZones;

				if(zone->flagSolidFluid==0)
				{
					flagViscoeff = zone->Flowsol->flagViscoeff;
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);

					domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
					flagrho = zone->Flowsol->flagRho;
					rho = zone->Flowsol->Rho;
					KE = zone->Flowsol->Scal[0];
					Omega = zone->Flowsol->Scal[1];
					a1 = zone->Flowsol->Tcoeff->a1;
					gammastar = zone->Flowsol->Tcoeff->gammastar;

					cellindex = 0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						Vorticity = *(zone->Flowsol->Vorticity+cellindex);
						nu = *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff)/(*(rho+flagrho*cellindex));
						y = *(zone->ZGeom->minWallDistance+cellindex);
						arg2 = max(2.0*sqrt(*(KE+cellindex))/(0.09* *(Omega+cellindex)*y) , 500.0*nu/(y*y* *(Omega+cellindex)));
						F2 = tanh(arg2*arg2);

						*(domainPtr+cellindex) = *(rho+flagrho*cellindex) * a1 * *(KE+cellindex)/max(a1 * *(Omega+cellindex),Vorticity*F2);
						
						if((*(domainPtr+cellindex)) < 0.0)
						{
							*(domainPtr+cellindex) = (1.0e-4* *(zone->Flowsol->Viscoeff+cellindex*flagViscoeff));
						}
						cellindex++;
					}

					numZbc=(zone)->Number_of_BC_Patches;

					for(zBc=0;zBc<numZbc;zBc++)
					{
						index = 0;
						dir = (zone->ZoneBc+zBc)->Direction;
						ptrRange=(zone->ZoneBc+zBc)->Point_Range;

						slabViscoeffPtr = (zone->ZoneBc+zBc)->ZSlab->Viscoeff;
						slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
						slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
						slabKE =  (zone->ZoneBc+zBc)->ZSlab->Scal[0];
						slabOmega =  (zone->ZoneBc+zBc)->ZSlab->Scal[1];

						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							{
								*(slabBcptr+index) = *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
								if( *(slabBcptr+index) < 0.0 )
								{
									*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
								}
								index++;
							}
						}
						if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								*(slabBcptr+index) = *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
								if( *(slabBcptr+index) < 0.0 )
								{
									*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
								}
								index++;
							}
						}
						if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								*(slabBcptr+index) = *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
								if( *(slabBcptr+index) < 0.0 )
								{
									*(slabBcptr+index) = 1.0e-4 * *(slabViscoeffPtr+index*flagViscoeff);
								}
								index++;
							}
						}
					}
				}
			}
		break;

		case 8:		//Spalart-Allmaras Model
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = base->CgnsZone+numZones;
				
				if(zone->flagSolidFluid==0)
				{
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);
					
					flagrho = zone->Flowsol->flagRho;
					flagViscoeff = zone->Flowsol->flagViscoeff;
					
					domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
					KE = zone->Flowsol->Scal[0];//Current time step Kinematic Eddy Viscosity
					//Epsilon = zone->PrevFlowsol->Scal[0];//Previous time step Kinematic Eddy Viscosity
					rho = zone->Flowsol->Rho;
					Rho = *(rho);
					mu = *(zone->Flowsol->Viscoeff);
					
					cellindex = 0;
					
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(flagrho == 1)
							Rho = *(rho+cellindex);
						if(flagViscoeff == 1)
							mu = *(zone->Flowsol->Viscoeff+cellindex);
						
						nu = mu/Rho;
						nu_t = *(KE+cellindex);
						
						X = nu_t/nu;//X = *(Epsilon+cellindex)/nu;
						y = X*X*X;
						fv1 = y/(y+357.911);//fv1 = y/(y+Cv1*Cv1*Cv1);//Cv1 = 7.1;
						
						*(domainPtr+cellindex) = Rho * nu_t * fv1;
						
						cellindex++;
					}
					
					numZbc = zone->Number_of_BC_Patches;
					
					for(zBc=0;zBc < numZbc;zBc++)
					{
						zBcPtr = zone->ZoneBc+zBc;
						dir = zBcPtr->Direction;
						
						slabBcptr = zBcPtr->ZSlab->TViscoeff;
						
						switch(dir)
						{
							case 1://Imax face
							case 2://Imin face
								totalCells = (NJ-1) * (NK-1);
								break;
							case 3://Jmax face
							case 4://Jmin face
								totalCells = (NI-1) * (NK-1);
								break;
							case 5://Kmax face
							case 6://Kmin face
								totalCells = (NI-1) * (NJ-1);
								break;
						}
		
						if(zBcPtr->Boundary_Tag == 22)
						{
							for(index=0; index < totalCells; index++)
							{
								*(slabBcptr+index) = 0.0;
							}
						}
						else
						{
							slabKE =  zBcPtr->ZSlab->Scal[0];
							//slabEpsilon =  zBcPtr->ZSlab->ScalOld[0];
							
							slabrho = zBcPtr->ZSlab->Rho;
							slabViscoeffPtr = zBcPtr->ZSlab->Viscoeff;
							
							for(index=0; index < totalCells; index++)
							{
								if(flagrho == 1)
									Rho = *(slabrho+index);
								if(flagViscoeff == 1)
									mu  = *(slabViscoeffPtr+index);
								
								nu = mu/Rho;
								nu_t = *(slabKE+index);
								
								X = nu_t/nu;//X = *(slabEpsilon+index)/nu;
								y = X*X*X;
								fv1 = y/(y+357.911);//fv1 = y/(y+Cv1*Cv1*Cv1);//Cv1 = 7.1;
								
								*(slabBcptr+index) = Rho * nu_t * fv1;
							}
						}
					}
				}
			}
		break;
		
// 		case 9:		//Baldwin-Lomax (Zero Equation) Model
// 			for(numZones=0;numZones<nZones;numZones++)
// 			{
// 				zone = base->CgnsZone+numZones;
// 
// 				if(zone->flagSolidFluid==0)
// 				{
// 					ptrSize = zone->Size;
// 					NI=*(ptrSize+0);
// 					NJ=*(ptrSize+1);
// 					NK=*(ptrSize+2);
// 
// 					flagViscoeff = zone->Flowsol->flagViscoeff;
// 					domainPtr = zone->Flowsol->Tcoeff->TViscoeff;
// 					flagrho = zone->Flowsol->flagRho;
// 					rho = zone->Flowsol->Rho;
// 					
// 					mu  = *(zone->Flowsol->Viscoeff);
// 					Rho = *rho;
// 					
// 					cellindex = 0;
// 					for(k=0;k<NK-1;k++)
// 					for(j=0;j<NJ-1;j++)
// 					for(i=0;i<NI-1;i++)
// 					{
// 						if(flagViscoeff == 1)
// 						{
// 							mu = *(zone->Flowsol->Viscoeff+cellindex);
// 						}
// 						if(flagrho == 1)
// 						{
// 							Rho = *(rho+cellindex);
// 						}
// 						
// 						y  = *(zone->ZGeom->minWallDistance+cellindex);
// 						lmix = MixingLength(y,Utau,mu,Rho);
// 						vorticity = *(zone->Flowsol->Vorticity+cellindex);
// 						
// 						mu_ti = lmix * lmix * vorticity;
// 						
// 						ymax = Fmax();
// 						Udif = Udif();
// 						
// 						F_wake = F_wake();
// 						F_Kleb = F_Kleb(y,ymax/C_Kleb);
// 						
// 						mu_to = alpha * Ccp * Fwake * Fkleb;
// 						
// 						*(domainPtr+cellindex) = InnerOROuter(mu_ti,mu_to,y,ym,cellindex) ;
// 						cellindex++;
// 					}
// 
// 					numZbc=(zone)->Number_of_BC_Patches;
// 
// 					for(zBc=0;zBc<numZbc;zBc++)
// 					{
// 						dir = (zone->ZoneBc+zBc)->Direction;
// 
// 						slabViscoeffPtr = (zone->ZoneBc+zBc)->ZSlab->Viscoeff;
// 						slabBcptr = (zone->ZoneBc+zBc)->ZSlab->TViscoeff;
// 						slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
// 
// 						index = 0;
// 						if(dir==1||dir==2)
// 						{
// 							for(k=0; k<NK-1; k++) 
// 							for(j=0; j<NJ-1; j++)
// 							{
// 								*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
// 								index++;
// 							}
// 						}
// 						if(dir==3||dir==4)
// 						{
// 							for(k=0; k<NK-1; k++) 
// 							for(i=0; i<NI-1; i++)
// 							{
// 								*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
// 								index++;
// 							}
// 						}
// 						if(dir==5||dir==6)
// 						{
// 							for(j=0; j<NJ-1; j++)
// 							for(i=0; i<NI-1; i++)
// 							{
// 								*(slabBcptr+index) = gammastar * *(slabrho+flagrho*index) * *(slabKE+index)/(*(slabOmega+index));
// 								index++;
// 							}
// 						}
// 					}
// 				}
// 			}
// 			break;

		default:
			printf("Improper Turbulence Model Tag in the function CalcTurbViscosity()\n");
			exit(748);
		break;
	}

	//Updating Turbulent Viscosity in the Interfaces of all FLUID Zones
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = base->CgnsZone+numZones;

		if(zone->flagSolidFluid==0)
		{
			numIface = zone->Number_of_Interfaces;

			for(iface=0;iface<numIface;iface++)
			{
				IfacePtr = zone->Iface+iface;
				ptrZoneNeighbor = (Zone *)FindZonePtr(IfacePtr->Adjoining_Zone_Name,base);

				if(ptrZoneNeighbor->flagSolidFluid == 0)
				{
					FillerInterfaceVariables(715,IfacePtr,ptrZoneNeighbor,(IfacePtr->Direction)-1);
				}
				else
				{
					//If the adjacent zone is SOLID, boundary condition has to be applied for turbulent viscosity
					//i.e., Turbulent Viscosity can be made equal to 'Zero' for all turbulence models except standard and RNG k-epsilon models.
					if(base->flagTurbModel==1 || base->flagTurbModel==3)	//Check the Interface Condition for "log version of standard k-epsilon model"
					{
// 						printf("Updating turbulent viscosity in Zone: %s \n", zone->zonename);
						HomogeneousNeumannIC(zone,IfacePtr,715);//NonOrthogonalHomogeneousNeumannIC(base,zone,IfacePtr,715);//HomogeneousNeumannIC(zone,IfacePtr,715);
// 						printf("Updation is over ...\n");
					}
					else
						HomogeneousDirichletIC(715,IfacePtr);
				}
				//Added by chandrasekhar--- To assign the mixture turbulent viscosity to phase turbulent viscosity
				if(base->flagMPFV==1)
				{
					dir=IfacePtr->Direction;
					ptrRange=IfacePtr->Point_Range;
					index = 0;
					if(dir==1||dir==2)
					{
						if(dir==1)
							i=NI-1;
						else
							i=1;
						for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
						for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
						{
							for(nophase=0;nophase<(base->nphase);nophase++)
							{	
								*(IfacePtr->ZSlab->TViscoeffofphase[nophase]+index)=*(IfacePtr->ZSlab->TViscoeff+index);//*(*(IfacePtr->ZSlab->RhoofPhase[nophase]+index*0)/(*(IfacePtr->ZSlab->Rho+index)));
							//printf("\nZone-%s$$Ifacedir-%d$$Tv[nophase]=%lf",zone->zonename,dir,*(IfacePtr->ZSlab->TViscoeffofphase[nophase]+index));
							}
							index++;
						}
					}
					if(dir==3||dir==4)
					{
						if(dir==3)
							j=NJ-1;
						else
							j=1;
						for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
						for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
						{
							for(nophase=0;nophase<(base->nphase);nophase++)
							{
								*(IfacePtr->ZSlab->TViscoeffofphase[nophase]+index)=*(IfacePtr->ZSlab->TViscoeff+index);//*(*(IfacePtr->ZSlab->RhoofPhase[nophase]+index*0)/(*(IfacePtr->ZSlab->Rho+index)));
							//printf("\nZone-%s$$Ifacedir-%d$$Tv[nophase]=%lf",zone->zonename,dir,*(IfacePtr->ZSlab->TViscoeffofphase[nophase]+index));
							}
							index++;
						}
					}
					if(dir==5||dir==6)
					{
						if(dir==5)
							k=NK-1;
						else
							k=1;
						for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
						for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
						{
							for(nophase=0;nophase<(base->nphase);nophase++)
							{
								*(IfacePtr->ZSlab->TViscoeffofphase[nophase]+index)=*(IfacePtr->ZSlab->TViscoeff+index);//*(*(IfacePtr->ZSlab->RhoofPhase[nophase]+index*0)/(*(IfacePtr->ZSlab->Rho+index)));
							//printf("\nZone-%s$$Ifacedir-%d$$Tv[nophase]=%lf",zone->zonename,dir,*(IfacePtr->ZSlab->TViscoeffofphase[nophase]+index));
							}
							index++;
						}
					}
				
				}
				//End --- chandrasekhar
			}
		}
	}

	return 0;
}


