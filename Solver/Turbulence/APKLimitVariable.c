#include "APKLimitVariable.h"

int LimitVariable(Base *base,int vartag)
{
	int nZones,numZones,i,j,k,NI,NJ,NK,count,*ptrSize,index;
	int LowerLimitSlabCount=0,UpperLimitSlabCount=0,LowerLimitCount=0,UpperLimitCount=0,totalCells;
	int flagRho,flagViscoeff,flagTurbLimiter;
	double upperlimit,lowerlimit,mu,Uref=10.0,*KE,*Epsilon,kmin,epsmin,kmax,Epsmax,rho,mu_t,vol,TProduction;
	Zone *zone;
	Interface *IfacePtr;
	ZoneBoundary *bcPtr;
	//FILE *fp;
	
// #ifdef FILE_OPEN
// {
// 	fp = fopen("Settings/TurbLimiter.dat", "r");
// 	fscanf(fp,"%d", &flagTurbLimiter);
// 	fscanf(fp,"%lf", &Uref);
// 	fclose(fp);
// }
// #endif	
	
	flagTurbLimiter = base->BaseParameters->Turb->flagLimKEps;

	if(flagTurbLimiter == 1)
	{
		kmin   = base->BaseParameters->Turb->kmin;
		epsmin = base->BaseParameters->Turb->epsmin;
		
		Uref   = base->BaseParameters->Turb->Uref;
		kmax   = 0.1*Uref*Uref;//0.2*Uref*Uref/2.0;
		//kmax = base->BaseParameters->Turb->kmax;
		
		//if(base->flagTurbModel == 1 || base->flagTurbModel == 3 || base->flagTurbModel == 4 || base->flagTurbModel == 5)
		if(base->flagTurbModel == 1 || base->flagTurbModel == 3)
		{
			nZones = base->Number_of_zones;
			switch(vartag)
			{
				case 600:
					lowerlimit = kmin;// * Uref*Uref;//1.0e-10;
					upperlimit = kmax;
					
					for(numZones=0;numZones<nZones;numZones++)
					{
						zone = base->CgnsZone+numZones;
						
						if(zone->flagSolidFluid==0)
						{
							ptrSize = zone->Size;
							NI=*(ptrSize+0);
							NJ=*(ptrSize+1);
							NK=*(ptrSize+2);
							
							KE = zone->Flowsol->Scal[0];
							Epsilon = zone->Flowsol->Scal[1];
							
							flagRho = zone->Flowsol->flagRho;
							flagViscoeff = zone->Flowsol->flagViscoeff;
							rho = *(zone->Flowsol->Rho);
							
							index = 0;
							for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
							{
								//Limiting turbulent kinetic energy using lower bound value
								if(*(KE+index) < lowerlimit)
								{
									//printf("TKE Lower Limit:%e %s NI=%d NJ=%d NK=%d k=%d j=%d i=%d TKE=%e\n", lowerlimit,zone->zonename,NI,NJ,NK,k,j,i,*(KE+index));
									*(KE+index) = lowerlimit;
									LowerLimitCount++;
								}
								//Limiting turbulent kinetic energy using upper bound value
								if(*(KE+index) > upperlimit)
								{
									//printf("TKE Upper Limit:%e %s NI=%d NJ=%d NK=%d k=%d j=%d i=%d TKE=%e\n", upperlimit,zone->zonename,NI,NJ,NK,k,j,i,*(KE+index));
									*(KE+index) = upperlimit;
									
									//Calculating the upper limit for epsilon
// 									if(flagRho == 1 || flagViscoeff == 1)
// 									{
// 										rho = *(zone->Flowsol->Rho+index*flagRho);
// 									}
									
									mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+index);
									TProduction = *(zone->Flowsol->Tcoeff->TProduction+index);
									vol = *(zone->ZGeom->Volume+index);
									
									Epsmax = mu_t*TProduction/vol;//rho*mu_t*TProduction/vol;

									Epsmax = 0.09 * rho*kmax*kmax/mu_t;
									if(*(Epsilon+index) < Epsmax)
										*(Epsilon+index) = Epsmax;
									
									UpperLimitCount++;
								}
								index++;
							}
			
							for(count=0; count < zone->Number_of_BC_Patches; count++)
							{
								bcPtr = zone->ZoneBc+count;
								KE    = bcPtr->ZSlab->Scal[0];
								Epsilon = bcPtr->ZSlab->Scal[1];
								
								totalCells = bcPtr->ZSlab->Total_Cells;
								
								for(index=0; index < totalCells; index++)
								{
									//Limiting turbulent kinetic energy using lower bound value
									if(*(KE+index) < lowerlimit)
									{
										//printf("TKE Lower Limit:%e %s Boundary Dir=%d slab index=%d TKE=%e\n", lowerlimit,zone->zonename,bcPtr->Direction,index,*(KE+index));
										*(KE+index) = fabs(*(KE+index));//lowerlimit;
										LowerLimitSlabCount++;
									}
									//Limiting turbulent kinetic energy using upper bound value
									if(*(KE+index) > upperlimit)
									{
										//printf("TKE Upper Limit:%e %s Boundary Dir=%d slab index=%d TKE=%e\n", upperlimit,zone->zonename,bcPtr->Direction,index,*(KE+index));
										*(KE+index) = upperlimit;
										
										//Calculating the upper limit for epsilon
// 										if(flagRho == 1 || flagViscoeff == 1)
// 										{
// 											rho = *(bcPtr->ZSlab->Rho+index*flagRho);
// 										}
										
										mu_t = *(bcPtr->ZSlab->TViscoeff+index);
										Epsmax = 0.09 * rho*kmax*kmax/mu_t;
										
										if(*(Epsilon+index) < Epsmax)
											*(Epsilon+index) = Epsmax;
										
										UpperLimitSlabCount++;
									}
								}
							}
							
							for(count=0; count < zone->Number_of_Interfaces; count++)
							{
								IfacePtr = zone->Iface+count;
								if(IfacePtr->ConjugateIfaceTag == 1)
								{
									KE    =  IfacePtr->ZSlab->Scal[0];
									Epsilon = IfacePtr->ZSlab->Scal[1];
									
									totalCells = IfacePtr->ZSlab->Total_Cells;
									
									for(index=0; index < totalCells; index++)
									{
										//Limiting turbulent kinetic energy using lower bound value
										if(*(KE+index) < lowerlimit)
										{
											//printf("TKE Lower Limit:%e %s Iface Dir=%d slab index=%d TKE=%e\n", lowerlimit,zone->zonename,bcPtr->Direction,index,*(KE+index));
											*(KE+index) = lowerlimit;//fabs(*(KE+index));//lowerlimit;
											LowerLimitSlabCount++;
										}
										//Limiting turbulent kinetic energy using upper bound value
										if(*(KE+index) > upperlimit)
										{
											//printf("TKE Upper Limit:%e %s Iface Dir=%d slab index=%d TKE=%e\n", upperlimit,zone->zonename,bcPtr->Direction,index,*(KE+index));
											*(KE+index) = upperlimit;
										
											//Calculating the upper limit for epsilon
// 											if(flagRho == 1 || flagViscoeff == 1)
// 											{
// 												rho = *(IfacePtr->ZSlab->Rho+index*flagRho);
// 											}
											
											mu_t = *(IfacePtr->ZSlab->TViscoeff+index);
											Epsmax = 0.09 * rho*kmax*kmax/mu_t;
											
											if(*(Epsilon+index) < Epsmax)
												*(Epsilon+index) = Epsmax;
												
											UpperLimitSlabCount++;
										}
									}
							
								}
							}
						}
					}
					if(UpperLimitCount!=0 || UpperLimitSlabCount!=0)
						printf(">>> Max. Turbulent Kinetic Energy is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", upperlimit,UpperLimitCount,UpperLimitSlabCount);
					if(LowerLimitCount!=0 || LowerLimitSlabCount!=0)
						printf(">>> Min. Turbulent Kinetic Energy is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", lowerlimit,LowerLimitCount,LowerLimitSlabCount);
					break;
				case 601:
					lowerlimit = epsmin;// * Uref*Uref*Uref*Uref/(1.0e-4/780.0);//1.0e-12;
					
					for(numZones=0;numZones<nZones;numZones++)
					{
						zone = base->CgnsZone+numZones;
						
						if(zone->flagSolidFluid==0)
						{
							ptrSize = zone->Size;
							NI=*(ptrSize+0);
							NJ=*(ptrSize+1);
							NK=*(ptrSize+2);
							
							Epsilon = zone->Flowsol->Scal[1];
							
							flagRho = zone->Flowsol->flagRho;
							flagViscoeff = zone->Flowsol->flagViscoeff;
							rho = *(zone->Flowsol->Rho);
							mu  = *(zone->Flowsol->Viscoeff);
							upperlimit = 900.0 * rho*kmax*kmax/mu; //0.09 * rho*kmax*kmax/(1.0e-4*mu);
							
							index = 0;
							for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
							{
								//Calculating the upper limit for epsilon
								if(flagRho == 1 || flagViscoeff == 1)
								{
									rho = *(zone->Flowsol->Rho+index*flagRho);
									mu = *(zone->Flowsol->Viscoeff+index*flagViscoeff);
									upperlimit = 900.0 * rho*kmax*kmax/mu; //0.09 * kmax*kmax/(1.0e-4*mu);
								}
								//Limiting turbulent kinetic energy using lower bound value
								if(*(Epsilon+index) < lowerlimit)
								{
									*(Epsilon+index) = lowerlimit;//fabs(*(Epsilon+index));//lowerlimit;
									LowerLimitCount++;
								}
								//Limiting turbulent kinetic energy using upper bound value
								if(*(Epsilon+index) > upperlimit)
								{
									*(Epsilon+index) = upperlimit;
									UpperLimitCount++;
								}
								index++;
							}
			
							for(count=0; count < zone->Number_of_BC_Patches; count++)
							{
								bcPtr   = zone->ZoneBc+count;
								Epsilon = bcPtr->ZSlab->Scal[1];
								
								totalCells = bcPtr->ZSlab->Total_Cells;
								
								for(index=0; index < totalCells; index++)
								{
									//Calculating the upper limit for epsilon
									if(flagRho == 1 || flagViscoeff == 1)
									{
										rho = *(bcPtr->ZSlab->Rho+index*flagRho);
										mu  = *(bcPtr->ZSlab->Viscoeff+index*flagViscoeff);
										upperlimit = 900.0 * rho*kmax*kmax/mu; //0.09 * kmax*kmax/(1.0e-4*mu);
									}
									//Limiting turbulent kinetic energy using lower bound value
									if(*(Epsilon+index) < lowerlimit)
									{
										*(Epsilon+index) = lowerlimit;//fabs(*(Epsilon+index));//lowerlimit;
										LowerLimitSlabCount++;
									}
									//Limiting turbulent kinetic energy using upper bound value
									if(*(Epsilon+index) > upperlimit)
									{
										*(Epsilon+index) = upperlimit;
										UpperLimitSlabCount++;
									}
								}
							}
							
							for(count=0; count < zone->Number_of_Interfaces; count++)
							{
								IfacePtr = zone->Iface+count;
								if(IfacePtr->ConjugateIfaceTag == 1)
								{
									Epsilon = IfacePtr->ZSlab->Scal[1];
									
									totalCells = IfacePtr->ZSlab->Total_Cells;
									
									for(index=0; index < totalCells; index++)
									{
										//Calculating the upper limit for epsilon
										if(flagRho == 1 || flagViscoeff == 1)
										{
											rho = *(IfacePtr->ZSlab->Rho+index*flagRho);
											mu  = *(IfacePtr->ZSlab->Viscoeff+index*flagViscoeff);
											upperlimit = 900.0 * rho*kmax*kmax/mu; //0.09 * kmax*kmax/(1.0e-4*mu);
										}
										//Limiting turbulent kinetic energy using lower bound value
										if(*(Epsilon+index) < lowerlimit)
										{
											*(Epsilon+index) = lowerlimit;//fabs(*(Epsilon+index));//lowerlimit;
											LowerLimitSlabCount++;
										}
										//Limiting turbulent kinetic energy using upper bound value
										if(*(Epsilon+index) > upperlimit)
										{
											*(Epsilon+index) = upperlimit;
											UpperLimitSlabCount++;
										}
									}
								}
							}
						}
					}
					if(UpperLimitCount!=0 || UpperLimitSlabCount!=0)
						printf(">>> Max. Epsilon is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", upperlimit,UpperLimitCount,UpperLimitSlabCount);
					if(LowerLimitCount!=0 || LowerLimitSlabCount!=0)
						printf(">>> Min. Epsilon is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", lowerlimit,LowerLimitCount,LowerLimitSlabCount);
					break;
				default:
					printf(">>> No limiter is applied to variable [%d] ...\n", vartag);
					return 0;
					break;
			}
		}
	}
	
	return 0;
}

//------------------------------------------------------
int LimitVariableZone(Base *base, Zone *zone, int vartag)
{
	int i,j,k,NI,NJ,NK,count,*ptrSize,*ptrRange,index,NIcell,NJcell,NKcell;
	int LowerLimitSlabCount=0,UpperLimitSlabCount=0,LowerLimitCount=0,UpperLimitCount=0,totalCells;
	double upperlimit,lowerlimit,Uref=100.0,*KE,*Epsilon,kmax;
	Interface *IfacePtr;
	ZoneBoundary *bcPtr;
	
	//NOTE: "Uref" should be obtained from user input.
	Uref=5.0;
	kmax = 0.2*Uref*Uref/2.0;
	
	if(base->flagTurbModel == 1 || base->flagTurbModel == 3)
	{
		switch(vartag)
		{
			case 600:
				lowerlimit = 1.0e-10;// * Uref*Uref;//1.0e-10;
				upperlimit = kmax;
				
				if(zone->flagSolidFluid==0)
				{
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);
					
					KE = zone->Flowsol->Scal[0];
					
					index = 0;
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
					{
						//Limiting turbulent kinetic energy using lower bound value
						if(*(KE+index) < 0.0)
						{
							*(KE+index) = lowerlimit;
							LowerLimitCount++;
						}
						//Limiting turbulent kinetic energy using upper bound value
						if(*(KE+index) > upperlimit)
						{
							*(KE+index) = upperlimit;
							UpperLimitCount++;
						}
						index++;
					}
	
					for(count=0; count < zone->Number_of_BC_Patches; count++)
					{
						bcPtr = zone->ZoneBc+count;
						KE    = bcPtr->ZSlab->Scal[0];
						
						ptrRange = bcPtr->Point_Range;
						NIcell = *(ptrRange+3) - *(ptrRange+0);
						NJcell = *(ptrRange+4) - *(ptrRange+1);
						NKcell = *(ptrRange+5) - *(ptrRange+2);
						
						switch(bcPtr->Direction)
						{
							case 1://Imax face
								case 2://Imin face
									totalCells = NJcell * NKcell;
									break;
									case 3://Jmax face
										case 4://Jmin face
										totalCells = NIcell * NKcell;
										break;
										case 5://Kmax face
										case 6://Kmin face
										totalCells = NIcell * NJcell;
										break;
						}
						
						for(index=0; index < totalCells; index++)
						{
							//Limiting turbulent kinetic energy using lower bound value
							if(*(KE+index) < 0.0)
							{
								*(KE+index) = lowerlimit;
								LowerLimitSlabCount++;
							}
							//Limiting turbulent kinetic energy using upper bound value
							if(*(KE+index) > upperlimit)
							{
								*(KE+index) = upperlimit;
								UpperLimitSlabCount++;
							}
						}
					}
					
					for(count=0; count < zone->Number_of_Interfaces; count++)
					{
						IfacePtr = zone->Iface+count;
						//if(IfacePtr->ConjugateIfaceTag == 1)
						{
							KE    =  IfacePtr->ZSlab->Scal[0];
							
							ptrRange = IfacePtr->Point_Range;
							NIcell = *(ptrRange+3) - *(ptrRange+0);
							NJcell = *(ptrRange+4) - *(ptrRange+1);
							NKcell = *(ptrRange+5) - *(ptrRange+2);
							
							switch(IfacePtr->Direction)
							{
								case 1://Imax face
									case 2://Imin face
										totalCells = NJcell * NKcell;
										break;
										case 3://Jmax face
										case 4://Jmin face
										totalCells = NIcell * NKcell;
										break;
										case 5://Kmax face
										case 6://Kmin face
										totalCells = NIcell * NJcell;
										break;
							}
							
							for(index=0; index < totalCells; index++)
							{
								//Limiting turbulent kinetic energy using lower bound value
								if(*(KE+index) < 0.0)
								{
									*(KE+index) = lowerlimit;
									LowerLimitSlabCount++;
								}
								//Limiting turbulent kinetic energy using upper bound value
								if(*(KE+index) > upperlimit)
								{
									*(KE+index) = upperlimit;
									UpperLimitSlabCount++;
								}
							}
					
						}
					}
				}
				if(UpperLimitCount!=0 || UpperLimitSlabCount!=0)
					printf(">>> Zone:%s Max. TKE is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", zone->zonename,upperlimit,UpperLimitCount,UpperLimitSlabCount);
				if(LowerLimitCount!=0 || LowerLimitSlabCount!=0)
					printf(">>> Zone:%s Min. TKE is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", zone->zonename,lowerlimit,LowerLimitCount,LowerLimitSlabCount);
				break;
			case 601:
				lowerlimit = 1.0e-12;// * Uref*Uref*Uref*Uref/(1.0e-4/780.0);//1.0e-12;
				upperlimit = 1.0e20;
				
				if(zone->flagSolidFluid==0)
				{
					ptrSize = zone->Size;
					NI=*(ptrSize+0);
					NJ=*(ptrSize+1);
					NK=*(ptrSize+2);
					
					Epsilon = zone->Flowsol->Scal[1];
					
					index = 0;
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
					{
						//Limiting turbulent kinetic energy using lower bound value
						if(*(Epsilon+index) < 0.0)
						{
							*(Epsilon+index) = lowerlimit;
							LowerLimitCount++;
						}
						//Limiting turbulent kinetic energy using upper bound value
						if(*(Epsilon+index) > upperlimit)
						{
							*(Epsilon+index) = upperlimit;
							UpperLimitCount++;
						}
						index++;
					}
	
					for(count=0; count < zone->Number_of_BC_Patches; count++)
					{
						bcPtr   = zone->ZoneBc+count;
						Epsilon = bcPtr->ZSlab->Scal[1];
						
						ptrRange = bcPtr->Point_Range;
						NIcell = *(ptrRange+3) - *(ptrRange+0);
						NJcell = *(ptrRange+4) - *(ptrRange+1);
						NKcell = *(ptrRange+5) - *(ptrRange+2);
						
						switch(bcPtr->Direction)
						{
							case 1://Imax face
								case 2://Imin face
									totalCells = NJcell * NKcell;
									break;
									case 3://Jmax face
										case 4://Jmin face
										totalCells = NIcell * NKcell;
										break;
										case 5://Kmax face
										case 6://Kmin face
										totalCells = NIcell * NJcell;
										break;
						}
						
						for(index=0; index < totalCells; index++)
						{
							//Limiting turbulent kinetic energy using lower bound value
							if(*(Epsilon+index) < 0.0)
							{
								*(Epsilon+index) = lowerlimit;
								LowerLimitSlabCount++;
							}
							//Limiting turbulent kinetic energy using upper bound value
							if(*(Epsilon+index) > upperlimit)
							{
								*(Epsilon+index) = upperlimit;
								UpperLimitSlabCount++;
							}
						}
					}
					
					for(count=0; count < zone->Number_of_Interfaces; count++)
					{
						IfacePtr = zone->Iface+count;
						//if(IfacePtr->ConjugateIfaceTag == 1)
						{
							Epsilon = IfacePtr->ZSlab->Scal[1];
							
							ptrRange = IfacePtr->Point_Range;
							NIcell = *(ptrRange+3) - *(ptrRange+0);
							NJcell = *(ptrRange+4) - *(ptrRange+1);
							NKcell = *(ptrRange+5) - *(ptrRange+2);
							
							switch(IfacePtr->Direction)
							{
								case 1://Imax face
									case 2://Imin face
										totalCells = NJcell * NKcell;
										break;
										case 3://Jmax face
										case 4://Jmin face
										totalCells = NIcell * NKcell;
										break;
										case 5://Kmax face
										case 6://Kmin face
										totalCells = NIcell * NJcell;
										break;
							}
							
							for(index=0; index < totalCells; index++)
							{
								//Limiting turbulent kinetic energy using lower bound value
								if(*(Epsilon+index) < 0.0)
								{
									*(Epsilon+index) = lowerlimit;
									LowerLimitSlabCount++;
								}
								//Limiting turbulent kinetic energy using upper bound value
								if(*(Epsilon+index) > upperlimit)
								{
									*(Epsilon+index) = upperlimit;
									UpperLimitSlabCount++;
								}
							}
						}
					}
				}
				if(UpperLimitCount!=0 || UpperLimitSlabCount!=0)
					printf(">>> Zone:%s Max. Epsilon is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", zone->zonename,upperlimit,UpperLimitCount,UpperLimitSlabCount);
				if(LowerLimitCount!=0 || LowerLimitSlabCount!=0)
					printf(">>> Zone:%s Min. Epsilon is limited to %e in %d Cell(s) and %d Slab(s) <<<\n", zone->zonename,lowerlimit,LowerLimitCount,LowerLimitSlabCount);
				break;
			default:
				printf(">>> No limiter is applied to variable [%d] ...\n", vartag);
				return 0;
				break;
		}
	}
	return 0;
}

