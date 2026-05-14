#include "APKprintslabs.h"
#define UREF 1.0
#define H 0.5
//H is half channel width
/**
This function, printslabs() prints the slab values(both at boundaries and interfaces) of a variable into a data file just to check
the values. If required, at a wall, 'Cfx', 'u_tau','k' etc. can be computed.
*/
int printslab(Base *base)
{
	int numzones,cntZ,domainGridIndexBottomLeft,domainGridIndexTopRight;
	int i,j,k,numIface,numZbc,dir,iface,zBc;
	int index,domainCellIndex,NI,NJ,NK;
	int *ptrRange,indices[3],flagRho=0,flagViscoeff=0;
	double Cf,uTau,TKE,Epsilon,theta,delta_star,delta,rho,mu;
	double Uref=UREF, X, Y, Z,Re_theta,Re_delta_star,Re_delta,Re_tau,Re_x,Uavg=UREF;
	Zone *zone;
	ZoneBoundary *ZBc;
	numzones= base->Number_of_zones;
	FILE *prslab,*prslab2;


	prslab= fopen("SOLUTIONS/printslab.dat","w");
	prslab2= fopen("SOLUTIONS/slabCfKEpUtau.dat","w");
	//printf("I am in printslab\n");
	if(base->BaseParameters->Turb->numTurbScalars >=2)
	{
		for(cntZ=0;cntZ<numzones;cntZ++)
		{
			zone = base->CgnsZone+cntZ;
			flagViscoeff = zone->Flowsol->flagViscoeff;
			flagRho = zone->Flowsol->flagRho;
			NI=*(zone->Size+0);
			NJ=*(zone->Size+1);
			NK=*(zone->Size+2);
			if(zone->flagSolidFluid == 0)
			{
				fprintf(prslab,"##### ZONE NAME :%s ########\n",zone->zonename);
				numIface=zone->Number_of_Interfaces;
				numZbc=zone->Number_of_BC_Patches;
				//printf("##### ZONE NAME :%s ########\n",zone->zonename);
				//printf("numIface=%d  numZBc=%d\n",numIface,numZbc);
				for(dir=1;dir<=6;dir++)
				{
					//printf("dir=%d\n",dir);
					for(iface=0;iface<numIface;iface++)
					{
						//printf("%d--\n",dir);
						if(((zone->Iface+iface)->Direction)==dir)
						{
							//printf("%d-\n",dir);
							ptrRange=(zone->Iface+iface)->Point_Range;

							fprintf(prslab,"Iface dir->%d\n",dir);
							fprintf(prslab,"~~~~~~~~~Interface Name: %s ~~~~~~~~~~\n",(((base)->CgnsZone+cntZ)->Iface+iface)->Interface_Name);
							//printf("Iface dir->%d\n",dir);
							//printf("~~~~~~~~~Interface Name: %s ~~~~~~~~~~\n",(((base)->CgnsZone+cntZ)->Iface+iface)->Interface_Name);

							if(dir==1||dir==2)
							{
								//printf("i m in direction 1\n");
								//printf("i ve crossed the flag\n");
								for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
								for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
								{
									index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
									TKE = *((zone->Iface+iface)->ZSlab->Scal[0]+index);
									Epsilon = *((zone->Iface+iface)->ZSlab->Scal[1]+index);

									if((base->flagLogKEps) == 1)
									{
										TKE = exp(TKE);
										Epsilon = exp(Epsilon);
									}

									fprintf(prslab,"\t%e\t%e\n",TKE,Epsilon);
									//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*((zone->Iface+iface)->ZSlab->P+index),*((zone->Iface+iface)->ZSlab->v+index),*((zone->Iface+iface)->ZSlab->w+index));
									//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,volume=%12.8lf\n",*((zone->Iface+iface)->ZSlab->T+index),i,j,k,*((zone->Iface+iface)->ZSlab->Volume+index));
								}
							}
							if(dir==3||dir==4)
							{
								//printf("i m in direction 3\n");
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									TKE = *((zone->Iface+iface)->ZSlab->Scal[0]+index);
									Epsilon = *((zone->Iface+iface)->ZSlab->Scal[1]+index);

									if((base->flagLogKEps) == 1)
									{
										TKE = exp(TKE);
										Epsilon = exp(Epsilon);
									}

									fprintf(prslab,"\t%e\t%e\n",TKE,Epsilon);

									//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*((zone->Iface+iface)->ZSlab->P+index),*((zone->Iface+iface)->ZSlab->v+index),*((zone->Iface+iface)->ZSlab->w+index));
									//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,volume=%12.8lf\n",*((zone->Iface+iface)->ZSlab->T+index),i,j,k,*((zone->Iface+iface)->ZSlab->Volume+index));
								}	
							}
							if(dir==5||dir==6)
							{
								//printf("i m in direction 5\n");
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									TKE = *((zone->Iface+iface)->ZSlab->Scal[0]+index);
									Epsilon = *((zone->Iface+iface)->ZSlab->Scal[1]+index);

									if((base->flagLogKEps) == 1)
									{
										TKE = exp(TKE);
										Epsilon = exp(Epsilon);
									}

									fprintf(prslab,"\t%e\t%e\n",TKE,Epsilon);

									//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*((zone->Iface+iface)->ZSlab->P+index),*((zone->Iface+iface)->ZSlab->v+index),*((zone->Iface+iface)->ZSlab->w+index));
									//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,volume=%12.8lf\n",*((zone->Iface+iface)->ZSlab->T+index),i,j,k,*((zone->Iface+iface)->ZSlab->Volume+index));
								}
							}
						}
					}
					for(zBc=0;zBc<numZbc;zBc++)
					{
						ZBc = zone->ZoneBc+zBc;
						if(ZBc->Direction==dir)
						{
							//printf("boundary dir=%d\n",dir);
							ptrRange=ZBc->Point_Range;
							//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
							fprintf(prslab,"ZoneBc dir->%d\n",dir);
							fprintf(prslab,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
							//printf("ZoneBc dir->%d\n",dir);
							//printf("~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
							if(dir==1||dir==2)
							{
								//printf("inside if\n");
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								{
									//printf("k=%d,j=%d\n",k,j);
									index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
									//printf("INDEX=%d\n",index);
									TKE = *(ZBc->ZSlab->Scal[0]+index);
									Epsilon = *(ZBc->ZSlab->Scal[1]+index);

									if((base->flagLogKEps) == 1)
									{
										TKE = exp(TKE);
										Epsilon = exp(Epsilon);
									}

									fprintf(prslab,"\t%e\t%e\n",TKE,Epsilon);
									//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->P+index),*(ZBc->ZSlab->v+index),*(ZBc->ZSlab->w+index));
									//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
								}
							}
							if(dir==3||dir==4)
							{

								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									TKE = *(ZBc->ZSlab->Scal[0]+index);
									Epsilon = *(ZBc->ZSlab->Scal[1]+index);

									if((base->flagLogKEps) == 1)
									{
										TKE = exp(TKE);
										Epsilon = exp(Epsilon);
									}

									fprintf(prslab,"\t%e\t%e\n",TKE,Epsilon);

									//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->P+index),*(ZBc->ZSlab->v+index),*(ZBc->ZSlab->w+index));
									//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
								}
							}
							if(dir==5||dir==6)
							{
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									TKE = *(ZBc->ZSlab->Scal[0]+index);
									Epsilon = *(ZBc->ZSlab->Scal[1]+index);

									if((base->flagLogKEps) == 1)
									{
										TKE = exp(TKE);
										Epsilon = exp(Epsilon);
									}

									fprintf(prslab,"\t%e\t%e\n",TKE,Epsilon);
									//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->P+index),*(ZBc->ZSlab->v+index),*(ZBc->ZSlab->w+index));
									//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
								}
							}
						}

						if(ZBc->Boundary_Tag==22)
						{
							if(base->flagTurbModel == 1 || base->flagTurbModel == 3)
							{
								if(ZBc->Direction==dir)
								{
									//printf("boundary dir=%d\n",dir);
									ptrRange=ZBc->Point_Range;
									//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
									fprintf(prslab2,"ZoneBc dir->%d\n",dir);
									fprintf(prslab2,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
									fprintf(prslab2,"\"Cfx\"\t \"k\" \t \"Epsilon\" \t  \"UTau\" \t\"index\"\n");
									//printf("ZoneBc dir->%d\n",dir);
									//printf("~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
									if(dir==1||dir==2)
									{
										//printf("inside if\n");
										for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
										for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
										{
											//printf("k=%d,j=%d\n",k,j);
											index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

											Cf = *(ZBc->ZSlab->Scal[0]+index)*0.6;
											TKE = *(ZBc->ZSlab->Scal[0]+index);
											Epsilon = *(ZBc->ZSlab->Scal[1]+index);
											uTau = *(ZBc->uTau+index);

											if((base->flagLogKEps) == 1)
											{
												TKE = exp(TKE);
												Epsilon = exp(Epsilon);
												Cf = 0.6*TKE;
											}
											//printf("INDEX=%d\n",index);
											//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(j-1)+(j-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(j-1)+(j-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(j-1)+(j-1)));
											fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%d\n",Cf,TKE,Epsilon,uTau,j);
											//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
											//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
										}
									}
									if(dir==3||dir==4)
									{

										for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
										for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
										{
											index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

											Cf = *(ZBc->ZSlab->Scal[0]+index)*0.6;
											TKE = *(ZBc->ZSlab->Scal[0]+index);
											Epsilon = *(ZBc->ZSlab->Scal[1]+index);
											uTau = *(ZBc->uTau+index);

											if((base->flagLogKEps) == 1)
											{
												TKE = exp(TKE);
												Epsilon = exp(Epsilon);
												Cf = 0.6*TKE;
											}

											//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));
											fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%d\n",Cf,TKE,Epsilon,uTau,i);
											//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
											//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
										}
									}
									if(dir==5||dir==6)
									{
										for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
										for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
										{
											index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

											Cf = *(ZBc->ZSlab->Scal[0]+index)*0.6;
											TKE = *(ZBc->ZSlab->Scal[0]+index);
											Epsilon = *(ZBc->ZSlab->Scal[1]+index);
											uTau = *(ZBc->uTau+index);


											if((base->flagLogKEps) == 1)
											{
												TKE = exp(TKE);
												Epsilon = exp(Epsilon);
												Cf = 0.6*TKE;
											}

											//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(j-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(j-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(j-1)*(i-1)+(i-1)));
											fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%d\n",Cf,TKE,Epsilon,uTau,i);
											//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
											//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
										}
									}
								}
							}
							else if(base->flagTurbModel == 2 || base->flagTurbModel == 6)
							{
								if(ZBc->Direction==dir)
								{
									//printf("boundary dir=%d\n",dir);
									ptrRange=ZBc->Point_Range;
									//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
									fprintf(prslab2,"ZoneBc dir->%d\n",dir);
									fprintf(prslab2,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\nU_ref=%lf\n",ZBc->BCName,Uref);
									fprintf(prslab2,"variables=\n\t\"X\"\t\t\"Y\"\t\t\"Z\"\t  \"Re_theta\"\t \"Re_delta*\" \t \"Re_delta\" \t  \"theta\"\t    \"delta*\" \t   \"delta\" \t   \"Cfx\"\t \"k\" \t\t  \"Omega\" \t  \"UTau\"\t\"index\" \"Re_tau\"\t\"Re_H\"\t\"Uavg\"\n");
									//printf("ZoneBc dir->%d\n",dir);
									//printf("~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
									if(dir==1||dir==2)
									{
										if(dir == 1)
											i=NI-1;
										else
											i=1;

										//printf("inside if\n");
										for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
										for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
										{
											//printf("k=%d,j=%d\n",k,j);
											index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

											TKE = *(ZBc->ZSlab->Scal[0]+index);
											Epsilon = *(ZBc->ZSlab->Scal[1]+index);
											uTau = *(ZBc->uTau+index);
											Cf = 2 * (uTau/Uref) * (uTau/Uref);
           										domainCellIndex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
											domainGridIndexBottomLeft = NI*NJ*(k-1)+NI*(j-1)+(i-1);
											domainGridIndexTopRight = NI*NJ*(k)+NI*(j)+(i-1);
											X = 0.5*(*(zone->CgnsGrid->X+domainGridIndexBottomLeft) + *(zone->CgnsGrid->X+domainGridIndexTopRight));
											Y = 0.5*(*(zone->CgnsGrid->Y+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Y+domainGridIndexTopRight));
											Z = 0.5*(*(zone->CgnsGrid->Z+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Z+domainGridIndexTopRight));
           										indices[0]=i-1;
											indices[1]=j-1;
											indices[2]=k-1;
											rho = *(ZBc->ZSlab->Rho+index*flagRho);
											mu  = *(ZBc->ZSlab->Viscoeff+index*flagViscoeff);

											//The function FindUavg() is written for turbulent channel flow case ...
											Uavg = FindUavg(zone,ZBc,indices);
											Cf = 2 * (uTau/Uavg) * (uTau/Uavg);
											Re_x = rho*Uavg*(2*H)/mu;
											//Re_x = rho*Uref*Y/mu;

											FindReTheta(zone,ZBc,indices, &theta, &delta_star, &delta);
											//printf("theta=%e\n",theta);
											Re_theta	= rho*Uref*theta/mu;
											Re_delta_star	= rho*Uref*delta_star/mu;
											Re_delta	= rho*Uref*delta/mu;
											Re_tau		= rho*uTau*H/mu;//not generalised and used for Channel flow configuration : H->Channel Half Width

											fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\t%12.8e\t%12.8lf\t%d\t%lf\t%e\t%12.8lf\n", X,Y,Z,Re_theta,Re_delta_star,Re_delta,theta,delta_star,delta,Cf,TKE,Epsilon,uTau,index,Re_tau,Re_x,Uavg);
											//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));

											//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
											//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
										}
									}
									if(dir==3||dir==4)
									{
										if(dir == 3)
											j=NJ-1;
										else
											j=1;

										for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
										for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
										{
											index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

											TKE = *(ZBc->ZSlab->Scal[0]+index);
											Epsilon = *(ZBc->ZSlab->Scal[1]+index);
											uTau = *(ZBc->uTau+index);
											Cf = 2 * (uTau/Uref) * (uTau/Uref);
											domainCellIndex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
											domainGridIndexBottomLeft = NI*NJ*(k-1)+NI*(j-1)+(i-1);
											domainGridIndexTopRight = NI*NJ*(k)+NI*(j-1)+(i);
											X = 0.5*(*(zone->CgnsGrid->X+domainGridIndexBottomLeft) + *(zone->CgnsGrid->X+domainGridIndexTopRight));
											Y = 0.5*(*(zone->CgnsGrid->Y+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Y+domainGridIndexTopRight));
											Z = 0.5*(*(zone->CgnsGrid->Z+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Z+domainGridIndexTopRight));
           										indices[0]=i-1;
											indices[1]=j-1;
											indices[2]=k-1;
											rho = *(ZBc->ZSlab->Rho+index*flagRho);
											mu  = *(ZBc->ZSlab->Viscoeff+index*flagViscoeff);

											//The function FindUavg() is written for turbulent channel flow case ...
											Uavg = FindUavg(zone,ZBc,indices);
											Cf = 2 * (uTau/Uavg) * (uTau/Uavg);
											Re_x = rho*Uavg*(2*H)/mu;
											//Re_x = rho*Uref*X/mu;

											FindReTheta(zone,ZBc,indices, &theta, &delta_star, &delta);
											Re_theta	= rho*Uref*theta/mu;
											Re_delta_star	= rho*Uref*delta_star/mu;
											Re_delta	= rho*Uref*delta/mu;
											Re_tau		= rho*uTau*H/mu;//not generalised and used for Channel flow configuration : H->Channel Half Width

											fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\t%12.8e\t%12.8lf\t%d\t%lf\t%e\t%12.8lf\n", X,Y,Z,Re_theta,Re_delta_star,Re_delta,theta,delta_star,delta,Cf,TKE,Epsilon,uTau,index,Re_tau,Re_x,Uavg);
											//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));
           										//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
											//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
										}
									}
									if(dir==5||dir==6)
									{
										if(dir == 5)
											k=NK-1;
										else
											k=1;
										for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
										for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
										{
											index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

											TKE = *(ZBc->ZSlab->Scal[0]+index);
											Epsilon = *(ZBc->ZSlab->Scal[1]+index);
											uTau = *(ZBc->uTau+index);
											Cf = 2 * (uTau/Uref) * (uTau/Uref);
											domainCellIndex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
											domainGridIndexBottomLeft = NI*NJ*(k-1)+NI*(j-1)+(i-1);
											domainGridIndexTopRight = NI*NJ*(k-1)+NI*(j)+(i);
											X = 0.5*(*(zone->CgnsGrid->X+domainGridIndexBottomLeft) + *(zone->CgnsGrid->X+domainGridIndexTopRight));
											Y = 0.5*(*(zone->CgnsGrid->Y+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Y+domainGridIndexTopRight));
											Z = 0.5*(*(zone->CgnsGrid->Z+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Z+domainGridIndexTopRight));
           										indices[0]=i-1;
											indices[1]=j-1;
											indices[2]=k-1;
											rho = *(ZBc->ZSlab->Rho+index*flagRho);
											mu  = *(ZBc->ZSlab->Viscoeff+index*flagViscoeff);

											//The function FindUavg() is written for turbulent channel flow case ...
											Uavg = FindUavg(zone,ZBc,indices);
											Cf = 2 * (uTau/Uavg) * (uTau/Uavg);
											Re_x = rho*Uavg*(2*H)/mu;
											//Re_x = rho*Uref*Z/mu;//verify this for a problem under consideration!!

											FindReTheta(zone,ZBc,indices, &theta, &delta_star, &delta);
											Re_theta	= rho*Uref*theta/mu;
											Re_delta_star	= rho*Uref*delta_star/mu;
											Re_delta	= rho*Uref*delta/mu;
											Re_tau		= rho*uTau*H/mu;//not generalised and used for Channel flow configuration : H->Channel Half Width

											fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\t%12.8e\t%12.8lf\t%d\t%lf\t%e\t%12.8lf\n", X,Y,Z,Re_theta,Re_delta_star,Re_delta,theta,delta_star,delta,Cf,TKE,Epsilon,uTau,index,Re_tau,Re_x,Uavg);
											//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));
           										//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
											//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
										}
									}
								}
							}

						}
					}
				}
			}
		}
	}
	else if(base->flagTurb==1)
	{
		for(cntZ=0;cntZ<numzones;cntZ++)
		{
			zone = base->CgnsZone+cntZ;
			flagViscoeff = zone->Flowsol->flagViscoeff;
			flagRho = zone->Flowsol->flagRho;
			NI=*(zone->Size+0);
			NJ=*(zone->Size+1);
			NK=*(zone->Size+2);
			if(zone->flagSolidFluid == 0)
			{
				fprintf(prslab,"##### ZONE NAME :%s ########\n",zone->zonename);
				numIface=zone->Number_of_Interfaces;
				numZbc=zone->Number_of_BC_Patches;
				
				for(zBc=0;zBc<numZbc;zBc++)
				{
					ZBc = zone->ZoneBc+zBc;

					if(ZBc->Boundary_Tag==22)
					{
						dir = ZBc->Direction;
						ptrRange=ZBc->Point_Range;
						//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
						fprintf(prslab2,"ZoneBc dir->%d\n",dir);
						fprintf(prslab2,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\nU_ref=%lf\n",ZBc->BCName,Uref);
						fprintf(prslab2,"variables=\n\t\"X\"\t\t\"Y\"\t\t\"Z\"\t  \"Re_theta\"\t \"Re_delta*\" \t \"Re_delta\" \t  \"theta\"\t    \"delta*\" \t   \"delta\" \t   \"Cfx\"\t  \"UTau\"\t\"index\" \"Re_tau\"\t\"Re_H\"\t\"Uavg\"\n");
						//printf("ZoneBc dir->%d\n",dir);
						//printf("~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
						if(dir==1||dir==2)
						{
							if(dir == 1)
								i=NI-1;
							else
								i=1;

							//printf("inside if\n");
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							{
								//printf("k=%d,j=%d\n",k,j);
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
								
								uTau = *(ZBc->uTau+index);
								Cf = 2 * (uTau/Uref) * (uTau/Uref);
								domainCellIndex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
								domainGridIndexBottomLeft = NI*NJ*(k-1)+NI*(j-1)+(i-1);
								domainGridIndexTopRight = NI*NJ*(k)+NI*(j)+(i-1);
								X = 0.5*(*(zone->CgnsGrid->X+domainGridIndexBottomLeft) + *(zone->CgnsGrid->X+domainGridIndexTopRight));
								Y = 0.5*(*(zone->CgnsGrid->Y+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Y+domainGridIndexTopRight));
								Z = 0.5*(*(zone->CgnsGrid->Z+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Z+domainGridIndexTopRight));
								indices[0]=i-1;
								indices[1]=j-1;
								indices[2]=k-1;
								rho = *(ZBc->ZSlab->Rho+index*flagRho);
								mu  = *(ZBc->ZSlab->Viscoeff+index*flagViscoeff);

								//The function FindUavg() is written for turbulent channel flow case ...
								Uavg = FindUavg(zone,ZBc,indices);
								Cf = 2 * (uTau/Uavg) * (uTau/Uavg);
								Re_x = rho*Uavg*(2*H)/mu;
								//Re_x = rho*Uref*Y/mu;

								FindReTheta(zone,ZBc,indices, &theta, &delta_star, &delta);
								//printf("theta=%e\n",theta);
								Re_theta	= rho*Uref*theta/mu;
								Re_delta_star	= rho*Uref*delta_star/mu;
								Re_delta	= rho*Uref*delta/mu;
								Re_tau		= rho*uTau*H/mu;//not generalised and used for Channel flow configuration : H->Channel Half Width

								fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%d\t%lf\t%lf\t%lf\n", X,Y,Z,Re_theta,Re_delta_star,Re_delta,theta,delta_star,delta,Cf,uTau,index,Re_tau,Re_x,Uavg);
								//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));

								//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
								//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
							}
						}
						if(dir==3||dir==4)
						{
							if(dir == 3)
								j=NJ-1;
							else
								j=1;

							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								
								uTau = *(ZBc->uTau+index);
								Cf = 2 * (uTau/Uref) * (uTau/Uref);
								domainCellIndex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
								domainGridIndexBottomLeft = NI*NJ*(k-1)+NI*(j-1)+(i-1);
								domainGridIndexTopRight = NI*NJ*(k)+NI*(j-1)+(i);
								X = 0.5*(*(zone->CgnsGrid->X+domainGridIndexBottomLeft) + *(zone->CgnsGrid->X+domainGridIndexTopRight));
								Y = 0.5*(*(zone->CgnsGrid->Y+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Y+domainGridIndexTopRight));
								Z = 0.5*(*(zone->CgnsGrid->Z+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Z+domainGridIndexTopRight));
								indices[0]=i-1;
								indices[1]=j-1;
								indices[2]=k-1;
								rho = *(ZBc->ZSlab->Rho+index*flagRho);
								mu  = *(ZBc->ZSlab->Viscoeff+index*flagViscoeff);

								//The function FindUavg() is written for turbulent channel flow case ...
								Uavg = FindUavg(zone,ZBc,indices);
								Cf = 2 * (uTau/Uavg) * (uTau/Uavg);
								Re_x = rho*Uavg*(2*H)/mu;
								//Re_x = rho*Uref*X/mu;

								FindReTheta(zone,ZBc,indices, &theta, &delta_star, &delta);
								Re_theta	= rho*Uref*theta/mu;
								Re_delta_star	= rho*Uref*delta_star/mu;
								Re_delta	= rho*Uref*delta/mu;
								Re_tau		= rho*uTau*H/mu;//not generalised and used for Channel flow configuration : H->Channel Half Width

								fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%d\t%lf\t%lf\t%lf\n", X,Y,Z,Re_theta,Re_delta_star,Re_delta,theta,delta_star,delta,Cf,uTau,index,Re_tau,Re_x,Uavg);
								//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));
								//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
								//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
							}
						}
						if(dir==5||dir==6)
						{
							if(dir == 5)
								k=NK-1;
							else
								k=1;
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								
								uTau = *(ZBc->uTau+index);
								Cf = 2 * (uTau/Uref) * (uTau/Uref);
								domainCellIndex = (NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+(i-1);
								domainGridIndexBottomLeft = NI*NJ*(k-1)+NI*(j-1)+(i-1);
								domainGridIndexTopRight = NI*NJ*(k-1)+NI*(j)+(i);
								X = 0.5*(*(zone->CgnsGrid->X+domainGridIndexBottomLeft) + *(zone->CgnsGrid->X+domainGridIndexTopRight));
								Y = 0.5*(*(zone->CgnsGrid->Y+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Y+domainGridIndexTopRight));
								Z = 0.5*(*(zone->CgnsGrid->Z+domainGridIndexBottomLeft) + *(zone->CgnsGrid->Z+domainGridIndexTopRight));
								indices[0]=i-1;
								indices[1]=j-1;
								indices[2]=k-1;
								rho = *(ZBc->ZSlab->Rho+index*flagRho);
								mu  = *(ZBc->ZSlab->Viscoeff+index*flagViscoeff);

								//The function FindUavg() is written for turbulent channel flow case ...
								Uavg = FindUavg(zone,ZBc,indices);
								Cf = 2 * (uTau/Uavg) * (uTau/Uavg);
								Re_x = rho*Uavg*(2*H)/mu;
								//Re_x = rho*Uref*Z/mu;//verify this for a problem under consideration!!

								FindReTheta(zone,ZBc,indices, &theta, &delta_star, &delta);
								Re_theta	= rho*Uref*theta/mu;
								Re_delta_star	= rho*Uref*delta_star/mu;
								Re_delta	= rho*Uref*delta/mu;
								Re_tau		= rho*uTau*H/mu;//not generalised and used for Channel flow configuration : H->Channel Half Width

								fprintf(prslab2,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%d\t%lf\t%lf\t%lf\n", X,Y,Z,Re_theta,Re_delta_star,Re_delta,theta,delta_star,delta,Cf,uTau,index,Re_tau,Re_x,Uavg);
								//fprintf(prslab,"\tu=%12.8lf\tv=%12.8lf\tw=%12.8lf\n",*(ZBc->ZSlab->u+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->v+(k-1)*(i-1)+(i-1)),*((zone->Iface+iface)->ZSlab->w+(k-1)*(i-1)+(i-1)));
								//fprintf(prslab,"%12.8lf\t %12.8lf\n",*(ZBc->uTau+index),*(ZBc->deltaP+index));
								//fprintf(prslab,"\tu=%12.8lf,\t%d,\t%d,\t%d,\tvolume=%12.8lf\n",*(ZBc->ZSlab->T+index),i,j,k,*(ZBc->ZSlab->Volume+index));
							}
						}
					}
				}
			}
		}
	}
	
	fclose(prslab);
	fclose(prslab2);
	return 0;
}

double FindReTheta(Zone *zone,ZoneBoundary *zBc,int indices[3], double *Theta, double *delta_star, double *delta)
{
	double u,L=0,dy,Uref;
	int NI,NJ,NK,i,j,k,maxIndex,minIndex,countBc,numZbc,cellIndex,tag=0;
	ZoneBoundary *NBc;

	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	i = indices[0];
	j = indices[1];
	k = indices[2];
	numZbc=zone->Number_of_BC_Patches;
	*Theta=0;
	*delta_star=0;
	*delta=0;

	if(zBc->Boundary_Tag == 22 && zone->Flowsol->flagTurb == 1)
	{
		switch(zBc->Direction)
		{
			case 1://Imax face
				minIndex = 0;
				for(countBc=0;countBc<numZbc;countBc++)
				{
					NBc = zone->ZoneBc+countBc;
					if(NBc->Direction == 2 && NBc->Boundary_Tag == 22)
					{
						minIndex = NI/2;
						countBc = numZbc;
					}
				}

				//Uref = *(zone->Flowsol->v+(NI-1)*(NJ-1)*k+(NI-1)*j+minIndex);
				Uref = UREF;
				for(i=NI-2;i>=minIndex;i--)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->v+cellIndex);
					dy = *(zone->CgnsGrid->X+(NI*NJ)*k+NI*j+(i+1)) - *(zone->CgnsGrid->X+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);
					u = u/Uref;

					L = L + dy;
					if(u >= 0.9999 && tag == 0)
					{
						*delta = L;
						tag = 1;
					}
					*Theta = *Theta + u*(1.0-u)*dy;
					*delta_star= *delta_star + (1.0-u)*dy;
				}
				*Theta = *Theta/L;
				*delta_star = *delta_star/L;
			break;
			case 2://Imin face
				maxIndex = NI-1;
				for(countBc=0;countBc<numZbc;countBc++)
				{
					NBc = zone->ZoneBc+countBc;
					if(NBc->Direction == 1 && NBc->Boundary_Tag == 22)
					{
						maxIndex = NI/2;
						countBc = numZbc;
					}
				}

				//Uref = *(zone->Flowsol->v+(NI-1)*(NJ-1)*k+(NI-1)*j+maxIndex);
				Uref = UREF;
				for(i=0;i<maxIndex;i++)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->v+cellIndex);
					dy = *(zone->CgnsGrid->X+(NI*NJ)*k+NI*j+(i+1)) - *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);
					u = u/Uref;

					L = L + dy;
					if(u >= 0.9999 && tag == 0)
					{
						*delta = L;
						tag = 1;
					}
					*Theta = *Theta + u*(1.0-u)*dy;
					*delta_star= *delta_star + (1.0-u)*dy;
				}
				*Theta = *Theta/L;
				*delta_star = *delta_star/L;
			break;
			case 3://Jmax face
				minIndex = 0;
				for(countBc=0;countBc<numZbc;countBc++)
				{
					NBc = zone->ZoneBc+countBc;
					if(NBc->Direction == 4 && NBc->Boundary_Tag == 22)
					{
						minIndex = NJ/2;
						countBc = numZbc;
					}
				}

				//Uref = *(zone->Flowsol->u+(NI-1)*(NJ-1)*k+(NI-1)*minIndex+i);
				Uref = UREF;
				for(j=NJ-2;j>=minIndex;j--)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->u+cellIndex);
					dy = *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j+1)+i) - *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j)+i);
					dy = fabs(dy);
					u = u/Uref;

					L = L + dy;
					if(u >= 0.9999 && tag == 0)
					{
						*delta = L;
						tag = 1;
					}
					*Theta = *Theta + u*(1.0-u)*dy;
					*delta_star= *delta_star + (1.0-u)*dy;
				}
				*Theta = *Theta/L;
				*delta_star = *delta_star/L;
			break;
			case 4://Jmin face
				maxIndex = NJ-1;
				for(countBc=0;countBc<numZbc;countBc++)
				{
					NBc = zone->ZoneBc+countBc;
					if(NBc->Direction == 3 && NBc->Boundary_Tag == 22)
					{
						maxIndex = NJ/2;
						countBc = numZbc;
					}
				}

				//Uref = *(zone->Flowsol->u+(NI-1)*(NJ-1)*k+(NI-1)*maxIndex+i);
				Uref = UREF;
				for(j=0;j<maxIndex;j++)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->u+cellIndex);
					dy = *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j+1)+i) - *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j)+i);
					dy = fabs(dy);
					u = u/Uref;

					L = L + dy;
					if(u >= 0.9999 && tag == 0)
					{
						*delta = L;
						tag = 1;
					}
					*Theta = *Theta + u*(1.0-u)*dy;
					*delta_star= *delta_star + (1.0-u)*dy;
					//printf("\ntheta=%lf\tu=%lf\tdy=%e\tUref=%lf\tL=%e\n", *Theta,u,dy,Uref,L);
				}
				*Theta = *Theta/L;
				*delta_star = *delta_star/L;
				//printf("\nx=%lf\ttheta=%lf\tu=%lf\tdy=%e\tUref=%lf\tL=%e\n", *(zone->CgnsGrid->X+(NI*NJ)*k+NI*0+i),*Theta,u,dy,Uref,L);
			break;
			case 5://Kmax face
				minIndex = 0;
				for(countBc=0;countBc<numZbc;countBc++)
				{
					NBc = zone->ZoneBc+countBc;
					if(NBc->Direction == 6 && NBc->Boundary_Tag == 22)
					{
						minIndex = NK/2;
						countBc = numZbc;
					}
				}

				//Uref = *(zone->Flowsol->u+(NI-1)*(NJ-1)*minIndex+(NI-1)*j+i);//u????
				Uref = UREF;
				for(k=NK-2;k>=minIndex;k--)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->u+cellIndex);//u????
					dy = *(zone->CgnsGrid->Z+(NI*NJ)*(k+1)+NI*j+i) - *(zone->CgnsGrid->Z+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);
					u = u/Uref;

					L = L + dy;
					if(u >= 0.9999 && tag == 0)
					{
						*delta = L;
						tag = 1;
					}
					*Theta = *Theta + u*(1.0-u)*dy;
					*delta_star= *delta_star + (1.0-u)*dy;
				}
				*Theta = *Theta/L;
				*delta_star = *delta_star/L;
			break;
			case 6://Kmin face
				maxIndex = NK-1;
				for(countBc=0;countBc<numZbc;countBc++)
				{
					NBc = zone->ZoneBc+countBc;
					if(NBc->Direction == 5 && NBc->Boundary_Tag == 22)
					{
						maxIndex = NK/2;
						countBc = numZbc;
					}
				}

				//Uref = *(zone->Flowsol->u+(NI-1)*(NJ-1)*maxIndex+(NI-1)*j+i);
				Uref = UREF;
				for(k=0;k<maxIndex;k++)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->u+cellIndex);
					dy = *(zone->CgnsGrid->Z+(NI*NJ)*(k+1)+NI*j+i) - *(zone->CgnsGrid->Z+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);
					u = u/Uref;

					L = L + dy;
					if(u >= 0.9999 && tag == 0)
					{
						*delta = L;
						tag = 1;
					}
					*Theta = *Theta + u*(1.0-u)*dy;
					*delta_star= *delta_star + (1.0-u)*dy;
				}
				*Theta = *Theta/L;
				*delta_star = *delta_star/L;
			break;
			default:
				printf("In FindReTheta(), Improper direction of the Boundary Patch .. dir=%d\n", zBc->Direction);
			break;
		}
	}

	return(Uref);
}


double FindUavg(Zone *zone,ZoneBoundary *zBc,int indices[3])
{
	double u=0.0,L=0.0,dy,Uavg=0.0;//Don't tamper with this value
	int NI,NJ,NK,i,j,k,maxIndex,minIndex,cellIndex;

	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	i = indices[0];
	j = indices[1];
	k = indices[2];

	if(zBc->Boundary_Tag == 22 && zone->Flowsol->flagTurb == 1)
	{
		switch(zBc->Direction)
		{
			case 1://Imax face
				minIndex = 0;

				for(i=NI-2;i>=minIndex;i--)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->v+cellIndex);
					dy = *(zone->CgnsGrid->X+(NI*NJ)*k+NI*j+(i+1)) - *(zone->CgnsGrid->X+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);

					L = L + dy;
					Uavg = Uavg + u*dy;
				}
				Uavg = Uavg/L;
			break;
			case 2://Imin face
				maxIndex = NI-1;

				for(i=0;i<maxIndex;i++)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->v+cellIndex);
					dy = *(zone->CgnsGrid->X+(NI*NJ)*k+NI*j+(i+1)) - *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);

					L = L + dy;
					Uavg = Uavg + u*dy;
				}
				Uavg = Uavg/L;
			break;
			case 3://Jmax face
				minIndex = 0;

				for(j=NJ-2;j>=minIndex;j--)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->u+cellIndex);
					dy = *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j+1)+i) - *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j)+i);
					dy = fabs(dy);

					L = L + dy;
					Uavg = Uavg + u*dy;
				}
				Uavg = Uavg/L;
			break;
			case 4://Jmin face
				maxIndex = NJ-1;

				for(j=0;j<maxIndex;j++)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					u = *(zone->Flowsol->u+cellIndex);
					dy = *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j+1)+i) - *(zone->CgnsGrid->Y+(NI*NJ)*k+NI*(j)+i);
					dy = fabs(dy);
					
					L = L + dy;
					Uavg = Uavg + u*dy;
					//printf("\nu=%lf\tdy=%e\tUavg=%lf\tL=%e\n", u,dy,Uavg,L);
				}
				Uavg = Uavg/L;
			break;
			case 5://Kmax face
				minIndex = 0;

				for(k=NK-2;k>=minIndex;k--)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					if((zone->Flowsol->w) != NULL)
						u = *(zone->Flowsol->w+cellIndex);//u????
					dy = *(zone->CgnsGrid->Z+(NI*NJ)*(k+1)+NI*j+i) - *(zone->CgnsGrid->Z+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);
					
					L = L + dy;
					Uavg = Uavg + u*dy;
					//printf("\nu=%lf\tdy=%e\tUavg=%lf\tL=%e\n", u,dy,Uavg,L);
				}
				Uavg = Uavg/L;

			break;
			case 6://Kmin face
				maxIndex = NK-1;

				for(k=0;k<maxIndex;k++)
				{
					cellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					if((zone->Flowsol->w) != NULL)
						u = *(zone->Flowsol->w+cellIndex);
					dy = *(zone->CgnsGrid->Z+(NI*NJ)*(k+1)+NI*j+i) - *(zone->CgnsGrid->Z+(NI*NJ)*k+NI*j+i);
					dy = fabs(dy);
					
					L = L + dy;
					Uavg = Uavg + u*dy;
					//printf("\nu=%lf\tdy=%e\tUavg=%lf\tL=%e\n", u,dy,Uavg,L);
				}
				Uavg = Uavg/L;
			break;
			default:
				printf("In FindReTheta(), Improper direction of the Boundary Patch .. dir=%d\n", zBc->Direction);
			break;
		}
	}

	Uavg = 1.0;//Temporarily written while running for flatplate boundary layer. Remove this in case of channel flow problem
	return (Uavg);
}

//FOR CHECKING SUSHIL
int printslabs(Base *base)
{
	int numzones,cntZ;
	int i,j,k,numIface,numZbc,dir,zBc;
	int index,NI,NJ,NK;
	int *ptrRange,flagRho=0,flagViscoeff=0,nophase,nop;
	//double rho,mu;
	double momentumx,momentumy,momentumz;
	double u,v,vartheta,ndensity,Cdensity;//adhi+shash
	Zone *zone;
	ZoneBoundary *ZBc;
	numzones= base->Number_of_zones;
	FILE *prslab;


	prslab= fopen("SOLUTIONS/printslab.dat","w");
	nop = base->nphase;
	if(base->flagMPFV == 1)
	{
		for(cntZ=0;cntZ<numzones;cntZ++)
		{
			zone = base->CgnsZone+cntZ;
			flagViscoeff = zone->Flowsol->flagViscoeff;
			flagRho = zone->Flowsol->flagRho;
			NI=*(zone->Size+0);
			NJ=*(zone->Size+1);
			NK=*(zone->Size+2);
			if(zone->flagSolidFluid == 0)
			{
				
				numIface=zone->Number_of_Interfaces;
				numZbc=zone->Number_of_BC_Patches;
				for(dir=1;dir<=4;dir++)
				{
					for(zBc=0;zBc<numZbc;zBc++)
					{
						ZBc = zone->ZoneBc+zBc;
						if(ZBc->Direction==dir)
						{
							//printf("boundary dir=%d\n",dir);
							ptrRange=ZBc->Point_Range;
							//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
							
							fprintf(prslab,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
							//printf("ZoneBc dir->%d\n",dir);
							//printf("~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
							if(dir==1||dir==2)
							{
								//printf("inside if\n");
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								{
									//printf("k=%d,j=%d\n",k,j);
									index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
									//printf("INDEX=%d\n",index);
									
									for(nophase=0;nophase<nop;nophase++)
									{
										momentumx = *(ZBc->ZSlab->Momentumofphasex[nophase]+index);
										momentumy = *(ZBc->ZSlab->Momentumofphasey[nophase]+index);
										momentumz = *(ZBc->ZSlab->Momentumofphasez[nophase]+index);
										fprintf(prslab,"\t%e\t%e\t%e\t",momentumx,momentumy,momentumz);
									}
									momentumx = *(ZBc->ZSlab->Momentumx+index);
									momentumy = *(ZBc->ZSlab->Momentumy+index);
									momentumz = *(ZBc->ZSlab->Momentumz+index);
									fprintf(prslab,"\t%e\t%e\t%e\n",momentumx,momentumy,momentumz);
									
								}
							}
							if(dir==3||dir==4)
							{

								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									
									
									for(nophase=0;nophase<nop;nophase++)
									{
										momentumx = *(ZBc->ZSlab->Momentumofphasex[nophase]+index);
										momentumy = *(ZBc->ZSlab->Momentumofphasey[nophase]+index);
										momentumz = *(ZBc->ZSlab->Momentumofphasez[nophase]+index);
										fprintf(prslab,"\t%e\t%e\t%e\t",momentumx,momentumy,momentumz);
									}
									momentumx = *(ZBc->ZSlab->Momentumx+index);
									momentumy = *(ZBc->ZSlab->Momentumy+index);
									momentumz = *(ZBc->ZSlab->Momentumz+index);
									fprintf(prslab,"\t%e\t%e\t%e\n",momentumx,momentumy,momentumz);
								}
							}
							if(dir==5||dir==6)
							{
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									for(nophase=0;nophase<nop;nophase++)
									{
										momentumx = *(ZBc->ZSlab->Momentumofphasex[nophase]+index);
										momentumy = *(ZBc->ZSlab->Momentumofphasey[nophase]+index);
										momentumz = *(ZBc->ZSlab->Momentumofphasez[nophase]+index);
										fprintf(prslab,"\t%e\t%e\t%e\t",momentumx,momentumy,momentumz);
									}
									momentumx = *(ZBc->ZSlab->Momentumx+index);
									momentumy = *(ZBc->ZSlab->Momentumy+index);
									momentumz = *(ZBc->ZSlab->Momentumz+index);
									fprintf(prslab,"\t%e\t%e\t%e\n",momentumx,momentumy,momentumz);
								}
							}
						}
					}
				}
				
			}
		}
	}
	//adhi+shash
	if((base->flagDroplet == 1)||(base->flagParticle == 1))
	{
		for(cntZ=0;cntZ<numzones;cntZ++)
		{
			zone = base->CgnsZone+cntZ;
			flagViscoeff = zone->Flowsol->flagViscoeffDisperse;
			flagRho = zone->Flowsol->flagRhoDisperse;
			NI=*(zone->Size+0);
			NJ=*(zone->Size+1);
			NK=*(zone->Size+2);
			if(zone->flagSolidFluid == 0)
			{
				numIface=zone->Number_of_Interfaces;
				numZbc=zone->Number_of_BC_Patches;
				for(dir=1;dir<=6;dir++)
				{
					for(zBc=0;zBc<numZbc;zBc++)
					{
						ZBc = zone->ZoneBc+zBc;
						if(ZBc->Direction==dir)
						{
							//printf("boundary dir=%d\n",dir);
							ptrRange=ZBc->Point_Range;
							//printf("\t\t %d %d %d %d %d %d \n",*(ptrRange+0),*(ptrRange+1),*(ptrRange+2),*(ptrRange+3),*(ptrRange+4),*(ptrRange+5));
							fprintf(prslab,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
							//printf("ZoneBc dir->%d\n",dir);
							//printf("~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",ZBc->BCName);
							fprintf(prslab,"Vx\tVy\tu\tv\tw\tVar\n");
							if(dir==1||dir==2)
							{
									//printf("inside if\n");
								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								{
									//printf("k=%d,j=%d\n",k,j);
									index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
									//printf("INDEX=%d\n",index);
									momentumx = *(ZBc->ZSlab->MomentumDispersex+index);
									momentumy = *(ZBc->ZSlab->MomentumDispersey+index);
									//momentumz = *(ZBc->ZSlab->MomentumDispersez+index);
									fprintf(prslab,"\t%e\t%e\t",momentumx,momentumy);
									u = *(ZBc->ZSlab->u+index);
									v = *(ZBc->ZSlab->v+index);
								//w = *(ZBc->ZSlab->w+index);
				
									fprintf(prslab,"%e\t%e\t",u,v);
									vartheta = *(ZBc->ZSlab->VarthetaDisperse+index);
									fprintf(prslab,"%e\t",vartheta);
									if(base->flagDroplet == 1)
									{
										ndensity  = *(ZBc->ZSlab->NumDensity+index);
										fprintf(prslab,"%e\t",ndensity);
				
										Cdensity  = *(ZBc->ZSlab->Rho+index);
										fprintf(prslab,"%e",Cdensity);
									}   
									fprintf(prslab,"\n");   
								}
							}
							if(dir==3||dir==4)
							{

								for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
				
				
									momentumx = *(ZBc->ZSlab->MomentumDispersex+index);
									momentumy = *(ZBc->ZSlab->MomentumDispersey+index);
									//momentumz = *(ZBc->ZSlab->MomentumDispersez+index);
									fprintf(prslab,"\t%e\t%e\t",momentumx,momentumy);
									u = *(ZBc->ZSlab->u+index);
									v = *(ZBc->ZSlab->v+index);
									//w = *(ZBc->ZSlab->w+index);
				
									fprintf(prslab,"%e\t%e\t",u,v);
									vartheta = *(ZBc->ZSlab->VarthetaDisperse+index);
									fprintf(prslab,"%e\t",vartheta);
									if(base->flagDroplet == 1)
									{
										ndensity  = *(ZBc->ZSlab->NumDensity+index);
										fprintf(prslab,"%e\n",ndensity);
										Cdensity  = *(ZBc->ZSlab->Rho+index);
										fprintf(prslab,"%e",Cdensity);
									}   
									fprintf(prslab,"\n");  
								}
							}
							if(dir==5||dir==6)
							{
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
									for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
								{
									index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
									momentumx = *(ZBc->ZSlab->MomentumDispersex+index);
									momentumy = *(ZBc->ZSlab->MomentumDispersey+index);
									//momentumz = *(ZBc->ZSlab->MomentumDispersez+index);
									fprintf(prslab,"\t%e\t%e\t",momentumx,momentumy);
									u = *(ZBc->ZSlab->u+index);
									v = *(ZBc->ZSlab->v+index);
									//w = *(ZBc->ZSlab->w+index);
				
									fprintf(prslab,"%e\t%e\t",u,v);
									vartheta = *(ZBc->ZSlab->VarthetaDisperse+index);
									fprintf(prslab,"%e\t",vartheta);
									if(base->flagDroplet == 1)
									{
										ndensity  = *(ZBc->ZSlab->NumDensity+index);
										fprintf(prslab,"%e\n",ndensity);
										Cdensity  = *(ZBc->ZSlab->Rho+index);
										fprintf(prslab,"%e",Cdensity);
									}   
									fprintf(prslab,"\n");
								}
							}
						}
					}
				}
		
			}
		}
	}
	
	fclose(prslab);
	
	return 0;
}


