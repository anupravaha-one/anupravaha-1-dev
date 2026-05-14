#include "APKmassflux2d.h"

/**
This function is used only once for calculating the massflux just before the first timestep using
the guess values of mass velocities. This is used for 2D geometries.

*/

int  massflux2d(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int numiface,numBc,face,iface,zBc,dir,Index,index;
	double *SurfArea,*F,*rho;
	int indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUmass,*domainVmass;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6];
	double *slabRho[6],*domainRho;
	Zone* zone;

	int NI,NJ,NK,*ptrSize;
	//double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,v_P,v_E,v_W,v_S,v_N;
	double vol_P,vol_E,vol_W,vol_N,vol_S;
	double u_p,u_e,u_w,u_s,u_n,v_p,v_e,v_w,v_s,v_n;
	double vol_p,vol_e,vol_w,vol_n,vol_s;
	double rho_P,rho_E,rho_W,rho_N,rho_S;
	double rhoue,rhove,rhouw,rhovw,rhoun,rhovn,rhous,rhovs;
	
	nZones=base->Number_of_zones;	
	for(numZones=0;numZones<nZones;numZones++)
	{
		/*Copying down u,v,w into umass,vmass,wmass ....domain as well as slabs*/
		ptrSize=((base)->CgnsZone+numZones)->Size;
		zone=(base->CgnsZone+numZones);
		if((zone->flagSolidFluid==0))
		{
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
				*((zone->Flowsol)->umass+indexP)=*((zone->Flowsol)->u+indexP);
				*((zone->Flowsol)->vmass+indexP)=*((zone->Flowsol)->v+indexP);

			}
			numiface=(zone)->Number_of_Interfaces;

			numBc=(zone)->Number_of_BC_Patches;

			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numiface;iface++)
				{
					if(((zone)->Iface+iface)->Direction==dir)
					{
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							if(dir==1||dir==2)
								index= k*(NJ-1)+j;
							if(dir==3||dir==4)
								index= k*(NI-1)+i;
							if(dir==5||dir==6)
								index= j*(NI-1)+i;
							*((((zone)->Iface+iface)->ZSlab)->umass+index) = *((((zone)->Iface+iface)->ZSlab)->u+index);
							*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);

						}

					}
				}
				for(zBc=0;zBc<numBc;zBc++)
				{
					if(((zone)->ZoneBc+zBc)->Direction==dir)
					{
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							if(dir==1||dir==2)
								index= k*(NJ-1)+j;
							if(dir==3||dir==4)
								index= k*(NI-1)+i; 
							if(dir==5||dir==6)
								index= j*(NI-1)+i;
							*((((zone)->ZoneBc+zBc)->ZSlab)->umass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->u+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->vmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->v+index);

						}
					}
				}
			}/*End of direction Loop*/
		}
	}/*End of Zone loop*/

	nZones=base->Number_of_zones;
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		if( (zone->flagSolidFluid==0) && ((zone->Flowsol)->flagRho==1) )
		{/*Done By Ashish for Conjugate Heat Transfer Case*/	
			ptrSize=zone->Size;

			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			F = zone->Flowsol->Massflux;

			SurfArea=zone->ZGeom->SurfaceArea;

			rho=zone->Flowsol->Rho;

			domainUmass=MakeDomainArray(zone,100);
			domainVmass=MakeDomainArray(zone,200);
			domainVolume=MakeDomainArray(zone,5);
			domainRho=MakeDomainArray(zone,710);/*Density*/

			MakeSlabPtrArray(zone,slabUmass,1);
			MakeSlabPtrArray(zone,slabVmass,2);
			MakeSlabPtrArray(zone,slabVolume,5);
			MakeSlabPtrArray(zone,slabRho,710);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
				indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
				indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
				indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
				indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;

				u_P=*(domainUmass+indexP);
				v_P=*(domainVmass+indexP);
				rho_P=*(domainRho+indexP);
				vol_P=*(domainVolume+indexP);

				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+(NJ-1)*k+j);
					v_E=*(slabVmass[0]+(NJ-1)*k+j);
					rho_E=*(slabRho[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+(NJ-1)*k+j);
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
					rho_E=*(domainRho+indexE);
					vol_E=*(domainVolume+indexE);
				}

				if(i==0)
				{
					u_W =*(slabUmass[1]+(NJ-1)*k+j);
					v_W =*(slabVmass[1]+(NJ-1)*k+j);
					rho_W =*(slabRho[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+(NJ-1)*k+j);
				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
					rho_W=*(domainRho+indexW);
					vol_W=*(domainVolume+indexW);
				}

				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+(NI-1)*k+i);
					v_N=*(slabVmass[2]+(NI-1)*k+i);
					rho_N=*(slabRho[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+(NI-1)*k+i);
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
					rho_N=*(domainRho+indexN);
					vol_N=*(domainVolume+indexN);
				}

				if(j==0)
				{
					u_S=*(slabUmass[3]+(NI-1)*k+i);
					v_S=*(slabVmass[3]+(NI-1)*k+i);
					rho_S=*(slabRho[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+(NI-1)*k+i);
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
					rho_S=*(domainRho+indexS);
					vol_S=*(domainVolume+indexS);
				}

				


				rhoue=((rho_E*u_E*vol_P+rho_P*u_P*vol_E)/(vol_P+vol_E));
				rhove=((rho_E*v_E*vol_P+rho_P*v_P*vol_E)/(vol_P+vol_E));

				rhouw=((rho_W*u_W*vol_P+rho_P*u_P*vol_W)/(vol_P+vol_W));
				rhovw=((rho_W*v_W*vol_P+rho_P*v_P*vol_W)/(vol_P+vol_W));

				rhoun=((rho_N*u_N*vol_P+rho_P*u_P*vol_N)/(vol_P+vol_N));
				rhovn=((rho_N*v_N*vol_P+rho_P*v_P*vol_N)/(vol_P+vol_N));

				rhous=((rho_S*u_S*vol_P+rho_P*u_P*vol_S)/(vol_P+vol_S));
				rhovs=((rho_S*v_S*vol_P+rho_P*v_P*vol_S)/(vol_P+vol_S));

				Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
				face = 0;  // west face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhouw*(*(SurfArea+Index+3*0+0))+rhovw*(*(SurfArea+Index+3*0+1)));
					//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 1;  // east face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoue*(*(SurfArea+Index+3*1+0))+rhove*(*(SurfArea+Index+3*1+1)));
					//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 2;  // south face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhous*(*(SurfArea+Index+3*2+0))+rhovs*(*(SurfArea+Index+3*2+1)));
					//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 3;  // north face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoun*(*(SurfArea+Index+3*3+0))+rhovn*(*(SurfArea+Index+3*3+1)));
					//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 4;  // bottom face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 5;  // top face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));;
			}/*End of looping over NI,NJ,NK*/
		}
		if( (zone->flagSolidFluid==0) && ((zone->Flowsol)->flagRho==0) )
		{/*Done By Ashish for Conjugate Heat Transfer Case*/	
			ptrSize=zone->Size;

			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			F = zone->Flowsol->Massflux;

			SurfArea=zone->ZGeom->SurfaceArea;

			rho=zone->Flowsol->Rho;

			domainUmass=MakeDomainArray(zone,100);
			domainVmass=MakeDomainArray(zone,200);
			domainVolume=MakeDomainArray(zone,5);
			domainRho=MakeDomainArray(zone,710);/*Density*/

			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			MakeSlabPtrArray(zone,slabVolume,5);
			MakeSlabPtrArray(zone,slabRho,710);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
				indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
				indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
				indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
				indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;

				u_p=*(domainUmass+indexP);
				v_p=*(domainVmass+indexP);
				vol_p=*(domainVolume+indexP);

				if(i==(NI-2))
				{
					u_e=*(slabUmass[0]+(NJ-1)*k+j);
					v_e=*(slabVmass[0]+(NJ-1)*k+j);
					vol_e=*(slabVolume[0]+(NJ-1)*k+j);
				}
				else
				{
					u_e=*(domainUmass+indexE);
					v_e=*(domainVmass+indexE);
					vol_e=*(domainVolume+indexE);
				}

				if(i==0)
				{
					u_w =*(slabUmass[1]+(NJ-1)*k+j);
					v_w =*(slabVmass[1]+(NJ-1)*k+j);
					vol_w =*(slabVolume[1]+(NJ-1)*k+j);
				}
				else
				{
					u_w=*(domainUmass+indexW);
					v_w=*(domainVmass+indexW);
					vol_w=*(domainVolume+indexW);
				}

				if(j==(NJ-2))
				{
					u_n=*(slabUmass[2]+(NI-1)*k+i);
					v_n=*(slabVmass[2]+(NI-1)*k+i);
					vol_n=*(slabVolume[2]+(NI-1)*k+i);
				}
				else
				{
					u_n=*(domainUmass+indexN);
					v_n=*(domainVmass+indexN);
					vol_n=*(domainVolume+indexN);
				}

				if(j==0)
				{
					u_s=*(slabUmass[3]+(NI-1)*k+i);
					v_s=*(slabVmass[3]+(NI-1)*k+i);
					vol_s=*(slabVolume[3]+(NI-1)*k+i);
				}
				else
				{
					u_s=*(domainUmass+indexS);
					v_s=*(domainVmass+indexS);
					vol_s=*(domainVolume+indexS);
				}
					
				rhoue=((u_e*vol_p+u_p*vol_e)/(vol_p+vol_e));
				rhouw=((u_w*vol_p+u_p*vol_w)/(vol_p+vol_w));
				rhous=((u_s*vol_p+u_p*vol_s)/(vol_p+vol_s));
				rhoun=((u_n*vol_p+u_p*vol_n)/(vol_p+vol_n));
				rhove=((v_e*vol_p+v_p*vol_e)/(vol_p+vol_e));
				rhovw=((v_w*vol_p+v_p*vol_w)/(vol_p+vol_w));
				rhovs=((v_s*vol_p+v_p*vol_s)/(vol_p+vol_s));
				rhovn=((v_n*vol_p+v_p*vol_n)/(vol_p+vol_n));
				
				
				//printf("%d %d %d rhoue=%lf\trhove=%lf\trhouw=%lf\trhovw=%lf\trhoun=%lf\trhovn=%lf\trhous=%lf\trhovs=%lf\n",i,j,k,rhoue,rhove,rhouw,rhovw,rhoun,rhovn,rhous,rhovs);
				Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
				//printf("\nI=%d\tJ=%d\tK=%d\t",i,j,k);
				face = 0;  // west face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(rhouw*(*(SurfArea+Index+3*0+0))+rhovw*(*(SurfArea+Index+3*0+1)));
					//printf("west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 1;  // east face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(rhoue*(*(SurfArea+Index+3*1+0))+rhove*(*(SurfArea+Index+3*1+1)));
					//printf("east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 2;  // south face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(rhous*(*(SurfArea+Index+3*2+0))+rhovs*(*(SurfArea+Index+3*2+1)));
					//printf("south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 3;  // north face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(rhoun*(*(SurfArea+Index+3*3+0))+rhovn*(*(SurfArea+Index+3*3+1)));
					//printf("north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 4;  // bottom face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					//printf("bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 5;  // top face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					//printf("top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			}/*End of looping over NI,NJ,NK*/
		}/*end of else */
	}/*End of looping over the zones */
	return 0;
}
