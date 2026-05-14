#include "APKmassfluxstar2d.h"
/**
This function is used for calculating the massflux from the current converged
values of massvelocities. This function is used for 3D domains.
See Eqn No(3.36) Project Report No. 1, Dec 2005.
*/

int massfluxstar2d(Base *base)
{
	//printf("I am in massfluxstar\n");
	int nZones,numZones;
	int k,j,i,face;
	int Index,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUmass,*domainVmass;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6];
	double *slabRho[6],*domainRho;

	double *SurfArea,*F,*rho;
	Zone* zone;
	nZones=(base->Number_of_zones);
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ue,ve,un,vn,us,vs;
	double u_P,u_E,u_W,u_S,u_N,v_P,v_E,v_W,v_S,v_N;
	double vol_P,vol_E,vol_W,vol_N,vol_S;
	double rho_P,rho_E,rho_W,rho_N,rho_S;
	double rhoue,rhove,rhouw,rhovw,rhoun,rhovn,rhous,rhovs;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		if( (zone->flagSolidFluid==0) && ((zone->Flowsol)->flagRho==1) )
		{
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
					//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			}/*End of looping over NI,NJ,NK*/
		}
		if( (zone->flagSolidFluid==0) && ((zone->Flowsol)->flagRho==0) )
		{
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
			
			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			MakeSlabPtrArray(zone,slabVolume,5);
			
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
				vol_P=*(domainVolume+indexP);

				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+(NJ-1)*k+j);
					v_E=*(slabVmass[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+(NJ-1)*k+j);
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
					vol_E=*(domainVolume+indexE);
				}

				if(i==0)
				{
					u_W =*(slabUmass[1]+(NJ-1)*k+j);
					v_W =*(slabVmass[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+(NJ-1)*k+j);
				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
					vol_W=*(domainVolume+indexW);
				}

				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+(NI-1)*k+i);
					v_N=*(slabVmass[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+(NI-1)*k+i);
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
					vol_N=*(domainVolume+indexN);
				}

				if(j==0)
				{
					u_S=*(slabUmass[3]+(NI-1)*k+i);
					v_S=*(slabVmass[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+(NI-1)*k+i);
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
					vol_S=*(domainVolume+indexS);
				}

				ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
				uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
				us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
				un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));

				ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
				vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
				vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
				vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));

				Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
				face = 0;  // west face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(uw*(*(SurfArea+Index+3*0+0))+vw*(*(SurfArea+Index+3*0+1)));
					//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 1;  // east face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(ue*(*(SurfArea+Index+3*1+0))+ve*(*(SurfArea+Index+3*1+1)));
					//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 2;  // south face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(us*(*(SurfArea+Index+3*2+0))+vs*(*(SurfArea+Index+3*2+1)));
					//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 3;  // north face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(un*(*(SurfArea+Index+3*3+0))+vn*(*(SurfArea+Index+3*3+1)));
					//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 4;  // bottom face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 5;  // top face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
					
			}/*End of looping over NI,NJ,NK*/
		}/*end of else */
	}/*End of looping over the zones */
	return 0;
}
