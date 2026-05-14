#include "APKMagneticflux.h"
/**
This function is used only once for calculating the magneticflux just before the first timestep using
the guess values of magnetic field. This is used for 3D geometries.
*/
int Magneticflux(Base *base)
{
	int nZones,numZones;
	int k,j,i,indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	int numiface,numBc,face,Index;
	//int iface,zBc,dir,Index;
	double *SurfArea,*F;
	double *domainVolume,*domainBx,*domainBy,*domainBz;//rajendra
	double *slabVolume[6],*slabBx[6],*slabBy[6],*slabBz[6];//rajendra
	Zone* zone;
	
	int NI,NJ,NK,*ptrSize;
	double Bxw,Byw,Bzw,Bxe,Bye,Bze,Bxn,Byn,Bzn,Bxs,Bys,Bzs,Bxt,Byt,Bzt,Bxb,Byb,Bzb;//rajendra
	double Bx_P,Bx_E,Bx_W,Bx_S,Bx_N,Bx_T,Bx_B,By_P,By_E,By_W,By_S,By_N,By_T,By_B,Bz_P,Bz_E,Bz_W,Bz_S,Bz_N,Bz_T,Bz_B;//r

	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	nZones=base->Number_of_zones;	
	for(numZones=0;numZones<nZones;numZones++)
	{
		
		ptrSize=((base)->CgnsZone+numZones)->Size;
		zone=(base->CgnsZone+numZones);
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		
		numiface=(zone)->Number_of_Interfaces;
		
		numBc=(zone)->Number_of_BC_Patches;
		
		
	}/*End of Zone loop*/

	nZones=base->Number_of_zones;
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->Magneticflux;

		SurfArea=zone->ZGeom->SurfaceArea;
		domainBx=MakeDomainArray(zone,56);
		domainBy=MakeDomainArray(zone,57);
		domainBz=MakeDomainArray(zone,58);
		domainVolume=MakeDomainArray(zone,5);
		
		MakeSlabPtrArray(zone,slabBx,56);
		MakeSlabPtrArray(zone,slabBy,57);
		MakeSlabPtrArray(zone,slabBz,58);
		
		MakeSlabPtrArray(zone,slabVolume,5);

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
// 			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			indexE=indexP+1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
			indexW=indexP-1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
			indexN=indexP+(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
			indexS=indexP-(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexT=indexP+(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=indexP-(NI-1)*(NJ-1);//(NJ-1)

			Bx_P=*(domainBx+indexP);
			By_P=*(domainBy+indexP);
			Bz_P=*(domainBz+indexP);
			vol_P=*(domainVolume+indexP);

			if(i==(NI-2))
			{
				Bx_E=*(slabBx[0]+(NJ-1)*k+j);
				By_E=*(slabBy[0]+(NJ-1)*k+j);
				Bz_E=*(slabBz[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				Bx_E=*(domainBx+indexE);
				By_E=*(domainBy+indexE);
				Bz_E=*(domainBz+indexE);
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				Bx_W =*(slabBx[1]+(NJ-1)*k+j);
				By_W =*(slabBy[1]+(NJ-1)*k+j);
				Bz_W =*(slabBz[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				Bx_W=*(domainBx+indexW);
				By_W=*(domainBy+indexW);
				Bz_W=*(domainBz+indexW);
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				Bx_N=*(slabBx[2]+(NI-1)*k+i);
				By_N=*(slabBy[2]+(NI-1)*k+i);
				Bz_N=*(slabBz[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				Bx_N=*(domainBx+indexN);
				By_N=*(domainBy+indexN);
				Bz_N=*(domainBz+indexN);
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				Bx_S=*(slabBx[3]+(NI-1)*k+i);
				By_S=*(slabBy[3]+(NI-1)*k+i);
				Bz_S=*(slabBz[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				Bx_S=*(domainBx+indexS);
				By_S=*(domainBy+indexS);
				Bz_S=*(domainBz+indexS);
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				Bx_T=*(slabBx[4]+(NI-1)*j+i);
				By_T=*(slabBy[4]+(NI-1)*j+i);
				Bz_T=*(slabBz[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				Bx_T=*(domainBx+indexT);
				By_T=*(domainBy+indexT);
				Bz_T=*(domainBz+indexT);
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				Bx_B=*(slabBx[5]+(NI-1)*j+i);
				By_B=*(slabBy[5]+(NI-1)*j+i);
				Bz_B=*(slabBz[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				Bx_B=*(domainBx+indexB);
				By_B=*(domainBy+indexB);
				Bz_B=*(domainBz+indexB);
				vol_B=*(domainVolume+indexB);
			}


			Bxe=((Bx_E*vol_P+Bx_P*vol_E)/(vol_P+vol_E));
			Bxw=((Bx_W*vol_P+Bx_P*vol_W)/(vol_P+vol_W));
			Bxs=((Bx_S*vol_P+Bx_P*vol_S)/(vol_P+vol_S));
			Bxn=((Bx_N*vol_P+Bx_P*vol_N)/(vol_P+vol_N));
			Bxt=((Bx_T*vol_P+Bx_P*vol_T)/(vol_P+vol_T));
			Bxb=((Bx_B*vol_P+Bx_P*vol_B)/(vol_P+vol_B));

			Bye=((By_E*vol_P+By_P*vol_E)/(vol_P+vol_E));
			Byw=((By_W*vol_P+By_P*vol_W)/(vol_P+vol_W));
			Bys=((By_S*vol_P+By_P*vol_S)/(vol_P+vol_S));
			Byn=((By_N*vol_P+By_P*vol_N)/(vol_P+vol_N));
			Byt=((By_T*vol_P+By_P*vol_T)/(vol_P+vol_T));
			Byb=((By_B*vol_P+By_P*vol_B)/(vol_P+vol_B));

			Bze=((Bz_E*vol_P+Bz_P*vol_E)/(vol_P+vol_E));
			Bzw=((Bz_W*vol_P+Bz_P*vol_W)/(vol_P+vol_W));
			Bzs=((Bz_S*vol_P+Bz_P*vol_S)/(vol_P+vol_S));
			Bzn=((Bz_N*vol_P+Bz_P*vol_N)/(vol_P+vol_N));
			Bzt=((Bz_T*vol_P+Bz_P*vol_T)/(vol_P+vol_T));
			Bzb=((Bz_B*vol_P+Bz_P*vol_B)/(vol_P+vol_B));

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxw*(*(SurfArea+Index+3*0+0))+Byw*(*(SurfArea+Index+3*0+1))+Bzw*(*(SurfArea+Index+3*0+2)));
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxe*(*(SurfArea+Index+3*1+0))+Bye*(*(SurfArea+Index+3*1+1))+Bze*(*(SurfArea+Index+3*1+2)));
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxs*(*(SurfArea+Index+3*2+0))+Bys*(*(SurfArea+Index+3*2+1))+Bzs*(*(SurfArea+Index+3*2+2)));
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 3;  // north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxn*(*(SurfArea+Index+3*3+0))+Byn*(*(SurfArea+Index+3*3+1))+Bzn*(*(SurfArea+Index+3*3+2)));
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 4;  // bottom face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxb*(*(SurfArea+Index+3*4+0))+Byb*(*(SurfArea+Index+3*4+1))+Bzb*(*(SurfArea+Index+3*4+2)));
				//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 5;  // top face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxt*(*(SurfArea+Index+3*5+0))+Byt*(*(SurfArea+Index+3*5+1))+Bzt*(*(SurfArea+Index+3*5+2)));
				//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	return 0;
}

/**
This function is used only once for calculating the Magneticflux just before the first timestep using
the guess values of Magneticfield. This is used for 2D geometries.
*/

int Magneticflux2d(Base *base)
{
	int nZones,numZones;
	int k,j,i,indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	int numiface,numBc,face,Index;
	//int iface,zBc,dir,Index;
	double *SurfArea,*F;
	double *domainVolume,*domainBx,*domainBy;//rajendra
	double *slabVolume[6],*slabBx[6],*slabBy[6];//rajendra
	Zone* zone;
	
	int NI,NJ,NK,*ptrSize;
	double Bxw,Byw,Bxe,Bye,Bxn,Byn,Bxs,Bys,Bxt,Byt,Bxb,Byb;//rajendra
	double Bx_P,Bx_E,Bx_W,Bx_S,Bx_N,Bx_T,Bx_B,By_P,By_E,By_W,By_S,By_N,By_T,By_B;//rajendra

	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	nZones=base->Number_of_zones;	
	for(numZones=0;numZones<nZones;numZones++)
	{
		
		ptrSize=((base)->CgnsZone+numZones)->Size;
		zone=(base->CgnsZone+numZones);
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		
		numiface=(zone)->Number_of_Interfaces;
		
		numBc=(zone)->Number_of_BC_Patches;
		
		
	}/*End of Zone loop*/
	nZones=base->Number_of_zones;
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->Magneticflux;

		SurfArea=zone->ZGeom->SurfaceArea;
		domainBx=MakeDomainArray(zone,56);
		domainBy=MakeDomainArray(zone,57);
		
		
		domainVolume=MakeDomainArray(zone,5);
		
		MakeSlabPtrArray(zone,slabBx,56);
		MakeSlabPtrArray(zone,slabBy,57);
		
		
		MakeSlabPtrArray(zone,slabVolume,5);

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
// 			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			indexE=indexP+1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
			indexW=indexP-1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
			indexN=indexP+(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
			indexS=indexP-(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexT=indexP+(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=indexP-(NI-1)*(NJ-1);//(NJ-1)

			Bx_P=*(domainBx+indexP);
			By_P=*(domainBy+indexP);
			
			vol_P=*(domainVolume+indexP);

			if(i==(NI-2))
			{
				Bx_E=*(slabBx[0]+(NJ-1)*k+j);
				By_E=*(slabBy[0]+(NJ-1)*k+j);
				
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				Bx_E=*(domainBx+indexE);
				By_E=*(domainBy+indexE);
				
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				Bx_W =*(slabBx[1]+(NJ-1)*k+j);
				By_W =*(slabBy[1]+(NJ-1)*k+j);
				
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				Bx_W=*(domainBx+indexW);
				By_W=*(domainBy+indexW);
				
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				Bx_N=*(slabBx[2]+(NI-1)*k+i);
				By_N=*(slabBy[2]+(NI-1)*k+i);
				
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				Bx_N=*(domainBx+indexN);
				By_N=*(domainBy+indexN);
				
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				Bx_S=*(slabBx[3]+(NI-1)*k+i);
				By_S=*(slabBy[3]+(NI-1)*k+i);
				
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				Bx_S=*(domainBx+indexS);
				By_S=*(domainBy+indexS);
				
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				Bx_T=*(slabBx[4]+(NI-1)*j+i);
				By_T=*(slabBy[4]+(NI-1)*j+i);
				
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				Bx_T=*(domainBx+indexT);
				By_T=*(domainBy+indexT);
				
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				Bx_B=*(slabBx[5]+(NI-1)*j+i);
				By_B=*(slabBy[5]+(NI-1)*j+i);
				
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				Bx_B=*(domainBx+indexB);
				By_B=*(domainBy+indexB);
				
				vol_B=*(domainVolume+indexB);
			}


			Bxe=((Bx_E*vol_P+Bx_P*vol_E)/(vol_P+vol_E));
			Bxw=((Bx_W*vol_P+Bx_P*vol_W)/(vol_P+vol_W));
			Bxs=((Bx_S*vol_P+Bx_P*vol_S)/(vol_P+vol_S));
			Bxn=((Bx_N*vol_P+Bx_P*vol_N)/(vol_P+vol_N));
			Bxt=((Bx_T*vol_P+Bx_P*vol_T)/(vol_P+vol_T));
			Bxb=((Bx_B*vol_P+Bx_P*vol_B)/(vol_P+vol_B));

			Bye=((By_E*vol_P+By_P*vol_E)/(vol_P+vol_E));
			Byw=((By_W*vol_P+By_P*vol_W)/(vol_P+vol_W));
			Bys=((By_S*vol_P+By_P*vol_S)/(vol_P+vol_S));
			Byn=((By_N*vol_P+By_P*vol_N)/(vol_P+vol_N));
			Byt=((By_T*vol_P+By_P*vol_T)/(vol_P+vol_T));
			Byb=((By_B*vol_P+By_P*vol_B)/(vol_P+vol_B));
			
			

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxw*(*(SurfArea+Index+3*0+0))+Byw*(*(SurfArea+Index+3*0+1)));
				
				
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxe*(*(SurfArea+Index+3*1+0))+Bye*(*(SurfArea+Index+3*1+1)));
				
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//printf("\n east=%e,surf=%e\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face),*(SurfArea+Index+3*1+1));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxs*(*(SurfArea+Index+3*2+0))+Bys*(*(SurfArea+Index+3*2+1)));
				
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 3;  // north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxn*(*(SurfArea+Index+3*3+0))+Byn*(*(SurfArea+Index+3*3+1)));
				
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 4;  // bottom face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxb*(*(SurfArea+Index+3*4+0))+Byb*(*(SurfArea+Index+3*4+1)));

				//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 5;  // top face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (Bxt*(*(SurfArea+Index+3*5+0))+Byt*(*(SurfArea+Index+3*5+1)));

				//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	return 0;
}
