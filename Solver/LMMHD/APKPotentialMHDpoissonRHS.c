#include "APKPotentialMHDpoissonRHS.h"
/**
This function is used for calculating the RHS of the potentialMHDpoisson equation from the current converged
values of velocities. This function is used for 3D domains.See Eqn No() Thesis Report of vetcha naveen.This function is similar to volumefluxstar.c
*/
int PotentialMHDpoissonRHS(Base *base)
{
	int nZones,numZones;
	int k,j,i,face;
	int Index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUcrossBX,*domainUcrossBY,*domainUcrossBZ;
	double *slabVolume[6],*slabUcrossBX[6],*slabUcrossBY[6],*slabUcrossBZ[6];
	double *SurfArea,*F;
	Zone* zone;

	nZones=(base->Number_of_zones);
       
	
	int NI,NJ,NK,*ptrSize;
	double ucrossbxw,ucrossbyw,ucrossbzw,ucrossbxe,ucrossbye,ucrossbze,ucrossbxn,ucrossbyn,ucrossbzn,ucrossbxs,ucrossbys,ucrossbzs,ucrossbxt,ucrossbyt,ucrossbzt,ucrossbxb,ucrossbyb,ucrossbzb;
	double ucrossbx_P,ucrossbx_E,ucrossbx_W,ucrossbx_S,ucrossbx_N,ucrossbx_T,ucrossbx_B,ucrossby_P,ucrossby_E,ucrossby_W,ucrossby_S,ucrossby_N,ucrossby_T,ucrossby_B,ucrossbz_P,ucrossbz_E,ucrossbz_W,ucrossbz_S,ucrossbz_N,ucrossbz_T,ucrossbz_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	//printf("\n I am in PotentialMHDpoissonRHS.c before the for loop\n");
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->PotentialMHDpoissonRHS;

		SurfArea=zone->ZGeom->SurfaceArea;

		domainUcrossBX=MakeDomainArray(zone,22);
		domainUcrossBY=MakeDomainArray(zone,23);
		domainUcrossBZ=MakeDomainArray(zone,24);
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray_com(zone,slabUcrossBX,22);
		MakeSlabPtrArray_com(zone,slabUcrossBY,23);
		MakeSlabPtrArray_com(zone,slabUcrossBZ,24);
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
			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			
			ucrossbx_P=*(domainUcrossBX+indexP);
			ucrossby_P=*(domainUcrossBY+indexP);
			ucrossbz_P=*(domainUcrossBZ+indexP);
			vol_P=*(domainVolume+indexP);
			//fprintf(fp,"\n u=%lf\t v=%lf\t w=%lf",ucrossbx_P,ucrossby_P,ucrossbz_P);

			if(i==(NI-2))
			{
				ucrossbx_E=*(slabUcrossBX[0]+(NJ-1)*k+j);
				ucrossby_E=*(slabUcrossBY[0]+(NJ-1)*k+j);
				ucrossbz_E=*(slabUcrossBZ[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				ucrossbx_E=*(domainUcrossBX+indexE);
				ucrossby_E=*(domainUcrossBY+indexE);
				ucrossbz_E=*(domainUcrossBZ+indexE);
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				ucrossbx_W =*(slabUcrossBX[1]+(NJ-1)*k+j);
				ucrossby_W =*(slabUcrossBY[1]+(NJ-1)*k+j);
				ucrossbz_W =*(slabUcrossBZ[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				ucrossbx_W=*(domainUcrossBX+indexW);
				ucrossby_W=*(domainUcrossBY+indexW);
				ucrossbz_W=*(domainUcrossBZ+indexW);
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				ucrossbx_N=*(slabUcrossBX[2]+(NI-1)*k+i);
				ucrossby_N=*(slabUcrossBY[2]+(NI-1)*k+i);
				ucrossbz_N=*(slabUcrossBZ[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				ucrossbx_N=*(domainUcrossBX+indexN);
				ucrossby_N=*(domainUcrossBY+indexN);
				ucrossbz_N=*(domainUcrossBZ+indexN);
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				ucrossbx_S=*(slabUcrossBX[3]+(NI-1)*k+i);
				ucrossby_S=*(slabUcrossBY[3]+(NI-1)*k+i);
				ucrossbz_S=*(slabUcrossBZ[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				ucrossbx_S=*(domainUcrossBX+indexS);
				ucrossby_S=*(domainUcrossBY+indexS);
				ucrossbz_S=*(domainUcrossBZ+indexS);
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				ucrossbx_T=*(slabUcrossBX[4]+(NI-1)*j+i);
				ucrossby_T=*(slabUcrossBY[4]+(NI-1)*j+i);
				ucrossbz_T=*(slabUcrossBZ[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				ucrossbx_T=*(domainUcrossBX+indexT);
				ucrossby_T=*(domainUcrossBY+indexT);
				ucrossbz_T=*(domainUcrossBZ+indexT);
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				ucrossbx_B=*(slabUcrossBX[5]+(NI-1)*j+i);
				ucrossby_B=*(slabUcrossBY[5]+(NI-1)*j+i);
				ucrossbz_B=*(slabUcrossBZ[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				ucrossbx_B=*(domainUcrossBX+indexB);
				ucrossby_B=*(domainUcrossBY+indexB);
				ucrossbz_B=*(domainUcrossBZ+indexB);
				vol_B=*(domainVolume+indexB);
			}


			ucrossbxe=((ucrossbx_E*vol_P+ucrossbx_P*vol_E)/(vol_P+vol_E));
			ucrossbxw=((ucrossbx_W*vol_P+ucrossbx_P*vol_W)/(vol_P+vol_W));
			ucrossbxs=((ucrossbx_S*vol_P+ucrossbx_P*vol_S)/(vol_P+vol_S));
			ucrossbxn=((ucrossbx_N*vol_P+ucrossbx_P*vol_N)/(vol_P+vol_N));
			ucrossbxt=((ucrossbx_T*vol_P+ucrossbx_P*vol_T)/(vol_P+vol_T));
			ucrossbxb=((ucrossbx_B*vol_P+ucrossbx_P*vol_B)/(vol_P+vol_B));

			ucrossbye=((ucrossby_E*vol_P+ucrossby_P*vol_E)/(vol_P+vol_E));
			ucrossbyw=((ucrossby_W*vol_P+ucrossby_P*vol_W)/(vol_P+vol_W));
			ucrossbys=((ucrossby_S*vol_P+ucrossby_P*vol_S)/(vol_P+vol_S));
			ucrossbyn=((ucrossby_N*vol_P+ucrossby_P*vol_N)/(vol_P+vol_N));
			ucrossbyt=((ucrossby_T*vol_P+ucrossby_P*vol_T)/(vol_P+vol_T));
			ucrossbyb=((ucrossby_B*vol_P+ucrossby_P*vol_B)/(vol_P+vol_B));

			ucrossbze=((ucrossbz_E*vol_P+ucrossbz_P*vol_E)/(vol_P+vol_E));
			ucrossbzw=((ucrossbz_W*vol_P+ucrossbz_P*vol_W)/(vol_P+vol_W));
			ucrossbzs=((ucrossbz_S*vol_P+ucrossbz_P*vol_S)/(vol_P+vol_S));
			ucrossbzn=((ucrossbz_N*vol_P+ucrossbz_P*vol_N)/(vol_P+vol_N));
			ucrossbzt=((ucrossbz_T*vol_P+ucrossbz_P*vol_T)/(vol_P+vol_T));
			ucrossbzb=((ucrossbz_B*vol_P+ucrossbz_P*vol_B)/(vol_P+vol_B));

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxw*(*(SurfArea+Index+3*0+0))+ucrossbyw*(*(SurfArea+Index+3*0+1))+ucrossbzw*(*(SurfArea+Index+3*0+2));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxe*(*(SurfArea+Index+3*1+0))+ucrossbye*(*(SurfArea+Index+3*1+1))+ucrossbze*(*(SurfArea+Index+3*1+2));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxs*(*(SurfArea+Index+3*2+0))+ucrossbys*(*(SurfArea+Index+3*2+1))+ucrossbzs*(*(SurfArea+Index+3*2+2));
			face = 3;  // north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxn*(*(SurfArea+Index+3*3+0))+ucrossbyn*(*(SurfArea+Index+3*3+1))+ucrossbzn*(*(SurfArea+Index+3*3+2));
			face = 4;  // bottom face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxb*(*(SurfArea+Index+3*4+0))+ucrossbyb*(*(SurfArea+Index+3*4+1))+ucrossbzb*(*(SurfArea+Index+3*4+2));
			face = 5;  // top face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxt*(*(SurfArea+Index+3*5+0))+ucrossbyt*(*(SurfArea+Index+3*5+1))+ucrossbzt*(*(SurfArea+Index+3*5+2));
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	//printf("\n I am in PotentialMHDpoissonRHS.c after the for loop\n");
	//fclose(fp);
	return 0;
}

/**
This function is used for calculating the RHS of the potentialMHDpoisson equation from the current converged
values of velocities. This function is used for 2D domains.See Eqn No() Thesis Report of vetcha naveen.This function is similar to volumefluxstar2d.c
*/
int PotentialMHDpoissonRHS2d(Base *base)
{
	int nZones,numZones;
	int k,j,i,face;
	int Index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUcrossBX,*domainUcrossBY,*domainUcrossBZ;
	double *slabVolume[6],*slabUcrossBX[6],*slabUcrossBY[6],*slabUcrossBZ[6];
	double *SurfArea,*F;
	Zone* zone;

	//  FILE *fp;
	//	fp= fopen("intermediate.dat","w");
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	double ucrossbxw,ucrossbyw,ucrossbzw,ucrossbxe,ucrossbye,ucrossbze,ucrossbxn,ucrossbyn,ucrossbzn,ucrossbxs,ucrossbys,ucrossbzs,ucrossbxt,ucrossbyt,ucrossbzt,ucrossbxb,ucrossbyb,ucrossbzb;
	double ucrossbx_P,ucrossbx_E,ucrossbx_W,ucrossbx_S,ucrossbx_N,ucrossbx_T,ucrossbx_B,ucrossby_P,ucrossby_E,ucrossby_W,ucrossby_S,ucrossby_N,ucrossby_T,ucrossby_B,ucrossbz_P,ucrossbz_E,ucrossbz_W,ucrossbz_S,ucrossbz_N,ucrossbz_T,ucrossbz_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->PotentialMHDpoissonRHS;

		SurfArea=zone->ZGeom->SurfaceArea;

		domainUcrossBX=MakeDomainArray(zone,22);
		domainUcrossBY=MakeDomainArray(zone,23);
		domainUcrossBZ=MakeDomainArray(zone,24);
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray_com(zone,slabUcrossBX,22);
		MakeSlabPtrArray_com(zone,slabUcrossBY,23);
		MakeSlabPtrArray_com(zone,slabUcrossBZ,24);
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
			indexB=indexP-(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
	
			ucrossbx_P=*(domainUcrossBX+indexP);
			ucrossby_P=*(domainUcrossBY+indexP);
			ucrossbz_P=*(domainUcrossBZ+indexP);


			vol_P=*(domainVolume+indexP);
			//		fprintf(fp,"%lf \t %lf \t %lf\n",ucrossbx_P,ucrossby_P,ucrossbz_P);

			if(i==(NI-2))
			{
				ucrossbx_E=*(slabUcrossBX[0]+(NJ-1)*k+j);
				ucrossby_E=*(slabUcrossBY[0]+(NJ-1)*k+j);
				ucrossbz_E=*(slabUcrossBZ[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				ucrossbx_E=*(domainUcrossBX+indexE);
				ucrossby_E=*(domainUcrossBY+indexE);
				ucrossbz_E=*(domainUcrossBZ+indexE);
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				ucrossbx_W =*(slabUcrossBX[1]+(NJ-1)*k+j);
				ucrossby_W =*(slabUcrossBY[1]+(NJ-1)*k+j);
				ucrossbz_W =*(slabUcrossBZ[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				ucrossbx_W=*(domainUcrossBX+indexW);
				ucrossby_W=*(domainUcrossBY+indexW);
				ucrossbz_W=*(domainUcrossBZ+indexW);
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				ucrossbx_N=*(slabUcrossBX[2]+(NI-1)*k+i);
				ucrossby_N=*(slabUcrossBY[2]+(NI-1)*k+i);
				ucrossbz_N=*(slabUcrossBZ[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				ucrossbx_N=*(domainUcrossBX+indexN);
				ucrossby_N=*(domainUcrossBY+indexN);
				ucrossbz_N=*(domainUcrossBZ+indexN);
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				ucrossbx_S=*(slabUcrossBX[3]+(NI-1)*k+i);
				ucrossby_S=*(slabUcrossBY[3]+(NI-1)*k+i);
				ucrossbz_S=*(slabUcrossBZ[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				ucrossbx_S=*(domainUcrossBX+indexS);
				ucrossby_S=*(domainUcrossBY+indexS);
				ucrossbz_S=*(domainUcrossBZ+indexS);
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				ucrossbx_T=*(slabUcrossBX[4]+(NI-1)*j+i);
				ucrossby_T=*(slabUcrossBY[4]+(NI-1)*j+i);
				ucrossbz_T=*(slabUcrossBZ[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
        		}
			else
			{
				ucrossbx_T=*(domainUcrossBX+indexT);
				ucrossby_T=*(domainUcrossBY+indexT);
				ucrossbz_T=*(domainUcrossBZ+indexT);
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
			 	ucrossbx_B=*(slabUcrossBX[5]+(NI-1)*j+i);
			  	ucrossby_B=*(slabUcrossBY[5]+(NI-1)*j+i);
			        ucrossbz_B=*(slabUcrossBZ[5]+(NI-1)*j+i);
			  	vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				ucrossbx_B=*(domainUcrossBX+indexB);
				ucrossby_B=*(domainUcrossBY+indexB);
				ucrossbz_B=*(domainUcrossBZ+indexB);
				vol_B=*(domainVolume+indexB);
			}


			ucrossbxe=((ucrossbx_E*vol_P+ucrossbx_P*vol_E)/(vol_P+vol_E));
			ucrossbxw=((ucrossbx_W*vol_P+ucrossbx_P*vol_W)/(vol_P+vol_W));
			ucrossbxs=((ucrossbx_S*vol_P+ucrossbx_P*vol_S)/(vol_P+vol_S));
			ucrossbxn=((ucrossbx_N*vol_P+ucrossbx_P*vol_N)/(vol_P+vol_N));
			ucrossbxt=((ucrossbx_T*vol_P+ucrossbx_P*vol_T)/(vol_P+vol_T));
			ucrossbxb=((ucrossbx_B*vol_P+ucrossbx_P*vol_B)/(vol_P+vol_B));

			ucrossbye=((ucrossby_E*vol_P+ucrossby_P*vol_E)/(vol_P+vol_E));
			ucrossbyw=((ucrossby_W*vol_P+ucrossby_P*vol_W)/(vol_P+vol_W));
			ucrossbys=((ucrossby_S*vol_P+ucrossby_P*vol_S)/(vol_P+vol_S));
			ucrossbyn=((ucrossby_N*vol_P+ucrossby_P*vol_N)/(vol_P+vol_N));
			ucrossbyt=((ucrossby_T*vol_P+ucrossby_P*vol_T)/(vol_P+vol_T));
			ucrossbyb=((ucrossby_B*vol_P+ucrossby_P*vol_B)/(vol_P+vol_B));
			
			ucrossbze=((ucrossbz_E*vol_P+ucrossbz_P*vol_E)/(vol_P+vol_E));
			ucrossbzw=((ucrossbz_W*vol_P+ucrossbz_P*vol_W)/(vol_P+vol_W));
			ucrossbzs=((ucrossbz_S*vol_P+ucrossbz_P*vol_S)/(vol_P+vol_S));
			ucrossbzn=((ucrossbz_N*vol_P+ucrossbz_P*vol_N)/(vol_P+vol_N));
			ucrossbzt=((ucrossbz_T*vol_P+ucrossbz_P*vol_T)/(vol_P+vol_T));
			ucrossbzb=((ucrossbz_B*vol_P+ucrossbz_P*vol_B)/(vol_P+vol_B));

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxw*(*(SurfArea+Index+3*0+0))+ucrossbyw*(*(SurfArea+Index+3*0+1))+ucrossbzw*(*(SurfArea+Index+3*0+2));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxe*(*(SurfArea+Index+3*1+0))+ucrossbye*(*(SurfArea+Index+3*1+1))+ucrossbze*(*(SurfArea+Index+3*1+2));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxs*(*(SurfArea+Index+3*2+0))+ucrossbys*(*(SurfArea+Index+3*2+1))+ucrossbzs*(*(SurfArea+Index+3*2+2));
			face = 3;  // north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxn*(*(SurfArea+Index+3*3+0))+ucrossbyn*(*(SurfArea+Index+3*3+1))+ucrossbzn*(*(SurfArea+Index+3*3+2));
			face = 4;  // bottom face
	                	  *(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxb*(*(SurfArea+Index+3*4+0))+ucrossbyb*(*(SurfArea+Index+3*4+1))+ucrossbzb*(*(SurfArea+Index+3*4+2));       
		       
			face = 5;  // top face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = ucrossbxt*(*(SurfArea+Index+3*5+0))+ucrossbyt*(*(SurfArea+Index+3*5+1))+ucrossbzt*(*(SurfArea+Index+3*5+2));

				
                 
		            	
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */

	//	 fclose(fp);
	return 0;
}

