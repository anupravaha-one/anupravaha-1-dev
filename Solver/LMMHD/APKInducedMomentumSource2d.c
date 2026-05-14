#include"APKInducedMomentumSource2d.h"
double InducedMomentumSource2d(Zone *zone,int NI,int NJ,int NK,int index[3],double grad[2])
{
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainBx,*domainBy;
	double *slabVolume[6],*slabBx[6],*slabBy[6];	
	double *F;
	//double source[6][3],sumSource;
	//Zone* zone;
	//nZones=(base->Number_of_zones);
	
	//int NI,NJ,NK,*ptrSize;
	//double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double Bx[6],By[6];
	double Bx_P,Bx_E,Bx_W,Bx_S,Bx_N,Bx_T,Bx_B,By_P,By_E,By_W,By_S,By_N,By_T,By_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	//for(numZones=0;numZones<nZones;numZones++)
	//{
		//zone=(base->CgnsZone+numZones);
		
		//ptrSize=zone->Size;

// 		NI=*(ptrSize+0);
// 		NJ=*(ptrSize+1);
// 		NK=*(ptrSize+2);
// 		
		F = zone->PrevFlowsol->Magneticflux;
		
		//SurfArea=zone->ZGeom->SurfaceArea;
		
		domainBx=MakeDomainArray(zone,56);
		domainBy=MakeDomainArray(zone,57);
		
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabBx,56);
		MakeSlabPtrArray(zone,slabBy,57);
		
		MakeSlabPtrArray(zone,slabVolume,5);

// 		for(k=0;k<NK-1;k++)
// 		for(j=0;j<NJ-1;j++)
// 		for(i=0;i<NI-1;i++)
// 		{	
			
	
			k=index[0],j=index[1],i=index[2];
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
			indexE=indexP+1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
			indexW=indexP-1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
			indexN=indexP+(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
			indexS=indexP-(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexT=indexP+(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=indexP-(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;

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

		
			Bx[1]=((Bx_E*vol_P+Bx_P*vol_E)/(vol_P+vol_E));
			Bx[0]=((Bx_W*vol_P+Bx_P*vol_W)/(vol_P+vol_W));
			Bx[2]=((Bx_S*vol_P+Bx_P*vol_S)/(vol_P+vol_S));
			Bx[3]=((Bx_N*vol_P+Bx_P*vol_N)/(vol_P+vol_N));
			Bx[5]=((Bx_T*vol_P+Bx_P*vol_T)/(vol_P+vol_T));
			Bx[4]=((Bx_B*vol_P+Bx_P*vol_B)/(vol_P+vol_B));

			By[1]=((By_E*vol_P+By_P*vol_E)/(vol_P+vol_E));
			By[0]=((By_W*vol_P+By_P*vol_W)/(vol_P+vol_W));
			By[2]=((By_S*vol_P+By_P*vol_S)/(vol_P+vol_S));
			By[3]=((By_N*vol_P+By_P*vol_N)/(vol_P+vol_N));
			By[5]=((By_T*vol_P+By_P*vol_T)/(vol_P+vol_T));
			By[4]=((By_B*vol_P+By_P*vol_B)/(vol_P+vol_B));

			

			grad[0]=((*(F+indexP*6+0))*Bx[0])+((*(F+indexP*6+1))*Bx[1])+((*(F+indexP*6+2))*Bx[2])+((*(F+indexP*6+3	))*Bx[3])+((*(F+indexP*6+4))*Bx[4])+((*(F+indexP*6+5))*Bx[5]);

			grad[1]=((*(F+indexP*6+0))*By[0])+((*(F+indexP*6+1))*By[1])+((*(F+indexP*6+2))*By[2])+((*(F+indexP*6+3	))*By[3])+((*(F+indexP*6+4))*By[4])+((*(F+indexP*6+5))*By[5]);

			

			
			
			
		
/*End of looping over NI,NJ,NK*/
		
	//}/*End of looping over the zones */
// 	printf("sourceee111e=%e\n",source);
	//printf("\nsource returned=%e\n",sumSource);
	return (0);
	
}
