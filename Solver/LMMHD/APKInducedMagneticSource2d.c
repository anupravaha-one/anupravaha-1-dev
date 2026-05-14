#include"APKInducedMagneticSource2d.h"
double InducedMagneticSource2d(Zone *zone,int NI,int NJ,int NK,int index[3],double grad[2])
{
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainU,*domainV;
	double *slabVolume[6],*slabU[6],*slabV[6];	
	double *F;
	//double source[6][3],sumSource;
	//Zone* zone;
	//nZones=(base->Number_of_zones);
	
	//int NI,NJ,NK,*ptrSize;
	//double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u[6],v[6];
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	//for(numZones=0;numZones<nZones;numZones++)
	//{
		//zone=(base->CgnsZone+numZones);
		
		//ptrSize=zone->Size;

// 		NI=*(ptrSize+0);
// 		NJ=*(ptrSize+1);
// 		NK=*(ptrSize+2);
// 		
		F = zone->Flowsol->Magneticflux;
		
		//SurfArea=zone->ZGeom->SurfaceArea;
		
		domainU=MakeDomainArray(zone,1);
		domainV=MakeDomainArray(zone,2);
		
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabU,1);
		MakeSlabPtrArray(zone,slabV,2);
		
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

			u_P=*(domainU+indexP);
			v_P=*(domainV+indexP);
			
			vol_P=*(domainVolume+indexP);
			
			
			if(i==(NI-2))
			{
				u_E=*(slabU[0]+(NJ-1)*k+j);
				v_E=*(slabV[0]+(NJ-1)*k+j);
				
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
				
			}
			else
			{
				u_E=*(domainU+indexE);
				v_E=*(domainV+indexE);
				
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				u_W =*(slabU[1]+(NJ-1)*k+j);
				v_W =*(slabV[1]+(NJ-1)*k+j);
				
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				u_W=*(domainU+indexW);
				v_W=*(domainV+indexW);
				
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				u_N=*(slabU[2]+(NI-1)*k+i);
				v_N=*(slabV[2]+(NI-1)*k+i);
				
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				u_N=*(domainU+indexN);
				v_N=*(domainV+indexN);
				
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				u_S=*(slabU[3]+(NI-1)*k+i);
				v_S=*(slabV[3]+(NI-1)*k+i);
				
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				u_S=*(domainU+indexS);
				v_S=*(domainV+indexS);
				
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				u_T=*(slabU[4]+(NI-1)*j+i);
				v_T=*(slabV[4]+(NI-1)*j+i);
				
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				u_T=*(domainU+indexT);
				v_T=*(domainV+indexT);
				
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				u_B=*(slabU[5]+(NI-1)*j+i);
				v_B=*(slabV[5]+(NI-1)*j+i);
				
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				u_B=*(domainU+indexB);
				v_B=*(domainV+indexB);
				
				vol_B=*(domainVolume+indexB);
			}

		
			u[1]=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
			u[0]=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
			u[2]=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
			u[3]=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));
			u[5]=((u_T*vol_P+u_P*vol_T)/(vol_P+vol_T));
			u[4]=((u_B*vol_P+u_P*vol_B)/(vol_P+vol_B));

			v[1]=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
			v[0]=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
			v[2]=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
			v[3]=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));
			v[5]=((v_T*vol_P+v_P*vol_T)/(vol_P+vol_T));
			v[4]=((v_B*vol_P+v_P*vol_B)/(vol_P+vol_B));

			

			grad[0]=((*(F+indexP*6+0))*u[0])+((*(F+indexP*6+1))*u[1])+((*(F+indexP*6+2))*u[2])+((*(F+indexP*6+3	))*u[3])+((*(F+indexP*6+4))*u[4])+((*(F+indexP*6+5))*u[5]);
			

			grad[1]=((*(F+indexP*6+0))*v[0])+((*(F+indexP*6+1))*v[1])+((*(F+indexP*6+2))*v[2])+((*(F+indexP*6+3	))*v[3])+((*(F+indexP*6+4))*v[4])+((*(F+indexP*6+5))*v[5]);

// 			printf("\n value=%e\n",grad[1]);
// 			exit(0);

			
			
			
		
/*End of looping over NI,NJ,NK*/
		
	//}/*End of looping over the zones */
// 	printf("sourceee111e=%e\n",source);
	//printf("\nsource returned=%e\n",sumSource);
	return (0);
	
}

