#include "APKCalcNusseltNo.h"
/**
The function CalcNusseltNo() is used to calculate the Nusselt no. at the zone boundaries. This function 
loops over all the zones and over all the zone boundaries and calculates the nusselt number. This function
can calculate Nusselt No like quantity for Scalars also. Hence a proper variableTag has to be supplied. This function is
not yet fully generalised.
*/

/*	Remember This function is for Calculating Nusselt No so send only vartag of Temperature 
	and Scalar. In the Datastructure there is space in both Interface and Boundary to store the 
	Average Nusselt No. of the respective face. also one file is created named
	SOLUTIONS/AverageNusseltNo.dat to check the value.
*/

// Test it with the Differentially heated 3D cube Test case

int CalcNusseltNo(Base *base,int vartag)
{
	FILE *s,*f,*v,*m;
	int iface,NI,NJ,NK,*ptrSize,nZone,INDEX,INDEX1;
	int i,j,k,count;
	char tag,name[33];
	void *ptr;
	double Temflux,Fdiff[6],WTemflux,T,P,C;
	double *slabVariable[6],*slabVolume[6],AreaSum;
	//double *deferCorrFlux,*upwindFlux,*diagConvect,*tempVar;
	double *vols,*beta,*variable,diffdia,*surf,Temfluxarea;
	double p,q,r,Area,*tempVar;
	int index[3];
	Zone *zone;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *FindZonePtr();
	s =fopen("SOLUTIONS/AverageNusseltNo.dat","a");
	if (s==NULL){ s =fopen("SOLUTIONS/AverageNusseltNo.dat","w"); }
	//v =fopen("SOLUTIONS/LocalNusseltNo.dat","a");
	//f =fopen("SOLUTIONS/WeitageNusseltNo.dat","w");
	//m =fopen("SOLUTIONS/FluxNusseltNo.dat","w");
	Temfluxarea=0.0;
	count = base->Number_of_zones;
	//for(nZone=0;nZone<base->Number_of_zones;nZone++)
	
	fprintf(s,"\n%e\t",base->BaseParameters->Time);
	      //{
		zone=(base->CgnsZone+1);
		//if(zone->flagSolidFluid==1)
		//{
			surf= zone->ZGeom->SurfaceArea;

			variable=MakeDomainArray(zone,vartag);//variableDomain stores variable pointer
			tempVar=MakeDomainArray(zone,vartag*10);//variableDomain stores variable pointer
			MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of u */
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			vols= zone->ZGeom->Volume;
			beta= zone->ZGeom->Beta;

			//fprintf(s,"\n~~~~~~~ZONENAME :%s ~~~~~~~~\n",(zone->zonename));
			//fprintf(f,"\n~~~~~~~ZONENAME :%s ~~~~~~~~\n",(zone->zonename));
			//fprintf(v,"\n~~~~~~~ZONENAME :%s ~~~~~~~~\n",(zone->zonename));
			//fprintf(m,"\n~~~~~~~ZONENAME :%s ~~~~~~~~\n",(zone->zonename));
			ptrSize=zone->Size;
			NI = *(ptrSize+0);
			NJ = *(ptrSize+1);
			NK = *(ptrSize+2);

			for(iface=0;iface<6;iface++)
			{
				ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
			
				if(tag=='b')
				{
					bcPtr=(ZoneBoundary *)ptr;
					
					switch(iface)
					{
						case 0:
						T=0.0;
					
						for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							{
								INDEX1 = (NJ-1)*k+j;
								T = T + 1.0/(*((bcPtr)->ZSlab->T+INDEX1));
                                
							}
						
							if(vartag==6)
							{
								(bcPtr->AvNusseltNoT) = (T)/((NJ-1)*(NK-1));
								fprintf(s,"%12.8e\t",(bcPtr->AvNusseltNoT));
							
							}
						
						break;
						case 1:
						T=0.0;
					
						for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							{
								INDEX1 = (NJ-1)*k+j;
								T = T + 1.0/(*((bcPtr)->ZSlab->T+INDEX1));
                                
							}
						
							if(vartag==6)
							{
								(bcPtr->AvNusseltNoT) = (T)/((NJ-1)*(NK-1));
								fprintf(s,"%12.8e\t",(bcPtr->AvNusseltNoT));
							
							}
						break;
						case 2:
						T=0.0;
					
						for(k=0;k<NK-1;k++)
						   for(i=0;i<NI-1;i++)
							{
								INDEX1 = (NI-1)*k+i;
								T = T + 1.0/(*((bcPtr)->ZSlab->T+INDEX1));
                                
							}
						
							if(vartag==6)
							{
								(bcPtr->AvNusseltNoT) = (T)/((NI-1)*(NK-1));
								fprintf(s,"%12.8e\t",(bcPtr->AvNusseltNoT));
							
							}
						break;
						case 3:
							T=0.0;
					
						for(k=0;k<NK-1;k++)
						   for(i=0;i<NI-1;i++)
							{
								INDEX1 = (NI-1)*k+i;
								T = T + 1.0/(*((bcPtr)->ZSlab->T+INDEX1));
                                
							}
						
							if(vartag==6)
							{
								(bcPtr->AvNusseltNoT) = (T)/((NI-1)*(NK-1));
								fprintf(s,"%12.8e\t",(bcPtr->AvNusseltNoT));
							
							}
						break;
						case 4:
							T=0.0;
					
						for(j=0;j<NJ-1;j++)
						   for(i=0;i<NI-1;i++)
							{
								INDEX1 = (NI-1)*j+i;
								T = T + 1.0/(*((bcPtr)->ZSlab->T+INDEX1));
                                
							}
						
							if(vartag==6)
							{
								(bcPtr->AvNusseltNoT) = (T)/((NI-1)*(NJ-1));
								fprintf(s,"%12.8e\t",(bcPtr->AvNusseltNoT));
							
							}
						break;
						case 5:
							T=0.0;
					
						for(j=0;j<NJ-1;j++)
						   for(i=0;i<NI-1;i++)
							{
								INDEX1 = (NI-1)*j+i;
								T = T + 1.0/(*((bcPtr)->ZSlab->T+INDEX1));
                                
							}
						
							if(vartag==6)
							{
								(bcPtr->AvNusseltNoT) = (T)/((NI-1)*(NJ-1));
								fprintf(s,"%12.8e\t",(bcPtr->AvNusseltNoT));
							
							}
						break;
					}
				}
			}
		//}
		
	//}
	fclose(s);
//	fclose(f);
//	fclose(v);
//	fclose(m);
	
	return (0);
}
