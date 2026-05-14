#include "APKShearStress.h"

/**
This function is under development and has to be tested.
*/

int ShearStress(Base *base,int vartag)
{
	FILE *s;
	int iface,NI,NJ,NK,*ptrSize,nZone,INDEX;
	int i,j,k,count,dir;
	char tag,name[33];
	void *ptr;
	double Temflux,*F,Fdiff[6];
	double *slabVariable[6],*slabVolume[6];
	double *vols,*beta,*variable,diffdia,*surf,Kcond,Temfluxarea;
	double p,q,r,Area,*tempVar,Temp;
	int index[3];
	Zone *zone;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	
	s = fopen("SOLUTIONS/BoundaryShearStres.dat","w");
	Temfluxarea=0.0;
	count = base->Number_of_zones;
	for(nZone=0;nZone<base->Number_of_zones;nZone++)
	{
		zone=(base->CgnsZone+nZone);
		surf= zone->ZGeom->SurfaceArea;
		
		variable=MakeDomainArray(zone,vartag);//variableDomain stores variable pointer
		tempVar=MakeDomainArray(zone,vartag*10);//variableDomain stores variable pointer
		
		fprintf(s,"\n~~~~~~~ZONENAME :%s ~~~~~~~~\n",(zone->zonename));
		ptrSize=zone->Size;
		F = zone->Flowsol->Massflux;
		NI = *(ptrSize+0);
		NJ = *(ptrSize+1);
		NK = *(ptrSize+2);
		//printf("I clear the check\n");
		for(iface=0;iface<6;iface++)
		{
			ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
			if(tag=='i')
			{
				facePtr=(Interface *)ptr;
			}
			else if(tag=='b')
			{
				MakeSlabPtrArray(zone,slabVariable,vartag);/*Called for making slabs of u */
				MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
				bcPtr=(ZoneBoundary *)ptr;
				dir=bcPtr->Direction;
				vols= zone->ZGeom->Volume;
				beta= zone->ZGeom->Beta;
				switch(dir)
				{
					case 1:
						Temflux=0.0;
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							index[0]=k,index[1]=j,index[2]=NI-2;
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							i=NI-2;
							INDEX = ((NI-1)*(NJ-1)*k+(NI-1)*j+i);
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							p=*(surf+(INDEX)*6*3+3*1+0);
							q=*(surf+(INDEX)*6*3+3*1+1);
							r=*(surf+(INDEX)*6*3+3*1+2);
							Area=sqrt(p*p+q*q+r*r);
							Temflux=Temflux +Kcond*Fdiff[1];
							Temp=(Kcond*Fdiff[1])/Area;
							if(vartag==6)
							{
								//*((bcPtr->ZSlab->GradT)+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradT)+((NJ-1)*k+j))=Temp;
								fprintf(s,"EAST Imax %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
							else
							{
								//*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradScal[vartag-600])+((NJ-1)*k+j))=Temp;
								fprintf(s,"EAST Imax %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
						}
						if(vartag==6)
						{
							(bcPtr->SumNorGradT)=Temflux;
							fprintf(s,"\nEAST Imax \tShearStress=%12.8e\n",(bcPtr->SumNorGradT));
						}
						else
						{
							(bcPtr->SumNorGradScal[vartag-600])=Temflux;
							fprintf(s,"\nEAST Imax \tShearStress=%12.8e\n",(bcPtr->SumNorGradScal[vartag-600]));
						}
					break;
					case 2:
						Temflux=0.0;
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							index[0]=k,index[1]=j,index[2]=0;
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							i=0;
							INDEX = ((NI-1)*(NJ-1)*k+(NI-1)*j+i);
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							p=*(surf+(INDEX)*6*3+3*0+0);
							q=*(surf+(INDEX)*6*3+3*0+1);
							r=*(surf+(INDEX)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
							Temflux=Temflux +Kcond*Fdiff[0];
							Temp=(Kcond*Fdiff[0])/Area;
							if(vartag==6)
							{
								//*((bcPtr->ZSlab->GradT)+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradT)+((NJ-1)*k+j))=Temp;
								fprintf(s,"WEST Imin %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
							else
							{
								//*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradScal[vartag-600])+((NJ-1)*k+j))=Temp;
								fprintf(s,"WEST Imin %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
						}
						if(vartag==6)
						{
							(bcPtr->SumNorGradT)=Temflux;
							fprintf(s,"\nWEST Imin \tShearStress=%12.8e\n",(bcPtr->SumNorGradT));
						}
						else
						{
							(bcPtr->SumNorGradScal[vartag-600])=Temflux;
							fprintf(s,"\nWEST Imin \tShearStress=%12.8e\n",(bcPtr->SumNorGradScal[vartag-600]));
						}
					break;
					case 3:
						Temflux=0.0;
						for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=k,index[1]=NJ-2,index[2]=i;
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							j=NJ-2;
							INDEX = ((NI-1)*(NJ-1)*k+(NI-1)*j+i);
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							p=*(surf+(INDEX)*6*3+3*3+0);
							q=*(surf+(INDEX)*6*3+3*3+1);
							r=*(surf+(INDEX)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							Temflux=Temflux +Kcond*Fdiff[3];
							Temp=(Kcond*Fdiff[3])/Area;
							if(vartag==6)
							{
								//*((bcPtr->ZSlab->GradT)+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradT)+((NI-1)*k+i))=Temp;
								fprintf(s,"NORTH Jmax %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
							else
							{
								//*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*k+i))=Temp;
								fprintf(s,"NORTH Jmax %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
						}
						if(vartag==6)
						{
							(bcPtr->SumNorGradT)=Temflux;
							fprintf(s,"\nNORTH Jmax \tShearStress=%12.8e\n",(bcPtr->SumNorGradT));
						}
						else
						{
							(bcPtr->SumNorGradScal[vartag-600])=Temflux;
							fprintf(s,"\nNORTH Jmax \tShearStress=%12.8e\n",(bcPtr->SumNorGradScal[vartag-600]));
						}
					break;
					case 4:
						Temflux=0.0;
						for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=k,index[1]=0,index[2]=i;
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							j=0;
							INDEX = ((NI-1)*(NJ-1)*k+(NI-1)*j+i);
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							p=*(surf+(INDEX)*6*3+3*2+0);
							q=*(surf+(INDEX)*6*3+3*2+1);
							r=*(surf+(INDEX)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							Temflux=Temflux +Kcond*Fdiff[2];
							Temp=(Kcond*Fdiff[2])/Area;
							if(vartag==6)
							{
								//*((bcPtr->ZSlab->GradT)+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradT)+((NI-1)*k+i))=Temp;
								fprintf(s,"SOUTH Jmin %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
							else
							{
								//*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*k+i))=Temp;
								fprintf(s,"SOUTH Jmin %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
						}
						if(vartag==6)
						{
							(bcPtr->SumNorGradT)=Temflux;
							fprintf(s,"\nSOUTH Jmin \tShearStress=%12.8e\n",(bcPtr->SumNorGradT));
						}
						else
						{
							(bcPtr->SumNorGradScal[vartag-600])=Temflux;
							fprintf(s,"\nSOUTH Jmin \tShearStress=%12.8e\n",(bcPtr->SumNorGradScal[vartag-600]));
						}
					break;
					case 5:
						Temflux=0.0;
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=NK-2,index[1]=j,index[2]=i;
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							k=NK-2;
							INDEX = ((NI-1)*(NJ-1)*k+(NI-1)*j+i);
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							p=*(surf+(INDEX)*6*3+3*5+0);
							q=*(surf+(INDEX)*6*3+3*5+1);
							r=*(surf+(INDEX)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							Temflux=Temflux +Kcond*Fdiff[5];
							Temp=(Kcond*Fdiff[5])/Area;
							if(vartag==6)
							{
								//*((bcPtr->ZSlab->GradT)+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradT)+((NI-1)*j+i))=Temp;
								fprintf(s,"TOP Kmax %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
							else
							{
								//*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*j+i))=Temp;
								fprintf(s,"TOP Kmax %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
						}
						if(vartag==6)
						{
							(bcPtr->SumNorGradT)=Temflux;
							fprintf(s,"\nTOP Kmax \tShearStress=%12.8e\n",(bcPtr->SumNorGradT));
						}
						else
						{
							(bcPtr->SumNorGradScal[vartag-600])=Temflux;
							fprintf(s,"\nTOP Kmax \tShearStress=%12.8e\n",(bcPtr->SumNorGradScal[vartag-600]));
						}
					break;
					case 6:
						Temflux=0.0;
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							index[0]=0,index[1]=j,index[2]=i;
							diffdia =Diffusion(index,variable,vols,ptrSize,slabVariable,slabVolume,beta,Fdiff);
							k=0;
							INDEX = ((NI-1)*(NJ-1)*k+(NI-1)*j+i);
							Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
							p=*(surf+(INDEX)*6*3+3*4+0);
							q=*(surf+(INDEX)*6*3+3*4+1);
							r=*(surf+(INDEX)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							Temflux=Temflux +Kcond*Fdiff[4];
							Temp=(Kcond*Fdiff[4])/Area;
							if(vartag==6)
							{
								//*((bcPtr->ZSlab->GradT)+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradT)+((NI-1)*j+i))=Temp;
								fprintf(s,"BOTTOM Kmin %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
							else
							{
								//*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*(NJ-1)*k+(NK-1)*j+i))=Temp;
								*((bcPtr->ZSlab->GradScal[vartag-600])+((NI-1)*j+i))=Temp;
								fprintf(s,"BOTTOM Kmin %d %d %d\tShearStress=%12.8e\n",i,j,k,Temp);
							}
						}
						if(vartag==6)
						{
							(bcPtr->SumNorGradT)=Temflux;
							fprintf(s,"\nBOTTOM Kmin \tShearStress=%12.8e\n",(bcPtr->SumNorGradT));
						}
						else
						{
							(bcPtr->SumNorGradScal[vartag-600])=Temflux;
							fprintf(s,"\nBOTTOM Kmin \tShearStress=%12.8e\n",(bcPtr->SumNorGradScal[vartag-600]));
						}
					break;
				}
			}
		}
	}
	fclose(s);
	return (0);
}
