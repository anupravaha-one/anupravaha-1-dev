#include "APKUpdateRadOldZSlabs.h"

/** Intensity value in slabs either boundary or Interface is updated at each Iteration by this function*/

int UpdateRadOldZSlabs(Base *base, Zone *zone)
	{
	 int i,j,k,l,m,ni,nj,nk,nt,np,iface,*ptrSize;
	 char tag,name[33];
	 void *ptr;
	 ZoneBoundary *bcPtr;
	 Interface *IfacePtr;
	 
	 ptrSize=zone->Size;
	 ni=*(ptrSize+0);
	 nj=*(ptrSize+1);
	 nk=*(ptrSize+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 	
		for(iface=0;iface<6;iface++)
			{
			 ptr=FindInterfaceBCPtrMod(zone, iface+1, name, &tag);
			 
			 if(tag=='b')
			 	{
				 bcPtr=(ZoneBoundary *)ptr;
				 
				 switch(iface+1)
				 {
				  case 1:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(j=0;j<nj-1;j++)
					if((zone->Flowsol->RadIntensity)!=NULL && (bcPtr->ZSlab->RadIntensity)!=NULL)
							*(bcPtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m)=*(bcPtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m);
				 break;
				 
				 case 2:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(j=0;j<nj-1;j++)
					if((zone->Flowsol->RadIntensity)!=NULL && (bcPtr->ZSlab->RadIntensity)!=NULL)
							*(bcPtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m)=*(bcPtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m);
				 break;
				 
				 case 3:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (bcPtr->ZSlab->RadIntensity)!=NULL)
							*(bcPtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m)=*(bcPtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m);
				 break;
				 
				 case 4:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (bcPtr->ZSlab->RadIntensity)!=NULL)
							*(bcPtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m)=*(bcPtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m);
				 break;
				 
				 case 5:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(j=0;j<nj-1;j++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (bcPtr->ZSlab->RadIntensity)!=NULL)
							*(bcPtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m)=*(bcPtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m);
				 break;
				 
				 case 6:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(j=0;j<nj-1;j++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (bcPtr->ZSlab->RadIntensity)!=NULL)
							*(bcPtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m)=*(bcPtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m);
				 break;
				 }
			 	}
			if(tag=='i')
				{
				 IfacePtr=(Interface *)ptr;
				 
				 switch(iface+1)
				 {
				  case 1:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(j=0;j<nj-1;j++)
					if((zone->Flowsol->RadIntensity)!=NULL && (IfacePtr->ZSlab->RadIntensity)!=NULL)
							*(IfacePtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m)=*(IfacePtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m);
				 break;
				 
				 case 2:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(j=0;j<nj-1;j++)
					if((zone->Flowsol->RadIntensity)!=NULL && (IfacePtr->ZSlab->RadIntensity)!=NULL)
							*(IfacePtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m)=*(IfacePtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m);
				 break;
				 
				 case 3:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (IfacePtr->ZSlab->RadIntensity)!=NULL)
							*(IfacePtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m)=*(IfacePtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m);
				 break;
				 
				 case 4:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(k=0;k<nk-1;k++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (IfacePtr->ZSlab->RadIntensity)!=NULL)
							*(IfacePtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m)=*(IfacePtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m);
				 break;
				 
				 case 5:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(j=0;j<nj-1;j++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (IfacePtr->ZSlab->RadIntensity)!=NULL)
							*(IfacePtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m)=*(IfacePtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m);
				 break;
				 
				 case 6:
				  for(l=0;l<nt;l++)
				  	for(m=0;m<np;m++)
						for(j=0;j<nj-1;j++)
							for(i=0;i<ni-1;i++)
					if((zone->Flowsol->RadIntensity)!=NULL && (IfacePtr->ZSlab->RadIntensity)!=NULL)
							*(IfacePtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m)=*(IfacePtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m);
				 break;
				 }
				}
			}
	 return 0;
	}
