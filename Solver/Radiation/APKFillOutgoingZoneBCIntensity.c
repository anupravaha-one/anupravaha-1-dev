#include "APKFillOutgoingZoneBCIntensity.h"

/**This function put the boundary cell value of outgoing Intensity to boundary patch*/

int FillOutgoingZoneBCIntensity(Base *base, Zone *z)
	{
	 int i,j,k,l,m,ni,nj,nk,nt,np,direction,face, cellindex;
	 double *Ds;
	 char tag,name[33];
	 void *ptr;
	 ZoneBoundary *bcPtr;
	 
	 ni=*(z->Size+0);
	 nj=*(z->Size+1);
	 nk=*(z->Size+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 Ds=z->ZGeom->DmF;
	 	 
	 for(direction=0;direction<6;direction++)
	 	{
		
		ptr=(void *)FindInterfaceBCPtrMod(z,direction+1,name,&tag);
		
		if(tag=='b')
			{ 
			bcPtr=(ZoneBoundary *)ptr;
			
			switch(direction+1)
			{
			case 1://East Zone face
			i=ni-2;
			face=1;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				for(k=0;k<nk-1;k++)
				for(j=0;j<nj-1;j++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>=0)
					*(bcPtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*nt*np+np*l+m);
					}
			break;
			
			case 2: //West Zone Face
			i=0;
			face=0;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				for(k=0;k<nk-1;k++)
				for(j=0;j<nj-1;j++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>=0)
					*(bcPtr->ZSlab->RadIntensity+((nj-1)*k+j)*nt*np+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*nt*np+np*l+m);
					}
			break;
			
			case 3://North Zone face
			j=nj-2;
			face=3;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				for(k=0;k<nk-1;k++)
				for(i=0;i<ni-1;i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>=0)
					*(bcPtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*nt*np+np*l+m);
					}
			break;
			
			case 4://South Zone face
			j=0;
			face=2;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				for(k=0;k<nk-1;k++)
				for(i=0;i<ni-1;i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>=0)
					*(bcPtr->ZSlab->RadIntensity+((ni-1)*k+i)*nt*np+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*nt*np+np*l+m);
					}
			break;
			
			case 5://Top Zone face
			k=nk-2;
			face=5;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				for(j=0;j<nj-1;j++)
				for(i=0;i<ni-1;i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>=0)
					*(bcPtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*nt*np+np*l+m);
					}
			break;
			
			case 6://Bottom Zone face
			k=0;
			face=4;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				for(j=0;j<nj-1;j++)
				for(i=0;i<ni-1;i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>=0)
					*(bcPtr->ZSlab->RadIntensity+((ni-1)*j+i)*nt*np+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*nt*np+np*l+m);
					}
			break;
			}
			}
		}
return 0;
	}
