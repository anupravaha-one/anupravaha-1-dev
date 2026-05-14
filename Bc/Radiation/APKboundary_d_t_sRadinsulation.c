#include "APKboundary_d_t_sRadinsulation.h"
/** This function sums all the energy coming to boundary and returns to domains isotropically*/

int boundary_d_t_sRadinsulation(Base *base , int direction, Zone *z,ZoneBoundary *zBc, double *dummy)
	{
	 int ni,nj,nk,nt,np,i,j,k,l,m,cellindex,face;
	 double *Ds,*Int,sum,PI,nRef,*AreaPointer,Ax,Ay,Az,Area;
	 
	 PI=4.0*atan(1.0);
	 ni=*(z->Size+0);
	 nj=*(z->Size+1);
	 nk=*(z->Size+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 Ds=z->ZGeom->DmF;
	 AreaPointer=z->ZGeom->SurfaceArea; 
	 //printf("%d\n",direction);
	 	switch(direction)
			{
			case 1://East Face
				i=ni-2;
				face=1;
				nRef=*(zBc->ZSlab->RefIndex);
				Int=zBc->ZSlab->RadIntensity;
				for(k=0;k<(nk-1);k++)
					for(j=0;j<(nj-1);j++)
						{
						cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						Ax=*(AreaPointer+cellindex*6*3+3*face+0);
						Ay=*(AreaPointer+cellindex*6*3+3*face+1);
						Az=*(AreaPointer+cellindex*6*3+3*face+2);
						Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>0)
							 sum=sum+*(Int+((nj-1)*k+j)*nt*np+np*l+m)*(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m))/Area;
						
						if(base->BaseParameters->flagRadiativeTemperature==1)
						{
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
							*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=sum/PI;
						}
						if(base->BaseParameters->flagRadiativeTemperature==0)
						{
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
							*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=sum/PI;
						}
							
						}
			break;
			case 2://west face
			i=0;
			face=0;
			nRef=*(zBc->ZSlab->RefIndex);
			Int=zBc->ZSlab->RadIntensity;
			for(k=0;k<(nk-1);k++)
				for(j=0;j<(nj-1);j++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					Ax=*(AreaPointer+cellindex*6*3+3*face+0);
					Ay=*(AreaPointer+cellindex*6*3+3*face+1);
					Az=*(AreaPointer+cellindex*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					sum=0.0;
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>0)
						 sum=sum+*(Int+((nj-1)*k+j)*nt*np+np*l+m)*(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m))/Area;
						
					if(base->BaseParameters->flagRadiativeTemperature==1)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=sum/PI;
						}
					if(base->BaseParameters->flagRadiativeTemperature==0)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=sum/PI;
						}
							
						}
			break;
			
			case 3://North face
			j=nj-2;
			face=3;
			nRef=*(zBc->ZSlab->RefIndex);
			Int=zBc->ZSlab->RadIntensity;
			for(k=0;k<(nk-1);k++)
				for(i=0;i<(ni-1);i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					Ax=*(AreaPointer+cellindex*6*3+3*face+0);
					Ay=*(AreaPointer+cellindex*6*3+3*face+1);
					Az=*(AreaPointer+cellindex*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					sum=0.0;
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>0)
						 sum=sum+*(Int+((ni-1)*k+i)*nt*np+np*l+m)*(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m))/Area;
						
					if(base->BaseParameters->flagRadiativeTemperature==1)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=sum/PI;
						}
					if(base->BaseParameters->flagRadiativeTemperature==0)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=sum/PI;
						}
							
						}
			break;
			case 4://South face
			j=0;
			face=2;
			nRef=*(zBc->ZSlab->RefIndex);
			Int=zBc->ZSlab->RadIntensity;
			for(k=0;k<(nk-1);k++)
				for(i=0;i<(ni-1);i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					Ax=*(AreaPointer+cellindex*6*3+3*face+0);
					Ay=*(AreaPointer+cellindex*6*3+3*face+1);
					Az=*(AreaPointer+cellindex*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					sum=0.0;
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>0)
						 sum=sum+*(Int+((ni-1)*k+i)*nt*np+np*l+m)*(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m))/Area;
						
					if(base->BaseParameters->flagRadiativeTemperature==1)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=sum/PI;
						}
					if(base->BaseParameters->flagRadiativeTemperature==0)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=sum/PI;
						}
							
						}
			break;
			case 5://Top Face
			k=nk-2;
			face=5;
			nRef=*(zBc->ZSlab->RefIndex);
			Int=zBc->ZSlab->RadIntensity;
			for(j=0;j<(nj-1);j++)
				for(i=0;i<(ni-1);i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					Ax=*(AreaPointer+cellindex*6*3+3*face+0);
					Ay=*(AreaPointer+cellindex*6*3+3*face+1);
					Az=*(AreaPointer+cellindex*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					sum=0.0;
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>0)
						 sum=sum+*(Int+((ni-1)*j+i)*nt*np+np*l+m)*(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m))/Area;
						
					if(base->BaseParameters->flagRadiativeTemperature==1)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=sum/PI;
						}
					if(base->BaseParameters->flagRadiativeTemperature==0)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=sum/PI;
						}
							
						}
			break;
			case 6://Bottom Face
			k=0;
			face=4;
			nRef=*(zBc->ZSlab->RefIndex);
			Int=zBc->ZSlab->RadIntensity;
			for(j=0;j<(nj-1);j++)
				for(i=0;i<(ni-1);i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					Ax=*(AreaPointer+cellindex*6*3+3*face+0);
					Ay=*(AreaPointer+cellindex*6*3+3*face+1);
					Az=*(AreaPointer+cellindex*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					sum=0.0;
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)>0)
						 sum=sum+*(Int+((ni-1)*j+i)*nt*np+np*l+m)*(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m))/Area;
						
					if(base->BaseParameters->flagRadiativeTemperature==1)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=sum/PI;
						}
					if(base->BaseParameters->flagRadiativeTemperature==0)
					{
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
						*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=sum/PI;
						}
							
					}
			break;
			}
	return 0;
	}
