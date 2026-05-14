#include "APKboundary_d_t_sInsulation.h"

/** This function is used to implement the both convective + radiative flux specified boundary condition in conjuction with the fuction boundary_d_t_sRobin()*/

int boundary_d_t_sInsulation(Base *base , int direction, Zone *z,ZoneBoundary *zBc, Dataset *dset, double *dummy)
	{
	 int ni,nj,nk,nt,np,i,j,k,l,m,indexP,face,*ptrsize,index[3];
	 double *Ds,eps,PI,nRef,*Temp,*AreaPointer,*vols,*value,*beta,Area,sigma=5.67e-8;
	 double Ax,Ay,Az;//*flowsol,*slabVelocity[6],*slabVolume[6],Fdiff[6],Diff[6],BetaCoeff;
	 
	 PI=4.0*atan(1.0);
	 ni=*(z->Size+0);
	 nj=*(z->Size+1);
	 nk=*(z->Size+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 Ds=z->ZGeom->DmF;
	 AreaPointer=z->ZGeom->SurfaceArea; 
	 ptrsize=z->Size;
			//MakeSlabPtrArray(z,slabVelocity,6);/*Called for making slabs of Temperature */
			//MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/
			//flowsol=z->Flowsol->T;
	 	switch(direction)
			{
			case 1://East Face
				i=ni-2;
				face=1;
				vols= z->ZGeom->Volume;
				eps=zBc->ZSlab->epsilon;
				nRef=*(zBc->ZSlab->RefIndex);
				value=(zBc->Data_Set+4)->Data_Array;
				beta= z->ZGeom->Beta;
				//printf("\ndirection=%d\tvalue=%e\n",direction,*value);
				Temp=zBc->ZSlab->T;
				for(k=0;k<(nk-1);k++)
					for(j=0;j<(nj-1);j++)
						{
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						//BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						Ax=*(AreaPointer+indexP*6*3+3*face+0);
						Ay=*(AreaPointer+indexP*6*3+3*face+1);
						Az=*(AreaPointer+indexP*6*3+3*face+2);
						Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
						index[0]=k,index[1]=j,index[2]=ni-2;
						//CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)<0)
							 {
							 *(dummy+((nj-1)*k+j)*nt*np+np*l+m)=sigma*pow(*(Temp+(nj-1)*k+j),4)/PI;
							//*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=(*(value+0)*Area-BetaCoeff*Fdiff[1])/(PI*Area);
							}
						}
			break;
			
			case 2://west face
			i=0;
			face=0;
			vols= z->ZGeom->Volume;
			eps=zBc->ZSlab->epsilon;
			nRef=*(zBc->ZSlab->RefIndex);
			value=(zBc->Data_Set+4)->Data_Array;
			beta= z->ZGeom->Beta;
			//printf("\ndirection=%d\tvalue=%e\n",direction,*value);
			Temp=zBc->ZSlab->T;
			for(k=0;k<(nk-1);k++)
				for(j=0;j<(nj-1);j++)
					{
					indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					//BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
					Ax=*(AreaPointer+indexP*6*3+3*face+0);
					Ay=*(AreaPointer+indexP*6*3+3*face+1);
					Az=*(AreaPointer+indexP*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					index[0]=k,index[1]=j,index[2]=0;
					//CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)<0)
						 {
						 *(dummy+((nj-1)*k+j)*nt*np+np*l+m)=sigma*pow(*(Temp+(nj-1)*k+j),4)/PI;
						 //*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=(*(value+0)*Area-BetaCoeff*Fdiff[0])/(PI*Area);
						 }
					}
			break;
			
			case 3://North face
			j=nj-2;
			face=3;
			vols= z->ZGeom->Volume;
			eps=zBc->ZSlab->epsilon;
			nRef=*(zBc->ZSlab->RefIndex);
			value=(zBc->Data_Set+4)->Data_Array;
			beta= z->ZGeom->Beta;
			//printf("\ndirection=%d\tvalue=%e\n",direction,*value);
			Temp=zBc->ZSlab->T;
			for(k=0;k<(nk-1);k++)
				for(i=0;i<(ni-1);i++)
					{
					indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					//BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
					Ax=*(AreaPointer+indexP*6*3+3*face+0);
					Ay=*(AreaPointer+indexP*6*3+3*face+1);
					Az=*(AreaPointer+indexP*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					index[0]=k,index[1]=nj-2,index[2]=i;
					//CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)<0)
						 {
						 *(dummy+((ni-1)*k+i)*nt*np+np*l+m)=sigma*pow(*(Temp+((ni-1)*k+i)),4)/PI;
						//*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=(*(value+0)*Area-BetaCoeff*Fdiff[3])/(PI*Area);
						}	
					}
			break;
			
			case 4://South face
			j=0;
			face=2;
			vols= z->ZGeom->Volume;
			eps=zBc->ZSlab->epsilon;
			nRef=*(zBc->ZSlab->RefIndex);
			value=(zBc->Data_Set+4)->Data_Array;
			beta= z->ZGeom->Beta;
			//printf("\ndirection=%d\tvalue=%e\n",direction,*value);
			Temp=zBc->ZSlab->T;
			for(k=0;k<(nk-1);k++)
				for(i=0;i<(ni-1);i++)
					{
					indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					//BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
					Ax=*(AreaPointer+indexP*6*3+3*face+0);
					Ay=*(AreaPointer+indexP*6*3+3*face+1);
					Az=*(AreaPointer+indexP*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					index[0]=k,index[1]=0,index[2]=i;
					//CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)<0)
						 {
						 *(dummy+((ni-1)*k+i)*nt*np+np*l+m)=sigma*pow(*(Temp+((ni-1)*k+i)),4)/PI;
						 //*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=(*(value+0)*Area-BetaCoeff*Fdiff[2])/(PI*Area);
						 }
					}
			break;
			case 5://Top Face
			k=nk-2;
			face=5;
			vols= z->ZGeom->Volume;
			eps=zBc->ZSlab->epsilon;
			nRef=*(zBc->ZSlab->RefIndex);
			value=(zBc->Data_Set+4)->Data_Array;
			beta= z->ZGeom->Beta;
			//printf("\ndirection=%d\tvalue=%e\n",direction,*value);
			Temp=zBc->ZSlab->T;
			for(j=0;j<(nj-1);j++)
				for(i=0;i<(ni-1);i++)
					{
					indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					//BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
					Ax=*(AreaPointer+indexP*6*3+3*face+0);
					Ay=*(AreaPointer+indexP*6*3+3*face+1);
					Az=*(AreaPointer+indexP*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					index[0]=nk-2,index[1]=j,index[2]=i;
					//CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)<0)
						 {
						 *(dummy+((ni-1)*j+i)*nt*np+np*l+m)=sigma*pow(*(Temp+((ni-1)*j+i)),4)/PI;
						 //*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=(*(value+0)*Area-BetaCoeff*Fdiff[5])/(PI*Area);
						 }
					}
			break;
			case 6://Bottom Face
			k=0;
			face=4;
			vols= z->ZGeom->Volume;
			eps=zBc->ZSlab->epsilon;
			nRef=*(zBc->ZSlab->RefIndex);
			value=(zBc->Data_Set+4)->Data_Array;
			beta= z->ZGeom->Beta;
			//printf("\ndirection=%d\tvalue=%e\n",direction,*value);
			Temp=zBc->ZSlab->T;
			for(j=0;j<(nj-1);j++)
				for(i=0;i<(ni-1);i++)
					{
					indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					//BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
					Ax=*(AreaPointer+indexP*6*3+3*face+0);
					Ay=*(AreaPointer+indexP*6*3+3*face+1);
					Az=*(AreaPointer+indexP*6*3+3*face+2);
					Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
					index[0]=0,index[1]=j,index[2]=i;
					//CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
					for(l=0;l<nt;l++)
						for(m=0;m<np;m++)
						 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)<0)
						 {
						 *(dummy+((ni-1)*j+i)*nt*np+np*l+m)=sigma*pow(*(Temp+((ni-1)*j+i)),4)/PI;
						 //*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=(*(value+0)*Area-BetaCoeff*Fdiff[4])/(PI*Area);
						 }
					}
			break;
			}
	return 0;
	}
