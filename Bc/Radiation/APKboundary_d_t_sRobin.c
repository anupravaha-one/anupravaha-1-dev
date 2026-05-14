#include "APKboundary_d_t_sRobin.h"
/** This function is used to implement the both convective and radiative flux specified boundary condition in conjuction with the fuction boundary_d_t_sInsulation()*/

int boundary_d_t_sRobin(Base *base, int direction, Zone *z, ZoneBoundary *zBc, Dataset *dset, double *dummy)
	{
	int i,j,k,l,m,ni,nj,nk,nt,np,*ptrsize,indexP,index[3],face;
	double p,q,r,Area,res,*vols,*beta,*surf,*slabVelocity[6],*slabVolume[6],Fdiff[6],sumN,sumO,*flowsol,*Ds;
	double diagonal,Diff[6],*slabPrevVar,Alpha,BetaCoeff,Gamma,sigma=5.67e-8,*varNew,*varOld,Residual,*Int,*SInt,sum;
	
	ni=*(z->Size+0);
	nj=*(z->Size+1);
	nk=*(z->Size+2); 
	nt=base->BaseParameters->ntheta;
	np=base->BaseParameters->nphi;
	ptrsize=z->Size;
			MakeSlabPtrArray(z,slabVelocity,6);/*Called for making slabs of Temperature */
			MakeSlabPtrArray(z,slabVolume,5);/*Called for making slabs of Volumes*/
			flowsol=z->Flowsol->T;
			Int=z->Flowsol->RadIntensity;
			SInt=zBc->ZSlab->RadIntensity;
			Ds=z->ZGeom->DmF;
		switch(direction)
		{
			case 1://east
				varOld=(double *)malloc((nj-1)*(nk-1)*sizeof(double));
				vols= z->ZGeom->Volume;
				Gamma=*(dset->Data_Array+0);
				Alpha=zBc->ZSlab->epsilon;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				slabPrevVar=zBc->ZSlab->TOld;
				varNew=slabPrevVar;
				
				do
				{
					slabVelocity[0]=varNew;
					for(k=0;k<nk-1;k++) 
					for(j=0;j<nj-1;j++)
						*(varOld+(nj-1)*k+j)=*(varNew+(nj-1)*k+j);
						//printf("Direction=%d\n",direction);
					for(k=0;k<nk-1;k++) 
					for(j=0;j<nj-1;j++)
					{
						i=ni-2;
						face=1;
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						p=*(surf+(indexP)*6*3+3*1+0);
						q=*(surf+(indexP)*6*3+3*1+1);
						r=*(surf+(indexP)*6*3+3*1+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=j,index[2]=ni-2;
						CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)>0)
							 sum=sum+*(SInt+((nj-1)*k+j)*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
							  //sum=sum+*(Int+indexP*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
						res=Gamma*Area-BetaCoeff*Fdiff[1]-Area*Alpha*sigma*pow(*(varOld+(nj-1)*k+j),4)+sum;//Diff;
						diagonal=BetaCoeff*Diff[1]+Alpha*Area*sigma*pow(*(varOld+(nj-1)*k+j),3);
						*(varNew+(nj-1)*k+j)=*(varOld+(nj-1)*k+j)+0.5*res/diagonal;
					}
					sumN=0.0;
					sumO=0.0;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
					sumN=sumN+*(varNew+(nj-1)*k+j);
					sumO=sumO+*(varOld+(nj-1)*k+j);
					}
					Residual=fabs(sumN-sumO)/fabs(sumN);
					//printf("Boundary Residual=%e\n",Residual);
				}while(Residual>1.0e-5);
				for(k=0;k<nk-1;k++) 
					for(j=0;j<nj-1;j++)
					*(dummy+((nj-1)*k+j))=*(varNew+(nj-1)*k+j);
					free(varOld);
		break;
			
			case 2://west
				varOld=(double *)malloc((nj-1)*(nk-1)*sizeof(double));
				vols= z->ZGeom->Volume;
				Gamma=*(dset->Data_Array+0);
				Alpha=zBc->ZSlab->epsilon;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				slabPrevVar=zBc->ZSlab->TOld;
				varNew=slabPrevVar;
				do
				{
					for(k=0;k<nk-1;k++) 
					for(j=0;j<nj-1;j++)
						*(varOld+(nj-1)*k+j)=*(varNew+(nj-1)*k+j);
						slabVelocity[1]=varNew;
						//printf("Direction=%d\n",direction);
					for(k=0;k<nk-1;k++) 
					for(j=0;j<nj-1;j++)
					{
						i=0;
						face=0;
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						p=*(surf+(indexP)*6*3+3*0+0);
						q=*(surf+(indexP)*6*3+3*0+1);
						r=*(surf+(indexP)*6*3+3*0+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=j,index[2]=0;
						CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)>0)
							 sum=sum+*(SInt+((nj-1)*k+j)*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
							 // sum=sum+*(Int+indexP*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
						res=Gamma*Area-BetaCoeff*Fdiff[0]-Area*Alpha*sigma*pow(*(varOld+(nj-1)*k+j),4)+sum;//Diff;
						diagonal=-BetaCoeff*Diff[0]+Alpha*Area*sigma*pow(*(varOld+(nj-1)*k+j),3);
						*(varNew+(nj-1)*k+j)=*(varOld+(nj-1)*k+j)+0.5*res/diagonal;
					}
					sumN=0.0;
					sumO=0.0;
					for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{
					sumN=sumN+*(varNew+(nj-1)*k+j);
					sumO=sumO+*(varOld+(nj-1)*k+j);
					}
					Residual=fabs(sumN-sumO)/fabs(sumN);
					//printf("Boundary Residual=%e\n",Residual);
				}while(Residual>1.0e-5);
				for(k=0;k<nk-1;k++) 
					for(j=0;j<nj-1;j++)
					*(dummy+((nj-1)*k+j))=*(varNew+(nj-1)*k+j);
					free(varOld);
				
				
			break;
			
			case 3://north
			//printf("Direction=%d\n",direction);
				varOld=(double *)malloc((ni-1)*(nk-1)*sizeof(double));
				vols= z->ZGeom->Volume;
				Gamma=*(dset->Data_Array+0);
				Alpha=zBc->ZSlab->epsilon;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				slabPrevVar=zBc->ZSlab->TOld;
				varNew=slabPrevVar;
				do
				{	slabVelocity[2]=varNew;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
						*(varOld+(ni-1)*k+i)=*(varNew+(ni-1)*k+i);
						
						
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						j=nj-2;
						face=3;
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						p=*(surf+(indexP)*6*3+3*3+0);
						q=*(surf+(indexP)*6*3+3*3+1);
						r=*(surf+(indexP)*6*3+3*3+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=nj-2,index[2]=i;
						CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)>0)
							  sum=sum+*(SInt+((ni-1)*k+i)*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
							  //sum=sum+*(Int+indexP*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
						res=Gamma*Area-BetaCoeff*Fdiff[3]-Area*Alpha*sigma*pow(*(varOld+(ni-1)*k+i),4)+sum;//Diff;
						diagonal=BetaCoeff*Diff[3]+Alpha*Area*sigma*pow(*(varOld+(ni-1)*k+i),3);
						*(varNew+(ni-1)*k+i)=*(varOld+(ni-1)*k+i)+0.5*res/diagonal;
					}
					sumN=0.0;
					sumO=0.0;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
					sumN=sumN+*(varNew+(ni-1)*k+i);
					sumO=sumO+*(varOld+(ni-1)*k+i);
					}
					Residual=fabs(sumN-sumO)/fabs(sumN);
					//printf("Boundary Residual=%e\n",Residual);
				}while(Residual>1.0e-5);
				for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					*(dummy+((ni-1)*k+i))=*(varNew+(ni-1)*k+i);
					free(varOld);
			break;
			
			case 4://south	
				varOld=(double *)malloc((ni-1)*(nk-1)*sizeof(double));
				vols= z->ZGeom->Volume;
				Gamma=*(dset->Data_Array+0);
				Alpha=zBc->ZSlab->epsilon;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				slabPrevVar=zBc->ZSlab->TOld;
				varNew=slabPrevVar;
				do
				{
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
						*(varOld+(ni-1)*k+i)=*(varNew+(ni-1)*k+i);
						slabVelocity[3]=varNew;
						//printf("Direction=%d\n",direction);
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
						j=0;
						face=2;
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						p=*(surf+(indexP)*6*3+3*2+0);
						q=*(surf+(indexP)*6*3+3*2+1);
						r=*(surf+(indexP)*6*3+3*2+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=k,index[1]=0,index[2]=i;
						CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)>0)
							 sum=sum+*(SInt+((ni-1)*k+i)*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
							 //sum=sum+*(Int+indexP*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
						res=Gamma*Area-BetaCoeff*Fdiff[2]-Area*Alpha*sigma*pow(*(varOld+(ni-1)*k+i),4)+sum;//Diff;
						diagonal=-BetaCoeff*Diff[2]+Alpha*Area*sigma*pow(*(varOld+(ni-1)*k+i),3);
						*(varNew+(ni-1)*k+i)=*(varOld+(ni-1)*k+i)+0.5*res/diagonal;
					}
					sumN=0.0;
					sumO=0.0;
					for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					{
					sumN=sumN+*(varNew+(ni-1)*k+i);
					sumO=sumO+*(varOld+(ni-1)*k+i);
					}
					Residual=fabs(sumN-sumO)/fabs(sumN);
					//printf("Boundary Residual=%e\n",Residual);
				}while(Residual>1.0e-5);
				for(k=0;k<nk-1;k++) 
					for(i=0;i<ni-1;i++)
					*(dummy+((ni-1)*k+i))=*(varNew+(ni-1)*k+i);
					free(varOld);
					
			break;
			
			case 5://top
				varOld=(double *)malloc((ni-1)*(nj-1)*sizeof(double));
				vols= z->ZGeom->Volume;
				Gamma=*(dset->Data_Array+0);
				Alpha=zBc->ZSlab->epsilon;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				slabPrevVar=zBc->ZSlab->TOld;
				varNew=slabPrevVar;
				do
				{
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
						*(varOld+(ni-1)*j+i)=*(varNew+(ni-1)*j+i);
						slabVelocity[4]=varNew;
						//printf("Direction=%d\n",direction);
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						k=nk-2;
						face=5;
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						p=*(surf+(indexP)*6*3+3*5+0);
						q=*(surf+(indexP)*6*3+3*5+1);
						r=*(surf+(indexP)*6*3+3*5+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=nk-2,index[1]=j,index[2]=i;
						CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)>0)
							  sum=sum+*(SInt+((ni-1)*j+i)*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
							  //sum=sum+*(Int+indexP*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
						res=Gamma*Area-BetaCoeff*Fdiff[5]-Area*Alpha*sigma*pow(*(varOld+(ni-1)*j+i),4)+sum;//Diff;
						diagonal=BetaCoeff*Diff[5]+Alpha*Area*sigma*pow(*(varOld+(ni-1)*j+i),3);
						*(varNew+(ni-1)*j+i)=*(varOld+(ni-1)*j+i)+0.5*res/diagonal;
					}
					sumN=0.0;
					sumO=0.0;
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
					sumN=sumN+*(varNew+(ni-1)*j+i);
					sumO=sumO+*(varOld+(ni-1)*j+i);
					}
					Residual=fabs(sumN-sumO)/fabs(sumN);
					//printf("Boundary Residual=%e\n",Residual);
				}while(Residual>1.0e-5);
				for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					*(dummy+((ni-1)*j+i))=*(varNew+(ni-1)*j+i);
					free(varOld);
				
			break;
			
			case 6://bottom
				varOld=(double *)malloc((ni-1)*(nj-1)*sizeof(double));
				vols= z->ZGeom->Volume;
				Gamma=*(dset->Data_Array+0);
				Alpha=zBc->ZSlab->epsilon;
				beta= z->ZGeom->Beta;
				surf= z->ZGeom->SurfaceArea;
				slabPrevVar=zBc->ZSlab->TOld;
				varNew=slabPrevVar;
				do
				{
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
						*(varOld+(ni-1)*j+i)=*(varNew+(ni-1)*j+i);
						slabVelocity[5]=varNew;
						//printf("Direction=%d\n",direction);
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
						k=0;
						face=4;
						indexP=(ni-1)*(nj-1)*k+(ni-1)*j+i;
						BetaCoeff=CalcDiffusionCoefCommon(base,z,6,indexP);
						p=*(surf+(indexP)*6*3+3*4+0);
						q=*(surf+(indexP)*6*3+3*4+1);
						r=*(surf+(indexP)*6*3+3*4+2);
						Area=sqrt(p*p+q*q+r*r);
						index[0]=0,index[1]=j,index[2]=i;
						CalGradRad(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff,Diff);
						sum=0.0;
						for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							 if(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m)>0)
							 sum=sum+*(SInt+((ni-1)*j+i)*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
							  //sum=sum+*(Int+indexP*nt*np+np*l+m)*(*(Ds+indexP*6*nt*np+np*nt*face+np*l+m));
						res=Gamma*Area-BetaCoeff*Fdiff[4]-Area*Alpha*sigma*pow(*(varOld+(ni-1)*j+i),4)+sum;//Diff;
						diagonal=-BetaCoeff*Diff[4]+Alpha*Area*sigma*pow(*(varOld+(ni-1)*j+i),3);
						*(varNew+(ni-1)*j+i)=*(varOld+(ni-1)*j+i)+0.5*res/diagonal;
					}
					sumN=0.0;
					sumO=0.0;
					for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{
					sumN=sumN+*(varNew+(ni-1)*j+i);
					sumO=sumO+*(varOld+(ni-1)*j+i);
					}
					Residual=fabs(sumN-sumO)/fabs(sumN);
					//printf("Boundary Residual=%e\n",Residual);
				}while(Residual>1.0e-5);
				for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					*(dummy+((ni-1)*j+i))=*(varNew+(ni-1)*j+i);
					free(varOld);
				
			break;
		}	
	return 0;
}
