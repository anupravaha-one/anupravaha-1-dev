#include "APKboundary_d_t_sRadsymmetry.h"

/** This function is used to apply the symmetric boundary condition. First it find the complementary direction and then equate the complementary
Intensity value to boundary condition*/

int boundary_d_t_sRadsymmetry(Base *base,int direction,Zone *z,double *dummy)
{
	//FILE *fp1,*fp2;
	 int i,j,k,l,m,li,mi,ni,nj,nk,nt,np,face,cellindex,lc,mc,lf,mf;
	 double *Ds, *NX, *NY, *NZ, *X, *Y, *Z, Xi, Yi, Zi,*Th,*Ph,dtheta,dphi,Pi,ldi,mdi,*CompInt;
	 double Omega1,Omega2,Omega3,Omega4,dthetaC,dphiC,dthetaF,dphiF,thetaC,thetaF,phiM;
	  
	 ni=*(z->Size+0);
	 nj=*(z->Size+1);
	 nk=*(z->Size+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 
	 Pi=4.0*atan(1);
	 
	 dtheta=Pi/nt;
	 dphi=2.0*Pi/np;
	 
	 Th=(double *)malloc(nt*sizeof(double));
	 Ph=(double *)malloc(np*sizeof(double));
	 
	 Ds=z->ZGeom->DmF;
	 
	 NX=z->ZGeom->NormalX;
	 NY=z->ZGeom->NormalY;
	 NZ=z->ZGeom->NormalZ;
	 
	 X=z->ZGeom->VectorX;
	 Y=z->ZGeom->VectorY;
	 Z=z->ZGeom->VectorZ;
	// printf("%d\n",direction);
		
		CompInt=(double *)malloc(nt*np*sizeof(double));
		
		//fp1=fopen("SOLUTIONS/RADIATION/Check.dat","w");
		//fp2=fopen("SOLUTIONS/RADIATION/SymCh.dat","w");
		
	switch(direction)
	{
		case 1://East Zone face
		i=ni-2;
		face=1;
		//printf("\n\tSym-direction=%d",direction);
		for(k=0;k<nk-1;k++)
		for(j=0;j<nj-1;j++)
			{
			cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
				{
				//printf("l=%d\tm=%d\t",l,m);
				for(li=0;li<nt;li++)
				for(mi=0;mi<np;mi++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*li+mi)>0)
				{
				Xi=*(X+np*li+mi)-(*(X+np*l+m)-2.0*(*(NX+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Yi=*(Y+np*li+mi)-(*(Y+np*l+m)-2.0*(*(NY+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Zi=*(Z+np*li+mi)-(*(Z+np*l+m)-2.0*(*(NZ+cellindex*6*nt*np+np*nt*face+np*l+m)));
				//printf("Xi=%e\tYi=%e\tZi=%e\n",Xi,Yi,Zi);
				if((fabs(Xi)<1.0e-8)&&(fabs(Yi)<1.0e-8)&&(fabs(Zi)<1.0e-8))
					{
					*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=*(dummy+((nj-1)*k+j)*nt*np+np*li+mi);
					//printf("li=%d\tmi=%d\n",li,mi);
					}
				}
				}
			}
			break;
			
		case 2: //West Zone Face
		i=0;
		face=0;
		//printf("\n\tSym-direction=%d",direction);
		for(k=0;k<nk-1;k++)
		for(j=0;j<nj-1;j++)
			{
			cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
				{
				//printf("l=%d\tm=%d\t",l,m);
				for(li=0;li<nt;li++)
				for(mi=0;mi<np;mi++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*li+mi)>0)
				{
				Xi=*(X+np*li+mi)-(*(X+np*l+m)-2.0*(*(NX+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Yi=*(Y+np*li+mi)-(*(Y+np*l+m)-2.0*(*(NY+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Zi=*(Z+np*li+mi)-(*(Z+np*l+m)-2.0*(*(NZ+cellindex*6*nt*np+np*nt*face+np*l+m)));
				//printf("Xi=%e\tYi=%e\tZi=%e\n",Xi,Yi,Zi);
				if((fabs(Xi)<1.0e-8)&&(fabs(Yi)<1.0e-8)&&(fabs(Zi)<1.0e-8))
					{
					*(dummy+((nj-1)*k+j)*nt*np+np*l+m)=*(dummy+((nj-1)*k+j)*nt*np+np*li+mi);
					//printf("li=%d\tmi=%d\n",li,mi);
					}
				}
				}
			}
			break;
			
		case 3://North Zone face
		j=nj-2;
		face=3;
		//printf("\n\tSym-direction=%d",direction);
		for(k=0;k<nk-1;k++)
		for(i=0;i<ni-1;i++)
			{
			cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
				{
				//printf("l=%d\tm=%d\t",l,m);
				for(li=0;li<nt;li++)
				for(mi=0;mi<np;mi++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*li+mi)>0)
				{
				Xi=*(X+np*li+mi)-(*(X+np*l+m)-2.0*(*(NX+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Yi=*(Y+np*li+mi)-(*(Y+np*l+m)-2.0*(*(NY+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Zi=*(Z+np*li+mi)-(*(Z+np*l+m)-2.0*(*(NZ+cellindex*6*nt*np+np*nt*face+np*l+m)));
				//printf("Xi=%e\tYi=%e\tZi=%e\n",fabs(Xi),fabs(Yi),fabs(Zi));
				if((fabs(Xi)<1.0e-8)&&(fabs(Yi)<1.0e-8)&&(fabs(Zi)<1.0e-8))
					{
					*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=*(dummy+((ni-1)*k+i)*nt*np+np*li+mi);
					//printf("li=%d\tmi=%d\n",li,mi);
					}
				}
				}
			}
			break;
			
		case 4://South Zone face
		j=0;
		face=2;
		//printf("\n\tSym-direction=%d",direction);
		for(k=0;k<nk-1;k++)
		for(i=0;i<ni-1;i++)
			{
			cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
				{
				//printf("l=%d\tm=%d\t",l,m);
				for(li=0;li<nt;li++)
				for(mi=0;mi<np;mi++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*li+mi)>0)
				{
				Xi=*(X+np*li+mi)-(*(X+np*l+m)-2.0*(*(NX+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Yi=*(Y+np*li+mi)-(*(Y+np*l+m)-2.0*(*(NY+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Zi=*(Z+np*li+mi)-(*(Z+np*l+m)-2.0*(*(NZ+cellindex*6*nt*np+np*nt*face+np*l+m)));
				//printf("Xi=%e\tYi=%e\tZi=%e\n",fabs(Xi),fabs(Yi),fabs(Zi));
				if((fabs(Xi)<1.0e-8)&&(fabs(Yi)<1.0e-8)&&(fabs(Zi)<1.0e-8))
					{
					*(dummy+((ni-1)*k+i)*nt*np+np*l+m)=*(dummy+((ni-1)*k+i)*nt*np+np*li+mi);
					//printf("li=%d\tmi=%d\n",li,mi);
					}
				}
				}
			}
			break;
			
		case 5://Top Zone face
		k=nk-2;
		face=5;
		//printf("\n\tSym-direction=%d",direction);
		for(j=0;j<nj-1;j++)
		for(i=0;i<ni-1;i++)
			{
			cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++) 
			if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
				{
				//printf("l=%d\tm=%d\t",l,m);
				//for(li=0;li<nt;li++)
				//for(mi=0;mi<np;mi++)
				//if(*(Ds+cellindex*6*nt*np+np*nt*face+np*li+mi)>0)
				//{
				Xi=(*(X+np*l+m)-2.0*(*(NX+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Yi=(*(Y+np*l+m)-2.0*(*(NY+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Zi=(*(Z+np*l+m)-2.0*(*(NZ+cellindex*6*nt*np+np*nt*face+np*l+m)));
				//printf("Xi=%e\tYi=%e\tZi=%e\n",Xi,Yi,Zi);
				//if((fabs(Xi)<1.0e-8)&&(fabs(Yi)<1.0e-8)&&(fabs(Zi)<1.0e-8))
					//{
					//*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=*(dummy+((ni-1)*j+i)*nt*np+np*li+mi);
					//printf("li=%d\tmi=%d\n",li,mi);
					//}
				//}
				*(Th+l)=acos(Zi);
				phiM=asin(Yi/sin(*(Th+l)));
								
				/*(Th+l)=acos(sin(0.1)*(*(Y+np*l+m))-cos(0.1)*(*(Z+np*l+m)));
				phiM=asin((cos(0.1)*(*(Y+np*l+m))+sin(0.1)*(*(Z+np*l+m)))/sin(*(Th+l)));*/
								
				//if((phiM>0)&&(fabs(sin(dtheta*(l+0.5))*cos(dphi*(m+0.5))-sin(*(Th+l))*cos(phiM))<1.0e-5))
				//else if((phiM<0)&&(fabs(sin(dtheta*(l+0.5))*cos(dphi*(m+0.5))-sin(*(Th+l))*cos(phiM))<1.0e-5))
				if((phiM>0)&&(fabs(Xi-sin(*(Th+l))*cos(phiM))<1.0e-8))
				*(Ph+m)=phiM;
				else if((phiM<0)&&(fabs(Xi-sin(*(Th+l))*cos(phiM))<1.0e-8))
				*(Ph+m)=2*Pi+phiM;
				else 
				*(Ph+m)=Pi-phiM;
				
				//fprintf(fp1,"\n%lf\t%lf\t%lf\t%lf",dtheta*(l+0.5)*180.0/Pi,dphi*(m+0.5)*180.0/Pi,*(Th+l)*180.0/Pi,*(Ph+m)*180.0/Pi);
				//printf("\n%lf\t%lf\t%lf\t%lf",dtheta*(l+0.5)*180.0/Pi,dphi*(m+0.5)*180.0/Pi,*(Th+l)*180.0/Pi,*(Ph+m)*180.0/Pi);
				
					ldi=*(Th+l)/dtheta;
					mdi=*(Ph+m)/dphi;
					lc=(int)ceil(ldi); mc=(int)ceil(mdi); lf=(int)floor(ldi); mf=(int)floor(mdi);
					if(mc==np)
					mc=0;
					
					if(lc==nt)
					{
					dphiC=fabs(dphi*(mc+0.5)-*(Ph+m));dphiF=fabs(*(Ph+m)-dphi*(mf+0.5)); 
					dthetaF=fabs(*(Th+l)-dphi*(lf+0.5));	thetaF=(*(Th+l)+dtheta*(lf+0.5))/2.0;
					Omega1=2.0*dphiF*sin(thetaF)*sin(dthetaF/2.0); Omega2=2.0*dphiC*sin(thetaF)*sin(dthetaF/2.0);
					if(dthetaF==0)
					*(CompInt+np*l+m)=*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf);
					else
					*(CompInt+np*l+m)=(Omega1*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mc))+Omega2*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf)))/(Omega1+Omega2);
					 }
					 else
					 {
					dphiC=fabs(dphi*(mc+0.5)-*(Ph+m));dphiF=fabs(*(Ph+m)-dphi*(mf+0.5));
					dthetaC=fabs(dtheta*(lc+0.5)-*(Th+l));dthetaF=fabs(*(Th+l)-dphi*(lf+0.5));
					thetaC=(*(Th+l)+dtheta*(lc+0.5))/2.0; thetaF=(*(Th+l)+dtheta*(lf+0.5))/2.0;
					Omega1=2.0*dphiF*sin(thetaF)*sin(dthetaF/2.0); Omega2=2.0*dphiC*sin(thetaF)*sin(dthetaF/2.0);
					Omega3=2.0*dphiF*sin(thetaC)*sin(dthetaC/2.0); Omega4=2.0*dphiC*sin(thetaC)*sin(dthetaC/2.0);

					*(CompInt+np*l+m)=(Omega1*(*(dummy+((ni-1)*j+i)*nt*np+np*lc+mc))+Omega4*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf))
							  +Omega3*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mc))+Omega2*(*(dummy+((ni-1)*j+i)*nt*np+np*lc+mf)))/(Omega1+Omega2+Omega3+Omega4);
					}
				
				}
				//fprintf(fp1,"\n\n");
				//printf("\n\n");
				
			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
					*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=*(CompInt+np*l+m);
					
			}
			break;
			
		case 6://Bottom Zone face
		k=0;
		face=4;
		//printf("\n\tSym-direction=%d",direction);
		for(j=0;j<nj-1;j++)
		for(i=0;i<ni-1;i++)
			{
			cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
			for(l=0;l<nt;l++)
			for(m=0;m<np;m++)
			if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
				{
				
				//*(Th+l)=acos(-sin(0.1)*(*(Y+np*l+m))-cos(0.1)*(*(Z+np*l+m)));
				//phiM=asin((cos(0.1)*(*(Y+np*l+m))-sin(0.1)*(*(Z+np*l+m)))/sin(*(Th+l)));
				
				Xi=(*(X+np*l+m)-2.0*(*(NX+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Yi=(*(Y+np*l+m)-2.0*(*(NY+cellindex*6*nt*np+np*nt*face+np*l+m)));
				Zi=(*(Z+np*l+m)-2.0*(*(NZ+cellindex*6*nt*np+np*nt*face+np*l+m)));
				
				*(Th+l)=acos(Zi);
				phiM=asin(Yi/sin(*(Th+l)));
				
				if((phiM>0)&&(fabs(Xi-sin(*(Th+l))*cos(phiM))<1.0e-8))
				*(Ph+m)=phiM;
				else if((phiM<0)&&(fabs(Xi-sin(*(Th+l))*cos(phiM))<1.0e-8))
				*(Ph+m)=2*Pi+phiM;
				else 
				*(Ph+m)=Pi-phiM;
				
				//if((phiM>0)&&(fabs(sin(dtheta*(l+0.5))*cos(dphi*(m+0.5))-sin(*(Th+l))*cos(phiM))<1.0e-5))
				//*(Ph+m)=phiM;
				//else if((phiM<0)&&(fabs(sin(dtheta*(l+0.5))*cos(dphi*(m+0.5))-sin(*(Th+l))*cos(phiM))<1.0e-5))
				//*(Ph+m)=2*Pi+phiM;
				//else 
				//*(Ph+m)=Pi-phiM;
				
				//fprintf(fp2,"\n%lf\t%lf\t%lf\t%lf",dtheta*(l+0.5)*180.0/Pi,dphi*(m+0.5)*180.0/Pi,*(Th+l)*180.0/Pi,*(Ph+m)*180.0/Pi);
				//printf("\n%lf\t%lf\t%lf\t%lf",dtheta*(l+0.5)*180.0/Pi,dphi*(m+0.5)*180.0/Pi,*(Th+l)*180.0/Pi,*(Ph+m)*180.0/Pi);
				
					ldi=*(Th+l)/dtheta;
					mdi=*(Ph+m)/dphi;
					lc=(int)ceil(ldi); mc=(int)ceil(mdi); lf=(int)floor(ldi); mf=(int)floor(mdi);
				
				
					if(mc==np)
					mc=0;
					
				//printf("\n%lf\t%lf\t%d\t%d\t%d\t%d",ldi,mdi,lf,lc,mf,mc);	
					if(lc==nt)
					{
					dphiC=fabs(dphi*(mc+0.5)-*(Ph+m));dphiF=fabs(*(Ph+m)-dphi*(mf+0.5)); 
					dthetaF=fabs(*(Th+l)-dphi*(lf+0.5));	thetaF=(*(Th+l)+dtheta*(lf+0.5))/2.0;
					Omega1=2.0*dphiF*sin(thetaF)*sin(dthetaF/2.0); Omega2=2.0*dphiC*sin(thetaF)*sin(dthetaF/2.0);
					if(dthetaF==0)
					*(CompInt+np*l+m)=*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf);
					else
					*(CompInt+np*l+m)=(Omega1*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mc))+Omega2*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf)))/(Omega1+Omega2);
					//printf("\n%lf\t%lf\t%lf\t%lf",dthetaF,dphiC,dphiF,thetaF);
					 //printf("\n%lf",*(CompInt+np*l+m));
					 }
					 else
					 {
					dphiC=fabs(dphi*(mc+0.5)-*(Ph+m));dphiF=fabs(*(Ph+m)-dphi*(mf+0.5));
					dthetaC=fabs(dtheta*(lc+0.5)-*(Th+l));dthetaF=fabs(*(Th+l)-dphi*(lf+0.5));
					thetaC=(*(Th+l)+dtheta*(lc+0.5))/2.0; thetaF=(*(Th+l)+dtheta*(lf+0.5))/2.0;
					Omega1=2.0*dphiF*sin(thetaF)*sin(dthetaF/2.0); Omega2=2.0*dphiC*sin(thetaF)*sin(dthetaF/2.0);
					Omega3=2.0*dphiF*sin(thetaC)*sin(dthetaC/2.0); Omega4=2.0*dphiC*sin(thetaC)*sin(dthetaC/2.0);

					*(CompInt+np*l+m)=(Omega1*(*(dummy+((ni-1)*j+i)*nt*np+np*lc+mc))+Omega4*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf))
							  +Omega3*(*(dummy+((ni-1)*j+i)*nt*np+np*lf+mc))+Omega2*(*(dummy+((ni-1)*j+i)*nt*np+np*lc+mf)))/(Omega1+Omega2+Omega3+Omega4);
							  //printf("\n%lf\t%lf\t%lf\t%lf",*(dummy+((ni-1)*j+i)*nt*np+np*lc+mc),*(dummy+((ni-1)*j+i)*nt*np+np*lf+mf),*(dummy+((ni-1)*j+i)*nt*np+np*lf+mc),*(dummy+((ni-1)*j+i)*nt*np+np*lc+mf));
							  //printf("\n%lf",*(CompInt+np*l+m));
					}
				
				}
				//fprintf(fp2,"\n\n");
				//printf("\n\n");
				
			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
				if(*(Ds+cellindex*6*nt*np+np*nt*face+np*l+m)<0)
					*(dummy+((ni-1)*j+i)*nt*np+np*l+m)=*(CompInt+np*l+m);
					
			}
			break;
	}
	free(Th);
	free(Ph);
	free(CompInt);
	//fclose(fp1);
	//fclose(fp2);
		//fflush(stdout);
return 0;
}
