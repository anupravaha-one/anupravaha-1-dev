#include "APKCalculateNormalVector.h"



int CalculateNormalVector(Base *CgnsBase, Zone *z, int index[3], int face, double *surfVector)
	{
	int i, j, k,  l, m, NI, NJ, NK, nt, np, cellindex;
	double theta,phi,dtheta,dphi,PI,*CalcX,*CalcY,*CalcZ,Ax,Ay,Az,Area;
	
	PI=4.0*atan(1.0);
	nt=CgnsBase->BaseParameters->ntheta;
	np=CgnsBase->BaseParameters->nphi;
	NI=*(z->Size+0);
	NJ=*(z->Size+1);
	NK=*(z->Size+2);
	
	CalcX=z->ZGeom->NormalX;
	CalcY=z->ZGeom->NormalY;
	CalcZ=z->ZGeom->NormalZ;
	i=index[0]; j=index[1]; k=index[2];
	dtheta=PI/nt;
	dphi=2.0*PI/np;
	
	switch(face)
		{
		case 0://west face
			for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
				for(m=0;m<np;m++)
				{
				cellindex=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				phi=dphi*(m+0.5);
				Ax=*(surfVector+0);
				Ay=*(surfVector+1);
				Az=*(surfVector+2);
				Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
				//printf("Aw=%e\t", Area);
				*(CalcX+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcY+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+1)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcZ+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+2)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				}
			}
		break;
		case 1://East face
			for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
				for(m=0;m<np;m++)
				{
				cellindex=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				phi=dphi*(m+0.5);
				Ax=*(surfVector+0);
				Ay=*(surfVector+1);
				Az=*(surfVector+2);
				Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
				//printf("Ae=%e\t", Area);
				*(CalcX+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcY+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+1)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcZ+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+2)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				}
			}
		break;
		case 2: //south face
			for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
				for(m=0;m<np;m++)
				{
				cellindex=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				phi=dphi*(m+0.5);
				Ax=*(surfVector+0);
				Ay=*(surfVector+1);
				Az=*(surfVector+2);
				Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
				//printf("As=%e\t", Area);
				*(CalcX+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcY+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+1)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcZ+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+2)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				}
			}
		break;
		case 3: //north face
			for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
				for(m=0;m<np;m++)
				{
				cellindex=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				phi=dphi*(m+0.5);
				Ax=*(surfVector+0);
				Ay=*(surfVector+1);
				Az=*(surfVector+2);
				Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
				//printf("An=%e\t", Area);
				*(CalcX+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcY+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+1)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcZ+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+2)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				}
			}
		break;
		case 4: //bottom face
			for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
				for(m=0;m<np;m++)
				{
				cellindex=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				phi=dphi*(m+0.5);
				Ax=*(surfVector+0);
				Ay=*(surfVector+1);
				Az=*(surfVector+2);
				Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
				//printf("Ab=%e\t", Area);
				*(CalcX+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcY+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+1)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcZ+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+2)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				}
			}
		break;
		case 5: //top face
			for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
				for(m=0;m<np;m++)
				{
				cellindex=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				phi=dphi*(m+0.5);
				Ax=*(surfVector+0);
				Ay=*(surfVector+1);
				Az=*(surfVector+2);
				Area=sqrt(Ax*Ax+Ay*Ay+Az*Az);
				//printf("At=%e\t", Area);
				*(CalcX+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcY+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+1)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				*(CalcZ+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+2)*(*(surfVector+0)*sin(theta)*cos(phi)+*(surfVector+1)*sin(theta)*sin(phi)+*(surfVector+2)*cos(theta))/(Area*Area);
				}
			}
		break;
		}
		
	return 0;
	}




int CalculateVector(Base *CgnsBase, double *CalcX, double *CalcY, double *CalcZ)
	{
	 int nt, np, l, m;
	 double theta, phi, dtheta, dphi, PI;
	 
	 PI=4.0*atan(1.0);
	 nt=CgnsBase->BaseParameters->ntheta;
	 np=CgnsBase->BaseParameters->nphi;
	 dtheta=PI/nt;
	 dphi=2.0*PI/np;
	 
	 	for(l=0;l<nt;l++)
			{
			theta=dtheta*(l+0.5);
			for(m=0;m<np;m++)
				{
				 phi=dphi*(m+0.5);
				 *(CalcX+np*l+m)=sin(theta)*cos(phi);
				 *(CalcY+np*l+m)=sin(theta)*sin(phi);
				 *(CalcZ+np*l+m)=cos(theta);
				 }
				 
			}
			
	return 0;
	}
