#include "APKCalculateDmF.h"

/** This function calculate the directional weight and control solid angle in the preprocessing stage*/

int CalculateDmF(Base *CgnsBase, Zone *z, int index[3], int face, double *surfVector)
	{
	int i, j, k,  l, m, NI, NJ, NK, nt, np, cellindex;
	double theta,phi,dtheta,dphi,PI,*CalcDmF;
	
	PI=4.0*atan(1.0);
	nt=CgnsBase->BaseParameters->ntheta;
	np=CgnsBase->BaseParameters->nphi;
	NI=*(z->Size+0);
	NJ=*(z->Size+1);
	NK=*(z->Size+2);
	
	CalcDmF=z->ZGeom->DmF;
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
				*(CalcDmF+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+1)*(sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+2)*(dphi*sin(theta)*cos(theta)*sin(dtheta));
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
				*(CalcDmF+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+1)*(sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+2)*(dphi*sin(theta)*cos(theta)*sin(dtheta));
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
				*(CalcDmF+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+1)*(sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+2)*(dphi*sin(theta)*cos(theta)*sin(dtheta));
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
				*(CalcDmF+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+1)*(sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+2)*(dphi*sin(theta)*cos(theta)*sin(dtheta));
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
				*(CalcDmF+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+1)*(sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+2)*(dphi*sin(theta)*cos(theta)*sin(dtheta));
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
				*(CalcDmF+cellindex*6*nt*np+np*nt*face+np*l+m)=*(surfVector+0)*(cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+1)*(sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta)))+*(surfVector+2)*(dphi*sin(theta)*cos(theta)*sin(dtheta));
				}
			}
		break;
		}
		
	return 0;
	}




int CalculateOmega(Base *CgnsBase, double *CalcOmega)
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
				 *(CalcOmega+np*l+m)=2.0*dphi*sin(theta)*sin(dtheta/2.0);
				 }
				 
			}
			
	return 0;
	}
