#include "APKPressureForce.h"

/**
This function loops over all boundaries and calculates the pressure force (p.dA) in x,y and z 
directions and stores them in bcPtr->PressForce.
*/

int PressureForce(Base *base)
{
	FILE *Pres;
	int iface,NI,NJ,NK,*ptrSize,nZone,numZone;
	int i,j,k,dir,INDEX;
	char tag,name[33];
	void *ptr;
	double Px,Py,Pz;
	double *Press,*surf;
	double p,q,r,Area;
	Zone *zone;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *FindZonePtr();
	Pres = fopen("SOLUTIONS/PressForce.dat","w");
	numZone=(base->Number_of_zones);
	for(nZone=0;nZone<numZone;nZone++)
	{
		
		zone=base->CgnsZone+nZone;
		fprintf(Pres,"\n~~~~~~~ZONENAME :%s ~~~~~~~~\n",(zone->zonename));
		ptrSize = zone->Size;
		NI = *(ptrSize+0);
		NJ = *(ptrSize+1);
		NK = *(ptrSize+2);
		printf("NI=%d\tNJ=%d\tNK=%d\n",NI,NJ,NK);
		//Press=zone->Flowsol->P;
		surf= zone->ZGeom->SurfaceArea;
		for(iface=0;iface<6;iface++)
		{

			ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
			if(tag=='i')
			{
				//printf("I am in Interface\n");
				facePtr=(Interface *)ptr;
			}
			else if(tag=='b')
			{
				//printf("I am in Boundary\n");
				bcPtr=(ZoneBoundary *)ptr;
				Press=bcPtr->ZSlab->P;
				dir=bcPtr->Direction;
				switch(dir)
				{
					case 1:
					Px=0.0;
					Py=0.0;
					Pz=0.0;
					//printf("I am in Case 1\n");
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							i=NI-2;
							INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
							//printf("INDEX Value %d\n",INDEX);
							p=*(surf+(INDEX)*6*3+3*1+0);
							q=*(surf+(INDEX)*6*3+3*1+1);
							r=*(surf+(INDEX)*6*3+3*1+2);
							//printf("\nEAST Imax\tPressX=%lf PressY=%lf PressZ=%lf\n",p,q,r);
							Area=sqrt(p*p+q*q+r*r);
							Px= Px+*(Press+((NJ-1)*k+j))*p;
							Py= Py+*(Press+((NJ-1)*k+j))*q;
							Pz= Pz+*(Press+((NJ-1)*k+j))*r;
						}
						(bcPtr->PressForce[0])=Px;
						(bcPtr->PressForce[1])=Py;
						(bcPtr->PressForce[2])=Pz;
						fprintf(Pres,"\nEAST Imax\tdir=%d\tPressX=%lf PressY=%lf PressZ=%lf\n",dir,(bcPtr->PressForce[0]),(bcPtr->PressForce[1]),(bcPtr->PressForce[2]));
					break;
					case 2:
					Px=0.0;
					Py=0.0;
					Pz=0.0;
					//printf("I am in Case 2\n");
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							i=0;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							p=*(surf+(INDEX)*6*3+3*0+0);
							q=*(surf+(INDEX)*6*3+3*0+1);
							r=*(surf+(INDEX)*6*3+3*0+2);
							Area=sqrt(p*p+q*q+r*r);
							Px= Px+*(Press+((NJ-1)*k+j))*p;
							Py= Py+*(Press+((NJ-1)*k+j))*q;
							Pz= Pz+*(Press+((NJ-1)*k+j))*r;
						}
						(bcPtr->PressForce[0])=Px;
						(bcPtr->PressForce[1])=Py;
						(bcPtr->PressForce[2])=Pz;
						fprintf(Pres,"\nWEST Imin\tdir=%d\tPressX=%lf PressY=%lf PressZ=%lf\n",dir,(bcPtr->PressForce[0]),(bcPtr->PressForce[1]),(bcPtr->PressForce[2]));
					break;
					case 3:
					Px=0.0;
					Py=0.0;
					Pz=0.0;
					//printf("I am in Case 3\n");
					for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							j=NJ-2;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							p=*(surf+(INDEX)*6*3+3*3+0);
							q=*(surf+(INDEX)*6*3+3*3+1);
							r=*(surf+(INDEX)*6*3+3*3+2);
							Area=sqrt(p*p+q*q+r*r);
							Px= Px+*(Press+((NI-1)*k+i))*p;
							Py= Py+*(Press+((NI-1)*k+i))*q;
							Pz= Pz+*(Press+((NI-1)*k+i))*r;
						}
						(bcPtr->PressForce[0])=Px;
						(bcPtr->PressForce[1])=Py;
						(bcPtr->PressForce[2])=Pz;
						fprintf(Pres,"\n\nNORTH Jmax\tdir=%d\tPressX=%lf PressY=%lf PressZ=%lf\n",dir,(bcPtr->PressForce[0]),(bcPtr->PressForce[1]),(bcPtr->PressForce[2]));
					break;
					case 4:
					Px=0.0;
					Py=0.0;
					Pz=0.0;
					//printf("I am in Case 4\n");
					for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							j=0;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							p=*(surf+(INDEX)*6*3+3*2+0);
							q=*(surf+(INDEX)*6*3+3*2+1);
							r=*(surf+(INDEX)*6*3+3*2+2);
							Area=sqrt(p*p+q*q+r*r);
							Px= Px+*(Press+((NI-1)*k+i))*p;
							Py= Py+*(Press+((NI-1)*k+i))*q;
							Pz= Pz+*(Press+((NI-1)*k+i))*r;
						}
						(bcPtr->PressForce[0])=Px;
						(bcPtr->PressForce[1])=Py;
						(bcPtr->PressForce[2])=Pz;
						fprintf(Pres,"\n\nSOUTH Jmin\tdir=%d\tPressX=%lf PressY=%lf PressZ=%lf\n",dir,(bcPtr->PressForce[0]),(bcPtr->PressForce[1]),(bcPtr->PressForce[2]));
					break;
					case 5:
					Px=0.0;
					Py=0.0;
					Pz=0.0;
					//printf("I am in Case 5\n");
					for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							k=NK-2;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							p=*(surf+(INDEX)*6*3+3*5+0);
							q=*(surf+(INDEX)*6*3+3*5+1);
							r=*(surf+(INDEX)*6*3+3*5+2);
							Area=sqrt(p*p+q*q+r*r);
							Px= Px+*(Press+((NI-1)*j+i))*p;
							Py= Py+*(Press+((NI-1)*j+i))*q;
							Pz= Pz+*(Press+((NI-1)*j+i))*r;
						}
						(bcPtr->PressForce[0])=Px;
						(bcPtr->PressForce[1])=Py;
						(bcPtr->PressForce[2])=Pz;
						fprintf(Pres,"\n\nTOP Kmax\tdir=%d\tPressX=%lf PressY=%lf PressZ=%lf\n",dir,(bcPtr->PressForce[0]),(bcPtr->PressForce[1]),(bcPtr->PressForce[2]));
					break;
					case 6:
					Px=0.0;
					Py=0.0;
					Pz=0.0;
					//printf("I am in Case 6\n");
					for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							k=0;
							INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
							p=*(surf+(INDEX)*6*3+3*4+0);
							q=*(surf+(INDEX)*6*3+3*4+1);
							r=*(surf+(INDEX)*6*3+3*4+2);
							Area=sqrt(p*p+q*q+r*r);
							Px= Px+*(Press+((NI-1)*j+i))*p;
							Py= Py+*(Press+((NI-1)*j+i))*q;
							Pz= Pz+*(Press+((NI-1)*j+i))*r;
						}
						(bcPtr->PressForce[0])=Px;
						(bcPtr->PressForce[1])=Py;
						(bcPtr->PressForce[2])=Pz;
						fprintf(Pres,"\n\nBOTTAM Kmin\tdir=%d\tPressX=%lf PressY=%lf PressZ=%lf\n",dir,(bcPtr->PressForce[0]),(bcPtr->PressForce[1]),(bcPtr->PressForce[2]));
					break;
				}
			}
		}
	}
	fclose(Pres);
	return (0);
}
