#include "APKFillNewRadsolToPrevRadsol.h"

/** This function updates the previous interation value Intensity with new calculated iteration value.*/

int FillNewRadsolToPrevRadsol(Base *base,Zone *z)
	{
	 int ni,nj,nk,nt,np,i,j,k,l,m,cellindex;
	 
	 ni=*(z->Size+0);
	 nj=*(z->Size+1);
	 nk=*(z->Size+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 
	 
	 for(l=0;l<nt;l++)
	 	for(m=0;m<np;m++)
			for(k=0;k<(nk-1);k++)
				for(j=0;j<(nj-1);j++)
					for(i=0;i<(ni-1);i++)
					{
					cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
					*(z->PrevFlowsol->RadIntensity+cellindex*np*nt+np*l+m)=*(z->Flowsol->RadIntensity+cellindex*np*nt+np*l+m);
					}
					
			//UpdateRadOldZSlabs(base, z);
	 return 0; 
	}
