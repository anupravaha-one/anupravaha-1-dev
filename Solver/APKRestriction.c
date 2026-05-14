#include "APKRestriction.h"

double Restriction(Base *base,Zone *zone,int g)
{
	int i,j,k,ic,jc,kc,NI,NJ,NK,NIC,NJC,NKC;
	double *A,*B;
	
	A=(zone->grid+g)->Gresidual;
	B=(zone->grid+g+1)->GRHS;
	
	NIC = *((zone->grid+g+1)->GSize);
	NJC = *((zone->grid+g+1)->GSize+1);
	NKC = *((zone->grid+g+1)->GSize+2);
	NI = *((zone->grid+g)->GSize);
	NJ = *((zone->grid+g)->GSize+1);
	NK = *((zone->grid+g)->GSize+2);
	
	if(base->BaseParameters->flag2d == 1)//2D case
	{
		if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
		{
//			printf("i am here in 1\n");
/*			kc = 0;
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}*/
//----------------------------------------------------------------------------------------------------------1-----------
			kc = 0;
			jc = NJC-2;
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			kc = 0;
			ic = NIC-2;
			for(jc=0;jc<NJC-2;jc++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i));
			}
			kc = 0;
			jc = NJC-2, ic = NIC-2;

			i = 2*ic;j = 2*jc;k = kc;
					
			*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=1*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

			kc = 0;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
//-------------------------------------------------------------------------------------------------------------2---------
/*			kc = 0;
			jc = NJC-2;
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			kc = 0;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}*/

		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
		{
//			printf("i am here in 2\n");
/*			kc = 0;
			ic = NIC-2;
			for(jc=0;jc<NJC-1;jc++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i));
			}
			kc = 0;
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}*/
//----------------------------------------------------------------------------------------------------------1----------
/*			kc = 0;
			jc = NJC-2;
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			kc = 0;
			ic = NIC-2;
			for(jc=0;jc<NJC-2;jc++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i));
			}
			kc = 0;
			jc = NJC-2, ic = NIC-2;

			i = 2*ic;j = 2*jc;k = kc;
					
			*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=1*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

			kc = 0;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}*/
//--------------------------------------------------------------------------------------------------------2-----------
			kc = 0;
			jc = NJC-2;
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			kc = 0;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
		{
//			printf("i am here in 3\n");
			kc = 0;
			jc = NJC-2;
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			kc = 0;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
		{
//			printf("i am here in 4\n");
			kc = 0;
			jc = NJC-2;
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			kc = 0;
			ic = NIC-2;
			for(jc=0;jc<NJC-2;jc++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i));
			}
			kc = 0;
			jc = NJC-2, ic = NIC-2;

			i = 2*ic;j = 2*jc;k = kc;
					
			*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=1*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

			kc = 0;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
					i = 2*ic;j = 2*jc;k = kc;
					
					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic)=4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
		}
	}
	else//3D case
	{
		if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
		{
			for(kc=0;kc<NKC-1;kc++)
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
					i = 2*ic;j = 2*jc;k = 2*kc;

					*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) = 8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												    +*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));

			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
		{
			ic = NIC-2;
			for(kc=0;kc<NKC-1;kc++)
			for(jc=0;jc<NJC-1;jc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;
				
				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i));
			}
			for(kc=0;kc<NKC-1;kc++)
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
		{
			jc = NJC-2;
			for(kc=0;kc<NKC-1;kc++)
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1));
			}
			for(kc=0;kc<NKC-1;kc++)
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
		{
			kc = NKC-2;
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
			for(kc=0;kc<NKC-2;kc++)
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
		{
			kc = NKC-2;
			jc = NJC-2;
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			
			jc = NJC-2;
			for(kc=0;kc<NKC-2;kc++)
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1));
			}
			kc = NKC-2;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
			for(kc=0;kc<NKC-2;kc++)
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-1;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
		{
			kc = NKC-2;
			ic = NIC-2;
			for(jc=0;jc<NJC-1;jc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i));
			}
			ic = NIC-2;
			for(kc=0;kc<NKC-2;kc++)
			for(jc=0;jc<NJC-1;jc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i));
			}
			kc = NKC-2;
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
			for(kc=0;kc<NKC-2;kc++)
			for(jc=0;jc<NJC-1;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
		{
			jc = NJC-2;
			ic = NIC-2;
			for(kc=0;kc<NKC-1;kc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i));
			}
			ic = NIC-2;
			for(kc=0;kc<NKC-1;kc++)
			for(jc=0;jc<NJC-2;jc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i));
			}
			jc = NJC-2;
			for(kc=0;kc<NKC-1;kc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1));
			}
			for(kc=0;kc<NKC-1;kc++)
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
		else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
		{
			kc = NKC-2;
			jc = NJC-1;
			ic = NIC-1;

			i = 2*ic;j = 2*jc;k = 2*kc;

			*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i);

			kc = NKC-2;ic = NIC-2;
			for(jc=0;jc<NJC-2;jc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i));
			}
			kc = NKC-2;jc = NJC-2;
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1));
			}
			ic = NIC-2;jc = NJC-2;
			for(kc=0;kc<NKC-2;kc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	2*0.5*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i));
			}
			kc = NKC-2;
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1));
			}
			jc = NJC-2;
			for(kc=0;kc<NKC-2;kc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1));
			}
			ic = NIC-2;
			for(kc=0;kc<NKC-2;kc++)
			for(jc=0;jc<NJC-2;jc++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	4*0.25*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i));
			}
			for(kc=0;kc<NKC-2;kc++)
			for(jc=0;jc<NJC-2;jc++)
			for(ic=0;ic<NIC-2;ic++)
			{
				i = 2*ic;j = 2*jc;k = 2*kc;

				*(B+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic) =	8*0.125*(*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1)
												+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i)+*(A+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1));
			}
		}
	}
	return 0;
}
