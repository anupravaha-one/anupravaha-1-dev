#include "APKProlongation2d.h"

double Prolongation2d(Base *base,Zone *zone,int gdi)
{
	int i,j,k,ic,jc,kc,NI,NJ,NK,NIC,NJC,NKC;
	double Vp,Ve,Vw,Vn,Vs,Vne,Vse,Vnw,Vsw;
	double *A,*B;
	
	A=(zone->grid+gdi)->V;
	B=(zone->grid+gdi-1)->E;

	NIC = *((zone->grid+gdi)->GSize);
	NJC = *((zone->grid+gdi)->GSize+1);
	NKC = *((zone->grid+gdi)->GSize+2);
	NI = *((zone->grid+gdi-1)->GSize);
	NJ = *((zone->grid+gdi-1)->GSize+1);
	NK = *((zone->grid+gdi-1)->GSize+2);
	
	if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
	{
//		printf("entered prolongation ieven jeven\n");
		kc = 0;k = 0;
		
		ic = 0;i = 0;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-2;j = NJ-3;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc=1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
		
		ic = NIC-2;i = NI-3;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-2;j = NJ-3;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
		
		jc=1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
		
		jc = 0;j = 0;
		ic=1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			ic++;
		}
		
		jc = NJC-2;j = NJ-3;
		ic=1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			ic++;
		}
		
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				ic++;
			}
			jc++;
		}
	}
	
	else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
	{
//		printf("entered prolongation iodd jeven\n");
		kc = 0;k = 0;
		
		ic = 0;i = 0;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-2;j = NJ-3;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc=1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
		
		ic = NIC-3;i = NI-4;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-2;j = NJ-3;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc=1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
		
		jc = 0;j = 0;
		ic=1;
		for(i=2;i<NI-4;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			ic++;
		}
		
		jc = NJC-2;j = NJ-3;
		ic=1;
		for(i=2;i<NI-4;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			ic++;
		}
		
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-4;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				ic++;
			}
			jc++;
		}
		
		ic = NIC-2;i = NI-2;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		
		jc = NJC-2;j = NJ-3;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			jc++;
		}
	}
	else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
	{
//		printf("entered prolongation ieven jodd\n");
		kc = 0;k = 0;
		
		ic = 0;i = 0;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-3;j = NJ-4;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc=1;
		for(j=2;j<NJ-4;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
		
		ic = NIC-2;i = NI-3;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-3;j = NJ-4;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc=1;
		for(j=2;j<NJ-4;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
		
		jc = 0;j = 0;
		ic=1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			ic++;
		}
		
		jc = NJC-3;j = NJ-4;
		ic=1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			ic++;
		}
		
		jc = 1;
		for(j=2;j<NJ-4;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				ic++;
			}
			jc++;
		}

		jc = NJC-2;j = NJ-2;
		
		ic = 0;i = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		
		ic = NIC-2;i = NI-3;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;

		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			ic++;
		}
	}
	else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
	{
//		printf("entered prolongation iodd jodd\n");
		kc = 0;k = 0;
		
		ic = 0;i = 0;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		
		jc = NJC-2;j = NJ-2;
		
		ic = 0;i = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		
		ic = 1;
		for(i=2;i<NI-4;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			ic++;
		}

		ic = NIC-2;i = NI-2;
		
		jc = 0;j = 0;
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		
		jc = 1;
		for(j=2;j<NJ-4;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			jc++;
		}

		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+NI-2) = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+NIC-2);

		jc = 1;
		for(j=2;j<NJ-4;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-4;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				ic++;
			}
			jc++;
		}
	}
	return 0;
}

