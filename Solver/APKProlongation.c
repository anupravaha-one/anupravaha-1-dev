#include "APKProlongation.h"

double Prolongation(Base *base,Zone *zone,int gdi)
{
	int i,j,k,ic,jc,kc,NI,NJ,NK,NIC,NJC,NKC;
	double Vp,Ve,Vw,Vn,Vs,Vne,Vse,Vnw,Vsw,Vt,Vb,Vte,Vtw,Vtn,Vts,Vtnw,Vtne,Vtsw,Vtse;
	double *A,*B;
	
	A=(zone->grid+gdi)->V;
	B=(zone->grid+gdi-1)->E;

	NIC = *((zone->grid+gdi)->GSize);
	NJC = *((zone->grid+gdi)->GSize+1);
	NKC = *((zone->grid+gdi)->GSize+2);
	NI = *((zone->grid+gdi-1)->GSize);
	NJ = *((zone->grid+gdi-1)->GSize+1);
	NK = *((zone->grid+gdi-1)->GSize+2);
	
	if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = Vt;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (Vt)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vtn)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3.5*Vte+3.5*Vtn)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (Vt)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vte)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vte)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (Vt)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-3;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3.5*Vte+3.5*Vtn)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-3;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-3;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-3;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vtn)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-3;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = Vt;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (Vt)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vte)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vte)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtw)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-2;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-2;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-2;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-2;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-2;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = Vt;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (Vt)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vtn)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vts)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-3;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-3;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-3;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-3;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+7*Vtn)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-3;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-3;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3.5*Vte+3.5*Vtn)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-3;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3.5*Ve+3.5*Vn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3.5*Vte+3.5*Vtn)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-3;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-3;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-3;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-3;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-3;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-3;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-3;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-3;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3.5*Ve+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3.5*Vte+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-3;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-3;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-3;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-3;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-3;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-3;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-3;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3.5*Vw+3.5*Vn)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3.5*Vtw+3.5*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-2;
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-2;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-2;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-2;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-3;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-3;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-3;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-2;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = Vt;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+7*Vte)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-3;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vtw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-3;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-2;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-2;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+7*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+7*Vtn)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-2;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-2;
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-3;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-3;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-2;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
	{
//----------------------------------------------------------------------------------------------(1)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
		Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
		Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
//----------------------------------------------------------------------------------------------(2)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(3)----------------
		kc = 0, k = 0;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
		Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
		
		Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
		Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
		Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
//----------------------------------------------------------------------------------------------(4)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
		Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
		Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
//----------------------------------------------------------------------------------------------(5)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(6)----------------
		kc = NKC-2, k = NK-2;
		jc = 0, j = 0;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
		Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
		Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
		Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
		Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
//----------------------------------------------------------------------------------------------(7)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
		Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
		Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
		Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
//----------------------------------------------------------------------------------------------(8)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(9)----------------
		kc = 0, k = 0;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		
		Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
		Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
		Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
		Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
		Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
		
		Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
		Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
		Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
//----------------------------------------------------------------------------------------------(10)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-2;
		ic = 0, i = 0;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
		Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
		Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
		Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
//----------------------------------------------------------------------------------------------(11)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-2;
		ic = 1;
		for(i=2;i<NI-2;i=i+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			
			ic++;
		}
//----------------------------------------------------------------------------------------------(12)----------------
		kc = NKC-2, k = NK-2;
		jc = NJC-2, j = NJ-2;
		ic = NIC-2, i = NI-2;
		
		Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
		Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
		Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
		Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
		Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
		Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
		
		Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
		Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
		Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
		
		Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
		Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
		Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
		
		Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
		
		*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
//----------------------------------------------------------------------------------------------(13)----------------
		ic = 0;i = 0;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(14)----------------
		ic = NIC-2;i = NI-2;
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(15)----------------
		ic = 0;i = 0;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3.5*Vw+3.5*Vs)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3.5*Vtw+3.5*Vts)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(16)----------------
		ic = NIC-2;i = NI-2;
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
			
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			
			kc++;
		}
//----------------------------------------------------------------------------------------------(17)----------------
		ic = 0;i = 0;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
			Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
			Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
			Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
			Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(18)----------------
		ic = 0;i = 0;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
			Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
			Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
			Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
			Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
			Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(19)----------------
		ic = NIC-2;i = NI-2;
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
			Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
			Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
			Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
			Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
			Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
			Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
			Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
			*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(20)----------------
		ic = NIC-2;i = NI-2;
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
			Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
			Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
			Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
			Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
			Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
			
			Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
			Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
			Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
			Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
			
			Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
			Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
			Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
			Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
			
			Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
			Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
			*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
			jc++;
		}
//----------------------------------------------------------------------------------------------(21)----------------
		ic = 0;i = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc+1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*kc+jc-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc+1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[1]+(NJC-1)*(kc+1)+jc-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(22)----------------
		ic = NIC-2;i = NI-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*kc+jc-1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[0]+(NJC-1)*(kc+1)+jc-1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				jc++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(23)----------------
		kc = 0;k = 0;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *((zone->grid+gdi)->SlabV[5]+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(24)----------------
		kc = NKC-2;k = NK-2;
		jc = 1;
		for(j=2;j<NJ-2;j=j+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic);
				Vtw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic-1);
				Vte = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic);
				
				Vtnw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[4]+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				
				ic++;
			}
			jc++;
		}
//----------------------------------------------------------------------------------------------(25)----------------
		jc = 0;j = 0;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic);
				Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
				Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
				Vsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic-1);
				Vse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*kc+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic);
				Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
				Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
				Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
				Vtsw = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic-1);
				Vtse = *((zone->grid+gdi)->SlabV[3]+(NIC-1)*(kc+1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(26)----------------
		jc = NJC-2;j = NJ-2;
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			ic = 1;
			for(i=2;i<NI-2;i=i+2)
			{
				Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
				
				Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
				Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
				Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
				Vn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic);
				Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
				Vnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic-1);
				Vne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*kc+ic+1);
				Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
				Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
				
				Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
				Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
				Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
				Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
				Vtn = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic);
				Vtnw = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic-1);
				Vtne = *((zone->grid+gdi)->SlabV[2]+(NIC-1)*(kc+1)+ic+1);
				Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
				Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
				
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
				*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
				*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
				
				ic++;
			}
			kc++;
		}
//----------------------------------------------------------------------------------------------(27)----------------
		kc = 1;
		for(k=2;k<NK-2;k=k+2)
		{
			jc = 1;
			for(j=2;j<NJ-2;j=j+2)
			{
				ic = 1;
				for(i=2;i<NI-2;i=i+2)
				{
					Vp = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic);
					
					Vw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic-1);
					Ve = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+ic+1);
					Vs = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic);
					Vn = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic);
					Vb = *(A+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+ic);
					Vnw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic-1);
					Vne = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+ic+1);
					Vsw = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic-1);
					Vse = *(A+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+ic+1);
					
					Vt = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic);
					Vtw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic-1);
					Vte = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+ic+1);
					Vts = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic);
					Vtn = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic);
					Vtnw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic-1);
					Vtne = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+ic+1);
					Vtsw = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic-1);
					Vtse = *(A+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+ic+1);
					
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (9*Vp+3*Vw+3*Vs+Vsw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1) = (9*Vp+3*Ve+3*Vs+Vse)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i) = (9*Vp+3*Vw+3*Vn+Vnw)/16;
					*(B+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1) = (9*Vp+3*Ve+3*Vn+Vne)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i) = (9*Vt+3*Vtw+3*Vts+Vtsw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1) = (9*Vt+3*Vte+3*Vts+Vtse)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i) = (9*Vt+3*Vtw+3*Vtn+Vtnw)/16;
					*(B+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i+1) = (9*Vt+3*Vte+3*Vtn+Vtne)/16;
					
					ic++;
				}
				jc++;
			}
			kc++;
		}
	}
	return 0;
}
