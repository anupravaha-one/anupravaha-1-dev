#include "APKDummyVariableProlongate.h"

int DummyVariableProlongate(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface;
	char tag,name[33];
	void *ptr;
	int i,j,k,ic,jc,kc,NI,NJ,NK,NIC,NJC,NKC;
	double Vp,Vw,Ve,Vb,Vt,Vs,Vn,Vtw,Vte,Vbw,Vbe,Vwn,Ven,Vws,Ves,Vts,Vtn,Vbs,Vbn;
	
	NI = *((zone->grid+gdi-1)->GSize+0);
	NJ = *((zone->grid+gdi-1)->GSize+1);
	NK = *((zone->grid+gdi-1)->GSize+2);
	
	NIC = *((zone->grid+gdi)->GSize+0);
	NJC = *((zone->grid+gdi)->GSize+1);
	NKC = *((zone->grid+gdi)->GSize+2);
	
	for(iface=0;iface<6;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			switch(iface)
			{
				case 0:
					if(((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
						
						kc = NKC-2;k = NK-3;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = Vp;
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-3;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-3;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
						
						kc = NKC-2;k = NK-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
						
						kc = NKC-2;k = NK-3;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-3;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-2;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
					else if(((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
						
						kc = NKC-2;k = NK-2;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vb);
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-2;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-3;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
						
						kc = NKC-2;k = NK-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
						
						kc = NKC-2;k = NK-2;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-2;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-2;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+NIC-2);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+NIC-2);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+NIC-2);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+NIC-2);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+NIC-2);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+NIC-2);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+NI-2) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+NI-2) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}

				break;
				case 1:
					if(((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
						
						kc = NKC-2;k = NK-3;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = Vp;
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-3;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-3;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
						
						kc = NKC-2;k = NK-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
						
						kc = NKC-2;k = NK-3;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vs);
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-3;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-2;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
					else if(((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
						
						kc = NKC-2;k = NK-2;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vb);
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-2;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-3;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
						
						kc = 0;k = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
						
						kc = NKC-2;k = NK-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
						Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
						
						kc = NKC-2;k = NK-2;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
						Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
						
						kc = 0;k = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							jc++;
						}

						kc = NKC-2;k = NK-2;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							jc++;
						}
						
						jc = 0;j = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
							Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
							kc++;
						}

						jc = NJC-2;j = NJ-2;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
							Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
							Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*jc+0);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*jc+0);
								Vts = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc-1)+0);
								Vtn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(jc+1)+0);
								Vbs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc-1)+0);
								Vbn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(jc+1)+0);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.25*(Vp+Vb+Vs+Vbs);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+0) = 0.25*(Vp+Vt+Vs+Vts);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.25*(Vp+Vb+Vn+Vbn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+0) = 0.25*(Vp+Vt+Vn+Vtn);
								jc++;
							}
							kc++;
						}
					}
				break;
				case 2:
					if(((NI-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-3;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
						
						kc = NKC-2;k = NK-3;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = Vp;
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-3;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-3;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-3;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
						
						kc = NKC-2;k = NK-3;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
						
						kc = NKC-2;k = NK-3;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-3;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-2;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-2;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
					else if(((NI-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-2;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						
						kc = NKC-2;k = NK-2;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vb);
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-2;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-3;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-3;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
						
						kc = NKC-2;k = NK-2;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						
						kc = NKC-2;k = NK-2;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-2;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-2;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-2;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*(NJC-2)+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*(NJC-2)+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(NJ-2)+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
				break;
				case 3:
					if(((NI-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-3;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
						
						kc = NKC-2;k = NK-3;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = Vp;
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-3;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-3;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-3;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
						
						kc = NKC-2;k = NK-3;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
						
						kc = NKC-2;k = NK-3;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vw);
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-3;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-2;
						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-3;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-2;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
					else if(((NI-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Vt);
						
						kc = NKC-2;k = NK-2;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						
						kc = NKC-2;k = NK-2;
						ic = NIC-2;i = NI-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Vb);
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-2;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-3;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.5*(Vp+Vt);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-3;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						kc = 0;k = 0;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
						
						kc = 0;k = 0;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
						
						kc = NKC-2;k = NK-2;
						ic = 0;i = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
						Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
						
						kc = NKC-2;k = NK-2;
						ic = NIC-2;i = NI-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
						Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
						
						kc = 0;k = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							ic++;
						}

						kc = NKC-2;k = NK-2;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							ic++;
						}
						
						ic = 0;i = 0;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
							Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(Vp+Vb);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.5*(Vp+Vt);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
							kc++;
						}

						ic = NIC-2;i = NI-2;
						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
							Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
							Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
							Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
							Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
							kc++;
						}

						kc = 1;
						for(k=2;k<NK-2;k=k+2)
						{
							ic = 1;
							for(i=2;i<NI-2;i=i+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1);
								Vb = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic);
								Vt = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic);
								Vtw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic-1);
								Vte = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc+1)+(NIC-1)*0+ic+1);
								Vbw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic-1);
								Vbe = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(kc-1)+(NIC-1)*0+ic+1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.25*(Vp+Vb+Vw+Vbw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i) = 0.25*(Vp+Vt+Vw+Vtw);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.25*(Vp+Vb+Ve+Vbe);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(k+1)+(NI-1)*0+i+1) = 0.25*(Vp+Vt+Ve+Vte);
								ic++;
							}
							kc++;
						}
					}
				break;
				case 4:
					if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
						
						ic = NIC-2;i = NI-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
						
						ic = NIC-2;i = NI-3;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = Vp;
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-3;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-3;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
						
						ic = NIC-2;i = NI-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						
						ic = NIC-2;i = NI-2;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-2;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-3;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
					else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						
						ic = NIC-2;i = NI-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
						
						ic = NIC-2;i = NI-3;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-3;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-2;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						
						ic = NIC-2;i = NI-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						
						ic = NIC-2;i = NI-2;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-2;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-2;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*(NKC-2)+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
				break;
				case 5:
					if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
						
						ic = NIC-2;i = NI-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
						
						ic = NIC-2;i = NI-3;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = Vp;
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-3;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-3;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
						
						ic = NIC-2;i = NI-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						
						ic = NIC-2;i = NI-2;
						jc = NJC-2;j = NJ-3;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-2;
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-3;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Ve);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-3;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
					else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						
						ic = NIC-2;i = NI-3;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
						
						ic = NIC-2;i = NI-3;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-3;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.5*(Vp+Vn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-2;
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
					{
						ic = 0;i = 0;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = Vp;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
						
						ic = 0;i = 0;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
						Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
						
						ic = NIC-2;i = NI-2;
						jc = 0;j = 0;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
						Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
						
						ic = NIC-2;i = NI-2;
						jc = NJC-2;j = NJ-2;
						
						Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
						Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
						Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
						Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
						
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
						
						ic = 0;i = 0;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vs);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.5*(Vp+Vn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							jc++;
						}

						ic = NIC-2;i = NI-2;
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							jc++;
						}
						
						jc = 0;j = 0;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
							Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.5*(Vp+Vw);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.5*(Vp+Ve);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
							ic++;
						}

						jc = NJC-2;j = NJ-2;
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
							Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
							Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
							Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
							Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
							Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
							
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
							ic++;
						}

						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							jc = 1;
							for(j=2;j<NJ-2;j=j+2)
							{
								Vp = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic);
								Vs = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic);
								Vn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic);
								Vw = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic-1);
								Ve = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*jc+ic+1);
								Ves = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic+1);
								Ven = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic+1);
								Vws = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc-1)+ic-1);
								Vwn = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*0+(NIC-1)*(jc+1)+ic-1);
								
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = 0.25*(Vp+Vw+Vs+Vws);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i) = 0.25*(Vp+Vw+Vn+Vwn);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*j+i+1) = 0.25*(Vp+Ve+Vs+Ves);
								*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*0+(NI-1)*(j+1)+i+1) = 0.25*(Vp+Ve+Vn+Ven);
								jc++;
							}
							ic++;
						}
					}
				break;
			}
		}
	}
	return 0;
}
