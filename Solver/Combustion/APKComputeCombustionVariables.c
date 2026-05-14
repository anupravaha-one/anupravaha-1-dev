#include "APKComputeCombustionVariables.h"
/**
Updation of different combustion vairables is done here.  including, fuel,Oxygen,Product, mass fracions, enthalpy, C_p etc are also updated in this function.
**/
int ComputeCombustionVariables(Base *CgnsBase)
{
	int i,j,k,m,NI,NJ,NK,zcount,*size,Spcount,totalcells,nzones;
	int Numspecies,index,numIface,iface,numZbc,zBc,dir,*ptrRange,flagrho;
	double Mix=0.0,Mixst=0.0,FMassFrac=0.0,OxMassFrac=0.0,IMassFrac=0.0,PMassFrac=0.0,enthalpy=0.0;
	double StConst=0.0,Temp=0.0,C1=0.0,C2=0.0,*nMoles,*SpMolWt,*Cpi,MolWtrecp,netMolWt,netCp,P,InertMFrac[10];
	int NIngases;

	double a,b,c,d,e,f,g,h,T,Tref,href,*Cpref,netCpref;                                                   //nitesh

	Zone *z;
	nzones = CgnsBase->Number_of_zones;
	NIngases = CgnsBase->BaseParameters->NIngases;
	z = CgnsBase->CgnsZone;

	StConst = *(CgnsBase->BaseParameters->Stoich);
        Numspecies = CgnsBase->BaseParameters->Nspecies;
	nMoles = CgnsBase->BaseParameters->Moles;
        SpMolWt = CgnsBase->BaseParameters->MolWt;
	Tref = *(CgnsBase->BaseParameters->InitialTemp);                                              //nitesh 
	href = *(CgnsBase->BaseParameters->enthalpyOxiInlet);                                         //nitesh
	Cpref = CgnsBase->BaseParameters->Cp_ref;                                                     //nitesh
	if(CgnsBase->BaseParameters->flagfueltype==0) //Single fuel case
        {
                C1 = (*(nMoles+Numspecies-2))*(*(SpMolWt+Numspecies-2));
                C2 = (*(nMoles+Numspecies-1))*(*(SpMolWt+Numspecies-1));
        }
        
        for(i=0;i<10;i++)//Initialization of InertMFrac[10] to zero
        {
        	InertMFrac[i] = 0.0;
        }

	if(CgnsBase->flagDiffKinetCont==1) //Diffusion controlled reaction
	{
        	Mixst = *(CgnsBase->BaseParameters->Zst);
        	for(zcount=0;zcount<nzones;zcount++)
        	{
			size = (z+zcount)->Size;
			NI =*(size+0);
			NJ = *(size+1);
			NK = *(size+2);

                        for(k=0;k<(NK-1);k++) // In the domain
			for(j=0;j<(NJ-1);j++)
			for(i=0;i<(NI-1);i++)
			{
                        	Mix = *(((z+zcount)->Flowsol)->MixFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i);

                                if(Mix >= 0.0 && Mix <= Mixst)
                                {
					*(((z+zcount)->Flowsol)->FuelMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = 0.0;
                                        FMassFrac = 0.0;
					*(((z+zcount)->Flowsol)->OxiMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                        OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                }
                                else if(Mix > Mixst && Mix <= 1.0)
                                {
					*(((z+zcount)->Flowsol)->FuelMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = ((Mix - Mixst)/(1.0 - Mixst))*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
					*(((z+zcount)->Flowsol)->OxiMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = 0.0;
                                        OxMassFrac = 0.0;
                                }

                                *(((z+zcount)->Flowsol)->enthalpy+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
                                enthalpy = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
                                IMassFrac = 0.0;
                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                {
                                *(((z+zcount)->Flowsol)->InMassFrac[m]+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                }
				*(((z+zcount)->Flowsol)->ProdMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);


                        }

                        numIface=(z+zcount)->Number_of_Interfaces;
			numZbc=(z+zcount)->Number_of_BC_Patches;

                        for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numIface;iface++)
				{
					if((((z+zcount)->Iface+iface)->Direction)==dir)
					{

						ptrRange=((z+zcount)->Iface+iface)->Point_Range;

						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
							for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

								Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);

                                                                if(Mix >= 0.0 && Mix <= Mixst)
                                				{
									*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index) = 0.0;
                                        				FMassFrac = 0.0;
									*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                        				OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                				}
                                				else if(Mix > Mixst && Mix <= 1.0)
                                				{
									*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index) = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        				FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
									*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = 0.0;
                                        				OxMassFrac = 0.0;
                                				}

                                				*((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
                                                                IMassFrac = 0.0;
                                                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                				*((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                				IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                }
								*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

							}
						}
						else if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);

                                                                if(Mix >= 0.0 && Mix <= Mixst)
                                				{
									*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index) = 0.0;
                                        				FMassFrac = 0.0;
									*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                        				OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                				}
                                				else if(Mix > Mixst && Mix <= 1.0)
                                				{
									*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index) = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        				FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
									*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = 0.0;
                                        				OxMassFrac = 0.0;
                                				}

                                				*((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
								IMassFrac = 0.0;
                                                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                				*((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                				IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                }
								*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
							}
						}
						else if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);

                                                                if(Mix >= 0.0 && Mix <= Mixst)
                                				{
									*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index) = 0.0;
                                        				FMassFrac = 0.0;
									*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                        				OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                				}
                                				else if(Mix > Mixst && Mix <= 1.0)
                                				{
									*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index) = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        				FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
									*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = 0.0;
                                        				OxMassFrac = 0.0;
                                				}

                                				*((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
								IMassFrac = 0.0;
                                                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                				*((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                				IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                }
								*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

							}
						}
					}
				}
				for(zBc=0;zBc<numZbc;zBc++)
				{
					if(((z+zcount)->ZoneBc+zBc)->Direction==dir)
					{


						ptrRange=((z+zcount)->ZoneBc+zBc)->Point_Range;

						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
								Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);

                                                                if(Mix >= 0.0 && Mix <= Mixst)
                                				{
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index) = 0.0;
                                        				FMassFrac = 0.0;
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                        				OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                				}
                                				else if(Mix > Mixst && Mix <= 1.0)
                                				{
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index) = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        				FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = 0.0;
                                        				OxMassFrac = 0.0;
                                				}

                                				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
                                				IMassFrac = 0.0;
                                                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                				IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                }
								*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

							}
						}
						else if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);

                                                                if(Mix >= 0.0 && Mix <= Mixst)
                                				{
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index) = 0.0;
                                        				FMassFrac = 0.0;
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                        				OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                				}
                                				else if(Mix > Mixst && Mix <= 1.0)
                                				{
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index) = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        				FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = 0.0;
                                        				OxMassFrac = 0.0;
                                				}

                                				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
                                				IMassFrac = 0.0;
                                                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                				IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                }
								*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
							}
						}
						else if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);

                                                                if(Mix >= 0.0 && Mix <= Mixst)
                                                                {
                                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index) = 0.0;
                                                                        FMassFrac = 0.0;
                                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                                                        OxMassFrac = (1.0 - Mix/Mixst)*(*(CgnsBase->BaseParameters->OxiMassFracInlet));
                                                                }
                                				else if(Mix > Mixst && Mix <= 1.0)
                                				{
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index) = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
                                        				FMassFrac = (Mix - Mixst)/(1.0 - Mixst)*(*(CgnsBase->BaseParameters->FuelMassFracInlet));
									*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = 0.0;
                                        				OxMassFrac = 0.0;
                                				}

                                				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index) = (*(CgnsBase->BaseParameters->enthalpyFuelInlet))*Mix + (1.0 - Mix)*(*(CgnsBase->BaseParameters->enthalpyOxiInlet));
                                				IMassFrac = 0.0;
                                                                for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                				IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                }
								*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

							}
						}
					}
				}
			}



        	}

	}

        else if(CgnsBase->flagDiffKinetCont==2) //Kinetically Controlled
        {

        	for(zcount=0;zcount<nzones;zcount++)
        	{
			size = (z+zcount)->Size;
			NI =*(size+0);
			NJ = *(size+1);
			NK = *(size+2);



                        for(k=0;k<(NK-1);k++)
			for(j=0;j<(NJ-1);j++)
			for(i=0;i<(NI-1);i++)
			{
        			if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                		{
					Mix = *(((z+zcount)->Flowsol)->MixFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
					FMassFrac = *(((z+zcount)->Flowsol)->FuelMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
                        		enthalpy =  *(((z+zcount)->Flowsol)->enthalpy+(NI-1)*(NJ-1)*k+(NI-1)*j+i);

                        		OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        		*(((z+zcount)->Flowsol)->OxiMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = OxMassFrac;
                        		IMassFrac = 0.0;
                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                        {
                                        *(((z+zcount)->Flowsol)->InMassFrac[m]+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                        }
                         		*(((z+zcount)->Flowsol)->ProdMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

                		}
                                /*else
                                {
					Mix = *(((z+zcount)->Flowsol)->MixFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
					FMassFrac = *(((z+zcount)->Flowsol)->FuelMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
                        		enthalpy =  *(((z+zcount)->Flowsol)->enthalpy+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
                                        IMassFrac = *(((z+zcount)->Flowsol)->InMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i);

                                        OxMassFrac = StConst*(FMassFrac - Mix);
                                        *(((z+zcount)->Flowsol)->OxiMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = OxMassFrac;
                                        *(((z+zcount)->Flowsol)->ProdMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
				}*/
			}

			numIface=(z+zcount)->Number_of_Interfaces;
			numZbc=(z+zcount)->Number_of_BC_Patches;

                        for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numIface;iface++)
				{
					if((((z+zcount)->Iface+iface)->Direction)==dir)
					{
						ptrRange=((z+zcount)->Iface+iface)->Point_Range;

						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
							for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

								if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                						{
									Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);
									FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index);

                        						OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        						*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                        						IMassFrac = 0.0;
                                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                        {
                                                                        *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        }
                         						*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

                						}
                                				/*else
                                				{
									Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);
									FMassFrac =*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index);
                                        				IMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac + index);


                                        				OxMassFrac = StConst*(FMassFrac - Mix);
                                        				*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                                        				*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				}*/

							}
						}
						else if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                                if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                						{
									Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);
									FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index);

                        						OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        						*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = OxMassFrac;

                        						IMassFrac = 0.0;
                                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                        {
                                                                        *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        }
                                                                	*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

                						}
                                				/*else
                                				{
									Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);
									FMassFrac =*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index);
                                        				IMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac + index);


                                        				OxMassFrac = StConst*(FMassFrac - Mix);
                                        				*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                                        				*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				}*/

							}
						}
						else if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
                                                                if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                						{
									Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);
									FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index);

                        						OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        						*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = OxMassFrac;

                        						IMassFrac = 0.0;
                                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                        {
                                                                        *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        }
                                                                	*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                						}
                                				/*else
                                				{
									Mix = *((((z+zcount)->Iface+iface)->ZSlab)->MixFrac + index);
									FMassFrac =*((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy + index);
                                        				IMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac + index);

                                        				OxMassFrac = StConst*(FMassFrac - Mix);
                                        				*((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                                        				*((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				}*/

							}
						}
					}
				}
				for(zBc=0;zBc<numZbc;zBc++)
				{
					if(((z+zcount)->ZoneBc+zBc)->Direction==dir)
					{

						ptrRange=((z+zcount)->ZoneBc+zBc)->Point_Range;

						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
								if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                						{
									Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);
									FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index);

                        						OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        						*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = OxMassFrac;

                        						IMassFrac = 0.0;
                                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                        {
                                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        }
                                                                	*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

                						}
                                				/*else
                                				{
									Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);
									FMassFrac =*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index);
                                        				IMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac + index);


                                        				OxMassFrac = StConst*(FMassFrac - Mix);
                                        				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                                        				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				}*/

							}
						}
						else if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                						{
									Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);
									FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index);

                        						OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        						*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = OxMassFrac;

                        						IMassFrac = 0.0;
                                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                        {
                                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        }
                                                                	*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

                						}
                                				/*else
                                				{
									Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);
									FMassFrac =*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index);
                                        				IMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac + index);


                                        				OxMassFrac = StConst*(FMassFrac - Mix);
                                        				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                                        				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				}*/

							}
						}
						else if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								if(CgnsBase->flagPreNonPre==2) //Non-Premixed
                						{
									Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);
									FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index);

                        						OxMassFrac = StConst*(FMassFrac - Mix*(*(CgnsBase->BaseParameters->FuelMassFracInlet)) - (Mix - 1.0)*(*(CgnsBase->BaseParameters->OxiMassFracInlet))/StConst);
                        						*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = OxMassFrac;

                        						IMassFrac = 0.0;
                                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                        {
                                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m] + index) = (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        IMassFrac = IMassFrac + (1.0 - Mix)*(*(CgnsBase->BaseParameters->InertgasOxiInlet[m])) + Mix*(*(CgnsBase->BaseParameters->InertgasFuelInlet[m]));
                                                                        }
                                                                	*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);

                						}
                                				/*else
                                				{
									Mix = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MixFrac + index);
									FMassFrac =*((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac + index);
                        						enthalpy =  *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy + index);
                                        				IMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac + index);


                                        				OxMassFrac = StConst*(FMassFrac - Mix);
                                        				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac + index) = OxMassFrac;
                                        				*((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac + index) = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                        				PMassFrac = 1.0 - (FMassFrac + OxMassFrac + IMassFrac);
                                				}*/


							}
						}
					}
				}
			}

		}

        }


	//Calculation of Cps of species in the domain and slabs
       	for(zcount=0;zcount<nzones;zcount++)
        {
		size = (z+zcount)->Size;
		NI =*(size+0);
		NJ = *(size+1);
		NK = *(size+2);

		totalcells = (NI-1)*(NJ-1)*(NK-1);

		for(Spcount=0;Spcount<Numspecies;Spcount++)
                for(k=0;k<(NK-1);k++)
		for(j=0;j<(NJ-1);j++)
		for(i=0;i<(NI-1);i++)
		{
                	index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			Temp = *(((z+zcount)->Flowsol)->T+index);
                        //*((((z+zcount)->Flowsol)->Tcoeff)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
			
			a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
			b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
			c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
			d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
			e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
			f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
			g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
			h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
			T=Temp;
			
			*((((z+zcount)->Flowsol)->Tcoeff)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0 + h*T*T*T*T*T*T*T/8.0;          //nitesh
			
		}




		numIface=(z+zcount)->Number_of_Interfaces;
		numZbc=(z+zcount)->Number_of_BC_Patches;

                for(dir=1;dir<=6;dir++)
                {
                        for(iface=0;iface<numIface;iface++)
                        {
                                if((((z+zcount)->Iface+iface)->Direction)==dir)
                                {

                                        ptrRange=((z+zcount)->Iface+iface)->Point_Range;


                                        if(dir==1||dir==2)
                                        {
                                        	totalcells = (*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+4)- *(ptrRange+1));
                                                for(Spcount=0;Spcount<Numspecies;Spcount++)
                                                for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
                                                for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

                                                        Temp = *((((z+zcount)->Iface+iface)->ZSlab)->T+index);
                                                        //*((((z+zcount)->Iface+iface)->ZSlab)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
							
							a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
							b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
							T=Temp;
							
							*((((z+zcount)->Iface+iface)->ZSlab)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0  + h*T*T*T*T*T*T*T/8.0;         //nitesh

                                                }
                                        }
                                        else if(dir==3||dir==4)
                                        {
                                        	totalcells=(*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(Spcount=0;Spcount<Numspecies;Spcount++)
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                        Temp = *((((z+zcount)->Iface+iface)->ZSlab)->T+index);
                                                        // *((((z+zcount)->Iface+iface)->ZSlab)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
						       
						       	a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
							b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
							T=Temp;
							
							*((((z+zcount)->Iface+iface)->ZSlab)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0 + h*T*T*T*T*T*T*T/8.0;          //nitesh


                                                }
                                        }
                                        else if(dir==5||dir==6)
                                        {
                                        	totalcells=(*(ptrRange+4)- *(ptrRange+1)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(Spcount=0;Spcount<Numspecies;Spcount++)
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
                                                        Temp = *((((z+zcount)->Iface+iface)->ZSlab)->T+index);
                                                        //*((((z+zcount)->Iface+iface)->ZSlab)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
							
							a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
							b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
							T=Temp;
							
							*((((z+zcount)->Iface+iface)->ZSlab)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0   + h*T*T*T*T*T*T*T/8.0;        //nitesh



                                                }
                                        }
                                }
                        }
                        for(zBc=0;zBc<numZbc;zBc++)
                        {
                                if(((z+zcount)->ZoneBc+zBc)->Direction==dir)
                                {

                                        ptrRange=((z+zcount)->ZoneBc+zBc)->Point_Range;

                                        if(dir==1||dir==2)
                                        {
                                        	totalcells = (*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+4)- *(ptrRange+1));
                                                for(Spcount=0;Spcount<Numspecies;Spcount++)
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

                                                        Temp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index);
                                                        //*((((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
							
							a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
							b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
							T=Temp;
							
							*((((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0    + h*T*T*T*T*T*T*T/8.0;       //nitesh


                                                }
                                        }
                                        else if(dir==3||dir==4)
                                        {
                                        	totalcells=(*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(Spcount=0;Spcount<Numspecies;Spcount++)
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
                                                        Temp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index);
                                                        //*((((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
							
							
							a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
							b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
							T=Temp;
							
							*((((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0     + h*T*T*T*T*T*T*T/8.0;      //nitesh


                                                }
                                        }
                                        else if(dir==5||dir==6)
                                        {
                                                totalcells=(*(ptrRange+4)- *(ptrRange+1)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(Spcount=0;Spcount<Numspecies;Spcount++)
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                        Temp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index);
                                                        //*((((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps + Spcount*totalcells + index) = *(CgnsBase->BaseParameters->CpCoeff+Spcount*3+0) + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+1))*Temp + (*(CgnsBase->BaseParameters->CpCoeff+Spcount*3+2))*Temp*Temp;
							
							a=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+0);
							b=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+1);
							c=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+2);
							d=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+3);
							e=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+4);
							f=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+5);
							g=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+6);
							h=*(CgnsBase->BaseParameters->CpCoeff+Spcount*8+7);
							T=Temp;
							
							*((((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps + Spcount*totalcells + index) = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f*T*T*T*T*T/6.0 + g*T*T*T*T*T*T/7.0    + h*T*T*T*T*T*T*T/8.0;        //nitesh


                                                }
                                        }
                                }
                        }
                }

	}

        //Calculation of net Cp ,Mol.Wt
	for(zcount=0;zcount<nzones;zcount++)
        {
		size = (z+zcount)->Size;
		NI =*(size+0);
		NJ = *(size+1);
		NK = *(size+2);

                totalcells = (NI-1)*(NJ-1)*(NK-1);
                Cpi = (((z+zcount)->Flowsol)->Tcoeff)->Cps;
                flagrho = ((z+zcount)->Flowsol)->flagRho;

                for(k=0;k<(NK-1);k++) //In the domain
		for(j=0;j<(NJ-1);j++)
		for(i=0;i<(NI-1);i++)
                {
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;

			FMassFrac = *(((z+zcount)->Flowsol)->FuelMassFrac+index);
			OxMassFrac = *(((z+zcount)->Flowsol)->OxiMassFrac+index);
                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                        {
			InertMFrac[m] = *(((z+zcount)->Flowsol)->InMassFrac[m]+index);
                        }
			PMassFrac = *(((z+zcount)->Flowsol)->ProdMassFrac+index);
                        Temp = *(((z+zcount)->Flowsol)->T+index);
                        P = *(((z+zcount)->Flowsol)->P+index);

			if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                        {
				netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
				
				for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                {
                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
				}
                        }
                        else //Mixture of fuels case
                        {
				netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
				
				for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                {
                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
				}
                        }
			*(((z+zcount)->Flowsol)->SpecificHeat+index) = netCp;
			if(CgnsBase->BaseParameters->flagfueltype == 0)
                        {
                        	MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                        	for(m=0;m<NIngases;m++)
                                {
					MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                }
			}
                        else //Mixture of fuels case
                        {
                       		MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                for(m=0;m<NIngases;m++)
                                {
					MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                }
                        }
                        *((((z+zcount)->Flowsol)->Tcoeff)->MolWteq + index) = 1.0/MolWtrecp;
                        netMolWt = 1.0/MolWtrecp;


                }

                numIface=(z+zcount)->Number_of_Interfaces;
                numZbc=(z+zcount)->Number_of_BC_Patches;

                for(dir=1;dir<=6;dir++) //In the slabs
                {
                        for(iface=0;iface<numIface;iface++)
                        {
                                if((((z+zcount)->Iface+iface)->Direction)==dir)
                                {

                                        ptrRange=((z+zcount)->Iface+iface)->Point_Range;
                                        Cpi = (((z+zcount)->Iface+iface)->ZSlab)->Cps;

                                        if(dir==1||dir==2)
                                        {

                                                totalcells = (*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+4)- *(ptrRange+1));
                                                for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
                                                for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

                                                        FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac+index);
                                                        OxMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac+index);
                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                        {
                                                        InertMFrac[m] = *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m]+index);
                                                        }
                                                        PMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac+index);
                                                        Temp = *((((z+zcount)->Iface+iface)->ZSlab)->T+index);
                                                        P = *((((z+zcount)->Iface+iface)->ZSlab)->P+index);

							if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        else
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        *((((z+zcount)->Iface+iface)->ZSlab)->SpecificHeat+index) = netCp;

                                                        if(CgnsBase->BaseParameters->flagfueltype == 0)
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        else
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        *((((z+zcount)->Iface+iface)->ZSlab)->MolWteq + index) = 1.0/MolWtrecp;
                                                        netMolWt = 1.0/MolWtrecp;


                                                }
                                        }
                                        else if(dir==3||dir==4)
                                        {

                                                totalcells=(*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                        FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac+index);
                                                        OxMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac+index);
                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                        {
                                                        InertMFrac[m] = *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m]+index);
                                                        }
                                                        PMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac+index);
                                                        Temp = *((((z+zcount)->Iface+iface)->ZSlab)->T+index);
							P = *((((z+zcount)->Iface+iface)->ZSlab)->P+index);
							if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        else
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        *((((z+zcount)->Iface+iface)->ZSlab)->SpecificHeat+index) = netCp;

                                                        if(CgnsBase->BaseParameters->flagfueltype == 0)
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        else
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        *((((z+zcount)->Iface+iface)->ZSlab)->MolWteq + index) = 1.0/MolWtrecp;
                                                        netMolWt = 1.0/MolWtrecp;


                                                }
                                        }
                                        else if(dir==5||dir==6)
                                        {
                                                totalcells=(*(ptrRange+4)- *(ptrRange+1)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                        FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac+index);
                                                        OxMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->OxiMassFrac+index);
                                                       	for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                        {
                                                        InertMFrac[m] = *((((z+zcount)->Iface+iface)->ZSlab)->InMassFrac[m]+index);
                                                        }
                                                        PMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->ProdMassFrac+index);
                                                        Temp = *((((z+zcount)->Iface+iface)->ZSlab)->T+index);
							P = *((((z+zcount)->Iface+iface)->ZSlab)->P+index);

							if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        else
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        *((((z+zcount)->Iface+iface)->ZSlab)->SpecificHeat+index) = netCp;

                                                        if(CgnsBase->BaseParameters->flagfueltype == 0)
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        else
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        *((((z+zcount)->Iface+iface)->ZSlab)->MolWteq + index) = 1.0/MolWtrecp;

                                                        netMolWt = 1.0/MolWtrecp;


                                                }
                                        }
                                }
                        }
                        for(zBc=0;zBc<numZbc;zBc++)
                        {
                                if(((z+zcount)->ZoneBc+zBc)->Direction==dir)
                                {

                                        ptrRange=((z+zcount)->ZoneBc+zBc)->Point_Range;
                                        Cpi = (((z+zcount)->ZoneBc+zBc)->ZSlab)->Cps;

                                        if(dir==1||dir==2)
                                        {
                                                totalcells = (*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+4)- *(ptrRange+1));
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));

                                                        FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac+index);
                                                        OxMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac+index);
                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                        {
                                                        InertMFrac[m] = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m]+index);
                                                        }
                                                        PMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac+index);
                                                        Temp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index);
                                                        P = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->P+index);

							if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        else
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->SpecificHeat+index) = netCp;

                                                        if(CgnsBase->BaseParameters->flagfueltype == 0)
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        else
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MolWteq + index) = 1.0/MolWtrecp;

                                                        netMolWt = 1.0/MolWtrecp;

                                                }
                                        }
                                        else if(dir==3||dir==4)
                                        {

                                                totalcells=(*(ptrRange+5)- *(ptrRange+2)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                        FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac+index);
                                                        OxMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac+index);
                                                        for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                        {
                                                        InertMFrac[m] = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m]+index);
                                                        }
                                                        PMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac+index);
                                                        Temp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index);
                                                        P = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->P+index);

							if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        else
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->SpecificHeat+index) = netCp;

                                                        if(CgnsBase->BaseParameters->flagfueltype == 0)
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        else
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MolWteq + index) = 1.0/MolWtrecp;

                                                        netMolWt = 1.0/MolWtrecp;

                                                }
                                        }
                                        else if(dir==5||dir==6)
                                        {

                                                totalcells=(*(ptrRange+4)- *(ptrRange+1)) * (*(ptrRange+3)- *(ptrRange+0));
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));

                                                        FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac+index);
                                                        OxMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->OxiMassFrac+index);
							for(m=0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                        {
                                                        InertMFrac[m] = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->InMassFrac[m]+index);
                                                        }
                                                        PMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->ProdMassFrac+index);
                                                        Temp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index);
                                                        P = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->P+index);

							if(CgnsBase->BaseParameters->flagfueltype == 0) //Single fuel
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*C1/(C1+C2)*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*C2/(C1+C2)*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        else
                                                        {
                                                                netCp = FMassFrac*(*(Cpi + (totalcells*0) + index)) + OxMassFrac*(*(Cpi+(totalcells*1) + index)) + PMassFrac*(*(Cpi + (totalcells*(2+NIngases))+ index)) + PMassFrac*(*(Cpi+(totalcells*(3+NIngases)) + index));
								
								for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                                                                {
                                                                netCp = netCp + InertMFrac[m]*(*(Cpi+(totalcells*(2+m)) + index));
								}
                                                        }
                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->SpecificHeat+index) = netCp;

                                                        if(CgnsBase->BaseParameters->flagfueltype == 0)
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + C1/(C1+C2)*PMassFrac/(*(SpMolWt+2+NIngases)) + C2/(C1+C2)*PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        else
                                                        {
                                                                MolWtrecp = FMassFrac/(*(SpMolWt+0)) + OxMassFrac/(*(SpMolWt+1)) + PMassFrac/(*(SpMolWt+2+NIngases)) + PMassFrac/(*(SpMolWt+3+NIngases));
                                                                for(m=0;m<NIngases;m++)
                                                                {
                                                                        MolWtrecp = MolWtrecp + InertMFrac[m]/(*(SpMolWt+2+m));
                                                                }
                                                        }
                                                        *((((z+zcount)->ZoneBc+zBc)->ZSlab)->MolWteq + index) = 1.0/MolWtrecp;

                                                        netMolWt = 1.0/MolWtrecp;


                                                }
                                        }
                                }
                        }
                }



        }

        //Calculating new temperature
	
	for(zcount=0;zcount<nzones;zcount++)
        {
		size = (z+zcount)->Size;
		NI =*(size+0);
		NJ = *(size+1);
		NK = *(size+2);

                totalcells = (NI-1)*(NJ-1)*(NK-1);

                for(k=0;k<(NK-1);k++) //In the domain
		for(j=0;j<(NJ-1);j++)
		for(i=0;i<(NI-1);i++)
                {
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			
			netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
			for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        {
                		netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
			}
			
			FMassFrac = *(((z+zcount)->Flowsol)->FuelMassFrac+index);

			netCp = *(((z+zcount)->Flowsol)->SpecificHeat+index);

                        enthalpy = *(((z+zcount)->Flowsol)->enthalpy+index);

                        //*(((z+zcount)->Flowsol)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
			
			*(((z+zcount)->Flowsol)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;      //nitesh
			
			//if(i==95 && j==5)
			//printf("temp=%lf \t ",*(((z+zcount)->Flowsol)->T+index));                                      //nitesh

                }

                numIface=(z+zcount)->Number_of_Interfaces;
                numZbc=(z+zcount)->Number_of_BC_Patches;

                for(dir=1;dir<=6;dir++) //In the slabs
                {
                        for(iface=0;iface<numIface;iface++)
                        {
                                if((((z+zcount)->Iface+iface)->Direction)==dir)
                                {

                                        ptrRange=((z+zcount)->Iface+iface)->Point_Range;


                                        if(dir==1||dir==2)
                                        {


                                                for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
                                                for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
							
							netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
							for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        				{
                						netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
							}
                                
							FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac+index);
                                                        
							netCp = *((((z+zcount)->Iface+iface)->ZSlab)->SpecificHeat+index);
                                                        
							enthalpy = *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy+index);

                                                        //*((((z+zcount)->Iface+iface)->ZSlab)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
							
							*((((z+zcount)->Iface+iface)->ZSlab)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;        //nitesh


                                                }
                                        }
                                        else if(dir==3||dir==4)
                                        {


                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							
							netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
							for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        				{
                						netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
							}
                                
							FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac+index);
                                                        
							netCp = *((((z+zcount)->Iface+iface)->ZSlab)->SpecificHeat+index);
                                                        
							enthalpy = *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy+index);

                                                        //*((((z+zcount)->Iface+iface)->ZSlab)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
							
							*((((z+zcount)->Iface+iface)->ZSlab)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;        //nitesh



                                                }
                                        }
                                        else if(dir==5||dir==6)
                                        {

                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							
							netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
							for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        				{
                						netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
							}
                                
							FMassFrac = *((((z+zcount)->Iface+iface)->ZSlab)->FuelMassFrac+index);
                                                        
							netCp = *((((z+zcount)->Iface+iface)->ZSlab)->SpecificHeat+index);
                                                        
							enthalpy = *((((z+zcount)->Iface+iface)->ZSlab)->enthalpy+index);

                                                        //*((((z+zcount)->Iface+iface)->ZSlab)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
							
							*((((z+zcount)->Iface+iface)->ZSlab)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;        //nitesh
							

                                                }
                                        }
                                }
                        }
                        for(zBc=0;zBc<numZbc;zBc++)
                        {
                                if(((z+zcount)->ZoneBc+zBc)->Direction==dir)
                                {

                                        ptrRange=((z+zcount)->ZoneBc+zBc)->Point_Range;

                                        if(dir==1||dir==2)
                                        {

                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
							
							netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
							for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        				{
                						netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
							}
                                
							FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac+index);
                                                        
							netCp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->SpecificHeat+index);
                                                        
							enthalpy = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy+index);

                                                        //*((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
														
							*((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;        //nitesh
                                                }
                                        }
                                        else if(dir==3||dir==4)
                                        {
                                                for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							
							netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
							for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        				{
                						netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
							}
                                
							FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac+index);
                                                        
							netCp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->SpecificHeat+index);
                                                        
							enthalpy = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy+index);

                                                        //*((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
														
							*((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;        //nitesh

                                                }
                                        }
                                        else if(dir==5||dir==6)
                                        {
                                                for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
                                                for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
                                                {
                                                        index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
							
							netCpref = *(CgnsBase->BaseParameters->OxiMassFracInlet) *(*(Cpref + 1));                                                                   //nitesh
                        
							for(m = 0;m<CgnsBase->BaseParameters->NIngases;m++)
                        				{
                						netCpref = netCpref + *(CgnsBase->BaseParameters->InertgasOxiInlet[m])*(*(Cpref +(2+m)));
							}
                                
							FMassFrac = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->FuelMassFrac+index);
                                                        
							netCp = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->SpecificHeat+index);
                                                        
							enthalpy = *((((z+zcount)->ZoneBc+zBc)->ZSlab)->enthalpy+index);

                                                        //*((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index) = (enthalpy - FMassFrac*(*(CgnsBase->BaseParameters->delH)))/netCp;
														
							*((((z+zcount)->ZoneBc+zBc)->ZSlab)->T+index) = (enthalpy - href - FMassFrac*(*(CgnsBase->BaseParameters->delH)) + netCpref * Tref )/netCp;        //nitesh

	                                        }
                                        }
                                }
                        }
                }
        }
return(0);

}


