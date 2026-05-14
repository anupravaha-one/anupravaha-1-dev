#include "APKWireCurrentDensity.h"
/**
This function WireCurrentDensity() is used for calculating the value of total
current reaching the plate (boundary tag = 24 or wall viscous isothermal) and
wire (boundary tag = 22 or wall viscous boundary). This function also calculates
the current density at each point on the grounded plate and prints it in ".dat"
file.

*/
int WireCurrentDensity(Base *base,int vartag)
{
	FILE *CharDens, *PlateCur;
	int iface,NI,NJ,NK,*ptrSize,nZone,INDEX;
	int i,j,k,count,dir;
	char tag,name[33];
	void *ptr;
	double Fdiff[6],Rdiff[6];
	double hprev, hlat, TotalCurrent, Chargeprev, Charge, centroidprev, PlateCurrent=0.0, WireCurrent = 0.0;
	double *slabPotential[6],*slabRho[6],*slabVolume[6];
	double *vols,*beta,*Potential,*Rho,diffdia,*surf,Kcond,Chargearea;
	double p,q,r,Area,DiffCoef;
	double *coordX,*coordY,*coordZ;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3], CoeffN[6];
	int index[3];
	Zone *zone;
	Interface *facePtr;
	ZoneBoundary *bcPtr; 
	Zone *FindZonePtr();
	//CharDens = fopen("Settings/ESPOutput/CurrentDensity.dat","w");
	PlateCur = fopen("Settings/ESPOutput/PlateCurrent.dat", "w");
	Chargearea=0.0;
	count = base->Number_of_zones;
	
	for(nZone=0;nZone<base->Number_of_zones;nZone++)
	{
		
		
		
		zone=(base->CgnsZone+nZone);
		surf= zone->ZGeom->SurfaceArea;
		
		coordX=zone->CgnsGrid->X;
		coordY=zone->CgnsGrid->Y; 
		coordZ=zone->CgnsGrid->Z; 
		
		
		
		Potential=MakeDomainArray(zone,8);//variableDomain stores Potential pointer
		Rho=MakeDomainArray(zone,9);
		
		ptrSize=zone->Size;
		//F = zone->Flowsol->Massflux;
		NI = *(ptrSize+0);
		NJ = *(ptrSize+1);
		NK = *(ptrSize+2);
		//printf("I clear the check\n");
		for(iface=0;iface<6;iface++)
		{
			ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
			if(tag=='i')
			{
				facePtr=(Interface *)ptr;
			}
			else if(tag=='b') 
			{
				bcPtr=(ZoneBoundary *)ptr;
				if((bcPtr->Boundary_Tag == 22)||(bcPtr->Boundary_Tag == 24))
				{
				
					MakeSlabPtrArray(zone,slabPotential,8);/*Called for making slabs of u */
					MakeSlabPtrArray(zone,slabRho,9);
					MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/

					dir=bcPtr->Direction;
					vols= zone->ZGeom->Volume; 
					beta= zone->ZGeom->Beta;
					switch(dir)
					{
						case 1://East Face
							Charge=0.0;TotalCurrent = 0.0;hlat = 0.0; hprev = 0.0; 
							Chargeprev = 0.0; centroid[0]=0.0; centroidprev = 0.0;
							CharDens = fopen("Settings/ESPOutput/CurrentDensity_East.dat","a");
							fprintf(CharDens,"VARIABLES = \"CoordinateX\"\n\"CoordinateY\"\n\"CoordinateZ\"\n\"Current Density\"\n\n");
							for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							{
								i=NI-2;
								index[0]=k,index[1]=j,index[2]=i;
								ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
								ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
								ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
								FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

								//Get the Diffusion of Charge Density at the boundary
								diffdia = Diffusion(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff);

								//Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
								diffdia = Diffusion(index,Rho,vols,ptrSize,slabRho,slabVolume,beta,Rdiff);

								INDEX = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
								
								Kcond = *(zone->Flowsol->MobIon)*(*((bcPtr->ZSlab)->RhoIon+((NJ-1)*k+j)));//
								DiffCoef = *(zone->Flowsol->ChargeDiffcoeff);
								p=*(surf+(INDEX)*6*3+3*1+0);
								q=*(surf+(INDEX)*6*3+3*1+1);
								r=*(surf+(INDEX)*6*3+3*1+2);
								Area=sqrt(p*p+q*q+r*r);
								
								Charge=-Kcond*Fdiff[1]/Area;
								//TotalCurrent = TotalCurrent - Kcond*Fdiff[1] - DiffCoef*Rdiff[1];
								TotalCurrent = TotalCurrent - Kcond*Fdiff[1];
								if(bcPtr->Boundary_Tag == 24)
									fprintf(CharDens,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\n",centroid[0],centroid[1],centroid[2],Charge);
								
							}
							
							
							if(bcPtr->Boundary_Tag == 22)
							{
								WireCurrent = WireCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nEast Wire = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On East Wire = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							else
							{
								PlateCurrent = PlateCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nEast Plate = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On East Plate = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							fclose(CharDens);
						break;
						case 2://West Face
							Charge=0.0;TotalCurrent = 0.0;hlat = 0.0; hprev = 0.0; 
							Chargeprev = 0.0; centroid[0]=0.0; centroidprev = 0.0;
							CharDens = fopen("Settings/ESPOutput/CurrentDensity_West.dat","a");
							fprintf(CharDens,"VARIABLES = \"CoordinateX\"\n\"CoordinateY\"\n\"CoordinateZ\"\n\"Current Density\"\n\n");
							for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							{
								i=0;
								index[0]=k,index[1]=j,index[2]=i;
								ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
								ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
								ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
								FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

								//Get the Diffusion of Charge Density at the boundary
								diffdia = Diffusion(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff);

								//Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
								diffdia = Diffusion(index,Rho,vols,ptrSize,slabRho,slabVolume,beta,Rdiff);

								INDEX = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
								
								Kcond = *(zone->Flowsol->MobIon)*(*((bcPtr->ZSlab)->RhoIon+((NJ-1)*k+j)));//
								DiffCoef = *(zone->Flowsol->ChargeDiffcoeff);
								p=*(surf+(INDEX)*6*3+3*1+0);
								q=*(surf+(INDEX)*6*3+3*1+1);
								r=*(surf+(INDEX)*6*3+3*1+2);
								Area=sqrt(p*p+q*q+r*r);
								
								Charge=Kcond*Fdiff[0]/Area;
								//TotalCurrent = TotalCurrent + Kcond*Fdiff[0] - DiffCoef*Rdiff[0];
								TotalCurrent = TotalCurrent + Kcond*Fdiff[0];
								if(bcPtr->Boundary_Tag == 24)
									fprintf(CharDens,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\n",centroid[0],centroid[1],centroid[2],Charge);

							}
							
							
							if(bcPtr->Boundary_Tag == 22)
							{
								WireCurrent = WireCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nWest Wire = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On West Wire = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							else
							{
								PlateCurrent = PlateCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nWest Plate = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On West Plate = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							fclose(CharDens);

						break;
						case 3://North Face
							Charge=0.0;TotalCurrent = 0.0;hlat = 0.0; hprev = 0.0; 
							Chargeprev = 0.0; centroid[0]=0.0; centroidprev = 0.0;
							CharDens = fopen("Settings/ESPOutput/CurrentDensity_North.dat","a");
							fprintf(CharDens,"VARIABLES = \"CoordinateX\"\n\"CoordinateY\"\n\"CoordinateZ\"\n\"Current Density\"\n\n");
							for(k=0;k<NK-1;k++)
							for(i=0;i<NI-1;i++)
							{
								j=NJ-2;
								index[0]=k,index[1]=NJ-2,index[2]=i;
								
								ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
								ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
								ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
								FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
								
								//Get the Diffusion of Charge Density at the boundary
								diffdia = Diffusion(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff);
								
								//Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
								diffdia = Diffusion(index,Rho,vols,ptrSize,slabRho,slabVolume,beta,Rdiff);
								
								INDEX = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
								
								Kcond = *(zone->Flowsol->MobIon)*(*((bcPtr->ZSlab)->RhoIon+((NI-1)*k+i)));//
								DiffCoef = *(zone->Flowsol->ChargeDiffcoeff);
								p=*(surf+(INDEX)*6*3+3*3+0);
								q=*(surf+(INDEX)*6*3+3*3+1);
								r=*(surf+(INDEX)*6*3+3*3+2);
								
								Area=sqrt(p*p+q*q+r*r);
								
								Charge=-Kcond*Fdiff[3]/Area;
								//TotalCurrent = TotalCurrent - Kcond*Fdiff[3] - DiffCoef*Rdiff[3];
								TotalCurrent = TotalCurrent - Kcond*Fdiff[3];
								if(bcPtr->Boundary_Tag == 24)
									fprintf(CharDens,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\n",centroid[0],centroid[1],centroid[2],Charge);
								

							}
							
							
							if(bcPtr->Boundary_Tag == 22)
							{
								WireCurrent = WireCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nNorth Wire = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On North Wire = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							else
							{
								PlateCurrent = PlateCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nNorth Plate = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On North Plate = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							fclose(CharDens);
						break;
						case 4://South Face
							Charge=0.0;TotalCurrent = 0.0;hlat = 0.0; hprev = 0.0; 
							Chargeprev = 0.0; centroid[0]=0.0; centroidprev = 0.0;
							CharDens = fopen("Settings/ESPOutput/CurrentDensity_South.dat","a");
							fprintf(CharDens,"VARIABLES = \"CoordinateX\"\n\"CoordinateY\"\n\"CoordinateZ\"\n\"Current Density\"\n\n");
							for(k=0;k<NK-1;k++)
							for(i=0;i<NI-1;i++)
							{
								j=0;
								index[0]=k,index[1]=j,index[2]=i;
								
								//ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
								//ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
								//ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
								//FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

								//Get the Diffusion of Potential at the boundary
								//diffdia = Diffusion(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff);
								  diffdia = Diffusion2d(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff,CoeffN);
								
								//Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
								//diffdia = Diffusion(index,Rho,vols,ptrSize,slabRho,slabVolume,beta,Rdiff);
								
								INDEX = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
								
								Kcond = *(zone->Flowsol->MobIon)*(*((bcPtr->ZSlab)->RhoIon+((NI-1)*k+i)));//
								//printf("\nIon Mobility = %e\tCharge = %e\n",*(zone->Flowsol->MobIon),(*((bcPtr->ZSlab)->RhoIon+((NI-1)*k+i))) );
								//DiffCoef = *(zone->Flowsol->ChargeDiffcoeff);
								p=*(surf+(INDEX)*6*3+3*2+0);
								q=*(surf+(INDEX)*6*3+3*2+1);
								r=*(surf+(INDEX)*6*3+3*2+2);
								Area=sqrt(p*p+q*q+r*r);
								
								Charge=Kcond*Fdiff[2]/Area;
								//TotalCurrent = TotalCurrent + Kcond*Fdiff[2] - DiffCoef*Rdiff[2];
								TotalCurrent = TotalCurrent + Kcond*Fdiff[2];
								if(bcPtr->Boundary_Tag == 24)
									fprintf(CharDens,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\n",centroid[0],centroid[1],centroid[2],Charge);
					

							}
							
							
							if(bcPtr->Boundary_Tag == 22)
							{
								WireCurrent = WireCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nSouth Wire = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On South Wire = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							else
							{
								PlateCurrent = PlateCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nSouth Plate = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On South Plate = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							fclose(CharDens);
						break;
						case 5://Top Face
							Charge=0.0;TotalCurrent = 0.0;hlat = 0.0; hprev = 0.0; 
							Chargeprev = 0.0; centroid[0]=0.0; centroidprev = 0.0;
							CharDens = fopen("Settings/ESPOutput/CurrentDensity_Top.dat","a");
							fprintf(CharDens,"VARIABLES = \"CoordinateX\"\n\"CoordinateY\"\n\"CoordinateZ\"\n\"Current Density\"\n\n");
							for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
							{
								k=NK-2;
								index[0]=k,index[1]=j,index[2]=i;
								ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
								ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
								ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
								FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

								//Get the Diffusion of Charge Density at the boundary
								diffdia = Diffusion(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff);

								//Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
								diffdia = Diffusion(index,Rho,vols,ptrSize,slabRho,slabVolume,beta,Rdiff);

								INDEX = (NI-1)*(NJ-1)*k+(NI-1)*j+i;

								Kcond = *(zone->Flowsol->MobIon)*(*((bcPtr->ZSlab)->RhoIon+((NI-1)*j+i)));//
								DiffCoef = *(zone->Flowsol->ChargeDiffcoeff);
								p=*(surf+(INDEX)*6*3+3*5+0);
								q=*(surf+(INDEX)*6*3+3*5+1);
								r=*(surf+(INDEX)*6*3+3*5+2);
								Area=sqrt(p*p+q*q+r*r);
								
								Charge=-Kcond*Fdiff[5]/Area;
								//TotalCurrent = TotalCurrent - Kcond*Fdiff[5] - DiffCoef*Rdiff[5];
								TotalCurrent = TotalCurrent - Kcond*Fdiff[5];
								if(bcPtr->Boundary_Tag == 24)
									fprintf(CharDens,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\n",centroid[0],centroid[1],centroid[2],Charge);
							}
							
							
							if(bcPtr->Boundary_Tag == 22)
							{
								WireCurrent = WireCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nTop Wire = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On Top Wire = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							else
							{
								PlateCurrent = PlateCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nTop Plate = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On Top Plate = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							fclose(CharDens);
						break;
						case 6:
							Charge=0.0;TotalCurrent = 0.0;hlat = 0.0; hprev = 0.0; 
							Chargeprev = 0.0; centroid[0]=0.0; centroidprev = 0.0;
							CharDens = fopen("Settings/ESPOutput/CurrentDensity_Bottom.dat","a");
							fprintf(CharDens,"VARIABLES = \"CoordinateX\"\n\"CoordinateY\"\n\"CoordinateZ\"\n\"Current Density\"\n\n");
							for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
							{
								k=0;
								index[0]=k,index[1]=j,index[2]=i;
								ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
								ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
								ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
								FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

								//Get the Diffusion of Charge Density at the boundary
								diffdia = Diffusion(index,Potential,vols,ptrSize,slabPotential,slabVolume,beta,Fdiff);

								//Kcond=CalcDiffusionCoef(base,zone,vartag,INDEX);
								diffdia = Diffusion(index,Rho,vols,ptrSize,slabRho,slabVolume,beta,Rdiff);

								INDEX = (NI-1)*(NJ-1)*k+(NI-1)*j+i;

								Kcond = *(zone->Flowsol->MobIon)*(*((bcPtr->ZSlab)->RhoIon+((NI-1)*j+i)));//
								DiffCoef = *(zone->Flowsol->ChargeDiffcoeff);
								p=*(surf+(INDEX)*6*3+3*4+0);
								q=*(surf+(INDEX)*6*3+3*4+1);
								r=*(surf+(INDEX)*6*3+3*4+2);
								Area=sqrt(p*p+q*q+r*r);
								
								Charge=Kcond*Fdiff[4]/Area;
								//TotalCurrent = TotalCurrent + Kcond*Fdiff[4] - DiffCoef*Rdiff[4];
								TotalCurrent = TotalCurrent + Kcond*Fdiff[4];
								if(bcPtr->Boundary_Tag == 24)
									fprintf(CharDens,"%12.8lf\t%12.8lf\t%12.8lf\t%12.8e\n",centroid[0],centroid[1],centroid[2],Charge);
							}
							
							
							if(bcPtr->Boundary_Tag == 22)
							{
								WireCurrent = WireCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nBottom Wire = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On Bottom Wire = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							else
							{
								PlateCurrent = PlateCurrent + fabs(TotalCurrent);
								fprintf(PlateCur,"\nBottom Plate = %12.8e", fabs(TotalCurrent));
								printf("\n\n**************************\nTotal Current On Bottom Plate = %12.8lf mA\n**************************\n", fabs(TotalCurrent*1.0e3));
							}
							fclose(CharDens);
						break;
					}//End of Switch
				
				}//End of Checking
			}
		}
	}
	printf("\nTotal Current on Plates = %12.8e Amp/m\nTotal Current on Wire = %12.8e Amp/m\n", PlateCurrent, WireCurrent);
	fprintf(PlateCur,"\n\n\n\n\t\t***************************\n\t\tTotal Current on All Plates = %12.8e\n\t\t***************************\n", fabs(PlateCurrent));
	fprintf(PlateCur,"\n\n\n\n\t\t***************************\n\t\tTotal Current on All Wires = %12.8e\n\t\t***************************\n", fabs(WireCurrent));
	fclose(PlateCur);
	return (0);
}
