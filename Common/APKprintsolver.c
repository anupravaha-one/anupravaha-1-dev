#include "APKprintsolver.h"
/**
This function,printsolver() prints the centroid values of the required variable into a data file, which goes into 'SOLUTIONS' 
folder. To this any other variable you wish to print into a data file can be added. 
*/

int printsolver(Base *CgnsBase)
{
	FILE *solu,*solv,*solw,*solT,*solP,*solS1,*solS2,*solS3,*vert,*horz,*Potential,*disperse;    
	solu=fopen("SOLUTIONS/solveru.dat","w");
	solv=fopen("SOLUTIONS/solverv.dat","w");
	solw=fopen("SOLUTIONS/solverw.dat","w");

	solT=fopen("SOLUTIONS/solverT.dat","w");

	solP=fopen("SOLUTIONS/solverP.dat","w");
	solS1=fopen("SOLUTIONS/solverS1.dat","w");
	solS2=fopen("SOLUTIONS/solverS2.dat","w");
	solS3=fopen("SOLUTIONS/solverS3.dat","w");
	vert=fopen("SOLUTIONS/verticalmidplane.dat","w");
	horz=fopen("SOLUTIONS/horizontalmidplane.dat","w");
	
	Potential=fopen("SOLUTIONS/ElectricPotential.dat","w");

	disperse=fopen("SOLUTIONS/disperse.dat","w");

	int NI,NJ,NK,nZones,numZones;
	int i,j,k,*ptrSize,num;
	double *coordX,*coordY,*coordZ,time;//,yPlus,utau;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3];

	int FindCellCentre();
	time = CgnsBase->BaseParameters->Time;
	num = time;	
	nZones=((CgnsBase)->Number_of_zones);
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		fprintf(solu,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(solv,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(solw,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(solT,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		
		fprintf(solP,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(solS1,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(solS2,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(solS3,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		fprintf(Potential,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		
		fprintf(disperse,"~~~~~~~ZONENAME :%s ~~~~~~~~Time =  %lf\n",((CgnsBase)->CgnsZone+numZones)->zonename,(CgnsBase->BaseParameters->Time));
 		fprintf(disperse,"~~~~~~~ZONENAME :%s ~~~~~~~~Time =  %lf\n",((CgnsBase)->CgnsZone+numZones)->zonename,time);

		ptrSize=((CgnsBase)->CgnsZone+numZones)->Size;
      
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		coordX=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->X;
		coordY=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->Y; 
		coordZ=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->Z; 
		
		/*fprintf(solS3,"\nX\tY\tZ\tMu_t\tMu_eff\tYplus  values at x=0.25\n");
		fprintf(vert,"\nX\tY\tZ\tMu_t\tMu_eff\tYplus  values at x=0.5\n");
		fprintf(horz,"\nX\tY\tZ\tMu_t\tMu_eff\tYplus  values at x=0.75\n");*/
//		if (numZones==0) 
		{
			if(CgnsBase->BaseParameters->flag2d==0) {
				fprintf(vert,"\n#Vertical Mid-plane Data at Zone No: %d \t Variables: X \t Y \t Z \t",numZones);
			} else {
				fprintf(vert,"\n#Vertical Mid-plane Data at Zone No: %d \t Variables: X \t Y \t",numZones);
			}
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL) fprintf(vert,"U\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->v!=NULL) fprintf(vert,"V\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->w!=NULL) fprintf(vert,"W\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->P!=NULL) fprintf(vert,"P\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->T!=NULL) fprintf(vert,"T\t");
			fprintf(vert,"\n");

			if(CgnsBase->BaseParameters->flag2d==0) {
				fprintf(horz,"\n#Horizontal Mid-plane Data at Zone No: %d \t Variables: X \t Y \t Z \t",numZones);
			} else {
				fprintf(horz,"\n#Horizontal Mid-plane Data at Zone No: %d \t Variables: X \t Y \t",numZones);
			}
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL) fprintf(horz,"U\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->v!=NULL) fprintf(horz,"V\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->w!=NULL) fprintf(horz,"W\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->P!=NULL) fprintf(horz,"P\t");
		        if((CgnsBase->CgnsZone+numZones)->Flowsol->T!=NULL) fprintf(horz,"T\t");
			fprintf(horz,"\n");
		}
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
			ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
			ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
			
			 
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
		
		
			
				if(((CgnsBase->CgnsZone+numZones)->Flowsol)->u!=NULL)
				{
					fprintf(solu,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
				}
				
				if(((CgnsBase->CgnsZone+numZones)->Flowsol)->v!=NULL)
				{
					fprintf(solv,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
				}
				
				if(((CgnsBase->CgnsZone+numZones)->Flowsol)->w!=NULL)
				{
					fprintf(solw,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->w+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
				}
				
				if(((CgnsBase->CgnsZone+numZones)->Flowsol)->P!=NULL)
				{
					fprintf(solP,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->P+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
				}
				if(((CgnsBase->CgnsZone+numZones)->Flowsol)->T!=NULL)
				{
					fprintf(solT,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->T+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
					//fprintf(solT,"%e\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->T+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
				}
				if((CgnsBase->CgnsZone+numZones)->Flowsol->Scal[0]!=NULL)
				{
					fprintf(solS1,"\n%12.8lf\t%12.8lf\t%12.8lf\t%.20e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Scal[0]+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
				}
				if((CgnsBase->CgnsZone+numZones)->Flowsol->Scal[1]!=NULL)
				{
					fprintf(solS2,"\n%12.8lf\t%12.8lf\t%12.8lf\t%.20e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Scal[1]+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
				}
				if((CgnsBase->CgnsZone+numZones)->Flowsol->Scal[2]!=NULL)
				{
					fprintf(solS3,"\n%12.8lf\t%12.8lf\t%12.8lf\t%.20e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Scal[2]+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
				}
//				if(numZones==0 && i==1*((NI-2)/2) && (CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL)
//				if(i==1*((NI-2)/2) && (CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL)
				if(i==(NI-2)/2 ) 
				{
					if(CgnsBase->BaseParameters->flag2d==0) {
						fprintf(vert,"%12.8lf\t%12.8lf\t%12.8lf\t",centroid[0],centroid[1],centroid[2]);
					} else {
						fprintf(vert,"%12.8lf\t%12.8lf\t",centroid[0],centroid[1]);
					}

				        if((CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL)
					fprintf(vert,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

				        if((CgnsBase->CgnsZone+numZones)->Flowsol->v!=NULL)
					fprintf(vert,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

				        if((CgnsBase->CgnsZone+numZones)->Flowsol->w!=NULL)
					fprintf(vert,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->w+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

				        if((CgnsBase->CgnsZone+numZones)->Flowsol->P!=NULL)
					fprintf(vert,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->P+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

				        if((CgnsBase->CgnsZone+numZones)->Flowsol->T!=NULL)
					fprintf(vert,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->T+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

					fprintf(vert,"\n");

				}
				
// 				if(numZones==0 && j==1*((NI-2)/2) && (CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL)
// 				if(j==1*((NJ-2)/2) && (CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL)
 				if(j==(NJ-2)/2) 
				{
					if(CgnsBase->BaseParameters->flag2d==0) {
						fprintf(horz,"%12.8lf\t%12.8lf\t%12.8lf\t",centroid[0],centroid[1],centroid[2]);
					} else {
						fprintf(horz,"%12.8lf\t%12.8lf\t",centroid[0],centroid[1]);
					}
 					if( (CgnsBase->CgnsZone+numZones)->Flowsol->u!=NULL)
					fprintf(horz,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

 					if( (CgnsBase->CgnsZone+numZones)->Flowsol->v!=NULL)
					fprintf(horz,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

 					if( (CgnsBase->CgnsZone+numZones)->Flowsol->w!=NULL)
					fprintf(horz,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->w+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

 					if( (CgnsBase->CgnsZone+numZones)->Flowsol->P!=NULL)
					fprintf(horz,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->P+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

 					if( (CgnsBase->CgnsZone+numZones)->Flowsol->T!=NULL)
					fprintf(horz,"%lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->T+(NI-1)*(NJ-1)*k+(NI-1)*j+i));

					fprintf(horz,"\n");

				}

			if((((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex!=NULL)&&((CgnsBase->flagDroplet == 1)))
			{

			fprintf(disperse,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->U_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersey+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->V_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->VarthetaDisperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->NumDensity+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->PrevFlowsol)->ParticleDia+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*(CgnsBase->BaseParameters->flagVariableDia)),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Rho+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*((CgnsBase->CgnsZone+numZones)->Flowsol->flagRho)),*((CgnsBase->CgnsZone+numZones)->Flowsol->Tdrop+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->T+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->DropFuelMassFrac+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
			}	
	
			if((((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex!=NULL)&&((CgnsBase->flagParticle == 1)))
			{
			fprintf(disperse,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->U_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersey+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->V_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->VarthetaDisperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
			}
			
			
			
			
			
			
			
			
			
			
			

                                 			
 	if((((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex!=NULL)&&((CgnsBase->flagDroplet == 1)))
	
	{

	fprintf(disperse,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->U_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersey+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->V_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->VarthetaDisperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->NumDensity+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->PrevFlowsol)->ParticleDia+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*(CgnsBase->BaseParameters->flagVariableDia)),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Rho+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*((CgnsBase->CgnsZone+numZones)->Flowsol->flagRho)));	
	}	
	
	if((((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex!=NULL)&&((CgnsBase->flagParticle == 1)))
	{
	fprintf(disperse,"\n%12.8lf\t%12.8lf\t%12.8lf\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersex+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->U_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->u+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->MomentumDispersey+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*((CgnsBase->CgnsZone+numZones)->Flowsol->V_Disperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->v+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->VarthetaDisperse+(NI-1)*(NJ-1)*k+(NI-1)*j+i));	
	}	
	
			

			//fprintf(Potential,"%12.8lf\t",*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Potential+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
			/*if((((CgnsBase->CgnsZone+numZones)->Flowsol)->Potential!=NULL)&&(i==1*((NI-2)/2)))
			{
				fprintf(Potential,"%12.8lf\t%12.8lf\t%12.8lf\t%lf\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Potential+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
			}*/
			
			//if((CgnsBase->CgnsZone+numZones)->Flowsol->Scal[2]!=NULL)
			/*if(numZones==0 && i==3*((NI-2)/4))
			{
				utau = *(((CgnsBase->CgnsZone+numZones)->ZoneBc)->uTau+i);
				yPlus = utau*centroid[1]/(*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Viscoeff));
				fprintf(solS3,"\n%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Tcoeff->TViscoeff+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Tcoeff->EffViscoeff+(NI-1)*(NJ-1)*k+(NI-1)*j+i),yPlus);
			}
			if(numZones==0 && i==1*((NI-2)/2))
			{
				utau = *(((CgnsBase->CgnsZone+numZones)->ZoneBc)->uTau+i);
				yPlus = utau*centroid[1]/(*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Viscoeff));
				fprintf(vert,"\n%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Tcoeff->TViscoeff+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Tcoeff->EffViscoeff+(NI-1)*(NJ-1)*k+(NI-1)*j+i),yPlus);
			}
			
			if(numZones==0 && i==1*((NI-2)/4))
			{
				utau = *(((CgnsBase->CgnsZone+numZones)->ZoneBc)->uTau+i);
				yPlus = utau*centroid[1]/(*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Viscoeff));
				fprintf(horz,"\n%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",centroid[0],centroid[1],centroid[2],*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Tcoeff->TViscoeff+(NI-1)*(NJ-1)*k+(NI-1)*j+i),*(((CgnsBase->CgnsZone+numZones)->Flowsol)->Tcoeff->EffViscoeff+(NI-1)*(NJ-1)*k+(NI-1)*j+i),yPlus);
			}*/
		}/*End of NI,NJ,NK loop*/
		
		
	
		
		
		
			
	
		
// 		Zone *zone;
// 		int direction,number,numiface,index,numBc;
// 		zone =(CgnsBase->CgnsZone+numZones); 
// 		numBc=zone->Number_of_BC_Patches;
// 		fprintf(solu,"\nNumber of boundaries=%d\n",numBc);
// 		for (index=0;index<numBc;index++)
// 		{
// 			fprintf(solu,"\nboundary number=%d\n",index);
// 			direction= (zone->ZoneBc+index)->Direction;
// 			if (direction==1 || direction==2)
// 			number = (NK-1)*(NJ-1);
// 			else if (direction==3 || direction==4)
// 			number = (NK-1)*(NI-1);
// 			else if (direction==5 || direction==6)
// 			number = (NI-1)*(NJ-1);
// 			fprintf(solu,"number=%d\n",number);
// 			for (i=0;i<number;i++)
// 			fprintf(solu,"%e\t",*(((zone->ZoneBc+index)->ZSlab)->u+i));
// 			//fprintf(fp,"%lf\n",g);
// 		}
// 		numiface=zone->Number_of_Interfaces;
// 		fprintf(solu,"\nNumber of Interfaces=%d\n",numiface);
// 		for (index=0;index<numiface;index++)
// 		{
// 			fprintf(solu,"\ninterface number=%d\n",index);
// 			direction= (zone->Iface+index)->Direction;
// 			if (direction==1 || direction==2)
// 			number = (NK-1)*(NJ-1);
// 			else if (direction==3 || direction==4)
// 			number = (NK-1)*(NI-1);
// 			else if (direction==5 || direction==6)
// 			number = (NI-1)*(NJ-1);
// 			fprintf(solu,"\nnumber=%d\n",number);
// 			for (i=0;i<number;i++)
// 			{
// 				if (i%(NI-1)==0)
// 				fprintf(solu,"\n");
// 				fprintf(solu,"%e\t",*(((zone->Iface+index)->ZSlab)->u+i));
// 			}
// 		}
		
		
	fprintf(Potential,"\n");
	
	}/*End of looping over Zones*/
	fclose(solu);
	fclose(solv);
	fclose(solw);
	fclose(solT);
	fclose(solP);
	fclose(solS1);
	fclose(solS2);
	fclose(solS3);
	fclose(vert);
	fclose(horz);
	fclose(Potential);
	fclose(disperse);
	
	return 0;
}

