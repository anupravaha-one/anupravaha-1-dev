#include "APKSolid_fraction.h"
//this function is written to calculate the massaveraged values
double Solid_fraction(Base *base)
{
        
	int numzones,cntZ=0,*ptrSize,NI,NJ,NK;
        int i,j,k,indexP,totalcells=0;
        Zone *zone;
        double *f,*u,*rho,*e,sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0,sum5=0.0,vol,total_volume=0.0,*newMassFlux;
	double test=0.0;
	char str[]="volumeavg.dat";
	FILE *fp;
        numzones= base->Number_of_zones;
	fp=fopen(str,"w");
	printf("file is  opened sucessfully and writing in print_enthalpy...\n");
        
	for(cntZ=0;cntZ<numzones;cntZ++)
        {
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
                NJ=*(ptrSize+1);
                NK=*(ptrSize+2);
		 
		newMassFlux = (double *)malloc(6*(NI-1)*(NJ-1)*(NK-1)*sizeof(double));
		CalcMassFlux(newMassFlux,zone);
		               
		f = zone->Flowsol->Rho;
		u = zone->Flowsol->u;
		e = zone->Flowsol->enthalpy;
		rho = zone->Flowsol->f;
		
		for(k=0;k<(NK-1);k++)
		for(i=0;i<(NI-1);i++)
		{
			total_volume=0.0;
			test=0.0;
			sum1=0.0;
			sum2=0.0;
			sum3=0.0;
			sum4=0.0;
			sum5=0.0;
			totalcells=0;
			for(j=0;j<(NJ-1);j++)
			{
				totalcells = totalcells+(NJ-1)*(NK-1);
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				vol = *(zone->ZGeom->Volume+indexP);
				test=test-(*(newMassFlux+indexP*6+0));
				sum1=sum1-(*(newMassFlux+indexP*6+0))*(*(u+indexP));
				sum2=sum2-(*(newMassFlux+indexP*6+0))*(*(rho+indexP));
				sum3=sum3-(*(newMassFlux+indexP*6+0))*(*(e+indexP));
				sum4=sum4-(*(newMassFlux+indexP*6+0))*(*(f+indexP));
				sum5=sum5-(*(newMassFlux+indexP*6+0))*((*(rho+indexP))*(*(f+indexP)/0.589));
				//total_volume = total_volume + vol;
			}
			//sum1 =sum1/test;
			sum2 =sum2/test;
			sum3 =sum3/test;
			sum4 =sum4/test;
			sum5 =sum5/test;
			sum1=sum1/test;
			if(k==0)
			fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,test,sum1,sum2,sum3,sum4,sum5);
			//exit(0);
			//fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\n",i,sum,sum1,sum2,sum3);
		}
		free(newMassFlux);
        }
	fclose(fp);
        
        return 0;
}





