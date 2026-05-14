#include "APKWriteWireCharge.h"

/**
This function is used for writing the wire charge density after certain
pre-determined number of timesteps as entered in WireChargeUpdateData.dat and
also writes the evolution of mean/average charge density on the wire.
*/

int WriteWireCharge(int direction,Zone *z,int size,double *dummy,int stp)
{
	int ni,nj,nk,i,j,k;
	FILE *MeanCharge,*WireCharge;
	double MnChar = 0.0;
	
		ni=*(z->Size+0);
		nj=*(z->Size+1);
		nk=*(z->Size+2); 

		
		switch(direction)
		{
			case 1://East
				WireCharge=fopen("Settings/ESPOutput/WireChargeDensity_East.dat","w");
				//fprintf(WireCharge, "\n*************   Time Step = %d   *************\n", stp);
				MeanCharge = fopen("Settings/ESPOutput/AverageWireCharge_East.dat","a");
				for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{

						fprintf(WireCharge, "%12.8e\n", *(dummy+(nj-1)*k+j));
						MnChar = MnChar + pow(*(dummy+(nj-1)*k+j),2);

					}
				MnChar = sqrt(MnChar/(nk-1)*(nj-1));
				fprintf(MeanCharge, "%12.8e\t%d\n", MnChar, stp);
				fclose(WireCharge);fclose(MeanCharge);
			break;
			
			case 2://West
				WireCharge=fopen("Settings/ESPOutput/WireChargeDensity_West.dat","w");
				//fprintf(WireCharge, "\n*************   Time Step = %d   *************\n", stp);
				MeanCharge = fopen("Settings/ESPOutput/AverageWireCharge_West.dat","a");
				for (k=0;k<nk-1;k++) 
					for (j=0;j<nj-1;j++)
					{

						fprintf(WireCharge, "%12.8e\n", *(dummy+(nj-1)*k+j));
						MnChar = MnChar + pow(*(dummy+(nj-1)*k+j),2);

					}
				MnChar = sqrt(MnChar/(nk-1)*(nj-1));
				fprintf(MeanCharge, "%12.8e\t%d\n", MnChar, stp);
				fclose(WireCharge);fclose(MeanCharge);
			break;

			case 3://North
				WireCharge=fopen("Settings/ESPOutput/WireChargeDensity_North.dat","w");
				//fprintf(WireCharge, "\n*************   Time Step = %d   *************\n", stp);
				MeanCharge = fopen("Settings/ESPOutput/AverageWireCharge_North.dat","a");
				for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{

						fprintf(WireCharge, "%12.8e\n", *(dummy+k*(ni-1)+i));
						MnChar = MnChar + pow(*(dummy+k*(ni-1)+i),2);

					}
				MnChar = sqrt(MnChar/(ni-1)*(nk-1));
				fprintf(MeanCharge, "%12.8e\t%d\n", MnChar, stp);
				fclose(WireCharge);fclose(MeanCharge);
			break;
			
			case 4://South
				WireCharge=fopen("Settings/ESPOutput/WireChargeDensity_South.dat","w");
				//fprintf(WireCharge, "\n*************   Time Step = %d   *************\n", stp);
				MeanCharge = fopen("Settings/ESPOutput/AverageWireCharge_South.dat","a");
				for (k=0;k<nk-1;k++) 
					for (i=0;i<ni-1;i++)
					{

						fprintf(WireCharge, "%12.8e\n", *(dummy+k*(ni-1)+i));
						MnChar = MnChar + pow(*(dummy+k*(ni-1)+i),2);

					}
				MnChar = sqrt(MnChar/(ni-1)*(nk-1));
				fprintf(MeanCharge, "%12.8e\t%d\n", MnChar, stp);
				fclose(WireCharge);fclose(MeanCharge);
			break;
			
			case 5:
				
				WireCharge=fopen("Settings/ESPOutput/WireChargeDensity_Top.dat","w");
				//fprintf(WireCharge, "\n*************   Time Step = %d   *************\n", stp);
				MeanCharge = fopen("Settings/ESPOutput/AverageWireCharge_Top.dat","a");
				for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{

						fprintf(WireCharge, "%12.8e\n", *(dummy+(ni-1)*j+i));
						MnChar = MnChar + pow(*(dummy+(ni-1)*j+i),2);

					}
				MnChar = sqrt(MnChar/(ni-1)*(nj-1));
				fprintf(MeanCharge, "%12.8e\t%d\n", MnChar, stp);
				fclose(WireCharge);fclose(MeanCharge);
			break;
			
			case 6:
				WireCharge=fopen("Settings/ESPOutput/WireChargeDensity_Bottom.dat","w");
				//fprintf(WireCharge, "\n*************   Time Step = %d   *************\n", stp);
				MeanCharge = fopen("Settings/ESPOutput/AverageWireCharge_Bottom.dat","a");
				for(j=0;j<nj-1;j++) 
					for(i=0;i<ni-1;i++)
					{

						fprintf(WireCharge, "%12.8e\n", *(dummy+(ni-1)*j+i));
						MnChar = MnChar + pow(*(dummy+(ni-1)*j+i),2);

					}
				MnChar = sqrt(MnChar/(ni-1)*(nj-1));
				fprintf(MeanCharge, "%12.8e\t%d\n", MnChar, stp);
				fclose(WireCharge);fclose(MeanCharge);
			break;


		}

	return (0);
}//End of main

