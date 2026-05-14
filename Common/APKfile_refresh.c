#include "APKfile_refresh.h"

int file_refresh(Root *CgnsRoot,char *inputfilename,char *outputfilename)
{
	int cntVariable;
	FILE *nscon;
	nscon=fopen("Settings/iwrite.dat","w");
	fprintf(nscon,"1");
	fclose(nscon);
	nscon=fopen("Settings/CountFileWrite.dat","w");
	fprintf(nscon,"5000");
	fclose(nscon);
	strcpy(CgnsRoot->inputgridname,inputfilename);
	strcpy(CgnsRoot->outputsolutionname,outputfilename);
	printf("\ninput filename is %s\n",CgnsRoot->inputgridname);
	printf("\noutput filename is %s\n",CgnsRoot->outputsolutionname);

 	nscon = fopen("SOLUTIONS/errEvolution.dat","w");
	if (nscon == NULL) {
	   system("mkdir SOLUTIONS");
 	   nscon = fopen("SOLUTIONS/errEvolution.dat","w");
	   if (nscon == NULL) { printf("Error opening File SOLUTIONS/errEvolution.dat Exiting.....\n");
	      exit(1);
           }
	}
// 	fprintf(nscon,"variables = \n\"#Time\"\t\t\"U rms\"\t\"V rms\"\t\"W rms\"\t\"nu_t rms\"\n");
//	fprintf(nscon,"variables = \n\"#Time\"\t\t\"U rms\"\t\"V rms\"\t\"W rms\"\t\"k rms\"\t\"Epsilon rms\"\n");
//	This fprintf is called now in Main
  	fclose(nscon);
	
	nscon=fopen("Settings/TurbRelaxation.dat","w");
	fprintf(nscon,"0.7");
	fclose(nscon);
	
	nscon=fopen("Settings/OmegaRelaxation.dat","w");
	fprintf(nscon,"0.7");
	fclose(nscon);

	nscon = fopen("Settings/PressureRelaxation.dat","w");
	fprintf(nscon,"1.0");
	fclose(nscon);
	
	nscon = fopen("Settings/VelocityRelaxation.dat","w");
	fprintf(nscon,"0.6");
	fclose(nscon);

	nscon=fopen("Settings/ChangeTimeStep.dat","w");
	fprintf(nscon,"0\n0.001");
	fprintf(nscon,"\n\n\n1. Change '0' in the first line to '1' for the solver to read the time step from the second line of this file\n2. Second line is the time step");
	fclose(nscon);
	
// 	nscon=fopen("Settings/TimeStep.dat","w");
// 	fprintf(nscon,"1e-3");
// 	fclose(nscon);

	nscon=fopen("Settings/TurbulenceStartTime.dat","w");
	fprintf(nscon,"0.0");
	fclose(nscon);

	nscon=fopen("Settings/DummyViscCofactor.dat","w");
	fprintf(nscon,"0.0");
	fclose(nscon);
	nscon=fopen("Settings/DummyViscocity.dat","w");
	fprintf(nscon,"0.0");
	fclose(nscon);


	nscon=fopen("Settings/printGamma.dat","w");
	fprintf(nscon,"0");
	fclose(nscon);
	
	nscon=fopen("Settings/Switch_Iter.dat","w");
	fprintf(nscon,"350");
	fclose(nscon);

	
	//system("rm SOLUTIONS/monitors/monitor*");
	//nscon=fopen("printEachTermKEeqns.dat","w");
	//fprintf(nscon,"0");
	//fclose(nscon);
	
	//fclose(fopen("K_Eps.dat","w"));
	
	//fclose(fopen("EachTermData.dat","w"));
	//fclose(fopen("KDetails.dat","w"));
	//fclose(fopen("EpsDetails.dat","w"));
	
	return 0;
}
