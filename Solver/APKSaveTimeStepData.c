#include "APKSaveTimeStepData.h"

/**
This function write the data to .cgns file depending upon the integer in CountFileWrite.dat.
If the integer is x, the data will be printed into the .cgns file after every x timesteps.
*/

int SaveTimeStepData(Root CgnsRoot,int cntTimeStep,double time)
{
	char *fileName;
	
	int fileTimeInter=5000;
	
#ifdef FILE_OPEN
{
	FILE *countWrite;

	countWrite=fopen("Settings/CountFileWrite.dat","r");
	
	fscanf(countWrite,"%d",&fileTimeInter);
	
	fclose(countWrite);
}
#endif
	
	if((cntTimeStep%fileTimeInter)==0)
	{
		fileName=FindFileNamePtr(CgnsRoot,time);
		printf("%s\n",fileName);
		
		if(CgnsRoot.CgnsBase->BaseParameters->tagVorticity==1)
			CalcVorticity(CgnsRoot.CgnsBase);
		if(CgnsRoot.CgnsBase->BaseParameters->TagEfield==1)
			ElectricFieldXYZ(CgnsRoot.CgnsBase);
	                 
	                  
                //WriteCgnsFile(fileName,CgnsRoot);
		
		printf("crossed the tag --2\n");
		printslab(CgnsRoot.CgnsBase);
		printsolver(CgnsRoot.CgnsBase);
	}
	return 0;
}


char *FindFileNamePtr(Root CgnsRoot,double time)
{
	char *str1,*str2,*strTemp;
	char *check,fileNameProblem[40],datafile[40];
	int i;
	
	strTemp=(char *)malloc(40*sizeof(char));
	str1=(char *)malloc(40*sizeof(char));
	str2=(char *)malloc(40*sizeof(char));

	strcpy(fileNameProblem,CgnsRoot.outputsolutionname);
	
	check= &fileNameProblem[0];
	i=0;
	while(*check!='.')
	{
		datafile[i]=*check;
		check++;
		i++;
	}
	*check=0;
	datafile[i]=0;
	strcat(datafile,".dat");
	sprintf(strTemp,"-%6.4lf.cgns",time);
	//str1=strcat(str1,str2);
	str1=strcat(fileNameProblem,strTemp);
	WriteCgnsFile(str1,CgnsRoot);
	Write_Solution(CgnsRoot.CgnsBase,datafile);
	return (str1);
}
