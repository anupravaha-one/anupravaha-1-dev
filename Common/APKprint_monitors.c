#include "APKprint_monitors.h"

int print_monitors(Base* base)
{
	int index,number,solvable_variables,*ptrSize,loop,inner,NI,NJ,NK;
	double *variable;
	Zone *zone;
	solvable_variables = base->total_var;
	char *filename,*temp;
	FILE *fp;
	filename = malloc(sizeof(char)*50);
	temp     = malloc(sizeof(char)*50);
	number = base->Monitors->number;
	
	for (loop=0;loop<number;loop++)
	{
		sprintf(temp,"%d.dat",loop+1);
		sprintf(filename,"SOLUTIONS/monitors/monitor");
		filename = strcat(filename,temp);
		//printf("file name is %s\n",filename);
		fp = fopen(filename,"a");
		zone = base->CgnsZone+(*(base->Monitors->zone_num+loop));
		ptrSize= zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		index = *(base->Monitors->i+loop)+(NI-1)*(*(base->Monitors->j+loop))+(NI-1)*(NJ-1)*(*(base->Monitors->k+loop));
		fprintf(fp,"%e\t",base->BaseParameters->Time);
		for(inner=0;inner<solvable_variables;inner++)
		{
			variable = MakeDomainArray(zone,base->VARTAG[inner]);
			fprintf(fp,"%e\t",*(variable+index));
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
	free(temp);
	free(filename);

	return 0;
}

