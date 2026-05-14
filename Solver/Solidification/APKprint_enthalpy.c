#include "APKprint_enthalpy.h"
int print_enthalpy(Base *CgnsBase)
{
	int NI,NJ,NK,no,numzones,total,i;
	int *size;
	int numBc,index,direction,number,numiface;
	char str[]="result_1.dat";
	FILE *fp;
	Zone *zone;
	numzones=CgnsBase->Number_of_zones;
	fp=fopen(str,"w");
	if (fp==NULL)
	{
	printf("file is not opened sucessfully so exiting in print_enthalpy...\n");
	exit(0);
	}
	
	
	for (no=0;no<numzones;no++)
	{
		zone = (CgnsBase->CgnsZone+no);
		size=((CgnsBase)->CgnsZone+no)->Size;
		//printf("No segmentation fault till now\n");
		NI = *(size);
		NJ = *(size+1);
		NK = *(size+2);
		total = (NI-1)*(NJ-1)*(NK-1);
		fprintf(fp,"\nIn enthalpy function %d %d %d\n",NI,NJ,NK);
		fprintf(fp,"\nTotal cells=%d\n",total);
		for (i=0;i<total;i++)
		{
			if (i%(NI-1)==0)
			fprintf(fp,"\n");
			fprintf(fp,"%e\t",*((((CgnsBase->CgnsZone+no))->Flowsol)->umass+i));
		}
		numBc=zone->Number_of_BC_Patches;
		fprintf(fp,"\nNumber of boundaries=%d\n",numBc);
		for (index=0;index<numBc;index++)
		{
			fprintf(fp,"\nboundary number=%d\n",index);
			direction= (zone->ZoneBc+index)->Direction;
			if (direction==1 || direction==2)
			number = (NK-1)*(NJ-1);
			else if (direction==3 || direction==4)
			number = (NK-1)*(NI-1);
			else if (direction==5 || direction==6)
			number = (NI-1)*(NJ-1);
			fprintf(fp,"number=%d\n",number);
			for (i=0;i<number;i++)
			fprintf(fp,"%e\t",*(((zone->ZoneBc+index)->ZSlab)->umass+i));
			//fprintf(fp,"%lf\n",g);
		}
		numiface=zone->Number_of_Interfaces;
		fprintf(fp,"\nNumber of Interfaces=%d\n",numiface);
		for (index=0;index<numiface;index++)
		{
			fprintf(fp,"\ninterface number=%d\n",index);
			direction= (zone->Iface+index)->Direction;
			if (direction==1 || direction==2)
			number = (NK-1)*(NJ-1);
			else if (direction==3 || direction==4)
			number = (NK-1)*(NI-1);
			else if (direction==5 || direction==6)
			number = (NI-1)*(NJ-1);
			fprintf(fp,"\nnumber=%d\n",number);
			for (i=0;i<number;i++)
			{
				if (i%(NI-1)==0)
				fprintf(fp,"\n");
				fprintf(fp,"%e\t",*(((zone->Iface+index)->ZSlab)->umass+i));
			}
		}
	}
	for (no=0;no<numzones;no++)
	{
		zone = (CgnsBase->CgnsZone+no);
		size=((CgnsBase)->CgnsZone+no)->Size;
		//printf("No segmentation fault till now\n");
		NI = *(size);
		NJ = *(size+1);
		NK = *(size+2);
		total = (NI-1)*(NJ-1)*(NK-1);
		fprintf(fp,"\nIn enthalpy function %d %d %d\n",NI,NJ,NK);
		fprintf(fp,"\nTotal cells=%d\n",total);
		for (i=0;i<total;i++)
		{
			if (i%(NI-1)==0)
			fprintf(fp,"\n");
			fprintf(fp,"%e\t",*((((CgnsBase->CgnsZone+no))->Flowsol)->vmass+i));
		}
		numBc=zone->Number_of_BC_Patches;
		fprintf(fp,"\nNumber of boundaries=%d\n",numBc);
		for (index=0;index<numBc;index++)
		{
			fprintf(fp,"\nboundary number=%d\n",index);
			direction= (zone->ZoneBc+index)->Direction;
			if (direction==1 || direction==2)
			number = (NK-1)*(NJ-1);
			else if (direction==3 || direction==4)
			number = (NK-1)*(NI-1);
			else if (direction==5 || direction==6)
			number = (NI-1)*(NJ-1);
			fprintf(fp,"number=%d\n",number);
			for (i=0;i<number;i++)
			fprintf(fp,"%e\t",*(((zone->ZoneBc+index)->ZSlab)->vmass+i));
			//fprintf(fp,"%lf\n",g);
		}
		numiface=zone->Number_of_Interfaces;
		fprintf(fp,"\nNumber of Interfaces=%d\n",numiface);
		for (index=0;index<numiface;index++)
		{
			fprintf(fp,"\ninterface number=%d\n",index);
			direction= (zone->Iface+index)->Direction;
			if (direction==1 || direction==2)
			number = (NK-1)*(NJ-1);
			else if (direction==3 || direction==4)
			number = (NK-1)*(NI-1);
			else if (direction==5 || direction==6)
			number = (NI-1)*(NJ-1);
			fprintf(fp,"\nnumber=%d\n",number);
			for (i=0;i<number;i++)
			{
				if (i%(NI-1)==0)
				fprintf(fp,"\n");
				fprintf(fp,"%e\t",*(((zone->Iface+index)->ZSlab)->vmass+i));
			}
		}
	}
	
	
	for (no=0;no<numzones;no++)
	{
		zone = (CgnsBase->CgnsZone+no);
		size=((CgnsBase)->CgnsZone+no)->Size;
		//printf("No segmentation fault till now\n");
		NI = *(size);
		NJ = *(size+1);
		NK = *(size+2);
		total = (NI-1)*(NJ-1)*(NK-1);
		fprintf(fp,"\nIn enthalpy function %d %d %d\n",NI,NJ,NK);
		fprintf(fp,"\nTotal cells=%d\n",total);
		for (i=0;i<total;i++)
		{
			if (i%(NI-1)==0)
			fprintf(fp,"\n");
			fprintf(fp,"%e\t",*((((CgnsBase->CgnsZone+no))->Flowsol)->P+i));
		}
		numBc=zone->Number_of_BC_Patches;
		fprintf(fp,"\nNumber of boundaries=%d\n",numBc);
		for (index=0;index<numBc;index++)
		{
			fprintf(fp,"\nboundary number=%d\n",index);
			direction= (zone->ZoneBc+index)->Direction;
			if (direction==1 || direction==2)
			number = (NK-1)*(NJ-1);
			else if (direction==3 || direction==4)
			number = (NK-1)*(NI-1);
			else if (direction==5 || direction==6)
			number = (NI-1)*(NJ-1);
			fprintf(fp,"number=%d\n",number);
			for (i=0;i<number;i++)
			fprintf(fp,"%e\t",*(((zone->ZoneBc+index)->ZSlab)->P+i));
			//fprintf(fp,"%lf\n",g);
		}
		numiface=zone->Number_of_Interfaces;
		fprintf(fp,"\nNumber of Interfaces=%d\n",numiface);
		for (index=0;index<numiface;index++)
		{
			fprintf(fp,"\ninterface number=%d\n",index);
			direction= (zone->Iface+index)->Direction;
			if (direction==1 || direction==2)
			number = (NK-1)*(NJ-1);
			else if (direction==3 || direction==4)
			number = (NK-1)*(NI-1);
			else if (direction==5 || direction==6)
			number = (NI-1)*(NJ-1);
			fprintf(fp,"\nnumber=%d\n",number);
			for (i=0;i<number;i++)
			{
				if (i%(NI-1)==0)
				fprintf(fp,"\n");
				fprintf(fp,"%e\t",*(((zone->Iface+index)->ZSlab)->P+i));
			}
		}
	}
	printf("Results are printed in %s\n",str);
	fclose(fp);
	
	
	return 0;
}
