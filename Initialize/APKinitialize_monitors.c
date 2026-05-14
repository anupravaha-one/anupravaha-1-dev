#include"APKinitialize_monitors.h"

int initialize_monitors(Base *base)
{
	int num,loop;
	int NI,NJ,NK;
	char *filename,*temp;
	FILE *fp,*var2,*var3;
	int inner,solvable_variables;
        
        var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
	solvable_variables = base->total_var;
	
	printf("\n\n\nMonitors are used for printing the solutions at a discrete points in the domain");
	printf("\nThis requires zone number and i,j,k indices of that cell as input to the solver");
	printf("\nThe solutions will be writen to bin/SOLUTIONS/monitors/ folder in tecplot data file format : ");
	printf("\nEnter the number of points you want to monitor : ");
	scanf("%d",&num);        
        fprintf(var2,"%d\n",num);        
        fprintf(var3,"%d \tEnter the number of points you want to monitor \n\n",num);

	base->Monitors = malloc(sizeof(monitors));
	base->Monitors->number = num;	
	filename = malloc(sizeof(char)*50);
	temp     = malloc(sizeof(char)*50);
	
	if(num>0)
	{
		base->Monitors->zone_num=malloc(sizeof(int)*num);
		base->Monitors->i	=malloc(sizeof(int)*num);
		base->Monitors->j	=malloc(sizeof(int)*num);
		base->Monitors->k	=malloc(sizeof(int)*num);
		printf("Use GridZ to findout the sequential number of the zone to be monitored\n");
		printf("Use Squential block traversal option and find out the number of the zone by counting each zone\n");
	}
	
	printf("The number of zones in this problem are %d\n",base->Number_of_zones);
	
	for (loop=0;loop<num;loop++)
	{
		
		if(base->Number_of_zones == 1) {
			*(base->Monitors->zone_num+loop) = 1;
			printf("Zone number for monitor %d\t is 1:",loop+1);
		} else {
			printf("Enter zone number for monitor %d\t:",loop+1);
			do{
			printf("Please only enter any number from  1----------> %d (including)\t\t",base->Number_of_zones);
			scanf("%d",(base->Monitors->zone_num+loop));
			}while(*(base->Monitors->zone_num+loop)<0 || *(base->Monitors->zone_num+loop)>base->Number_of_zones);        
                	fprintf(var2,"%d\n",*(base->Monitors->zone_num+loop));        
                	fprintf(var3,"%d \tEnter zone number for monitor %d \n\n",*(base->Monitors->zone_num+loop),loop+1);
			//printf("Going down\n");
		}
		*(base->Monitors->zone_num+loop)=*(base->Monitors->zone_num+loop)-1;

		NI =*((base->CgnsZone+(*(base->Monitors->zone_num+loop)))->Size);
		NJ =*((base->CgnsZone+(*(base->Monitors->zone_num+loop)))->Size+1);
		NK =*((base->CgnsZone+(*(base->Monitors->zone_num+loop)))->Size+2);
			
		if(base->BaseParameters->flag2d==1)
		{
			do{
				printf("Please only enter any number for\ni from 0----------> %d\t\t",NI-2);
				scanf("%d",(base->Monitors->i+loop));
			}while(((*(base->Monitors->i+loop)>NI-2 || *(base->Monitors->i+loop)<0)));        
                        fprintf(var2,"%d\n",*(base->Monitors->i+loop));        
                        fprintf(var3,"%d \t \n\n",*(base->Monitors->i+loop));
			
			do{
				printf("Please only enter any number for\nj from 0----------> %d\t\t",NJ-2);
				scanf("%d",(base->Monitors->j+loop));
			}while(((*(base->Monitors->j+loop)>NJ-2 || *(base->Monitors->j+loop)<0)));        
                        fprintf(var2,"%d\n",*(base->Monitors->j+loop));        
                        fprintf(var3,"%d \t \n\n",*(base->Monitors->j+loop));
			
			*(base->Monitors->k+loop)=0;
		}
		else
		{	
			do{
				printf("Please only enter any number for\ni from 0----------> %d\t\t",NI-2);
				scanf("%d",(base->Monitors->i+loop));
			}while(((*(base->Monitors->i+loop)>NI-2 || *(base->Monitors->i+loop)<0)));        
                        fprintf(var2,"%d\n",*(base->Monitors->i+loop));        
                        fprintf(var3,"%d \t \n\n",*(base->Monitors->i+loop));
			
			do{
				printf("Please only enter any number for\nj from 0----------> %d\t\t",NJ-2);
				scanf("%d",(base->Monitors->j+loop));
			}while(((*(base->Monitors->j+loop)>NJ-2 || *(base->Monitors->j+loop)<0)));        
                        fprintf(var2,"%d\n",*(base->Monitors->j+loop));                        
                        fprintf(var3,"%d \t \n\n",*(base->Monitors->j+loop));
			
			do{
				printf("Please only enter any number for\nk from 0----------> %d\t\t",NK-2);
				scanf("%d",(base->Monitors->k+loop));
			}while(((*(base->Monitors->k+loop)>NK-2 || *(base->Monitors->k+loop)<0)));        
                        fprintf(var2,"%d\n",*(base->Monitors->k+loop));                        
                        fprintf(var3,"%d \t \n\n",*(base->Monitors->k+loop));
			
		}
	}
		
	for (loop=0;loop<num;loop++)
	{
		sprintf(temp,"%d.dat",loop+1);
		sprintf(filename,"SOLUTIONS/monitors/monitor");
		filename = strcat(filename,temp);
		//printf("file name is %s\n",filename);
		fp = fopen(filename,"a");
		if (fp==NULL) {
		   printf("SOLUTIONS/monitors directory does not seem to exist....trying to create....");
		   system("mkdir SOLUTIONS/monitors");
		   fp = fopen(filename,"a");
		   if (fp==NULL) {printf("\n Can not open file SOLUTIONS/monitors/monitor* to write....Exiting......\n");
		     exit(1);
		   } else {printf("Created\n"); }
		}
		fprintf(fp,"VARIABLES = \n");
		fprintf(fp,"\"Time\"\t\t");
		for(inner=0;inner<solvable_variables;inner++)
			fprintf(fp,"\"%s\"\t",find_name(base,inner+1));
		
		fprintf(fp,"\n");
		fclose(fp);
	}
	                
        fclose(var2);
        fclose(var3);
	free(temp);
	free(filename);
	
	
	return 0;
}

