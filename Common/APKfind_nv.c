#include "APKfind_nv.h"
/**
This function, find_nv() is required to get the tags of the variables that are solved.
*/

int find_nv(Base *base,int vartag)
{
	int total,i;
	total = base->total_var;
	for (i=0;i<total;i++)
	{	
		if (base->VARTAG[i]==vartag)
		return(i);
	}
	printf ("Wrong Vartag %d supplied, exiting from find_nv\n",vartag);	
	exit(0);
}
