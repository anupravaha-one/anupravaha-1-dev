#include "APKFind_Direction.h"
#include "APKDirection.h"
/**
This function, Find_Direction() assigns direction values to the faces of a zone. That face can be an interface or boundary.
*/

int Find_Direction(int *direction,int *p1)
{

	int diff[3],i;
	
	for(i=0;i<3;i++)
		diff[i] = *(p1+i) - *(p1+i+3);

	if (diff[0]==0)
	{
		if(*p1 == 1) *direction = West;
		else	     *direction = East;
	}
	
	 else if (diff[1]==0)
	{
		if(*(p1+1) == 1) *direction = South;
		else	         *direction = North;
	}
	
	else if (diff[2]==0)
	{
		if(*(p1+2) == 1) *direction = Bottom;
		else	         *direction = Top;
	}
	else
		*direction = No_Direction;
	
	/*printf("Direction of Interface is:%d\n",*direction);
	fprintf(k,"Direction of Interface is:%d\n\n",*direction);*/
	return 0;
}

	
