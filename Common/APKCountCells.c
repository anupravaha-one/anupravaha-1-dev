#include "APKCountCells.h"
/**
This function calculates the number of Cells in a particular Zone. It takes Zone pointer and returns Number of 
Cells in that Zone.  
*/
int CountCells(Zone *zone)
{
	int totalCells;
	int *size;
	size = zone->Size;
	totalCells=(*(size+0)-1)*(*(size+1)-1)*(*(size+2)-1);
	
	return (totalCells);
}
