/* ******************************************************************************************** *
*                                     File : Misc.h                                             *
*  ******************************************************************************************** *
*  File Type   : C Header File                                                                  *
*  Description : This Header file contains function declarations(prototypes) which may be       *
*                needed to be included in different modules.Instead of declaring the same       *
*                function in each module wherever they are needed,they are declared once in this* 
*                header file,so that this header file is included wherever the following function*
*                prototypes are needed.                                                         *
*  ******************************************************************************************** */


typedef char WARNING;
typedef int INDEX;


WARNING Memory_Warning[]      = "Memory Allocation failed. Not enough memory available";
WARNING No_DataBase_Warning[] = "Invalid CGNS File. No Databases found";
WARNING No_Zone_Warning[]     = "Invalid CGNS File. No Zones found for one or more bases";
WARNING File_Warning[]        = "Requested file not found";
WARNING Direction_Warning[]   = "Direction Code -1 encountered";

