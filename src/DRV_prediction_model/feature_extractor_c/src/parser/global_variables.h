#ifndef __GLOBAL_VARIABLES_H__
#define __GLOBAL_VARIABLES_H__
#include <string>
using namespace std;

#define TECH 15 // 7 for ASAP7, 15 for NANGATE15

/////////// Nangate 15 ////////////
#if TECH == 15
	//Cost
	#define GCELL_SIZE 0.768 

	//Design Rules
	#define POLY_PITCH 0.064f
	#define M1_PITCH 0.064f
	#define M2_PITCH 0.064f
	#define M2_FIRSTLAST 0.064f
	#define CELL_HEIGHT 0.768f
	#define NUM_TRACK 12
	#define VIA_OVERHANG 0.09f //0.014(Via length/2)+0.032(Overhang)+0.044(EOL)

	//Layer Name
	#define METAL1 "M1"
	#define METAL2 "MINT1"
	#define METAL3 "MINT2"
	#define METAL4 "MINT3"
	#define METAL5 "MINT4"
	#define METAL6 "MINT5"
	#define VIA1 "V1_0"
	#define VIA2 "VINT1_0"
	#define VIA3 "VINT2_0"
	#define VIA4 "VINT3_0"
	#define VIA5 "VINT4_0"


///////////// ASAP 7 ///////////////
//Cost
#elif TECH == 7

	#define GCELL_SIZE 1.08 
	
	#define POLY_PITCH 0.072f
	#define M1_PITCH 0.072f // Originally 0.144 but the cells can be located on the half point of the grid
	#define M2_PITCH 0.144f
	#define M2_SPACE 0.072f
	#define M2_FIRSTLAST 0.180f
	#define M3_PITCH 0.144f

	#define CELL_HEIGHT 1.08f
	#define NUM_TRACK 7
	#define VIA_OVERHANG 0.056f //0.014(Via length/2)+0.032(Overhang)+0.044(EOL)

	//Layer Name
	#define METAL1 "M1"
	#define METAL2 "M2"
	#define METAL3 "M3"
	#define METAL4 "M4"
	#define METAL5 "M5"
	#define METAL6 "M6"
	#define VIA1 "V1"
	#define VIA2 "V2"
	#define VIA3 "V3"
	#define VIA4 "V4"
	#define VIA5 "V5"
#endif


#define PRECISION 0.0001

#define VDD "VDD"
#define VSS "VSS"

//P&R
#define OFFSET_X 5
#define OFFSET_Y 5

//Other Variables
#define MINF -10000000.0f
#define MAXF 10000000.0f
#define MINI -10000000
#define MAXI 10000000
#define HORIZONTAL true
#define VERTICAL false


#endif
