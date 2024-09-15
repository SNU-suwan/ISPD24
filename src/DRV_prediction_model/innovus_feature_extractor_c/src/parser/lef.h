#ifndef __LEF_H__
#define __LEF_H__

//#include "../global.h"
#include "default_classes.h"
#include "structures.h"
using namespace std;

class cell;
class macro;

class Lef{
public:
	Lef(string);
	~Lef();
public:
	string lefDirS;
	map<string,cell*> cellsMSCp;
	map<string,macro*> macrosMSMp;
public:
	void parse();
	void build_metal_layer_map();
	void build_metal_routing_track_map();
	void build_metal_routing_direction_map();
	void build_clock_units();
	void rectToPoly(vector<point2F*>&);
};

#endif
