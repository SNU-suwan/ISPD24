#ifndef __DEF_H__
#define __DEF_H__

#include "structures.h"
using namespace std;

class Def{
public:
	Def(string,Lef&,bool);
	~Def();
public:
	bool is_cadence;
	int scaleI;
	float versionF;
	point2F die_sizeP22[2]; // {x,y}, {x,y}
	string defDirS;
	string outDefDirS;

	vector<row*> rowsVp;
	
	Lef& lef;
	map<string,IOpin*> IOpinsMSpp;
	map<string,comp*> compsMSCp;
	map<string,comp_macro*> comp_macrosMSCp;
	map<string,net*> netsMSNp;
	map<float,vector<comp*>> compsMFVCp;

public:
	void parse();
	void write_def(string);
};

char * jumpNTimes(char *token, int n);
#endif
