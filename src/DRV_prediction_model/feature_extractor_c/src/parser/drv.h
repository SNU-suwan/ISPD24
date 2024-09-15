#ifndef __DRV_H__
#define __DRV_H__

#include "default_classes.h"
#include "structures.h"
using namespace std;

class Drv{
public:
	Drv();
	Drv(string);
	~Drv();
public:
	string drvDirS;
	vector<drv> drvV;
public:
	void parse();
};

#endif
