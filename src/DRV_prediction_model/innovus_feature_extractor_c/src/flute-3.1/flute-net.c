#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include "flute.h"


int main(int argc, char* argv[])
{
	int d=0;
    int x[MAXD], y[MAXD];
    Tree flutetree;
    int flutewl;
    
	int it=0;
	for(int i=1;i<argc;i+=2){
		x[it]=atof(argv[i])*1000;
		y[it]=atof(argv[i+1])*1000;
		it++;
	}
	d=it;

	readLUT();

    flutetree = flute(d, x, y, ACCURACY);
    //printf("FLUTE wirelength = %d\n", flutetree.length);

    flutewl = flute_wl(d, x, y, ACCURACY);
    //printf("FLUTE wirelength (without RSMT construction) = %d\n", flutewl);

	printtree(flutetree);
	plottree(flutetree);
}
