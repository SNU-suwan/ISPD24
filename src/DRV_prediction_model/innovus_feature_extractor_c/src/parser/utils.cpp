#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <map>
#include <unordered_map>
#include <ctime>
#include <sys/time.h>
#include <sys/stat.h>

#include "structures.h"
#include "default_classes.h"
#include "global_variables.h"
#include "utils.h"
#include "../global.h"
#include "../BBox.h"

using namespace std;


//map<string,int> METAL_LAYER_IDX;
void begin_timer(struct timeval* begin){	
	gettimeofday(begin,NULL);
}
double end_timer(struct timeval* begin, struct timeval* end){
	gettimeofday(end,NULL);
	double delta=((end->tv_sec-begin->tv_sec)*1000000u + end->tv_usec-begin->tv_usec)/1.e6;
	return delta;
}
string get_direction(const point2F& locP2, const point2F& targetP2){

	return targetP2.xF-locP2.xF>0?"R":targetP2.xF-locP2.xF<0?"L":targetP2.yF-locP2.yF>0?"U":"D";
}

vector<string> split( string s )
{
	vector<string> comp;
	stringstream Tmp;
	Tmp << s;
	string stmp;
	while ( Tmp >> stmp )
		comp.push_back(stmp);

	return comp;
}

point2I convert_coords_to_idx(point2F coordsPF){
	int xI=(coordsPF.xF-OFFSET_X)/GCELL_SIZE;
	int yI=(coordsPF.yF-OFFSET_Y)/GCELL_SIZE;
	point2I idxPI(xI,yI);
	return idxPI;
}

point2F convert_idx_to_coords(point2I idxPI){
	float xF=idxPI.xI*GCELL_SIZE+OFFSET_X;
	float yF=idxPI.yI*GCELL_SIZE+OFFSET_Y;
	point2F idxPF(xF,yF);
	return idxPF;
}

bool file_exist_test(const string& nameS){
	struct stat buffer;
	return(stat(nameS.c_str(),&buffer)==0);
}

float distance_sq(const point2F &point1, const point2F &point2){
	float dist=(point1.xF-point2.xF)*(point1.xF-point2.xF)+(point1.yF-point2.yF)*(point1.yF-point2.yF);
	return dist;
}

string get_origin_name(const string &cell_nameS){
	string origin_name=cell_nameS;
	size_t num_=count(cell_nameS.begin(),cell_nameS.end(),'_');
	if(num_==2){
		size_t idx_=cell_nameS.find("_");
		idx_=cell_nameS.find("_",idx_+1);
		origin_name=cell_nameS.substr(0,idx_);
	}
	return origin_name;
}

/*
point2I round_coords_to_idx(point2F coordsPF){
	float xF=(coordsPF.xF-OFFSET_X)/GCELL_SIZE;
	float yF=(coordsPF.yF-OFFSET_Y)/GCELL_SIZE;
	int xF_roundI=(int)round(xF);
	int yF_roundI=(int)round(yF);
	point2I idxPI(xF_roundI,yF_roundI);
	return idxPI;
}
*/

void debug_test(comp* & compp){
	cell*& cellp=compp->cellCp;
	debug_test(cellp);
}
void debug_test(cell* & cellp){
	//map<string,pin*> pinM=cellp->pinsMSPp;
	unordered_map<string,pin*> pinM=cellp->pinsMSPp;
	for(auto it=pinM.begin();it!=pinM.end();it++){
		pin *&pinp=it->second;
		vector<ap*> apsVp=pinp->apsVp;
		vector<ap*> debugap;
		
		for(auto it2=pinp->apsMIVAp.begin();it2!=pinp->apsMIVAp.end();it2++){
			vector<ap*> test =it2->second;
			for(size_t i=0;i<test.size();i++){
				debugap.push_back(test[i]);
			}
		}

		for(size_t i=0;i<apsVp.size();i++){
			bool isfound=false;
			for(size_t i2=0;i2<debugap.size();i2++){
				if(apsVp[i]==debugap[i2]) isfound=true;
			}
			assert(isfound==true);
		}

		assert(apsVp.size()==debugap.size());

	}
	
}

bool is_lines_intersect(const float* point1, const float* point2){
	float s1=min(point1[0],point1[1]);
	float f1=max(point1[0],point1[1]);
	float s2=min(point2[0],point2[1]);
	float f2=max(point2[0],point2[1]);

	bool cond1=(s2<s1 && s1<f2);
	bool cond2=(s2<f1 && f1<f2);
	bool cond3=(s1<s2 && s2<f1);
	bool cond4=(s1<f2 && f2<f1);
	if(cond1||cond2||cond3||cond4){
		return true;
	}
	return false;
}

BBox<double> pin_to_bbox(pin* pinp){
	Point<double> lb(MAXF,MAXF);
	Point<double> rt(MINF,MINF);
	for(size_t it=0;it<pinp->segmentsVp.size();it++){
		segment* segmentp=pinp->segmentsVp[it];
		for(size_t it2=0;it2<segmentp->pointsVp.size();it2++){
			point2F* coord=segmentp->pointsVp[it2];
			
			lb.x=coord->xF<lb.x?coord->xF:lb.x;
			lb.y=coord->yF<lb.y?coord->yF:lb.y;
			
			rt.x=coord->xF>rt.x?coord->xF:rt.x;
			rt.y=coord->yF>rt.y?coord->yF:rt.y;
		}
	}
	return BBox<double>(lb,rt);
}

map<int,string> build_metal_idx_layer(){
	map<int,string> metal_idx_layer;
	for(auto it=METAL_LAYER_IDX.begin();it!=METAL_LAYER_IDX.end();it++){
		metal_idx_layer.insert(make_pair(it->second,it->first));
	}
	return metal_idx_layer;
}
