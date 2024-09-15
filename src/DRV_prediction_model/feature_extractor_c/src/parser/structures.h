#ifndef __STRUCTURES__
#define __STRUCTURES__
#include <unordered_map>
#include <map>
#include <vector>
#include "../BBox.h"
#include "lef.h"
#include "global_variables.h"
#include "default_classes.h"

using namespace std;

class Lef;
class Def;

//lef
class instance;
class cell;
class macro;
class pin;
class segment;
class ap;

//def
class comp_instance;
class comp;
class comp_macro;
class row;
class IOpin;
class net;
cell* copyCell(cell*);
macro* copyMacro(macro*);

//drv
class drv;

// FLUTE
class graph;
class vertex;
class edge;
class gcell_grid;
class gcell;
class gedge;
class steiner_container;


enum vertex_type{
	//typeI
	CELL,
	STEINER,
	DUMMY
};

enum pin_type{
	//pinTypeI
	STEINER_PIN,
	PIN,
	IOPIN,
	DUMMY_PIN
};


class instance{
public:
	instance(string,point2F,unordered_map<string,pin*>);
	~instance();
public:
	string nameS;
	string dirS;
	point2F *originP2p;
	point2F sizeP2;
	unordered_map<string,pin*> pinsMSPp; 
public:
	void setPos(const point2F&,const point2F&);
	void setDir(const string&);
	void flipVer();
	void flipHor();
	void update_pin_center();
	void print_pins(bool is_print_ap=false);

};

class cell: public instance{
public:
	cell(string,point2F,unordered_map<string,pin*>);
	~cell();	
public:
	comp * compp;
public:
};

class macro: public instance{
public:
	macro(string,point2F,unordered_map<string,pin*>);
	~macro();	
public:
	comp_macro * compp;

public:
};

class pin{
public:
	pin(string);
	~pin();
public:
	string nameS;
	string dir;
	point2F centerP2; //centroid of the pin
	point2F center_ingridP2;
	double avg_ap_x_ingrid;
	vector<segment*> segmentsVp;
	vector<ap*> apsVp;
	unordered_map<int,vector<ap*>> apsMIVAp;
	net *netp;
	cell *cellp;
	macro *macrop;

	int targetTypeI;
	void* target;
	point2F* targetP2;
	
public:
	void createAps();
	void createApMap();
	void destroySegments();
	void destroyAps();
	void setPos(const point2F&,const point2F&);
	void flipVer(const point2F&,const point2F&);
	void flipHor(const point2F&,const point2F&);
};

class segment{
public:
	segment(string,vector<point2F*>);
	~segment();
public:
	string layerS;
	vector<point2F*> pointsVp;
	void setPos(const point2F&,const point2F&);
	void flipVer(const point2F&,const point2F&);
	void flipHor(const point2F&,const point2F&);
};

class ap{
public:
	ap(string,point2F);
	~ap();
public:
	string layerS;
	point2F pointP2;
	pin* pinp;
	int trackI;
public:
	void setPos(const point2F& locP2, const point2F&);
	void flipVer(const point2F&,const point2F&);
	void flipHor(const point2F&,const point2F&);
};

class comp_instance{
public:
	comp_instance(string,string,point2F);
public:
	string nameS;
	string dirS;
	point2F locP2;
public:
};


class comp: public comp_instance{
public:
	comp(string,string,point2F,cell*);
	~comp();
public:
	cell* cellCp;
public:
	bool operator < (const comp&);
	void update_pin_center();
	void setPos(const point2F&);
};

class comp_macro: public comp_instance{
public:
	comp_macro(string,string,point2F,macro*);
	~comp_macro();
public:
	macro* macroMp;
public:
	void update_pin_center();
	void setPos(const point2F&);
};

class row{
public:
	row(float,float);
	~row();
public:
	float rowHeightF;
	float rowX;
};


class IOpin{
public:
	IOpin(string,string,point2F,point2F*);
	~IOpin();
public:
	string nameS;
	string layerS;
	point2F centerP2;
	point2F sizeP22[2];
	net* netp;
};

class net{
public:
	net(string);
	~net();
public:
	string nameS;
	unordered_map<string,string> cellPinMSS; //be careful of PIN
	graph *steinerG;
};

class drv{
public:
	drv(int,string,string,vector<string>,BBox<double>);
	~drv();
public:
	int idI;
	string typeS;
	string layerS;
	vector<string> relatedNetsVS;
	BBox<double> bbox;
};

// for FLUTE
class graph{
public:
	graph(const vector<vertex*> &, vector<edge> &);
	~graph();
public:
	map<int,int> verticiesMII;//vertex idx to number of edges
	map<int,size_t> verIdxIdxMII;//vertex idx to verticiesppp idx
	vertex *** verticiesppp;
	vector<edge*> edgesVEp;
	int verticies_size;
public:
	//void addEdge(const edge&);
	void add_vertex(vertex*,vector<vertex*>);
	void delete_edge(edge*);
	void change_edge(int,int,int);
	void move_vertex(const int&,const point2F&);

	void delete_edge_and_add_vertex(edge*,vertex*,vector<vertex*>&);
	void print();
	
};

class vertex{
public:
	vertex();
	~vertex();
public:
	int idxI; // 1,2,3 for pin / 10001 10002 10003 for steiner
	int typeI; // CELL STEINER DUMMY
	int pinTypeI; // PIN IOPIN
	pin* pinp;
	IOpin* IOpinp;
	point2F centerP2;
public:
	bool operator == (const vertex&);
	vertex operator = (const vertex&);

};

class edge{
public:
	edge(bool,string,point2F*);
	edge(point2F, point2F);
	edge(vertex*, vertex*);
	edge();
	~edge();
public:
	bool direction;
	string layerS; //layer or direction
	point2F endPointsP22[2];
	vertex *endVerticies[2];
	int num_apsI;
	ap* source_app;
public:
	void clear();
	//bool operator < (const edge&) const;
	//edge operator = (const edge&);

};

// grid info

class gcell{
public:
	gcell(point2F);
	~gcell();
public:
	point2F originPF; //left bottom coordinates
	gedge *upgE;
	gedge *downgE;
	gedge *leftgE;
	gedge *rightgE;
	float demand;
};

class gedge{
public:
	gedge(point2I*,point2F*);
	~gedge();
public:
	float demandF;
	point2I endPointsPI2[2];
	point2F endPointsPF2[2];
};

class gcell_grid{
public:
	gcell_grid(Def&);
	~gcell_grid();
public:
	Def *def;
	int num_grids_x;
	int num_grids_y;
	map<point2I,gcell*> gcellsMPIGp;
public:
	void assign_net_info();
	void insert_graph_info_into_gcell_grid(graph*,bool is_remove=false);
	edge search_min_demand_edge(point2F&, point2F&,int);
	float search_demand_for_line(const point2F&, const point2F&);

	void print(int,int);
};

class steiner_container{
public:
	steiner_container(vertex*, const graph*, int);
	steiner_container();
	~steiner_container();
public:
	vertex* steiner_vertexp;
	edge sliding_edge[2];//[0] for horizontal, [1] for vertical 
	/*
	vertex* min_xvp;
	vertex* max_xvp;
	vertex* min_yvp;
	vertex* max_yvp;
	*/
	vector<vertex*> adjacent_verticies;
	vector<vertex*> adjacent_steiner;
	vector<steiner_container> adjacent_steiner_container;
	int vertex_degree;
	int vertex_idx;
	edge extended_sliding_range[2];
	//map<int,edge*> extended_sliding_edgeM[2]; //0 for x, 1 for y
	map<int,edge> extended_sliding_edgeM[2]; //0 for x, 1 for y
public:
	void destroy();
	void extend_sliding_edge(const steiner_container&,int);
};


#endif
