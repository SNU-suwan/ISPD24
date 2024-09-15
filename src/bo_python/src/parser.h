#ifndef PLACE_PARSER_H__
#define PLACE_PARSER_H__

#include <map>
#include <vector>
#include <tuple>
#include "BBox.h"
#include "cell.h"

using namespace std;

extern vector<string> CLOCK_CELLS;
extern vector<string> CLOCK_CELLS_NAME;
extern string INVERTER_NAME;

class Lef{
public:
	Lef(){};
	~Lef();
public:
	string lef_dir;
	int scale = 1000;
	std::map<std::string,standard_cell*> standard_cell_map;
public:
	void Parse();
};

class Def{
public:
	Def(){
		is_cadence=false;
	};
	~Def();

public:
	std::string def_dir;
	bool is_cadence;
	bool is_cts_design;
	
	int scale;
	int row_height;
	int cpp;

	std::string version;
	BBox<int> die_size;
	Point<int> offset;
    int chip_x;

	std::map<std::string,placed_cell*> placed_cells_map;
	std::map<int,std::vector<placed_cell*>> placed_cells_row_map;
	std::map<std::string,net*> net_map;

	std::map<std::string,tuple<string,string>> cell_prefix_postfix_map;

public:
	int GetScale();

	void Parse(Lef& lef);
	void WriteDef(string out_def_dir);
	void ValidationTest();
};

class ParserWrapper{
public:
	ParserWrapper(std::string lef_dir_, std::string def_dir_, bool is_cts_design_, bool is_cadence_)
		: _lef_dir(lef_dir_), _def_dir(def_dir_), is_cts_design(is_cts_design_) {
			lef.lef_dir=lef_dir_;
			def.def_dir=def_dir_;
			def.is_cts_design=is_cts_design_;
			def.is_cadence=is_cadence_;
		};
	~ParserWrapper();

	Lef lef;
	Def def;
	bool is_cts_design;
	void ParseLefDef();
private:	
	std::string _lef_dir;
	std::string _def_dir;
};

char * JumpNTimes(char *token, int n);
#endif
