#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <tuple>

#include "parser.h"
#include "cell.h"
#include "utils.h"

using namespace std;


char * JumpNTimes(char *token, int n){
	char * jumped_token;
	for (int i=0;i<n;i++){
		token=strtok(NULL," ");
	}
	jumped_token=token;
	return jumped_token;
}

Lef::~Lef(){
	for(auto it=standard_cell_map.begin();it!=standard_cell_map.end();it++){
		delete it->second;
	}
}

void Lef::Parse(){
	ifstream ifs;
	ifs.open(lef_dir);

	int cell_width;
	int cell_height;
	bool is_clock_cell;
	string standard_cell_name;
	
	if(ifs.is_open()){
		char lineC[1024];
		while(!ifs.eof()){
			ifs.getline(lineC,1024);
			char *token=strtok(lineC," ");
			while (token != NULL){
				string tmp_str(token);
				if (tmp_str=="MACRO"){
					token=strtok(NULL," ");
					string tmp_str(token);
					standard_cell_name=tmp_str;
				}
				else if(tmp_str=="SIZE"){
					token=strtok(NULL," ");
					float width=atof(token);
					token=strtok(NULL," ");
					token=strtok(NULL," ");
					float height=atof(token);
					
					cell_width=width*scale;
					cell_height=height*scale;
				}
				else if (tmp_str=="END"){
					token=strtok(NULL," ");
					if (token!=NULL){
						string tmp_str(token);
						if (tmp_str==standard_cell_name){
							standard_cell *new_standard_cell = new standard_cell{
								standard_cell_name, 
								cell_width, 
								cell_height
								};
							standard_cell_map.insert(make_pair(standard_cell_name,new_standard_cell));
						}
					}
				}
				token=strtok(NULL," ");
			}
		}
	}
	ifs.close();
}

Def::~Def(){
	for(auto it=placed_cells_map.begin();it!=placed_cells_map.end();it++){
		delete it->second;
	}
	placed_cells_map.clear();
}

void Def::Parse(Lef &lef){
	ifstream ifs;
	ifs.open(def_dir);

	char line[1024];
	bool is_write_components;
	bool is_offset_written=false;
	bool is_row_height_written=false;
	bool is_write_nets = false;
	bool is_write_net = false;

	std::string net_name;

	if(ifs.is_open()){
		while(!ifs.eof()){
			ifs.getline(line,1024);
			char *token=strtok(line," ");
			while(token !=NULL){
				string tmp_str(token);
				if (tmp_str=="UNITS"){
					token=JumpNTimes(token,3);
					scale=(int)atof(token);
				}
				else if (tmp_str=="VERSION"){
					token=strtok(NULL," ");
					version=string(token);
				}
				else if(tmp_str=="DIEAREA"){
					if(!is_cadence){
						token=JumpNTimes(token,2);
						int die_x_lb=atoi(token);
						token=JumpNTimes(token,1);
						int die_y_lb=atoi(token);
						token=JumpNTimes(token,4);
						int die_y_rt=atoi(token);
						token=JumpNTimes(token,3);
						int die_x_rt=atoi(token);
						die_size={die_x_lb,die_y_lb,die_x_rt,die_y_rt};
					}
					else{
						token=JumpNTimes(token,2);
						int die_x_lb=atoi(token);
						token=JumpNTimes(token,1);
						int die_y_lb=atoi(token);
						token=JumpNTimes(token,3);
						int die_x_rt=atoi(token);
						token=JumpNTimes(token,1);
						int die_y_rt=atoi(token);
						die_size={die_x_lb,die_y_lb,die_x_rt,die_y_rt};
					}
				}
				else if(tmp_str=="ROW"){
					if (!is_offset_written && !is_row_height_written){
						token=JumpNTimes(token,3);
						int offset_x=atoi(token);
						token=JumpNTimes(token,1);
						int offset_y=atoi(token);

						offset={offset_x, offset_y};
						is_offset_written=true;

						token=JumpNTimes(token,3);
                        int step = atoi(token);
						token=JumpNTimes(token,4);
						cpp=atoi(token);

                        chip_x=step*cpp;
					}
					else if(is_offset_written && !is_row_height_written){
						tmp_str = string(token);
						token=JumpNTimes(token,4);
						int second_row_position=atoi(token);
						row_height=second_row_position-offset.y;
						is_row_height_written=true;
					}
				}
				else if(tmp_str=="COMPONENTS"){
					is_write_components=true;
				}
				else if(tmp_str=="NETS"){
					is_write_nets=true;
				}
				else if(tmp_str=="-"){
					if (is_write_components){
						std::string prefix = "- ";
						std::string postfix = "";

						token=strtok(NULL," ");
						std::string cell_name(token);
						prefix += cell_name + " ";
						
						token=strtok(NULL," ");
						std::string cell_type(token);
						prefix += cell_type + " ";

						bool is_clock_cell = false;
					
						token=strtok(NULL," ");
						std::string tmp(token);
						prefix += tmp + " ";
					
						while (tmp!="("){
							token=strtok(NULL," ");
							tmp.assign(token);
							prefix += tmp + " ";
						}

						token=strtok(NULL," ");
						int x=atoi(token);
						token=strtok(NULL," ");
						int y=atoi(token);
						
						token=JumpNTimes(token,2);
						std::string direction=tmp.assign(token);
						placed_cell *new_placed_cell=new placed_cell {
							x,
							y,
							-1,
							false,
							cell_name,
							direction,
							lef.standard_cell_map[cell_type],
							is_clock_cell
							};
					
						placed_cells_map.insert(make_pair(cell_name,new_placed_cell));
						placed_cells_row_map[y].push_back(new_placed_cell);

						token=strtok(NULL," ");
						while(token !=NULL){
							tmp.assign(token);
							postfix += " " + tmp;
							token=strtok(NULL," ");
						}
						if (token == NULL && tmp != ";"){
							postfix += " ;";
						}
						std::tuple<std::string, std::string> prefix_postfix = {prefix, postfix};
						cell_prefix_postfix_map.insert(make_pair(cell_name, prefix_postfix));
					}
					if(is_write_nets){
						token=strtok(NULL," ");
						net_name.assign(token);
						net *new_net = new net(net_name);
						net_map.insert(make_pair(net_name, new_net));
						is_write_net=true;	
					}
				}
				else if(tmp_str == "("){
					if(is_write_net){
						token=strtok(NULL," ");
						std::string comp_name(token);
						token=strtok(NULL," ");
						std::string pin_name(token);
						tuple<std::string, std::string> cell_pin_tuple = make_tuple(comp_name, pin_name);
						net_map.at(net_name)->cell_pin_tuple_vec.push_back(cell_pin_tuple);
					}
				}
				else if(tmp_str == "+"){
					if(is_write_net){
						token=strtok(NULL," ");
						std::string tmp_str(token);
						if(tmp_str == "ROUTED" && is_cts_design){
							for(int i = 0; i < net_map.at(net_name)->cell_pin_tuple_vec.size(); i++){
								std::string &comp_name = get<0>(net_map.at(net_name)->cell_pin_tuple_vec[i]);
								if(placed_cells_map.find(comp_name) != placed_cells_map.end()){
									placed_cells_map.at(comp_name)->is_cts_cell = true;
								}
							}
							is_write_nets = false;
						}
						else if(tmp_str == "USE"){
							is_write_net = false;
						}
					}
				}
				else if(tmp_str=="END"){
					token=strtok(NULL," ");
					if (token!=NULL){
						tmp_str.assign(token);
						if (tmp_str=="COMPONENTS") is_write_components=false;
					}
				}
				token=strtok(NULL," ");
			}
		}
	}
	ifs.close();
	
	//Sort comps according to its x coordinate.
	for(auto it=placed_cells_row_map.begin();it!=placed_cells_row_map.end();it++){
		vector<placed_cell*>& placed_cells_vector=it->second;
		sort(placed_cells_vector.begin(),placed_cells_vector.end(),
				[](placed_cell* lcomp, placed_cell* rcomp) ->bool{
					return lcomp->x < rcomp->x;
				});
		for(size_t i_c = 0 ; i_c < placed_cells_vector.size(); i_c++){
			placed_cells_vector[i_c]->idx = i_c;
		}
	}
}

int Def::GetScale(){
	ifstream ifs;
	ifs.open(def_dir);

	char line[1024];

	if(ifs.is_open()){
		while(!ifs.eof()){
			ifs.getline(line,1024);
			char *token=strtok(line," ");
			while(token !=NULL){
				string tmp_str(token);
				if (tmp_str=="UNITS"){
					token=JumpNTimes(token,3);
					scale=(int)atof(token);
				}
				token=strtok(NULL," ");
			}
		}
	}
	ifs.close();
	return scale;
}

void Def::ValidationTest(){
	for(auto it = placed_cells_row_map.begin(); it != placed_cells_row_map.end(); it++){
		vector<placed_cell*> placed_cells = it->second;
		for(size_t i = 0; i < placed_cells.size() -1 ; i++){
			placed_cell* current_cell = placed_cells[i];
			placed_cell* next_cell = placed_cells[i+1];
			if(current_cell->x + current_cell->cell->cell_width > next_cell->x){
				cout<<"[ERROR] "
					<< current_cell->name <<" overlaps " << next_cell->name 
					<<" current right :" << current_cell->x + current_cell->cell->cell_width 
					<<" next left :" << next_cell->x
					<<endl;
			}
		}
	}
}

void Def::WriteDef(string out_def_dir){
	cout<<"Write def to "<<out_def_dir<<endl;
	ifstream ifs;
	ifs.open(def_dir);
	ofstream ofs;
	ofs.open(out_def_dir);
	
	char line[1024];
	string str_line;
	bool is_write=true;
	bool is_write_comps=false;
	if(ifs.is_open()){
		while (!ifs.eof()){
			ifs.getline(line,1024);
			string str_line=line;
			string tmp_str;
			char *token=strtok(line," ");
			if(token !=NULL){
				tmp_str.assign(token);
				if(tmp_str=="COMPONENTS"){
					token=strtok(NULL," ");
					tmp_str.assign(token);
					ofs<<"COMPONENTS "+tmp_str+" ;"<<endl;
					is_write_comps=true;
					is_write=false;
				}
				else if(tmp_str=="END"){
					token=strtok(NULL," ");
					tmp_str.assign(token);
					if(tmp_str=="COMPONENTS"){
						is_write_comps=false;
						is_write=true;
					}
				}
				if(is_write_comps){
					for(auto it=placed_cells_row_map.begin();it!=placed_cells_row_map.end();it++){
						vector<placed_cell*>& placed_cells_vector=it->second;
						for(size_t i=0;i<placed_cells_vector.size();i++){
							placed_cell* placed_cell_pointer=placed_cells_vector[i];

							std::string direction = placed_cell_pointer->original_direction;
							if(placed_cell_pointer->is_flipped){
									if(placed_cell_pointer->original_direction == "S") direction = "FS";
									else if(placed_cell_pointer->original_direction == "N") direction = "FN";
									else if(placed_cell_pointer->original_direction == "FS") direction = "S";
									else if(placed_cell_pointer->original_direction == "FN") direction = "N";
							}

							string write_comp = get<0>(cell_prefix_postfix_map[placed_cell_pointer->name])
								+ to_string(placed_cell_pointer->x)+" "
								+ to_string(placed_cell_pointer->y)+" ) "
								+ direction
								+ get<1>(cell_prefix_postfix_map[placed_cell_pointer->name]);
								
							ofs<<write_comp;
							ofs<<endl;
						}
					}
					is_write_comps=false;
				}
			}
			if(is_write){
				ofs<<str_line;
				ofs<<endl;
			}
		}
	}
	ifs.close();
	ofs.close();
}

ParserWrapper::~ParserWrapper(){};

void ParserWrapper::ParseLefDef(){
	std::cout<<"Lef file : "<< _lef_dir << endl;
	std::cout<<"Def file : "<< _def_dir << endl;
	lef.scale=def.GetScale();
	lef.Parse();
	def.Parse(lef);
	std::cout<<"Parsing Lef, Def file has done!"<<std::endl;
}
