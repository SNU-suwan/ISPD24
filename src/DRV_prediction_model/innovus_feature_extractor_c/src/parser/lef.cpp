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

#include "global_variables.h"
#include "default_classes.h"
#include "structures.h"
#include "utils.h"
#include "lef.h"
#include "../global.h"

map<string,int> METAL_LAYER_IDX;
//map<string,int> METAL_ROUTING_TRACK;
//map<string,string> METAL_ROUTING_DIRECTION;
map<int,int> METAL_ROUTING_TRACK;
map<int,string> METAL_ROUTING_DIRECTION;
vector<string> CLOCK_UNITS;

using namespace std;

Lef::Lef(string _lefDirS){
	lefDirS=_lefDirS;
	build_metal_layer_map();
	build_metal_routing_track_map();
	build_metal_routing_direction_map();
	build_clock_units();
}
Lef::~Lef(){
	for(auto it=cellsMSCp.begin();it!=cellsMSCp.end();it++){
		delete it->second;
	}
}
void Lef::build_metal_layer_map(){
	METAL_LAYER_IDX.insert(make_pair(METAL1,1));
	METAL_LAYER_IDX.insert(make_pair(METAL2,2));
	METAL_LAYER_IDX.insert(make_pair(METAL3,3));
	METAL_LAYER_IDX.insert(make_pair(METAL4,4));
	METAL_LAYER_IDX.insert(make_pair(METAL5,5));
	METAL_LAYER_IDX.insert(make_pair(METAL6,6));
}

void Lef::build_metal_routing_track_map(){
	METAL_ROUTING_TRACK.insert(make_pair(1,M1_TRACKS));
	METAL_ROUTING_TRACK.insert(make_pair(2,M2_TRACKS));
	METAL_ROUTING_TRACK.insert(make_pair(3,M3_TRACKS));
	METAL_ROUTING_TRACK.insert(make_pair(4,M4_TRACKS));
	METAL_ROUTING_TRACK.insert(make_pair(5,M5_TRACKS));
	METAL_ROUTING_TRACK.insert(make_pair(6,M6_TRACKS));
}

void Lef::build_metal_routing_direction_map(){
	METAL_ROUTING_DIRECTION.insert(make_pair(1,M1_ROUTING_DIRECTION));
	METAL_ROUTING_DIRECTION.insert(make_pair(2,M2_ROUTING_DIRECTION));
	METAL_ROUTING_DIRECTION.insert(make_pair(3,M3_ROUTING_DIRECTION));
	METAL_ROUTING_DIRECTION.insert(make_pair(4,M4_ROUTING_DIRECTION));
	METAL_ROUTING_DIRECTION.insert(make_pair(5,M5_ROUTING_DIRECTION));
	METAL_ROUTING_DIRECTION.insert(make_pair(6,M6_ROUTING_DIRECTION));
}

void Lef::build_clock_units(){
	//One can include more items.
	CLOCK_UNITS.push_back(BUF1);
	CLOCK_UNITS.push_back(BUF2);
	CLOCK_UNITS.push_back(INV1);
	CLOCK_UNITS.push_back(DFF1);
	CLOCK_UNITS.push_back(NAND1);
	CLOCK_UNITS.push_back(AOI1);
	CLOCK_UNITS.push_back(OAI1);
}

void Lef::parse(){

	ifstream ifs;
	ifs.open(lefDirS);

	bool isRect=false;
	bool isSkipPin=false;
	bool isMacro=false;

	string cellNameS;
	string segLayerS;
	string pinNameS;
	
	point2F sizeP2;
	
	unordered_map<string,pin*> pinsMSPp;


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
					cellNameS=tmp_str;
				}
				else if(tmp_str=="CLASS"){
					token=strtok(NULL," ");
					string tmp_str(token);
					if(tmp_str=="BLOCK"){
						isMacro=true;
					}
					else{
						isMacro=false;
					}
				}
				else if(tmp_str=="SIZE"){
					token=strtok(NULL," ");
					float xF=atof(token);
					token=strtok(NULL," ");
					token=strtok(NULL," ");
					float yF=atof(token);
					sizeP2={xF,yF};
				}
				else if(tmp_str=="PIN"){
					token=strtok(NULL," ");
					string tmp_str(token);
					if (tmp_str!=VDD && tmp_str!=VSS){
						pinNameS=tmp_str;
						pin* newPinPp=new pin(pinNameS);
						pinsMSPp.insert(make_pair(pinNameS,newPinPp));
						isSkipPin=false;
					}
					else{
						isSkipPin=true;
					}
				}
				else if(tmp_str=="OBS"){
					//Skip the OBS
					isSkipPin=true;
					break;
				}
				else if(tmp_str=="DIRECTION" && !isSkipPin){
					token = strtok(NULL, " ");
					pinsMSPp[pinNameS]->dir = token;
					break;
				}
				else if(tmp_str=="LAYER" && !isSkipPin){
					token=strtok(NULL," ");
					string tmp_str(token);
					segLayerS=tmp_str;
				}
				else if ((tmp_str=="RECT" || tmp_str=="POLYGON") && !isSkipPin){
					if (tmp_str=="RECT")	
						isRect=true;
					else  isRect=false;
						
					
					vector<point2F> tmpVP2;
					point2F tmpP2;
					int idx=-1;
					while(*token != ';'){
						idx+=1;
						token=strtok(NULL," ");
						if (*token != ';'){
							if (idx%2==0){
								tmpP2.xF=atof(token);
							}
							else{
								tmpP2.yF=atof(token);
								tmpVP2.push_back(tmpP2);
							}
						}
					}
					
					vector<point2F*> segCoordsVP2p;
					for(size_t itv=0;itv<tmpVP2.size();itv++){
						point2F* coordP2p=new point2F;
						*coordP2p=tmpVP2[itv];
						segCoordsVP2p.push_back(coordP2p);
					}
					if(isRect) 
						rectToPoly(segCoordsVP2p);
					segment *newSegment=new segment(segLayerS,segCoordsVP2p);
					pinsMSPp[pinNameS]->segmentsVp.push_back(newSegment);
					}
				else if (tmp_str=="END"){
					token=strtok(NULL," ");
					if (token!=NULL){
						string tmp_str(token);
						if (tmp_str==cellNameS){

							if(!isMacro){
								cell *newCell=new cell(cellNameS,sizeP2,pinsMSPp);
								cellsMSCp.insert(make_pair(cellNameS,newCell));
								//string origin_name=get_origin_name(cellNameS);

							}
							else{
								macro *newMacro=new macro(cellNameS,sizeP2,pinsMSPp);
								macrosMSMp.insert(make_pair(cellNameS,newMacro));
								//string origin_name=get_origin_name(cellNameS);
							}
							

							
							pinsMSPp.clear();

						}
						else if (tmp_str==pinNameS && !isSkipPin){
							if(!isMacro){
								pinsMSPp[pinNameS]->createAps();
								pinsMSPp[pinNameS]->createApMap();
							}
						}
					}
				}
				token=strtok(NULL," ");
				}
		}
	}
	ifs.close();
}

void Lef::rectToPoly(vector<point2F*>& rectV){
	//Convert rect to poly
	point2F ldP2=*rectV[0];
	point2F ruP2=*rectV[1];
	point2F rdP2={ruP2.xF,ldP2.yF};
	point2F luP2={ldP2.xF,ruP2.yF};
	
	for(size_t it=0;it<rectV.size();it++) delete rectV[it];
	rectV.clear();

	point2F *ldP2p=new point2F;
	point2F *ruP2p=new point2F;
	point2F *rdP2p=new point2F;
	point2F *luP2p=new point2F;

	ldP2p->xF=ldP2.xF;
	ldP2p->yF=ldP2.yF;
	ruP2p->xF=ruP2.xF;
	ruP2p->yF=ruP2.yF;
	rdP2p->xF=rdP2.xF;
	rdP2p->yF=rdP2.yF;
	luP2p->xF=luP2.xF;
	luP2p->yF=luP2.yF;
	
	rectV.push_back(ldP2p);
	rectV.push_back(rdP2p);
	rectV.push_back(ruP2p);
	rectV.push_back(luP2p);
}


