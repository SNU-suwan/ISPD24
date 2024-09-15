#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "global_variables.h"
#include "default_classes.h"
#include "structures.h"
#include "def.h"
#include "utils.h"
#include "../global.h"

using namespace std;

string bottom_most_routing_layer;
string upper_most_routing_layer;
//map<string,int> METAL_LAYER_IDX;

Def::Def(string _defDirS, Lef& _lef,bool _is_cadence):lef{_lef}{
	IOpinsMSpp.clear();
	compsMSCp.clear();
	netsMSNp.clear();
	compsMFVCp.clear();
	comp_macrosMSCp.clear();
	rowsVp.clear();

	defDirS=_defDirS;
	//lef=_lef;
	is_cadence=_is_cadence;

}

Def::~Def(){
	for(size_t it=0;it<rowsVp.size();it++){
		delete rowsVp[it];
	}
	for(auto it=IOpinsMSpp.begin();it!=IOpinsMSpp.end();it++){
		delete it->second;
	}
	// Problem point //
	for(auto it=compsMSCp.begin();it!=compsMSCp.end();it++){
		delete it->second;
	}
	/////////////////////
	for(auto it=netsMSNp.begin();it!=netsMSNp.end();it++){
		delete it->second;
	}
	rowsVp.clear();
	IOpinsMSpp.clear();
	compsMSCp.clear();
	netsMSNp.clear();
	compsMFVCp.clear();

}
void Def::parse(){
	ifstream ifs;
	ifs.open(defDirS);

	bool isWriteComps=false;
	bool isWritePins=false;
	bool isWriteNets=false;
	//bool isWriteSITE=false;
	int upper_most_idx=-10000;
	int bottom_most_idx=10000;
	char line[1024];
	map<int,string> metal_idx_layer=build_metal_idx_layer();

	if(ifs.is_open()){
		while(!ifs.eof()){
			ifs.getline(line,1024);
			char *token=strtok(line," ");
			while(token !=NULL){
				string tmp_str(token);

				if (tmp_str=="UNITS"){
					token=jumpNTimes(token,3);
					scaleI=(int)atof(token);
				}
				else if (tmp_str=="VERSION"){
					token=strtok(NULL," ");
					versionF=atof(token);
					//assert (("Only 5.8 version of DEF file is compatible.", versionF=5.8));
				}
				else if(tmp_str=="DIEAREA"){
					if(!is_cadence){
						token=jumpNTimes(token,2);
						float die_xlb=atof(token)/scaleI;
						token=jumpNTimes(token,1);
						float die_ylb=atof(token)/scaleI;
						token=jumpNTimes(token,4);
						float die_yrt=atof(token)/scaleI;
						token=jumpNTimes(token,3);
						float die_xrt=atof(token)/scaleI;
						point2F die_lb(die_xlb,die_ylb);
						point2F die_rt(die_xrt,die_yrt);
						die_sizeP22[0]=die_lb;
						die_sizeP22[1]=die_rt;
					}
					else{
						token=jumpNTimes(token,2);
						float die_xlb=atof(token)/scaleI;
						token=jumpNTimes(token,1);
						float die_ylb=atof(token)/scaleI;
						token=jumpNTimes(token,3);
						float die_xrt=atof(token)/scaleI;
						token=jumpNTimes(token,1);
						float die_yrt=atof(token)/scaleI;
						point2F die_lb(die_xlb,die_ylb);
						point2F die_rt(die_xrt,die_yrt);
						die_sizeP22[0]=die_lb;
						die_sizeP22[1]=die_rt;

					}
				}
				else if(tmp_str=="ROW"){
					token=jumpNTimes(token,3);
					float rowXF=atof(token)/scaleI;

					token=jumpNTimes(token,1);
					float rowYF=atof(token)/scaleI;

					row *newRow=new row(rowYF,rowXF);
					rowsVp.push_back(newRow);
					break;
				}

				else if(tmp_str=="COMPONENTS"){
					isWriteComps=true;
				}
				else if(tmp_str=="PINS"){
					isWritePins=true;
				}
				else if(tmp_str=="NETS"){
					isWriteNets=true;
				}
				/*
				else if(tmp_str=="SITE"){
					isWriteSITE=true;
				}
				*/
				else if(tmp_str=="-"){
					if (isWriteComps){
						token=strtok(NULL," ");
						string compNameS(token);
						
						token=strtok(NULL," ");
						string cellNameS(token);

						bool isMacro=false;
						for(auto it_macro=lef.macrosMSMp.begin();it_macro!=lef.macrosMSMp.end();it_macro++){
							if (cellNameS==it_macro->first){
								isMacro=true;
								break;
							}
						}

						token=strtok(NULL," ");
						string tmp(token);
						while (tmp!="PLACED" and tmp!="FIXED"){
							token=strtok(NULL," ");
							tmp.assign(token);
						}
				
						token=jumpNTimes(token,2);
						float xF=atof(token)/scaleI;
						token=strtok(NULL," ");
						float yF=atof(token)/scaleI;
						point2F locP2={xF,yF};
						//cout << cellNameS << endl;
						
						token=jumpNTimes(token,2);
						string dirS(token);
						if(!isMacro){
							cell *newCell=copyCell(lef.cellsMSCp[cellNameS]);
							comp *newComp=new comp(compNameS,dirS,locP2,newCell);
							compsMFVCp[yF].push_back(newComp);
							compsMSCp.insert(make_pair(compNameS,newComp));
							/*if (compNameS == "P1/P1/State_reg\\[1\\]"){
								cout << compNameS << endl;
								auto pin_list = compsMSCp[compNameS]->cellCp->pinsMSPp;
								for (auto pin_info : pin_list){
									auto pin_name = pin_info.first;
									auto single_pin = compsMSCp[compNameS]->cellCp->pinsMSPp[pin_name];

	
									cout << pin_name << " ";
									int pin_cx = std::round(single_pin->centerP2.xF * 1000);
									cout << pin_cx << " ";
            						int pin_cy = std::round(single_pin->centerP2.yF * 1000);
									cout << pin_cy << endl;
								}
							}*/
						}
						else{
							macro *newMacro=copyMacro(lef.macrosMSMp[cellNameS]);
							comp_macro *newCompMacro=new comp_macro(compNameS,dirS,locP2,newMacro);
							comp_macrosMSCp.insert(make_pair(compNameS,newCompMacro));

						}
					}
					else if (isWritePins){
						token=strtok(NULL, " ");
						tmp_str.assign(token);
						if(tmp_str=="VDD" || tmp_str=="VSS"){
							break;
						}
						string nameS=tmp_str;
						
						ifs.getline(line,1024);
						token=strtok(line," ");
						token=jumpNTimes(token,2);
						string layerS(token);

						token=jumpNTimes(token,2);
						float xLDF=atof(token)/scaleI;
						token=strtok(NULL," ");
						float yLDF=atof(token)/scaleI;

						token=jumpNTimes(token,3);
						float xRUF=atof(token)/scaleI;
						token=strtok(NULL," ");
						float yRUF=atof(token)/scaleI;

						ifs.getline(line,1024);
						token=strtok(line," ");
						token=jumpNTimes(token,3);
						tmp_str.assign(token);
						float xF=stof(tmp_str)/scaleI;

						token=strtok(NULL," ");
						tmp_str.assign(token);
						float yF=stof(tmp_str)/scaleI;
						
						point2F centerP2={xF,yF};
						point2F sizeP22[2];
						sizeP22[0]={xLDF,yLDF};
						sizeP22[1]={xRUF,yRUF};

						IOpin* newIOPin=new IOpin(nameS,layerS,centerP2,sizeP22);
						IOpinsMSpp.insert(make_pair(nameS,newIOPin));
					}
					else if (isWriteNets){
						vector<tuple<point2I,point2I,int,bool>> metal_vec;
						token=strtok(NULL," ");
						tmp_str.assign(token);
						string nameS=tmp_str;
						bool is_print = false;
						net* newNet=new net(nameS);

						if(!is_cadence){
							bool isReadNetsDone=false;
							bool isReadSegment=false;
							bool isNetContainsPin=false;
							bool isClockUnits=false;
							while(!isReadNetsDone){
								ifs.getline(line,4096);
								token=strtok(line," ");
								if (token==NULL){
									continue;
								}
								else{
									tmp_str.assign(token);
									if(tmp_str=="+"){
																				
										isReadSegment=true;
										token=strtok(NULL," ");
										tmp_str.assign(token);
										if(tmp_str=="USE"){
											netsMSNp.insert(make_pair(nameS,newNet));
											isReadNetsDone=true;
											isReadSegment=false;
										}
									}
									else if(tmp_str=="("){
										token=strtok(NULL," ");
										string cellNameS(token);
										token=strtok(NULL," ");
										string pinNameS(token);
										newNet->cellPinMSS.insert(make_pair(cellNameS,pinNameS));
										if(cellNameS!="PIN"){

											bool isMacro=false;
											for(auto it_macro=lef.macrosMSMp.begin();it_macro!=lef.macrosMSMp.end();it_macro++){
												if (cellNameS==it_macro->first){
													isMacro=true;
													break;
												}
											}
											if(!isMacro){
												compsMSCp[cellNameS]->cellCp->pinsMSPp[pinNameS]->netp=newNet;
												for(size_t it_B=0;it_B<CLOCK_UNITS.size();it_B++){
													if(compsMSCp[cellNameS]->cellCp->nameS.rfind(CLOCK_UNITS[it_B],0)==0){
														isClockUnits=true;
														break;
													}
												}
											}
											else{
												comp_macrosMSCp[cellNameS]->macroMp->pinsMSPp[pinNameS]->netp=newNet;
											}
										}
										else{
											IOpinsMSpp[pinNameS]->netp=newNet;
											isNetContainsPin=true;
										}
									}
									//if(isReadSegment && !isNetContainsPin && !isClockUnits){
									if(isReadSegment && !isNetContainsPin){
										token=strtok(NULL," ");
										tmp_str.assign(token);

										int cur_metal_idx=METAL_LAYER_IDX[tmp_str];

										token=jumpNTimes(token,5);
										if(token!=NULL){
											tmp_str.assign(token);
											if(tmp_str=="("){
												if(cur_metal_idx<bottom_most_idx){
													bottom_most_idx=cur_metal_idx;
													bottom_most_routing_layer=metal_idx_layer[bottom_most_idx];

												}
												//bottom_most_routing_layer=METAL2;
												/*
												if (bottom_most_routing_layer==METAL1){
													cout<<"here"<<endl;
												}
												*/
												
												if(cur_metal_idx>upper_most_idx){
													upper_most_idx=cur_metal_idx;
													upper_most_routing_layer=metal_idx_layer[upper_most_idx];
													
													
													/////////////////// delete here ///////////////////
													if (upper_most_idx>3 && upper_most_idx<5){
														upper_most_routing_layer=metal_idx_layer[3];
														//cout<<"here1 "<<nameS<<endl;
													}
													else if (upper_most_idx>5){
														upper_most_routing_layer=metal_idx_layer[5];
														//cout<<"here2 "<<nameS<<endl;
													}

													//////////////////////////////////////////////////
													
												}
											}
										}
									}
								}
							}
						}
						else{
							bool isReadNetsDone=false;
							bool isReadSegment=false;
							bool isNetContainsPin=false;
							bool isClockUnits=false;
							while(!isReadNetsDone){
								ifs.getline(line,4096);
								token=strtok(line," ");

								if (token==NULL){
									continue;
								}
								else{

									tmp_str.assign(token);
									if(tmp_str==";"){
										//netsMSNp.insert(make_pair(nameS,newNet));
										//cout << nameS << endl;
										if (is_print) {
											for (int i = 0; i < metal_vec.size(); i++){
												point2I p1 = get<0>(metal_vec[i]);
												point2I p2 = get<1>(metal_vec[i]);
												int x1 = p1.xI, y1 = p1.yI, x2 = p2.xI, y2 = p2.yI;
												int metal_len = get<2>(metal_vec[i]);
												bool is_vertical = get<3>(metal_vec[i]);
												
												BBox<int> metal_bbox1, metal_bbox2;
												int lx, ly, rx, ry;
												if (is_vertical){
													lx = x1 - 144;
													rx = x2 + 144;
													ly = y1 - 224;
													ry = y2 + 224;

													metal_bbox1.lb.x=lx;metal_bbox1.rt.x=rx;
													metal_bbox1.lb.y=ly;metal_bbox1.rt.y=ly+224;

													metal_bbox2.lb.x=lx;metal_bbox2.rt.x=rx;
													metal_bbox2.lb.y=ry-224;metal_bbox2.rt.y=ry;
													metal_len = metal_bbox2.rt.y - metal_bbox1.lb.y;
												}
												else{
													lx = x1 - 224;
													rx = x2 + 224;
													ly = y1 - 144;
													ry = y2 + 144;

													metal_bbox1.lb.x=lx;metal_bbox1.rt.x=lx+224;
													metal_bbox1.lb.y=ly;metal_bbox1.rt.y=ry;

													metal_bbox2.lb.x=rx-224;metal_bbox2.rt.x=rx;
													metal_bbox2.lb.y=ly;metal_bbox2.rt.y=ry;
													metal_len = metal_bbox2.rt.x - metal_bbox1.lb.x;
												}
												auto& pin_list = netsMSNp[nameS]->cellPinMSS;
												for(auto& pin_info : pin_list){
													string comp_name = pin_info.first;
													string pin_name = pin_info.second;
													if (comp_name == "PIN") continue;
													
													auto& single_pin = compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name];
													if (single_pin->net_direction != "default") continue;
													//cout << comp_name << " " << pin_name << " ";
													for(auto& seg : single_pin->segmentsVp){
														int seg_lx = (int)((seg->pointsVp)[0]->xF*4000);
														int seg_ly = (int)((seg->pointsVp)[0]->yF*4000);
														int seg_rx = (int)((seg->pointsVp)[2]->xF*4000);
														int seg_ry = (int)((seg->pointsVp)[2]->yF*4000);
														BBox<int> seg_bbox(seg_lx,seg_ly,seg_rx,seg_ry);
														if (metal_bbox1.overlap(seg_bbox) && ((is_vertical && seg->layerS == "M2") || (!is_vertical && seg->layerS == "M1"))){
															//cout << metal_bbox1 << endl;
															if (metal_len > 1200) {
																//if (is_vertical && seg->layerS == "M2") single_pin->net_direction = "UP";
																//else if (!is_vertical && seg->layerS == "M1") single_pin->net_direction = "RIGHT";
																if (is_vertical) single_pin->net_direction = "UP";
																else if (!is_vertical) single_pin->net_direction = "RIGHT";

																break;
															}
															else{
																point2I target_p1 = p1, target_p2 = p2;
																point2I before_p1(-1,-1), before_p2(-1,-1);
																bool is_iteration = true;
																while(is_iteration){
																	for (int j = 0; j < metal_vec.size(); j++){
																		point2I p1_2 = get<0>(metal_vec[j]);
																		point2I p2_2 = get<1>(metal_vec[j]);
																		//cout << "(" << p1_2.xI << "," << p1_2.yI << ") (" << p2_2.xI << "," << p2_2.yI << ") " << metal_vec.size() << endl;
																		int metal_len_2 = get<2>(metal_vec[j]);
																		bool is_vertical_2 = get<3>(metal_vec[j]);
																		if (i == j) continue;
																		if (p1_2 == before_p1 && p2_2 == before_p2) continue;
																		else if (p1_2 == target_p1 || p1_2 == target_p2){
																			if (metal_len_2 < 1200) {before_p1 = target_p1; before_p2 = target_p2; target_p1 = p1_2; target_p2 = p2_2;}
																			else{
																				if (is_vertical_2) single_pin->net_direction = "UP";
																				else single_pin->net_direction = "RIGHT";
																				is_iteration = false;
																				break;
																			}
																		}
																		else if (p2_2 == target_p1 || p2_2 == target_p2){
																			if (metal_len_2 < 1200) {before_p1 = target_p1; before_p2 = target_p2; target_p1 = p1_2; target_p2 = p2_2;}
																			else{
																				if (is_vertical_2) single_pin->net_direction = "DOWN";
																				else single_pin->net_direction = "LEFT";
																				is_iteration = false;
																				break;
																			}
																		}
																	}
																}
															}
														}
														else if (metal_bbox2.overlap(seg_bbox) && ((is_vertical && seg->layerS == "M2") || (!is_vertical && seg->layerS == "M1"))){
															//cout << metal_bbox2 << endl;
															if (metal_len > 1200) {
																if (is_vertical) single_pin->net_direction = "DOWN";
																else if (!is_vertical) single_pin->net_direction = "LEFT";
																break;
															}
															point2I target_p1 = p1, target_p2 = p2;
															point2I before_p1(-1,-1), before_p2(-1,-1);
															bool is_iteration = true;
															while(is_iteration){
																for (int j = 0; j < metal_vec.size(); j++){
																	point2I p1_2 = get<0>(metal_vec[j]);
																	point2I p2_2 = get<1>(metal_vec[j]);
																	//cout << "(" << p1_2.xI << "," << p1_2.yI << ") (" << p2_2.xI << "," << p2_2.yI << ") " << " ";
																	//cout << "(" << target_p1.xI << "," << target_p1.yI << ") (" << target_p2.xI << "," << target_p2.yI << ") " << endl;
																	int metal_len_2 = get<2>(metal_vec[j]);
																	bool is_vertical_2 = get<3>(metal_vec[j]);
																	if (i == j) continue;
																	if (p1_2 == before_p1 && p2_2 == before_p2) continue;
																	else if (p1_2 == target_p1 || p1_2 == target_p2){
																		if (metal_len_2 < 1200) {before_p1 = target_p1; before_p2 = target_p2; target_p1 = p1_2; target_p2 = p2_2;}
																		else{
																			if (is_vertical_2) single_pin->net_direction = "UP";
																			else single_pin->net_direction = "RIGHT";
																			is_iteration = false;
																			break;
																		}
																	}
																	else if (p2_2 == target_p1 || p2_2 == target_p2){
																		if (metal_len_2 < 1200) {before_p1 = target_p1; before_p2 = target_p2; target_p1 = p1_2; target_p2 = p2_2;}
																		else{
																			if (is_vertical_2) single_pin->net_direction = "DOWN";
																			else single_pin->net_direction = "LEFT";
																			is_iteration = false;
																			break;
																		}
																	}
																}
															}
														}
														//BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
														//cout << seg_bbox << endl;
													}
													//cout << single_pin->net_direction << endl;
												}
											}
										}
										netsMSNp.insert(make_pair(nameS,newNet));
										//if (nameS != "reset") isReadSegment=true;
										//token=strtok(NULL," ");
										//tmp_str.assign(token);

										metal_vec.clear();
										isReadNetsDone=true;
										isReadSegment=false;
										break;
									}/*
									else if(tmp_str=="+"){
										netsMSNp.insert(make_pair(nameS,newNet));
										if (nameS != "reset") isReadSegment=true;
										token=strtok(NULL," ");
										tmp_str.assign(token);
										*/
										/*if(tmp_str=="USE" | tmp_str=="SOURCE" | tmp_str==";"){
											netsMSNp.insert(make_pair(nameS,newNet));
											isReadNetsDone=true;
											isReadSegment=false;
										}*/
									//}
									else if(tmp_str=="("){
										while(token){
											token=strtok(NULL," ");
											string cellNameS(token);
											token=strtok(NULL," ");
										
											string pinNameS(token);
											newNet->cellPinMSS.insert(make_pair(cellNameS,pinNameS));

											if(cellNameS!="PIN"){
												compsMSCp[cellNameS]->cellCp->pinsMSPp[pinNameS]->netp=newNet;
											}
											else{
												IOpinsMSpp[pinNameS]->netp=newNet;
												isNetContainsPin=true;
											}
											token=jumpNTimes(token,2);
										}
									}
									if(isReadSegment){
										if(!is_print) continue;
										int x1,x2,y1,y2;
										int lx,ly,rx,ry;
										int metal_len;
										bool store_direction=false;
										bool is_vertical;
										token=strtok(NULL," ");
										while(token){
											tmp_str.assign(token);	
											if (tmp_str == "("){
												token = strtok(NULL," ");
												x1 = (int)atof(token);
												token = strtok(NULL," ");
												y1 = (int)atof(token);
												
												token = strtok(NULL," ");
												tmp_str.assign(token);
												if (tmp_str == "0") token=jumpNTimes(token,2);
												else token = strtok(NULL," ");
												tmp_str.assign(token);
												//token=jumpNTimes(token,2);
												//tmp_str.assign(token);

												if (tmp_str == "("){
													token = strtok(NULL," ");
													tmp_str.assign(token);
													if (tmp_str == "*") x2 = x1;
													else x2 = (int)atof(token);
													token = strtok(NULL," ");
													tmp_str.assign(token);
													if (tmp_str == "*") y2 = y1;
													else y2 = (int)atof(token);
													store_direction=true;

													if (x1 == x2){
														lx = x1 - 144;
														rx = x2 + 144;
														ly = y1 - 224;
														ry = y2 + 224;
														metal_len = ry - ly;
														is_vertical = true;
														metal_vec.push_back(make_tuple(point2I(x1,y1),point2I(x2,y2),metal_len,is_vertical));
													}
													if (y1 == y2){
														lx = x1 - 224;
														rx = x2 + 224;
														ly = y1 - 144;
														ry = y2 + 144;
														metal_len = rx - lx;
														is_vertical = false;
														metal_vec.push_back(make_tuple(point2I(x1,y1),point2I(x2,y2),metal_len,is_vertical));
													}
												}
												else break;
											}
											else token = strtok(NULL," ");
										}/*
										if (store_direction){
											BBox<int> metal_bbox1, metal_bbox2;
											int metal_len;
											if (is_vertical){
												metal_bbox1.lb.x=lx;metal_bbox1.rt.x=rx;
												metal_bbox1.lb.y=ly;metal_bbox1.rt.y=ly+224;

												metal_bbox2.lb.x=lx;metal_bbox2.rt.x=rx;
												metal_bbox2.lb.y=ry-224;metal_bbox2.rt.y=ry;
												metal_len = metal_bbox2.rt.y - metal_bbox1.lb.y;
											}
											else{
												metal_bbox1.lb.x=lx;metal_bbox1.rt.x=lx+224;
												metal_bbox1.lb.y=ly;metal_bbox1.rt.y=ry;

												metal_bbox2.lb.x=rx-224;metal_bbox2.rt.x=rx;
												metal_bbox2.lb.y=ly;metal_bbox2.rt.y=ry;
												metal_len = metal_bbox2.rt.x - metal_bbox1.lb.x;
											}
											auto& pin_list = netsMSNp[nameS]->cellPinMSS;
											for(auto& pin_info : pin_list){
												string comp_name = pin_info.first;
												string pin_name = pin_info.second;
												//cout << comp_name << " " << pin_name << endl;
												if (comp_name == "PIN") continue;
												auto& single_pin = compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name];
												for(auto& seg : single_pin->segmentsVp){
													int seg_lx = (int)((seg->pointsVp)[0]->xF*4000);
													int seg_ly = (int)((seg->pointsVp)[0]->yF*4000);
													int seg_rx = (int)((seg->pointsVp)[2]->xF*4000);
													int seg_ry = (int)((seg->pointsVp)[2]->yF*4000);
													BBox<int> seg_bbox(seg_lx,seg_ly,seg_rx,seg_ry);
													if (metal_bbox1.overlap(seg_bbox)){
														//cout << comp_name << " " << pin_name << " ";
														if (is_vertical && seg->layerS == "M2") single_pin->net_direction = "UP";
														else if (!is_vertical && seg->layerS == "M1") single_pin->net_direction = "RIGHT";
														if (metal_len < 1200) cout << single_pin->net_direction << endl;
														break;
													}
													else if (metal_bbox2.overlap(seg_bbox)){
														//cout << comp_name << " " << pin_name << " ";
														if (is_vertical && seg->layerS == "M2") single_pin->net_direction = "DOWN";
														else if (!is_vertical && seg->layerS == "M1") single_pin->net_direction = "LEFT";
														if (metal_len < 1200) cout << single_pin->net_direction << endl;
														break;
													}
													//BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
													//cout << seg_bbox << endl;
												}
											}
										}*/
										//net segment읽어서 direction 나타내기
										//if(!is_print) continue;
										/*int x1,x2,y1,y2;
										int lx,ly,rx,ry;
										
										bool store_direction=false;
										bool is_vertical;
										token=strtok(NULL," ");
										while(token){
											tmp_str.assign(token);	
											if (tmp_str == "("){
												token = strtok(NULL," ");
												x1 = (int)atof(token);
												token = strtok(NULL," ");
												y1 = (int)atof(token);
												token=jumpNTimes(token,2);
												tmp_str.assign(token);
												if (tmp_str == "("){
													token = strtok(NULL," ");
													tmp_str.assign(token);
													if (tmp_str == "*") x2 = x1;
													else x2 = (int)atof(token);
													token = strtok(NULL," ");
													tmp_str.assign(token);
													if (tmp_str == "*") y2 = y1;
													else y2 = (int)atof(token);
													store_direction=true;

													if (x1 == x2){
														lx = x1 - 144;
														rx = x2 + 144;
														ly = y1 - 224;
														ry = y2 + 224;
														is_vertical = true;
													}
													if (y1 == y2){
														lx = x1 - 224;
														rx = x2 + 224;
														ly = y1 - 144;
														ry = y2 + 144;
														is_vertical = false;
													}
												}
												else break;
											}
											else token = strtok(NULL," ");
										}
										if (store_direction){
											BBox<int> metal_bbox1, metal_bbox2;
											int metal_len;
											if (is_vertical){
												metal_bbox1.lb.x=lx;metal_bbox1.rt.x=rx;
												metal_bbox1.lb.y=ly;metal_bbox1.rt.y=ly+224;

												metal_bbox2.lb.x=lx;metal_bbox2.rt.x=rx;
												metal_bbox2.lb.y=ry-224;metal_bbox2.rt.y=ry;
												metal_len = metal_bbox2.rt.y - metal_bbox1.lb.y;
											}
											else{
												metal_bbox1.lb.x=lx;metal_bbox1.rt.x=lx+224;
												metal_bbox1.lb.y=ly;metal_bbox1.rt.y=ry;

												metal_bbox2.lb.x=rx-224;metal_bbox2.rt.x=rx;
												metal_bbox2.lb.y=ly;metal_bbox2.rt.y=ry;
												metal_len = metal_bbox2.rt.x - metal_bbox1.lb.x;
											}
											auto& pin_list = netsMSNp[nameS]->cellPinMSS;
											for(auto& pin_info : pin_list){
												string comp_name = pin_info.first;
												string pin_name = pin_info.second;
												//cout << comp_name << " " << pin_name << endl;
												if (comp_name == "PIN") continue;
												auto& single_pin = compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name];
												for(auto& seg : single_pin->segmentsVp){
													int seg_lx = (int)((seg->pointsVp)[0]->xF*4000);
													int seg_ly = (int)((seg->pointsVp)[0]->yF*4000);
													int seg_rx = (int)((seg->pointsVp)[2]->xF*4000);
													int seg_ry = (int)((seg->pointsVp)[2]->yF*4000);
													BBox<int> seg_bbox(seg_lx,seg_ly,seg_rx,seg_ry);
													if (metal_bbox1.overlap(seg_bbox)){
														//cout << comp_name << " " << pin_name << " ";
														if (is_vertical && seg->layerS == "M2") single_pin->net_direction = "UP";
														else if (!is_vertical && seg->layerS == "M1") single_pin->net_direction = "RIGHT";
														if (metal_len < 1200) cout << single_pin->net_direction << endl;
														break;
													}
													else if (metal_bbox2.overlap(seg_bbox)){
														//cout << comp_name << " " << pin_name << " ";
														if (is_vertical && seg->layerS == "M2") single_pin->net_direction = "DOWN";
														else if (!is_vertical && seg->layerS == "M1") single_pin->net_direction = "LEFT";
														if (metal_len < 1200) cout << single_pin->net_direction << endl;
														break;
													}
													//BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
													//cout << seg_bbox << endl;
												}
											}
										}*/
									}
								}
							}
						}
					}
				}
				/*
				else if(isWriteSITE){
					tmp_str.assign(token);
					if (tmp_str=="SIZE"){
						token=jumpNTimes(token,3);
						tmp_str.assign(token);
						SITE_Y=stod(tmp_str);
					}
					
				}
				*/
				else if(tmp_str=="END"){
					token=strtok(NULL," ");
					if (token!=NULL){
						tmp_str.assign(token);
						if (tmp_str=="COMPONENTS") isWriteComps=false;
						else if (tmp_str=="PINS") isWritePins=false;
						else if (tmp_str=="NETS") isWriteNets=false;
					}
				}
				token=strtok(NULL," ");
			}
		}
	}
	ifs.close();
	
	//Sort comps according to its x coordinate.
	map<float,vector<comp*>>::iterator it;
	for(it=compsMFVCp.begin();it!=compsMFVCp.end();it++){
		vector<comp*>& compVCp=it->second;
		sort(compVCp.begin(),compVCp.end(),
				[](comp* lcomp, comp* rcomp) ->bool{
					return lcomp->locP2.xF < rcomp->locP2.xF;
				});
	}
	//Sort rows
	sort(rowsVp.begin(),rowsVp.end(),	
			[](row* lrow, row* rrow) -> bool{
				return lrow->rowHeightF < rrow->rowHeightF;
			});
		
}

char * jumpNTimes(char *token, int n){
	char * jumpedToken;
	for (int i=0;i<n;i++){
		token=strtok(NULL," ");
	}
	jumpedToken=token;
	return jumpedToken;
}

void Def::write_def(string outDefDirS){
	ifstream ifs;
	ifs.open(defDirS);
	ofstream ofs;
	ofs.open(outDefDirS);
	
	char line[1024];
	string str_line;
	bool is_write_comps=false;

	if(ifs.is_open()){
		while (!ifs.eof()){
			ifs.getline(line,1024);
			string str_line=line;
			char *token=strtok(line," ");
			if(!is_write_comps){
				ofs<<str_line;
				ofs<<endl;
			}
			while(token !=NULL){
				string tmp_str(token);
				if(tmp_str!="COMPONENTS" && tmp_str!="END"){
					break;
				}
				else if(tmp_str=="COMPONENTS"){
					is_write_comps=true;
				}
				else if(tmp_str=="END"){
					token=strtok(NULL," ");
					tmp_str.assign(token);
					if(tmp_str=="COMPONENTS"){
						ofs<<str_line;
						ofs<<endl;
						is_write_comps=false;
					}
					else{
						break;
					}
				}
				if(is_write_comps){
					for(auto it=compsMFVCp.begin();it!=compsMFVCp.end();it++){
						vector<comp*>&compsVCp=it->second;
						for(size_t i=0;i<compsVCp.size();i++){
							comp* compp=compsVCp[i];
							string write_comp;

							if(is_cadence){
								int scaled_x=compp->locP2.xF*1000;
								int scaled_y=compp->locP2.yF*1000;
								int unscaled_s=scaleI/1000;

								int comp_x=scaled_x*unscaled_s;
								int comp_y=scaled_y*unscaled_s;

								write_comp=" - "+compp->nameS+" "+compp->cellCp->nameS+" + PLACED ( "+to_string(comp_x)+" "+to_string(comp_y)+" ) "+compp->cellCp->dirS+" ;";
							}
							else{
							
								write_comp=" - "+compp->nameS+" "+compp->cellCp->nameS+" + PLACED ( "+to_string((int)(compp->locP2.xF*scaleI))+" "+to_string((int)(compp->locP2.yF*scaleI))+" ) "+compp->cellCp->dirS+" ;";
							}
							ofs<<write_comp;
							ofs<<endl;
						}
						
					}
				}
				token=strtok(NULL," ");
			}
		}
	}
	ifs.close();
	ofs.close();
}


