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

						token=jumpNTimes(token,2);
						string dirS(token);
						if(!isMacro){
							cell *newCell=copyCell(lef.cellsMSCp[cellNameS]);
							comp *newComp=new comp(compNameS,dirS,locP2,newCell);
							compsMFVCp[yF].push_back(newComp);
							compsMSCp.insert(make_pair(compNameS,newComp));
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
						token=strtok(NULL," ");
						tmp_str.assign(token);
						string nameS=tmp_str;
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
							while(1){
								ifs.getline(line,1024);
								string testline=line;
								token=strtok(line," ");
								tmp_str.assign(token);
								if(tmp_str==";"){
									netsMSNp.insert(make_pair(nameS,newNet));
									break;
								}
								else{
									while(1){
										if (token==NULL){
											break;
										}
										else{
											tmp_str.assign(token);
											if(tmp_str=="("){
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
												}
											}
										}
										token=strtok(NULL," ");
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


