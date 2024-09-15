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
#include "drv.h"

using namespace std;
Drv::Drv(){
}

Drv::Drv(string _drvDirS){
	drvDirS=_drvDirS;
}
Drv::~Drv(){
}

void Drv::parse(){

	ifstream ifs;
	ifs.open(drvDirS);

	unordered_map<string,pin*> pinsMSPp;
	
	char delimit[]=" \t";
	bool is_starting_writing=false;
	if(ifs.is_open()){
		char lineC[1024];
		while(!ifs.eof()){
			ifs.getline(lineC,1024);
			char *token=strtok(lineC,delimit);
			while (token != NULL){
				string tmp_str(token);
				if (!is_starting_writing && tmp_str.find("Num",0)==0){
					is_starting_writing=true;
					break;
				}
				else if(!is_starting_writing){
					break;
				}
				else if(is_starting_writing){
					token=strtok(NULL,delimit);
					string tmp_str(token);
					int idI=stoi(tmp_str);

					token=strtok(NULL,delimit);
					tmp_str.assign(token);
					string typeS="";
					while (tmp_str.find("M",0)!=0 and tmp_str.find("V",0)!=0){
						typeS+=tmp_str+" ";
						token=strtok(NULL,delimit);
						tmp_str.assign(token);
					}
					typeS=typeS.substr(0,typeS.size()-1);
					
					string layerS=tmp_str;

					token=strtok(NULL,delimit); //(50)
					token=strtok(NULL,delimit); //Net:
					tmp_str.assign(token);

					vector<string> relatedNetsVS;
					if(tmp_str=="Net:"){
						token=strtok(NULL,delimit); 
						tmp_str.assign(token);
						relatedNetsVS.push_back(tmp_str);
						token=strtok(NULL,delimit); 
						tmp_str.assign(token); //Layer
					}
					else{
						while(tmp_str!="Layer:"){
							if(tmp_str.find("Net",0)!=0){
								relatedNetsVS.push_back(tmp_str);
							}
							token=strtok(NULL,delimit); 
							tmp_str.assign(token);
						}
					}

					token=strtok(NULL,delimit);//M1
					token=strtok(NULL,delimit);//(50)
					
					token=strtok(NULL,delimit);//(24.3130
					tmp_str.assign(token); 
					double lx=stod(tmp_str.substr(1,tmp_str.size()-1));
					token=strtok(NULL,delimit);//24.3130)
					tmp_str.assign(token); 
					double ly=stod(tmp_str.substr(0,tmp_str.size()-1));

					token=strtok(NULL,delimit);// (24.3150
					tmp_str.assign(token); 
					double rx=stod(tmp_str.substr(1,tmp_str.size()-1));
					token=strtok(NULL,delimit);// 24.3150)
					tmp_str.assign(token); 
					double ry=stod(tmp_str.substr(0,tmp_str.size()-1));
					BBox<double> bbox(lx,ly,rx,ry);
					drv new_drv(idI,typeS,layerS,relatedNetsVS,bbox);
					drvV.push_back(new_drv);
					break;

				}

				token=strtok(NULL,delimit);
			}
		}
	}
	ifs.close();
}



