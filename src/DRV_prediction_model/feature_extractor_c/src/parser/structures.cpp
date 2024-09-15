#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "lef.h"
#include "def.h"
#include "default_classes.h"
#include "global_variables.h"
#include "structures.h"
#include "utils.h"

using namespace std;

instance::instance(string _nameS, point2F _sizeP2, unordered_map<string,pin*> _pinsMSPp){
	nameS=_nameS;
	dirS="N";
	sizeP2=_sizeP2;
	pinsMSPp=_pinsMSPp;
}
instance::~instance(){
	for (auto &it:pinsMSPp){
		pin *pinp=it.second;
		delete pinp;
	}
	pinsMSPp.clear();
}

void instance::setPos(const point2F& locP2,const point2F& originP2){
	//Move to locP2 from originP2
	map<string,pin*>::iterator it;
	//for(it=pinsMSPp.begin();it!=pinsMSPp.end();it++)
	for(auto &it:pinsMSPp){
		//pin *pinp=it->second;
		pin *pinp=it.second;
		pinp->setPos(locP2,originP2);
	}
}
void instance::setDir(const string& _dirS){
	if(dirS!=_dirS){
		if(dirS.length()==_dirS.length()){
			//N-S, FN-FS
			flipVer();
			flipHor();
		}
		else if(dirS.back()==_dirS.back()){
			//S-FS, N-FN
			flipHor();
		}
		else{
			//N-FS, S-FN
			flipVer();
		}
	}
	assert (dirS==_dirS);
}
void instance::flipVer(){
	if(dirS=="N")dirS="FS";
	else if(dirS=="S")dirS="FN";
	else if(dirS=="FN")dirS="S";
	else if(dirS=="FS")dirS="N";

	map<string,pin*>::iterator it;
	//for(it=pinsMSPp.begin();it!=pinsMSPp.end();it++)
	for(auto &it:pinsMSPp){
		//pin *pinp=it->second;
		pin *pinp=it.second;
		pinp->flipVer(sizeP2,*originP2p);
	}
}

void instance::flipHor(){
	if(dirS=="N")dirS="FN";
	else if(dirS=="S")dirS="FS";
	else if(dirS=="FN")dirS="N";
	else if(dirS=="FS")dirS="S";

	//map<string,pin*>::iterator it;
	//for(it=pinsMSPp.begin();it!=pinsMSPp.end();it++){
	for(auto &it:pinsMSPp){
		//pin *pinp=it->second;
		pin *pinp=it.second;
		pinp->flipHor(sizeP2,*originP2p);
	}
}

void instance::update_pin_center(){
	for(auto it=pinsMSPp.begin();it!=pinsMSPp.end();it++){
		pin*& pinp=it->second;
		vector<ap*>&apsVp=pinp->apsVp;
		point2F centerP2={0,0};
		for(size_t i=0;i<apsVp.size();i++){
			point2F &pointP2=apsVp[i]->pointP2;
			centerP2+=pointP2;
		}
		if(pinp->nameS!="SN" && pinp->nameS!="RN" && pinp->nameS!="SE"){assert(apsVp.size()!=0);}
		if(apsVp.size()!=0){
			centerP2/=apsVp.size();
		}
		pinp->centerP2=centerP2;
	}
}
void instance::print_pins(bool is_print_ap){
	for(auto it=pinsMSPp.begin();it!=pinsMSPp.end();it++){
		pin* pinp=it->second;
		cout<<"Pin name:"<<it->first<<" "<<endl;
		assert(it->first==pinp->nameS);
		cout<<"  Center: "<<pinp->centerP2.xF<<" "<<pinp->centerP2.yF<<endl;
		if(is_print_ap){
			vector<ap*> &apsVp=pinp->apsVp;
			cout<<"  AP info"<<endl;
			for(size_t i=0;i<apsVp.size();i++){
				ap* app=apsVp[i];
				cout<<"  AP"<<i<<endl;
				cout<<"    Track :"<<app->trackI<<endl;
				cout<<"    Point :"<<app->pointP2.xF<<" "<<app->pointP2.yF<<endl;
				cout<<"    Layer :"<<app->layerS<<endl;
			}
		}
	}
	cout<<endl;
}

//cell
cell::cell(string _nameS, point2F _sizeP2, unordered_map<string,pin*> _pinsMSPp):instance(_nameS,_sizeP2,_pinsMSPp){
	compp=NULL;
	for(auto& it:pinsMSPp){
		it.second->cellp=this;
	}
}

cell::~cell(){
}

//macro
macro::macro(string _nameS, point2F _sizeP2, unordered_map<string,pin*> _pinsMSPp):instance(_nameS,_sizeP2,_pinsMSPp){
	compp=NULL;
	for(auto& it:pinsMSPp){
		it.second->macrop=this;
	}
}
macro::~macro(){
}

cell* copyCell(cell* cellp){
	//map<string,pin*> pinsMSPp;
	unordered_map<string,pin*> pinsMSPp;
	for(auto it=cellp->pinsMSPp.begin();it!=cellp->pinsMSPp.end();it++){
		pin* pinp=it->second;
		pin* newPin=new pin(*pinp);
		
		vector<ap*> apsVp;
		for(size_t it=0;it<pinp->apsVp.size();it++){
			ap* newAp=new ap(*pinp->apsVp[it]);
			newAp->pinp=newPin;
			apsVp.push_back(newAp);
		}

		vector<segment*> segmentsVp;
		for(size_t it=0;it<pinp->segmentsVp.size();it++){
			vector<point2F*> pointsVp;
			for(size_t it2=0;it2<pinp->segmentsVp[it]->pointsVp.size();it2++){
				point2F *newPoint2F=new point2F(*pinp->segmentsVp[it]->pointsVp[it2]);
				pointsVp.push_back(newPoint2F);
			}
			segment* newSegment=new segment(*pinp->segmentsVp[it]);
			newSegment->pointsVp=pointsVp;
			segmentsVp.push_back(newSegment);
		}

		newPin->segmentsVp=segmentsVp;
		newPin->apsVp=apsVp;
		newPin->createApMap();
		//pinsMSPp[it->first]=newPin;
		pinsMSPp.insert(make_pair(it->first,newPin));
	}
	cell *newCell=new cell(*cellp);
	newCell->pinsMSPp=pinsMSPp;
	//for(auto it=newCell->pinsMSPp.begin();it!=newCell->pinsMSPp.end();it++)
	for(auto &it:newCell->pinsMSPp){
		//it->second->cellp=newCell;
		it.second->cellp=newCell;
	}
	return newCell;
}


macro* copyMacro(macro* macrop){
	unordered_map<string,pin*> pinsMSPp;
	for(auto it=macrop->pinsMSPp.begin();it!=macrop->pinsMSPp.end();it++){
		pin* pinp=it->second;
		pin* newPin=new pin(*pinp);
		
		vector<segment*> segmentsVp;
		for(size_t it=0;it<pinp->segmentsVp.size();it++){
			vector<point2F*> pointsVp;
			for(size_t it2=0;it2<pinp->segmentsVp[it]->pointsVp.size();it2++){
				point2F *newPoint2F=new point2F(*pinp->segmentsVp[it]->pointsVp[it2]);
				pointsVp.push_back(newPoint2F);
			}
			segment* newSegment=new segment(*pinp->segmentsVp[it]);
			newSegment->pointsVp=pointsVp;
			segmentsVp.push_back(newSegment);
		}

		newPin->segmentsVp=segmentsVp;
		pinsMSPp.insert(make_pair(it->first,newPin));
	}
	macro *newMacro=new macro(*macrop);
	newMacro->pinsMSPp=pinsMSPp;
	for(auto &it:newMacro->pinsMSPp){
		it.second->macrop=newMacro;
	}
	return newMacro;
}


//pin
pin::pin(string _nameS){
	nameS=_nameS;
	segmentsVp.clear();
	apsVp.clear();
	netp=NULL;
	cellp=NULL;
	macrop=NULL;
	centerP2={0,0};
	center_ingridP2={0,0};
	avg_ap_x_ingrid = 0;

	
	targetTypeI = -1;
	target = NULL;
	targetP2 = nullptr;
}
pin::~pin(){
	for(size_t i=0;i<segmentsVp.size();i++){
		delete segmentsVp[i];
	}
	for(size_t i=0;i<apsVp.size();i++){
		delete apsVp[i];
	}
	segmentsVp.clear();
	apsVp.clear();
}

void pin::createAps(){
	for(size_t it=0;it<segmentsVp.size();it++){
		bool direction=HORIZONTAL;
		bool prevState=direction;
		float minXF=MAXF;
		float maxXF=MINF;
		float minYF=MAXF;
		float maxYF=MINF;
		float prevXF;
		float prevYF;
		string layerS=segmentsVp[it]->layerS;
		vector<point2F*>& pointsVp=segmentsVp[it]->pointsVp;
		vector<edge> edgesV;

		//convert segments into collection of straight lines
		for (size_t it2=0;it2<pointsVp.size();it2++){
			point2F &currP2=*pointsVp[it2];
			point2F &nextP2=*pointsVp[(it2+1)%(pointsVp.size())];
			point2F &nnextP2=*pointsVp[(it2+2)%(pointsVp.size())];

			if(currP2.xF==nextP2.xF) direction=VERTICAL;
			else direction=HORIZONTAL;
			if(direction==HORIZONTAL){
				//horizontal
				if (currP2.xF<minXF) minXF=currP2.xF;
				if (currP2.xF>maxXF) maxXF=currP2.xF;
				if (nextP2.xF<minXF) minXF=nextP2.xF;
				if (nextP2.xF>maxXF) maxXF=nextP2.xF;
			}
			else{
				//vertical	
				if (currP2.yF<minYF) minYF=currP2.yF;
				if (currP2.yF>maxYF) maxYF=currP2.yF;
				if (nextP2.yF<minYF) minYF=nextP2.yF;
				if (nextP2.yF>maxYF) maxYF=nextP2.yF;
			}

			prevState=direction;
			prevXF=currP2.xF;
			prevYF=currP2.yF;
			
			//Save the current edge if the direction changes
			if(nextP2.xF==nnextP2.xF) direction=VERTICAL;
			else direction=HORIZONTAL;
			if(prevState!=direction){ // If the direction changes, finish.
				if(prevState==HORIZONTAL){
					//horizontal
					point2F xMinP2={minXF,prevYF};
					point2F xMaxP2={maxXF,prevYF};
					point2F xP22[2]={xMinP2, xMaxP2};
					edge tmpEdge(prevState,segmentsVp[it]->layerS,xP22);
					edgesV.push_back(tmpEdge);
				}
				else{
					//vertical	
					point2F yMinP2={prevXF,minYF};
					point2F yMaxP2={prevXF,maxYF};
					point2F yP22[2]={yMinP2, yMaxP2};
					edge tmpEdge(prevState,segmentsVp[it]->layerS,yP22);
					edgesV.push_back(tmpEdge);
				}
				minXF=MAXF;
				maxXF=MINF;
				minYF=MAXF;
				maxYF=MINF;
			}

		}

		//Divide edges into the collection of grids
		map<float,vector<point2F>> sameYMFVP2;
		map<float,vector<point2F>> sameXMFVP2;
		for(size_t it_e=0;it_e<edgesV.size();it_e++){
			edge &curEdge=edgesV[it_e];
			//Update center
			centerP2+=curEdge.endPointsP22[0];
			centerP2+=curEdge.endPointsP22[1];

			string layerS=curEdge.layerS;
			if (curEdge.direction==HORIZONTAL){
				point2F startP2=curEdge.endPointsP22[0];
				point2F endP2=curEdge.endPointsP22[1];
				assert (startP2.yF==endP2.yF);
				int endGridI=endP2.xF/M1_PITCH;
				int startGridI=startP2.xF/M1_PITCH;

				for(int i=startGridI+1;i<endGridI+1;i++){
					point2F tmpP2={(float)(i*M1_PITCH),startP2.yF};
					tmpP2.containerS=layerS;
					sameXMFVP2[(float)(i*M1_PITCH)].push_back(tmpP2);
				}
			}
			else{
				point2F startP2=curEdge.endPointsP22[0];
				point2F endP2=curEdge.endPointsP22[1];
				assert (startP2.xF==endP2.xF);
				
				int endGridI=(TECH==7)?((endP2.yF-M2_FIRSTLAST)/M2_PITCH):endP2.yF/M2_PITCH;
				int startGridI=(TECH==7)?((startP2.yF-M2_FIRSTLAST)/M2_PITCH):startP2.yF/M2_PITCH;
				startGridI=(TECH==7 && startP2.yF<M2_FIRSTLAST)?-1:startGridI;


				for(int i=startGridI+1;i<endGridI+1;i++){
					point2F tmpP2={startP2.xF,(float)(i*M2_PITCH)};
					if(TECH==7){
						tmpP2.yF+=M2_FIRSTLAST;
					}
					
					tmpP2.containerS=layerS;
					
					if(TECH==7){
						sameYMFVP2[(float)(i*M2_PITCH)+M2_FIRSTLAST].push_back(tmpP2);
					}
					else{
						sameYMFVP2[(float)(i*M2_PITCH)].push_back(tmpP2);
					}
				}
			}
		}
		centerP2/=(2*edgesV.size()); // because each vertex is considered twice.
	
		//Create lines according to grid information
		map<float,vector<point2F>>::iterator itm;
		vector<edge> horLinesVE;
		vector<edge> verLinesVE;
		for(itm=sameYMFVP2.begin();itm!=sameYMFVP2.end();itm++){
			vector<point2F> &pointVP2=itm->second;
			sort(pointVP2.begin(),pointVP2.end());
			for (size_t itv=0;itv<pointVP2.size();itv+=2){
				point2F currP2=pointVP2[itv];
				point2F nextP2=pointVP2[itv+1];
				point2F pointP22[2]={currP2,nextP2};
				edge tmpEdge(HORIZONTAL,currP2.containerS,pointP22);
				horLinesVE.push_back(tmpEdge);
			}
		}
		if(TECH==7) {
			for(size_t ith=0;ith<horLinesVE.size();ith++){
				point2F startHorP2=horLinesVE[ith].endPointsP22[0];
				point2F endHorP2=horLinesVE[ith].endPointsP22[1];
				
				float start_x=startHorP2.xF+M1_PITCH/2;
				while(start_x<endHorP2.xF){
					point2F intersectionP2={start_x,startHorP2.yF};
					ap *newAp=new ap(layerS,intersectionP2);
					newAp->pinp=this;
					apsVp.push_back(newAp);
					start_x+=M1_PITCH;
				}
			}
		}
		else{
			for(itm=sameXMFVP2.begin();itm!=sameXMFVP2.end();itm++){
				vector<point2F> &pointVP2=itm->second;
				sort(pointVP2.begin(),pointVP2.end());
				for (size_t itv=0;itv<pointVP2.size();itv+=2){
					point2F currP2=pointVP2[itv];
					point2F nextP2=pointVP2[itv+1];
					point2F pointP22[2]={currP2,nextP2};
					edge tmpEdge(VERTICAL,currP2.containerS,pointP22);
					verLinesVE.push_back(tmpEdge);
				}
			}
			//Find intersections
			for(size_t ith=0;ith<horLinesVE.size();ith++){
				point2F startHorP2=horLinesVE[ith].endPointsP22[0];
				point2F endHorP2=horLinesVE[ith].endPointsP22[1];
				if(verLinesVE.size()>0){
					for(size_t itv=0;itv<verLinesVE.size();itv++){
						point2F startVerP2=verLinesVE[itv].endPointsP22[0];
						point2F endVerP2=verLinesVE[itv].endPointsP22[1];
						if(startHorP2.containerS==startVerP2.containerS){ //Check same layer
							if(startHorP2.yF>startVerP2.yF && startHorP2.yF<endVerP2.yF){
								if(startVerP2.xF>startHorP2.xF && startVerP2.xF<endHorP2.xF){
									point2F intersectionP2={startVerP2.xF,startHorP2.yF};
									ap *newAp=new ap(layerS,intersectionP2);
									newAp->pinp=this;
									apsVp.push_back(newAp);
								}
							}
						}
					}
				}
				else{
					//If gridless, set the midpoint as an access point
					point2F intersectionP2={(startHorP2.xF+endHorP2.xF)/2,startHorP2.yF};
					ap *newAp=new ap(layerS,intersectionP2);
					newAp->pinp=this;
					apsVp.push_back(newAp);
				}
			}
		}
	}
	//destroy segments
	//destroySegments();
	
	//delete same ap
	vector<size_t> sameApIndexVI;
	for(size_t it=0;it<apsVp.size();it++){
		for(size_t it2=it+1;it2<apsVp.size();it2++){
			if(apsVp[it]->pointP2.equalY(apsVp[it2]->pointP2) && abs(apsVp[it]->pointP2.xF-apsVp[it2]->pointP2.xF)<0.05){
				sameApIndexVI.push_back((size_t)it);
			}
		}
	}
	for(int it=sameApIndexVI.size()-1;it>=0;it--){
		delete apsVp[sameApIndexVI[it]];
		apsVp.erase(apsVp.begin()+sameApIndexVI[it]);
	}

}

void pin::createApMap(){
	apsMIVAp.clear();
	int countI=0;
	float offset=(cellp==NULL)?0:(cellp->compp==NULL)?0:cellp->compp->locP2.yF;
	for(size_t it=0;it<apsVp.size();it++){
		ap* app=apsVp[it];
		float& ap_yF=app->pointP2.yF;
		for(int i=0;i<NUM_TRACK+1;i++){
			float cur_trackF=i*M2_PITCH;
			if(TECH==7){
				if (i!=0){
					cur_trackF+=M2_FIRSTLAST;
					cur_trackF-=M2_PITCH;
				}
			}
			if(abs(ap_yF-cur_trackF-offset)<0.01){
				apsMIVAp[i].push_back(app);
				app->trackI=i;
				countI++;
				break;
			}
		}
	}
	assert(countI==apsVp.size());
}

void pin::destroyAps(){
	for(size_t it=0;it<apsVp.size();it++){
		delete apsVp[it];
	}
	apsVp.clear();
}

void pin::destroySegments(){
	for(size_t it=0;it<segmentsVp.size();it++){
		delete segmentsVp[it];
	}
	segmentsVp.clear();
}
void pin::setPos(const point2F& locP2,const point2F& originP2){
	for(size_t it=0;it<apsVp.size();it++){
		apsVp[it]->setPos(locP2,originP2);
	}
	for(size_t it=0;it<segmentsVp.size();it++){
		segmentsVp[it]->setPos(locP2,originP2);
	}
	centerP2-=originP2;
	centerP2+=locP2;
}
void pin::flipVer(const point2F& sizeP2, const point2F& originP2){
	//Flip AP
	for(size_t it=0;it<apsVp.size();it++){
		apsVp[it]->flipVer(sizeP2,originP2);
	}
	//Flip segment
	for(size_t it=0;it<segmentsVp.size();it++){
		segmentsVp[it]->flipVer(sizeP2,originP2);
	}
	//update ap map
	for(int i=0;i<NUM_TRACK/2;i++){
		if(apsMIVAp.count(i)>0){
			vector<ap*> tmpVAp=apsMIVAp[i];
			apsMIVAp[i]=apsMIVAp[NUM_TRACK-i];
			apsMIVAp[NUM_TRACK-i]=tmpVAp;

			if(apsMIVAp[i].size()==0){
				apsMIVAp.erase(i);
			}
			for(size_t i2=0;i2<tmpVAp.size();i2++){
				ap* app=tmpVAp[i2];
				app->trackI=NUM_TRACK-i;
			}
		}
	}
	//update center
	centerP2.yF=originP2.yF*2+sizeP2.yF-centerP2.yF;
}
void pin::flipHor(const point2F& sizeP2, const point2F& originP2){
	for(size_t it=0;it<apsVp.size();it++){
		apsVp[it]->flipHor(sizeP2,originP2);
	}
	for(size_t it=0;it<segmentsVp.size();it++){
		segmentsVp[it]->flipHor(sizeP2,originP2);
	}
	centerP2.xF=originP2.xF*2+sizeP2.xF-centerP2.xF;
}

//segment
segment::segment(string _layerS,vector<point2F*> _pointsVp){
	layerS=_layerS;
	pointsVp=_pointsVp;
}
segment::~segment(){
	for(size_t i=0;i<pointsVp.size();i++){
		delete pointsVp[i];
	}
	pointsVp.clear();
}
void segment::setPos(const point2F& locP2,const point2F& originP2){
	for(size_t i=0;i<pointsVp.size();i++){
		point2F * pointP2P=pointsVp[i];
		*pointP2P-=originP2;
		*pointP2P+=locP2;
	}
}
void segment::flipVer(const point2F& sizeP2, const point2F& originP2){
	for(size_t i=0;i<pointsVp.size();i++){
		point2F * pointP2P=pointsVp[i];
		pointP2P->yF =originP2.yF*2+sizeP2.yF-pointP2P->yF;
	}
}
void segment::flipHor(const point2F& sizeP2, const point2F& originP2){
	for(size_t i=0;i<pointsVp.size();i++){
		point2F * pointP2P=pointsVp[i];
		pointP2P->xF =originP2.xF*2+sizeP2.xF-pointP2P->xF;
	}
}

//ap
ap::ap(string _layerS,point2F _pointP2){
	layerS=_layerS;
	pointP2=_pointP2;
	pinp=NULL;
	trackI=0;
}
ap::~ap(){
	
}
void ap::setPos(const point2F& locP2,const point2F& originP2){
	pointP2-=originP2;
	pointP2+=locP2;
}
void ap::flipVer(const point2F& sizeP2, const point2F& originP2){
	pointP2.yF=originP2.yF*2+sizeP2.yF-pointP2.yF;
}
void ap::flipHor(const point2F& sizeP2, const point2F& originP2){
	pointP2.xF=originP2.xF*2+sizeP2.xF-pointP2.xF;
}


//comp_instance
comp_instance::comp_instance(string _nameS,string _dirS, point2F _locP2){
	nameS=_nameS;
	dirS=_dirS;
	locP2=_locP2;
}


//comp
comp::comp(string _nameS,string _dirS, point2F _locP2, cell* _cellCp):comp_instance(_nameS,_dirS,_locP2){

	point2F zero(0,0);	
	cellCp=_cellCp;
	cellCp->setPos(_locP2,zero);
	
	cellCp->originP2p=&locP2;
	cellCp->setDir(dirS);
	cellCp->compp=this;
	
	update_pin_center();
	
}
comp::~comp(){
	delete cellCp;
}

//comp_macro
comp_macro::comp_macro(string _nameS,string _dirS, point2F _locP2, macro* _macroMp):comp_instance(_nameS,_dirS,_locP2){

	point2F zero(0,0);	
	macroMp=_macroMp;
	macroMp->setPos(_locP2,zero);
	
	macroMp->originP2p=&locP2;
	macroMp->setDir(dirS);
	macroMp->compp=this;
	
	
}
comp_macro::~comp_macro(){
	delete macroMp;
}
bool comp::operator <(const comp& _comp){
	return locP2.xF<_comp.locP2.xF;
}
void comp::update_pin_center(){
	cellCp->update_pin_center();
}
void comp::setPos(const point2F &_locP2){
	//Move comp from locP2 to _locP2
	cellCp->setPos(_locP2,locP2);
	locP2=_locP2;
}


void comp_macro::setPos(const point2F &_locP2){
	//Move comp from locP2 to _locP2
	macroMp->setPos(_locP2,locP2);
	locP2=_locP2;
}
//row
row::row(float _rowHeightF,float _rowX){
	rowHeightF=_rowHeightF;
	rowX=_rowX;
}
row::~row(){
}

//IOpin
IOpin::IOpin(string _nameS, string _layerS, point2F _centerP2, point2F *_sizeP22){
	nameS=_nameS;
	layerS=_layerS;
	centerP2=_centerP2;
	sizeP22[0]=_sizeP22[0];
	sizeP22[1]=_sizeP22[1];
}
IOpin::~IOpin(){
}

//net
net::net(string _nameS){
	nameS=_nameS;
	cellPinMSS.clear();
	steinerG = nullptr;
}
net::~net(){
	if (steinerG != nullptr) {
		delete steinerG;
	}
}

//drv
drv::drv(int _idI, string _typeS, string _layerS, vector<string> _relatedNetsVS, BBox<double> _bbox){
	idI=_idI;
	typeS=_typeS;
	layerS=_layerS;
	relatedNetsVS=_relatedNetsVS;
	bbox=_bbox;
}
drv::~drv(){
}


// FOR FLUTE
//graph
graph::graph(const vector<vertex*> &verticiesVVp, vector<edge> &edgesVE){
	//num_verticiesI=verticiesVVp.size();
	/*
	for(size_t i=0;i<edgesVE.size();i++){
		edge * edgep=new edge;
		*edgep=edgesVE[i];
		edgesVEp.push_back(edgep);
	}
	*/
	/*
	//Sort edge vector by the index of the first end-vertex.
	sort(edgesVE.begin(),edgesVE.end(),
			[](const edge & ledge, const edge & redge) -> bool{
				return ledge.endVerticies[0]->idxI < redge.endVerticies[0]->idxI;
			});
	*/

	//Remove redundant edges
    vector<int> same_idx;
    for(int i1=0;i1<edgesVE.size()-1;i1++){
        int idx1_0=edgesVE[i1].endVerticies[0]->idxI;
        int idx1_1=edgesVE[i1].endVerticies[1]->idxI;
        for(int i2=i1+1;i2<edgesVE.size();i2++){
            int idx2_0=edgesVE[i2].endVerticies[0]->idxI;
            int idx2_1=edgesVE[i2].endVerticies[1]->idxI;

            if((idx1_0==idx2_0 && idx1_1==idx2_1) || (idx1_1==idx2_0 && idx1_0==idx2_1)){
                same_idx.push_back(i2);
            }    
        }    
    }    
    sort(same_idx.begin(), same_idx.end());
    same_idx.erase(unique(same_idx.begin(),same_idx.end()),same_idx.end());
    reverse(same_idx.begin(),same_idx.end());
    for(int i=0;i<same_idx.size();i++){
        edgesVE.erase(edgesVE.begin()+same_idx[i]);
    }
	    
	//Find the number of edges for each vertex
	for(size_t i=0;i<edgesVE.size();i++){
		if(verticiesMII.count(edgesVE[i].endVerticies[0]->idxI)==0){
			verticiesMII.insert(make_pair(edgesVE[i].endVerticies[0]->idxI,0));
		}
		if(verticiesMII.count(edgesVE[i].endVerticies[1]->idxI)==0){
			verticiesMII.insert(make_pair(edgesVE[i].endVerticies[1]->idxI,0));
		}
		/*
		if(debugA==true && (edgesVE[i].endVerticies[0]->idxI==10016 ||  edgesVE[i].endVerticies[1]->idxI==10016)){
			cout<<verticiesMII[10016]<<endl;
			cout<<edgesVE[i].endVerticies[0]->idxI<<endl;
			cout<<edgesVE[i].endVerticies[1]->idxI<<endl;
			verticiesMII[edgesVE[i].endVerticies[0]->idxI]+=1;
			verticiesMII[edgesVE[i].endVerticies[1]->idxI]+=1;
			cout<<verticiesMII[10016]<<endl;
			cout<<endl;
		}
		else{
			verticiesMII[edgesVE[i].endVerticies[0]->idxI]+=1;
			verticiesMII[edgesVE[i].endVerticies[1]->idxI]+=1;
		}
		*/
		verticiesMII[edgesVE[i].endVerticies[0]->idxI]+=1;
		verticiesMII[edgesVE[i].endVerticies[1]->idxI]+=1;
		edge * edgep=new edge;
		*edgep=edgesVE[i];
		edgesVEp.push_back(edgep);
		//assert(verticiesMII[edgesVE[i].endVerticies[0]->idxI]<5);
		assert(verticiesMII[edgesVE[i].endVerticies[0]->idxI]<10);
		
	}
	
	//Assign space for each edge
	verticies_size=verticiesVVp.size();
	verticiesppp=new vertex**[verticiesVVp.size()];
	for(size_t i=0;i<verticiesVVp.size();i++){
		int size=verticiesMII[verticiesVVp[i]->idxI];
		verticiesppp[i] = new vertex*[size+1];
	}

	//Initialize graph
	//map<int,size_t> verIdxIdxMII;//vertex idx to verticiesppp idx
	for(size_t i=0;i<verticiesVVp.size();i++){
		verticiesppp[i][0]=verticiesVVp[i];
		verIdxIdxMII.insert(make_pair(verticiesVVp[i]->idxI,i));
	}
	
	//Build the adjacency graph
	map<int,int> idxMII;
	for(size_t i=0;i<edgesVE.size();i++){
		vertex *startVp=edgesVE[i].endVerticies[0];
		vertex *endVp=edgesVE[i].endVerticies[1];
		int &startIdxI=startVp->idxI;
		int &endIdxI=endVp->idxI;

		verticiesppp[verIdxIdxMII[startIdxI]][idxMII[startIdxI]+1]=endVp;
		verticiesppp[verIdxIdxMII[endIdxI]][idxMII[endIdxI]+1]=startVp;
		idxMII[startIdxI]++;
		idxMII[endIdxI]++;
	
	}
	
}
graph::~graph(){
	/*
	int i=0;
	for(auto it=verticiesMII.begin();it!=verticiesMII.end();it++){
		for(size_t it2=0;it2<it->second;it2++){
			delete verticiesppp[i][it2];
		}
		delete verticiesppp[i];
		i++;
	}
	delete verticiesppp;
	*/

}

/*
void graph::addEdge(const edge& _edge){
	//Edge should have appropriate idx value.
	
	

}
*/

void graph::delete_edge_and_add_vertex(edge* edgep, vertex* vertexp, vector<vertex*>& connected_verticies){
	
	//delete edge
	vertex* vertex1p=edgep->endVerticies[0];
	vertex* vertex2p=edgep->endVerticies[1];

	int vertex1_idx=vertex1p->idxI;
	int vertex2_idx=vertex2p->idxI;

	//verticiesMII[vertex1_idx]--;
	//verticiesMII[vertex2_idx]--;

	int idx=-1;
	for(size_t i=0;i<edgesVEp.size();i++){
		if (edgesVEp[i]==edgep) idx=(int)i;
	}
	assert(idx>=0);
	edgesVEp.erase(edgesVEp.begin()+idx);

	//update verticiesMII and verIdxIdxMII
	int vertex_idx=vertexp->idxI;
	verticiesMII.insert(make_pair(vertex_idx,connected_verticies.size()));
	verIdxIdxMII.insert(make_pair(vertex_idx,(size_t)verticies_size++));

	//add edges
	for(size_t i=0;i<connected_verticies.size();i++){
		//verticiesMII[connected_verticies[i]->idxI]++;
		edge * edgep=new edge(vertexp,connected_verticies[i]);
		edgesVEp.push_back(edgep);
	}

	//size has been increased.
	vertex*** new_verticiesppp=new vertex**[verticies_size];
	for(auto it=verticiesMII.begin();it!=verticiesMII.end();it++){
		int vertex_idx_it=it->first;
		int vertex_degree=it->second;
		new_verticiesppp[verIdxIdxMII[vertex_idx_it]] = new vertex*[vertex_degree+1]; // self, others (so degree +1)
		
		bool is_deleted_edge_vertex1=vertex_idx_it==vertex1_idx || vertex_idx_it==vertex2_idx;
		int control_factor=0;
		if(it->first != vertex_idx){
			new_verticiesppp[verIdxIdxMII[vertex_idx_it]][0]=verticiesppp[verIdxIdxMII[vertex_idx_it]][0]; //assign self vertex
			for(int i=0;i<vertex_degree;i++){
				bool is_connected_vertex=false;
				for(size_t i2=0;i2<connected_verticies.size();i2++){
					if(connected_verticies[i2]->idxI==vertex_idx_it){
						is_connected_vertex=true;
					}
				}
				
				if(i!=vertex_degree-1 || !is_connected_vertex){
					vertex* vertex2p=verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1];
					int vertex_idx_it2=vertex2p->idxI;
					bool is_deleted_edge_vertex2=vertex_idx_it2==vertex1_idx || vertex_idx_it2==vertex2_idx;
					if(is_deleted_edge_vertex1 && is_deleted_edge_vertex2){
						control_factor=1;
					}
					new_verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1]=verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1+control_factor];
				}
				else{
					new_verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1]=vertexp;
					
				}
			}
		}
		else{
			new_verticiesppp[verIdxIdxMII[vertex_idx]][0]=vertexp;
			for(int i=0;i<connected_verticies.size();i++){
				new_verticiesppp[verIdxIdxMII[vertex_idx]][i+1]=connected_verticies[i];
			}
		}
	}

	

	//delete verticiesppp
	for(auto it=verIdxIdxMII.begin();it!=verIdxIdxMII.end();it++){
		if(it->first != vertex_idx){
			delete[] verticiesppp[it->second];
		}
	}
	delete[] verticiesppp;

	verticiesppp=new_verticiesppp;
}

void graph::delete_edge(edge* edgep){
	vertex* vertex1p=edgep->endVerticies[0];
	vertex* vertex2p=edgep->endVerticies[1];

	int vertex1_idx=vertex1p->idxI;
	int vertex2_idx=vertex2p->idxI;

	verticiesMII[vertex1_idx]--;
	verticiesMII[vertex2_idx]--;

	int idx=-1;
	for(size_t i=0;i<edgesVEp.size();i++){
		if (edgesVEp[i]==edgep) idx=(int)i;
	}
	assert(idx>=0);
	edgesVEp.erase(edgesVEp.begin()+idx);
}

void graph::add_vertex(vertex* vertexp, vector<vertex*> connected_verticies){
	int vertex_idx=vertexp->idxI;
	verticiesMII.insert(make_pair(vertex_idx,connected_verticies.size()));
	verIdxIdxMII.insert(make_pair(vertex_idx,(size_t)verticies_size++));

	//add edges
	for(size_t i=0;i<connected_verticies.size();i++){
		verticiesMII[connected_verticies[i]->idxI]++;
		//new_verticiesppp[verIdxIdxMII[connected_verticies[i]->idxI]][verticiesMII[connected_verticies[i]->idxI]]=vertexp;
		edge * edgep=new edge(vertexp,connected_verticies[i]);
		edgesVEp.push_back(edgep);
	}

	//size has been increased.
	vertex*** new_verticiesppp=new vertex**[verticies_size];
	for(auto it=verticiesMII.begin();it!=verticiesMII.end();it++){
		int vertex_idx_it=it->first;
		int vertex_degree=it->second;
		new_verticiesppp[verIdxIdxMII[vertex_idx_it]] = new vertex*[vertex_degree+1];
		if(it->first != vertex_idx){
			new_verticiesppp[verIdxIdxMII[vertex_idx_it]][0]=verticiesppp[verIdxIdxMII[vertex_idx_it]][0];
			for(int i=0;i<vertex_degree;i++){
				bool is_connected_vertex=false;
				for(size_t i2=0;i2<connected_verticies.size();i2++){
					if(connected_verticies[i2]->idxI==vertex_idx_it){
						is_connected_vertex=true;
					}
				}
				
				if(i!=vertex_degree-1 || !is_connected_vertex){
					new_verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1]=verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1];
				}
				else{
					new_verticiesppp[verIdxIdxMII[vertex_idx_it]][i+1]=vertexp;
					
				}
			}
		}
		else{
			new_verticiesppp[verIdxIdxMII[vertex_idx]][0]=vertexp;
			for(int i=0;i<connected_verticies.size();i++){
				new_verticiesppp[verIdxIdxMII[vertex_idx]][i+1]=connected_verticies[i];
			}
		}
	}

	

	//delete verticiesppp
	for(auto it=verIdxIdxMII.begin();it!=verIdxIdxMII.end();it++){
		if(it->first != vertex_idx){
			delete[] verticiesppp[it->second];
		}
	}
	delete[] verticiesppp;

	verticiesppp=new_verticiesppp;

	
}
void graph::change_edge(int idx1, int idx2, int reference_idx){
	//Change edge from idx1 - reference_idx to idx2 - reference_idx
	bool skip=false;
	for(size_t i=0;i<edgesVEp.size();i++){
		//If there already exists edge vetween idx2-ref then pass
		int edge_idx0=edgesVEp[i]->endVerticies[0]->idxI;
		int edge_idx1=edgesVEp[i]->endVerticies[1]->idxI;
		bool is_edge_idx0_idx2=(edge_idx0==idx2 && edge_idx1==reference_idx);
		bool is_edge_idx1_idx2=(edge_idx0==reference_idx && edge_idx1==idx2);
		if(is_edge_idx0_idx2 || is_edge_idx1_idx2){
			skip=true;
		}
		
	}

	if(!skip){
		verticiesMII[idx1]--;
		verticiesMII[idx2]++;
		vertex* vertex_idx2=verticiesppp[verIdxIdxMII[idx2]][0];
		bool edge_checker=false;
		for(size_t i=0;i<edgesVEp.size();i++){
			int edge_idx0=edgesVEp[i]->endVerticies[0]->idxI;
			int edge_idx1=edgesVEp[i]->endVerticies[1]->idxI;
			bool is_edge_idx0_idx1=(edge_idx0==idx1 && edge_idx1==reference_idx);
			bool is_edge_idx1_idx1=(edge_idx0==reference_idx && edge_idx1==idx1);
			if(is_edge_idx0_idx1 || is_edge_idx1_idx1){
				edge_checker=true;
				if(is_edge_idx0_idx1){
					edgesVEp[i]->endVerticies[0]=vertex_idx2;
					edgesVEp[i]->endPointsP22[0]=vertex_idx2->centerP2;
				}
				else{
					edgesVEp[i]->endVerticies[1]=vertex_idx2;
					edgesVEp[i]->endPointsP22[1]=vertex_idx2->centerP2;
				}
				break;
			}
		}
		assert(edge_checker); // Otherwise, no edge found that connects idx1 and reference_idx
		
		//add new vertex
		vertex** new_verticiespp=new vertex*[verticiesMII[idx2]];
		for(int i=0;i<verticiesMII[idx2];i++){
			new_verticiespp[i]=verticiesppp[verIdxIdxMII[idx2]][i];
		}
		new_verticiespp[verticiesMII[idx2]]=verticiesppp[verIdxIdxMII[reference_idx]][0];
		delete [] verticiesppp[verIdxIdxMII[idx2]];
		verticiesppp[verIdxIdxMII[idx2]]=new_verticiespp;
	}
}

void graph::move_vertex(const int& vertex_idx, const point2F& new_coord){
	size_t verticiesppp_idx=verIdxIdxMII[vertex_idx];
	vertex* vertexp=verticiesppp[verticiesppp_idx][0];
	vertexp->centerP2=new_coord;
}

void graph::print(){
	for(auto it=verticiesMII.begin();it!=verticiesMII.end();it++){
		//cout<<it->first<<" ";
		for(size_t i2=0;i2<it->second+1;i2++){
			cout<<verticiesppp[verIdxIdxMII[it->first]][i2]->idxI<<" "<<verticiesppp[verIdxIdxMII[it->first]][i2]->centerP2.xF<<" "<<verticiesppp[verIdxIdxMII[it->first]][i2]->centerP2.yF <<"/ ";
		}
		cout<<endl;
	}
}

//vertex
vertex::vertex(){
	idxI=0;
	typeI=0;
	pinTypeI=0;
	pinp=NULL;
	IOpinp=NULL;
	centerP2={0,0};
}
vertex::~vertex(){
}
bool vertex::operator == (const vertex& _vertex){
	if(idxI==_vertex.idxI && typeI==_vertex.typeI && centerP2==_vertex.centerP2){
		return true;
	}
	return false;
}
vertex vertex::operator = (const vertex& _vertex){
	idxI=_vertex.idxI;
	typeI=_vertex.typeI;
	centerP2=_vertex.centerP2;
	return *this;
}

//edge
edge::edge(bool _direction,string _layerS,point2F* _pointP22){
	direction=_direction;
	layerS=_layerS;
	endPointsP22[0]=_pointP22[0];
	endPointsP22[1]=_pointP22[1];
	endVerticies[0]=NULL;
	endVerticies[1]=NULL;
	num_apsI=-1;
	source_app=NULL;
}
edge::edge(){
	direction=false;
	endVerticies[0]=NULL;
	endVerticies[1]=NULL;
	num_apsI=-1;
	source_app=NULL;
}
/*
edge::edge(point2F* _pointP22){
	endPointsP22[0]=_pointP22[0];
	endPointsP22[1]=_pointP22[1];
}
*/
edge::edge(point2F _p1, point2F _p2){
	endPointsP22[0]=_p1;
	endPointsP22[1]=_p2;
	endVerticies[0]=NULL;
	endVerticies[1]=NULL;
	num_apsI=-1;
	source_app=NULL;
}
edge::edge(vertex* vp1, vertex* vp2){
	endPointsP22[0]=vp1->centerP2;
	endPointsP22[1]=vp2->centerP2;
	endVerticies[0]=vp1;
	endVerticies[1]=vp2;
	num_apsI=-1;
	source_app=NULL;
}
edge::~edge(){
}
void edge::clear(){
	//vertex dump1V;
	//vertex dump2V;
	point2F dump1P2={0,0};
	point2F dump2P2={0,0};

	direction=false;
	layerS.clear();

	endPointsP22[0]=dump1P2;
	endPointsP22[1]=dump2P2;
	/*
	delete endVerticies[0];
	delete endVerticies[1];
	*/
	/*
	endVerticies[0]=dump1V;
	endVerticies[1]=dump2V;
	*/
}
/*
edge edge::operator = (const edge& _edge){
	direction=_edge.direction;
	layerS=_edge.layerS;
	endPointsP22[0]=_edge.endPointsP22[0];
	endPointsP22[1]=_edge.endPointsP22[1];
	endVerticies[0]=_edge.endVerticies[0];
	endVerticies[1]=_edge.endVerticies[1];
	return *this;
}
*/

//gcell
gcell::gcell(point2F _originPF){
	originPF=_originPF;
	upgE=nullptr;
	downgE=nullptr;
	rightgE=nullptr;
	leftgE=nullptr;
	demand=0;
}
gcell::~gcell(){
	delete upgE;
	delete downgE;
	delete leftgE;
	delete rightgE;
	upgE=nullptr;
	downgE=nullptr;
	rightgE=nullptr;
	leftgE=nullptr;
}


//gedge
gedge::gedge(point2I* _endPointsPI2, point2F* _endPointsPF2){
	demandF=0.0;
	endPointsPI2[0]=_endPointsPI2[0];
	endPointsPI2[1]=_endPointsPI2[1];
	endPointsPF2[0]=_endPointsPF2[0];
	endPointsPF2[1]=_endPointsPF2[1];
}
gedge::~gedge(){
}

//gcell_grid
gcell_grid::gcell_grid(Def &def){
	this->def=&def;
	num_grids_x=((def.die_sizeP22[1].xF*1000)-OFFSET_X*2*1000)/(GCELL_SIZE*1000);
	num_grids_y=((def.die_sizeP22[1].yF*1000)-OFFSET_Y*2*1000)/(GCELL_SIZE*1000);
	for(int i_x=0;i_x<num_grids_x;i_x++){
		for(int i_y=0;i_y<num_grids_y;i_y++){
			point2I idxPI(i_x,i_y);
			point2I idx_upPI(i_x,i_y+1);
			point2I idx_rightPI(i_x+1,i_y);
			point2I idx_uprightPI(i_x+1,i_y+1);
			point2I idx_leftPI2[2]={idxPI,idx_upPI};
			point2I idx_downPI2[2]={idxPI,idx_rightPI};

			point2F coordsPF=convert_idx_to_coords(idxPI);
			point2F coords_upPF=convert_idx_to_coords(idx_upPI);
			point2F coords_rightPF=convert_idx_to_coords(idx_rightPI);
			point2F coords_uprightPF=convert_idx_to_coords(idx_uprightPI);
			point2F coords_leftPF2[2]={coordsPF,coords_upPF};
			point2F coords_downPF2[2]={coordsPF,coords_rightPF};

			gcell* gcellp=new gcell(coordsPF);
			gedge* lgedge=new gedge(idx_leftPI2,coords_leftPF2);
			gedge* dgedge=new gedge(idx_downPI2,coords_downPF2);
			gcellp->leftgE=lgedge;
			gcellp->downgE=dgedge;
			
			if(i_x>=1){
				point2I idx_leftPI(i_x-1,i_y);
				gcell* gcell_leftp=gcellsMPIGp[idx_leftPI];
				gcell_leftp->rightgE=lgedge;
			}
			if(i_y>=1){
				point2I idx_downPI(i_x,i_y-1);
				gcell* gcell_downp=gcellsMPIGp[idx_downPI];
				gcell_downp->upgE=dgedge;
			}
			if(i_x==num_grids_x-1){
				point2I idx_rightPI2[2]={idx_rightPI,idx_uprightPI};
				point2F coords_rightPF2[2]={coords_rightPF,coords_uprightPF};
				gedge* redge=new gedge(idx_rightPI2,coords_rightPF2);
				gcellp->rightgE=redge;
			}
			if(i_y==num_grids_y-1){
				point2I idx_upPI2[2]={idx_upPI,idx_uprightPI};
				point2F coords_upPF2[2]={coords_upPF,coords_uprightPF};
				gedge* uedge=new gedge(idx_upPI2,coords_upPF2);
				gcellp->upgE=uedge;
			}
			gcellsMPIGp.insert(make_pair(idxPI,gcellp));
		}
	}
	
	assign_net_info();
}
gcell_grid::~gcell_grid(){
	for(auto it=gcellsMPIGp.begin();it!=gcellsMPIGp.end();it++){
		gcell* gcellp=it->second;
		delete gcellp;
		gcellp=NULL;
	}
}

void gcell_grid::assign_net_info(){
	for(auto it_net=def->netsMSNp.begin();it_net!=def->netsMSNp.end();it_net++){
		net* netp=it_net->second;
		if(netp->steinerG!=NULL){
			graph* netgp=netp->steinerG;
			insert_graph_info_into_gcell_grid(netgp,false);
		}
	}
}
void gcell_grid::insert_graph_info_into_gcell_grid(graph* netgp, bool is_remove){
	vector<edge*> &edgesVEp=netgp->edgesVEp;
	for(size_t it_edge=0;it_edge<edgesVEp.size();it_edge++){
		edge* edgep=edgesVEp[it_edge];
		point2F start_pointPF=edgep->endPointsP22[0];
		point2F end_pointPF=edgep->endPointsP22[1];
		
		point2I start_pointPI=convert_coords_to_idx(start_pointPF);
		point2I end_pointPI=convert_coords_to_idx(end_pointPF);

		start_pointPI.xI=min(max(start_pointPI.xI,0),num_grids_x-1);
		start_pointPI.yI=min(max(start_pointPI.yI,0),num_grids_y-1);
		end_pointPI.xI=min(max(end_pointPI.xI,0),num_grids_x-1);
		end_pointPI.yI=min(max(end_pointPI.yI,0),num_grids_y-1);
		
		int start_x=min(start_pointPI.xI,end_pointPI.xI);
		int end_x=max(start_pointPI.xI,end_pointPI.xI);
		int start_y=min(start_pointPI.yI,end_pointPI.yI);
		int end_y=max(start_pointPI.yI,end_pointPI.yI);

		int x_diff=end_x-start_x;
		int y_diff=end_y-start_y;
		float weightF=(!is_remove)?1:-1;
		if(x_diff==0 && y_diff!=0){
			for(int it_y=0;it_y<y_diff;it_y++){
				point2I cur_coordPI={start_x,start_y+it_y};
				gcellsMPIGp[cur_coordPI]->upgE->demandF+=weightF;
			}
		}
		else if (x_diff!=0 && y_diff==0){
			for(int it_x=0;it_x<x_diff;it_x++){
				point2I cur_coordPI={start_x+it_x,start_y};
				gcellsMPIGp[cur_coordPI]->rightgE->demandF+=weightF;
			}
		}
		else{
			//x_diff !=0 && y_diff !=0
			weightF=(!is_remove)?0.5:-0.5;
			for(int it_y=0;it_y<y_diff;it_y++){
				point2I cur_startPI={start_x,start_y+it_y};
				point2I cur_endPI={end_x,start_y+it_y};
				gcellsMPIGp[cur_startPI]->upgE->demandF+=weightF;
				gcellsMPIGp[cur_endPI]->upgE->demandF+=weightF;
			}
			for(int it_x=0;it_x<x_diff;it_x++){
				point2I cur_startPI={start_x+it_x,start_y};
				point2I cur_endPI={start_x+it_x,end_y};
				gcellsMPIGp[cur_startPI]->rightgE->demandF+=weightF;
				gcellsMPIGp[cur_endPI]->rightgE->demandF+=weightF;
			}
		}
	}
}
float gcell_grid::search_demand_for_line(const point2F& start_point, const point2F& end_point){
	assert(start_point.xF==end_point.xF || start_point.yF==end_point.yF);
	point2I start_pointPI=convert_coords_to_idx(start_point);
	point2I end_pointPI=convert_coords_to_idx(end_point);
	
	float total_demand=0;
	if(start_pointPI.yI==end_pointPI.yI){
		for(int i_x=start_pointPI.xI;i_x<end_pointPI.xI;i_x++){
			point2I tmp_point={i_x,start_pointPI.yI};
			if(gcellsMPIGp.count(tmp_point)>0){
				gcell *gcellp=gcellsMPIGp[tmp_point];
				total_demand+=gcellp->rightgE->demandF;
			}
		}
	}
	else{
		for(int i_y=start_pointPI.yI;i_y<end_pointPI.yI;i_y++){
			point2I tmp_point={start_pointPI.xI,i_y};
			if(gcellsMPIGp.count(tmp_point)>0){
				gcell *gcellp=gcellsMPIGp[tmp_point];
				total_demand+=gcellp->upgE->demandF;
			}
		}
	}

	return total_demand;
}
edge gcell_grid::search_min_demand_edge(point2F& start_point, point2F& end_point, int xy){

	point2I start_pointPI=convert_coords_to_idx(start_point);
	point2I end_pointPI=convert_coords_to_idx(end_point);
	
	point2I min_startPI;
	point2I min_endPI;

	if(xy==0){
		//X
		float min_demand=MAXF;
		int min_point=0;
		for(int i1=start_pointPI.xI;i1<=end_pointPI.xI;i1++){
			float line_demand=0;
			for(int i2=start_pointPI.yI;i2<end_pointPI.yI;i2++){
				point2I tmp_point={i1,i2};
				if(gcellsMPIGp.count(tmp_point)!=0){
					gcell* gcellp=gcellsMPIGp[tmp_point];
					line_demand+=gcellp->upgE->demandF;
				}
			}
			if(line_demand<min_demand){
				min_demand=line_demand;
				min_point=i1;
			}
		}
		min_startPI={min_point,start_pointPI.yI};
		min_endPI={min_point,end_pointPI.yI};
	}
	else{
		//Y
		float min_demand=MAXF;
		int min_point=0;
		for(int i1=start_pointPI.yI;i1<=end_pointPI.yI;i1++){
			float line_demand=0;
			for(int i2=start_pointPI.xI;i2<end_pointPI.xI;i2++){
				point2I tmp_point={i2,i1};
				if(gcellsMPIGp.count(tmp_point)!=0){
					gcell* gcellp=gcellsMPIGp[tmp_point];
					line_demand+=gcellp->rightgE->demandF;
				}
			}
			if(line_demand<min_demand){
				min_demand=line_demand;
				min_point=i1;
			}
		}
		min_startPI={start_pointPI.xI,min_point};
		min_endPI={end_pointPI.xI,min_point};
	}
	
	point2F min_startF=convert_idx_to_coords(min_startPI);
	point2F min_endF=convert_idx_to_coords(min_endPI);

	if(xy==0){
		min_startF.xF=min(max(start_point.xF,min_startF.xF),end_point.xF);
		min_endF.xF=min_startF.xF;
		min_startF.yF=start_point.yF;
		min_endF.yF=end_point.yF;
	}
	else{
		min_startF.yF=min(max(start_point.yF,min_startF.yF),end_point.yF);
		min_endF.yF=min_startF.yF;
		min_startF.xF=start_point.xF;
		min_endF.xF=end_point.xF;
	
	}
	edge min_demand_edge(min_startF,min_endF);
	return min_demand_edge;
	
}

void gcell_grid::print(int x,int y){
	point2I tmpPI={x,y};
	gcell* gcellp=gcellsMPIGp[tmpPI];
	cout<<"Gcell Origin:"<<gcellp->originPF.xF<<","<<gcellp->originPF.yF<<endl;
	cout<<"upgE demand:"<<gcellp->upgE->demandF<<endl;
	cout<<"     I :("<<gcellp->upgE->endPointsPI2[0].xI<<","<<gcellp->upgE->endPointsPI2[0].yI<<"), ("<<gcellp->upgE->endPointsPI2[1].xI<<","<<gcellp->upgE->endPointsPI2[1].yI<<")"<<endl;
	cout<<"     F :("<<gcellp->upgE->endPointsPF2[0].xF<<","<<gcellp->upgE->endPointsPF2[0].yF<<"), ("<<gcellp->upgE->endPointsPF2[1].xF<<","<<gcellp->upgE->endPointsPF2[1].yF<<")"<<endl;


	cout<<"downgE demand:"<<gcellp->downgE->demandF<<endl;
	cout<<"     I :("<<gcellp->downgE->endPointsPI2[0].xI<<","<<gcellp->downgE->endPointsPI2[0].yI<<"), ("<<gcellp->downgE->endPointsPI2[1].xI<<","<<gcellp->downgE->endPointsPI2[1].yI<<")"<<endl;
	cout<<"     F :("<<gcellp->downgE->endPointsPF2[0].xF<<","<<gcellp->downgE->endPointsPF2[0].yF<<"), ("<<gcellp->downgE->endPointsPF2[1].xF<<","<<gcellp->downgE->endPointsPF2[1].yF<<")"<<endl;


	cout<<"leftgE demand:"<<gcellp->leftgE->demandF<<endl;
	cout<<"     I :("<<gcellp->leftgE->endPointsPI2[0].xI<<","<<gcellp->leftgE->endPointsPI2[0].yI<<"), ("<<gcellp->leftgE->endPointsPI2[1].xI<<","<<gcellp->leftgE->endPointsPI2[1].yI<<")"<<endl;
	cout<<"     F :("<<gcellp->leftgE->endPointsPF2[0].xF<<","<<gcellp->leftgE->endPointsPF2[0].yF<<"), ("<<gcellp->leftgE->endPointsPF2[1].xF<<","<<gcellp->leftgE->endPointsPF2[1].yF<<")"<<endl;


	cout<<"rightgE demand:"<<gcellp->rightgE->demandF<<endl;
	cout<<"     I :("<<gcellp->rightgE->endPointsPI2[0].xI<<","<<gcellp->rightgE->endPointsPI2[0].yI<<"), ("<<gcellp->rightgE->endPointsPI2[1].xI<<","<<gcellp->rightgE->endPointsPI2[1].yI<<")"<<endl;
	cout<<"     F :("<<gcellp->rightgE->endPointsPF2[0].xF<<","<<gcellp->rightgE->endPointsPF2[0].yF<<"), ("<<gcellp->rightgE->endPointsPF2[1].xF<<","<<gcellp->rightgE->endPointsPF2[1].yF<<")"<<endl;
}

//steiner_container
steiner_container::steiner_container(vertex* _steiner_vertexp, const graph* netgp, int _vertex_idx){
	assert(_steiner_vertexp->typeI==STEINER);
	steiner_vertexp=_steiner_vertexp;
	extended_sliding_edgeM[0].clear();
	extended_sliding_edgeM[1].clear();
	

	vertex *** verticiesppp=netgp->verticiesppp;
	map<int,int> verticiesMII=netgp->verticiesMII; //vertex idx to number of edges
	map<int,size_t> verIdxIdxMII=netgp->verIdxIdxMII; // vertex idx to verticiesppp idx

	vertex_idx=_vertex_idx;
	vertex_degree=verticiesMII.find(vertex_idx)->second;
	size_t verticiesppp_idx=verIdxIdxMII.find(vertex_idx)->second;

	//update adjacent steiner, vertex
	for (int i=1;i<vertex_degree+1;i++){
		vertex* vertexp=verticiesppp[verticiesppp_idx][i];
		if(vertexp->typeI==STEINER){
			adjacent_steiner.push_back(vertexp);
		}
		else if (vertexp->typeI==CELL){
			adjacent_verticies.push_back(vertexp);
		}
	}
	
	//update sliding range
	float max_x=MINF;
	float min_x=MAXF;
	float max_y=MINF;
	float min_y=MAXF;
	for(size_t i_av=0;i_av!=adjacent_verticies.size();i_av++){
		vertex* vp=adjacent_verticies[i_av];
		if(vp->centerP2.xF>max_x){
			max_x=vp->centerP2.xF;
			//max_xvp=vp;	
		}
		if(vp->centerP2.xF<min_x){
			min_x=vp->centerP2.xF;
			//min_xvp=vp;
		}
		if(vp->centerP2.yF>max_y){
			max_y=vp->centerP2.yF;
			//max_yvp=vp;
		}
		if(vp->centerP2.yF<min_y){
			min_y=vp->centerP2.yF;
			//min_yvp=vp;
		}
	}
	for(size_t i_sp=0;i_sp!=adjacent_steiner.size();i_sp++){
		vertex* vp=adjacent_steiner[i_sp];
		if(vp->centerP2.xF>max_x){
			max_x=vp->centerP2.xF;
			//max_xvp=vp;
		}
		if(vp->centerP2.xF<min_x){
			min_x=vp->centerP2.xF;
			//min_xvp=vp;
		}
		if(vp->centerP2.yF>max_y){
			max_y=vp->centerP2.yF;
			//max_yvp=vp;
		}
		if(vp->centerP2.yF<min_y){
			min_y=vp->centerP2.yF;
			//min_yvp=vp;
		}
	}
	point2F min_xP={min_x,steiner_vertexp->centerP2.yF};
	point2F max_xP={max_x,steiner_vertexp->centerP2.yF};
	point2F min_yP={steiner_vertexp->centerP2.xF,min_y};
	point2F max_yP={steiner_vertexp->centerP2.xF,max_y};
	
	sliding_edge[0]=edge(min_xP,max_xP);
	sliding_edge[1]=edge(min_yP,max_yP);
	extended_sliding_range[0]=sliding_edge[0];
	extended_sliding_range[1]=sliding_edge[1];
}
steiner_container::steiner_container(){
}
steiner_container::~steiner_container(){
}
void steiner_container::destroy(){
	/*
	if(sliding_edge!=NULL) delete sliding_edge;
	sliding_edge=NULL;
	*/
}
void steiner_container::extend_sliding_edge(const steiner_container &sc,int xy){
	assert(xy<2);
	int sc_vertex_idx=sc.vertex_idx;
	edge sc_vertex_sliding_edge=sc.sliding_edge[xy];
	//extended_sliding_edgeM[xy].insert(make_pair(sc_vertex_idx,&sc_vertex_sliding_edge));
	extended_sliding_edgeM[xy].insert(make_pair(sc_vertex_idx,sc_vertex_sliding_edge));

	//min max check
	if(xy==0){
		//X
		float sc_start_x=sc.sliding_edge->endPointsP22[0].xF;
		float sc_end_x=sc.sliding_edge->endPointsP22[1].xF;
		assert(sc_start_x<=sc_end_x);
		if(extended_sliding_range[0].endPointsP22[0].xF>sc_start_x){
			extended_sliding_range[0].endPointsP22[0].xF=sc_start_x;
		}
		if(extended_sliding_range[0].endPointsP22[1].xF<sc_end_x){
			extended_sliding_range[0].endPointsP22[1].xF=sc_end_x;
		}
	}
	else{
		//Y
		float sc_start_y=sc.sliding_edge->endPointsP22[0].yF;
		float sc_end_y=sc.sliding_edge->endPointsP22[1].yF;
		assert(sc_start_y<=sc_end_y);
		if(extended_sliding_range[1].endPointsP22[0].yF>sc_start_y){
			extended_sliding_range[1].endPointsP22[0].yF=sc_start_y;
		}
		if(extended_sliding_range[1].endPointsP22[1].yF<sc_end_y){
			extended_sliding_range[1].endPointsP22[1].yF=sc_end_y;
		}
	}

}
/*
void steiner_container::get_vertex_info(const vertex ***&verticiesppp,const size_t& verticiesppp_idx,const int& vertex_degree, bool is_updating_adj_steiner=false, vector<int>* done_steiner_idx=NULL){
	//update adjacent steiner, vertex
	for (int i=1;i<vertex_degree+1;i++){
		vertex* vertexp=verticiesppp[verticiesppp_idx][i];
		if(vertexp->typeI==STEINER && is_updating_adj_steiner){
			adjacent_steiner.push_back(vertexp);
			done_steiner_idx->push_back(vertexp->idxI);
		}
		else if (vertexp->typeI==CELL){
			adjacent_verticies.push_back(vertexp);
		}
	}
	
	//update sliding range
	float max_x=MINF;
	float min_x=MAXF;
	float max_y=MINF;
	float min_y=MAXF;
	for(size_t i_av=0;i_av!=adjacent_verticies.size();i_av++){
		vertex* vp=adjacent_verticies[i_av];
		if(vp->centerP2.xF>max_x){
			max_x=vp->centerP2.xF;
		}
		if(vp->centerP2.xF<min_x){
			min_x=vp->centerP2.xF;
		}
		if(vp->centerP2.yF>max_y){
			max_y=vp->centerP2.yF;
		}
		if(vp->centerP2.yF<min_y){
			min_y=vp->centerP2.yF;
		}
	}
	point2F max_xy={max_x,max_y};
	point2F min_xy={min_x,min_y};
	
	sliding_range=new edge(min_xy,max_xy);

}
*/