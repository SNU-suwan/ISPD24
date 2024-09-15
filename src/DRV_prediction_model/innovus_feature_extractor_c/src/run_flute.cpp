#include "run_flute.h"


void initialize_net_v2(Def &def){
	readLUT();
	for(auto it=def.netsMSNp.begin();it!=def.netsMSNp.end();it++){
		net* netp=it->second;
		unordered_map<string,string> cellPinMSS=netp->cellPinMSS;
		
		vector<point2F> center_comp_VP2;
		for(auto &it2:cellPinMSS){
			string cell_name_S=it2.first;
			string pin_name_S=it2.second;
			if(cell_name_S=="PIN"){
				IOpin* cur_IOpin_p=def.IOpinsMSpp[pin_name_S];
				point2F center_P2=cur_IOpin_p->centerP2;
				/*
				for(size_t i_row=0;i_row<rowsVp.size()-1;i_row++){
					if (rowsVp[i_row]->rowHeightF<=center_P2.yF && center_P2.yF<rowsVp[i_row+1]->rowHeightF){
						
					}
				}
				*/
				center_comp_VP2.push_back(center_P2);
			}
			else{
				comp* cur_comp_p=def.compsMSCp[cell_name_S];
				point2F center_P2={cur_comp_p->cellCp->pinsMSPp[pin_name_S]->centerP2.xF,cur_comp_p->locP2.yF+cur_comp_p->cellCp->sizeP2.yF/2};
				center_comp_VP2.push_back(center_P2);
			}	
		}

		//run FLUTE
		/*
		string flute_output_s;
		if(center_comp_VP2.size()>2 && !is_small_flute){
			flute_output_s=exe_flute(flute_dir_s,center_comp_VP2,ofs,netp->nameS);
		}

		if(center_comp_VP2.size()<NET_NUM_THRESHOLD && center_comp_VP2.size()>2){ // FLUTE cannot make large nets or direct net.
			//save the results if needed
			if(is_small_flute){
				flute_output_s=exe_flute(flute_dir_s,center_comp_VP2,ofs,netp->nameS);
			}
			graph *netGraphp=parse_flute(def,flute_output_s,netp->nameS);
			
			netp->steinerG=netGraphp;
		}
		*/
		if (center_comp_VP2.size()>2 && center_comp_VP2.size()<NET_NUM_THRESHOLD){
			graph* netGraphp=run_flute_for_each_net(def,center_comp_VP2,netp->nameS);
			netp->steinerG=netGraphp;
		}
		else if (center_comp_VP2.size()==2){
			graph *netGraphp=direct_connection(def,netp);
			netp->steinerG=netGraphp;
		}
	}
}

graph *run_flute_for_each_net(Def &def, vector<point2F> center_comp_VP2, string netnameS){

	int vec_size=center_comp_VP2.size();
	int x[vec_size], y[vec_size];

	for (int i=0;i<vec_size;i++){
		x[i]=center_comp_VP2[i].xF*1000;
		y[i]=center_comp_VP2[i].yF*1000;
	}
	
	Tree flutetree=flute(vec_size,x,y,ACCURACY);

	int i;
	string data="";
	for (i=0;i<flutetree.deg;i++){
		data+=" "+to_string(i)+" :   x = "+to_string(flutetree.branch[i].x)\
						 +" y = "+to_string(flutetree.branch[i].y)\
						+"  e = "+to_string(flutetree.branch[i].n)+"\n";
	}
	for (i=flutetree.deg;i<2*flutetree.deg-2;i++){
		data+=" s"+to_string(i)+" :   x = "+to_string(flutetree.branch[i].x)\
						 +" y = "+to_string(flutetree.branch[i].y)\
						+"  e = "+to_string(flutetree.branch[i].n)+"\n";
	}
	data+="\n";
	for (i=0; i<2*flutetree.deg-2; i++) {
        data+=to_string(flutetree.branch[i].x)+" "+to_string(flutetree.branch[i].y)+"\n";
		data+=to_string(flutetree.branch[flutetree.branch[i].n].x)+" "+to_string(flutetree.branch[flutetree.branch[i].n].y)+"\n";
		/*
		printf("%d %d\n", flutetree.branch[i].x, flutetree.branch[i].y);
        printf("%d %d\n\n", flutetree.branch[flutetree.branch[i].n].x,
               flutetree.branch[flutetree.branch[i].n].y);
		*/
    }

	string delimiter="\n";
	string token;
	size_t last=0;
	size_t next=0;
	bool is_write_tree_done=false;
	vector<vertex*> verticiesVVp;
	vector<point2F> edgeVP2;
	vector<edge> edgesVE;
	edge edgeE;
	int idx_vertex=-1;
	
	while((next=data.find(delimiter,last))!=string::npos){
		token=data.substr(last,next-last);
		last=next+1;
		
		vector<string> info=split(token);
		
		if(!is_write_tree_done && info.size()!=0){
			//Collect verticies
			idx_vertex++; //if idx_vertex is over 9, the idx changes. (ex. from 6[0] :[1] x[2]..  s10:[0] x[1] ...)
			vertex *newVp=new vertex;
			point2F centerP2;
			if(info[0].rfind("s",0)==0){
				//newVp->idxI=stoi(info[0].substr(1,info[0].size()-1))+STEINER_VERTEX;
				/*
				if(idx_vertex<10)
					newVp->idxI=stoi(info[0].substr(1,info[0].size()-1))+STEINER_VERTEX;
				else
					newVp->idxI=stoi(info[0].substr(1,info[0].size()-2))+STEINER_VERTEX;
				*/
				newVp->idxI=stoi(info[0].substr(1,info[0].size()-1))+STEINER_VERTEX;
				
				newVp->typeI=STEINER;
				newVp->pinTypeI=STEINER_PIN;
			}
			else{
				//newVp->idxI=stoi(info[0]);
				/*
				if(idx_vertex<10)
					newVp->idxI=stoi(info[0]);
				else
					newVp->idxI=stoi(info[0].substr(0,info[0].size()-1));
				*/
				newVp->idxI=stoi(info[0]);
				
				newVp->typeI=CELL;
			}
			centerP2.xF=stof(info[4])/1000;
			centerP2.yF=stof(info[7])/1000;
			/*
			if(idx_vertex<10){
				centerP2.xF=stof(info[4])/1000;
				centerP2.yF=stof(info[7])/1000;
			}
			else{
				centerP2.xF=stof(info[3])/1000;
				centerP2.yF=stof(info[6])/1000;
			}
			*/
			newVp->centerP2=centerP2;
			verticiesVVp.push_back(newVp);
		}
		else if(is_write_tree_done && info.size()!=0){
			//After collect verticies, collect edge
			float xF=stof(info[0])/1000;
			float yF=stof(info[1])/1000;
			point2F centerP2={xF,yF};
			edgeVP2.push_back(centerP2);
			if(edgeVP2.size()==2){
				edgeE.endPointsP22[0]=edgeVP2[0];
				edgeE.endPointsP22[1]=edgeVP2[1];
				//edgesVE.push_back(edgeE);
				
				if(edgeE.endPointsP22[0]!=edgeE.endPointsP22[1]){
					//same vertex test 1
					edgesVE.push_back(edgeE);
				}
				
				edgeVP2.clear();
				edgeE.clear();
			}
		}
		if(info.size()==0 && !is_write_tree_done){
			//Writing tree has finished
			is_write_tree_done=true;

			//Delete duplicate elements 
			vector<size_t> dup_idxVs;
			for(size_t i=0;i<verticiesVVp.size();i++){
				for(size_t i2=i+1;i2<verticiesVVp.size();i2++){
					if(verticiesVVp[i]->centerP2==verticiesVVp[i2]->centerP2){
						//Insert element only when there is no the same element.
						bool isNotIn=true;
						for(size_t iv=0;iv<dup_idxVs.size();iv++){
							if(dup_idxVs[iv]==i2) isNotIn=false;
						}
						if(isNotIn)
							dup_idxVs.push_back(i2);
					}
				}
			}
			//Sort values in a descending order
			std::sort(dup_idxVs.begin(),dup_idxVs.end(),greater<size_t>());
			//reverse(dup_idxVs.begin(),dup_idxVs.end());
			for(size_t i=0;i<dup_idxVs.size();i++){
				delete verticiesVVp[dup_idxVs[i]];
				verticiesVVp.erase(verticiesVVp.begin()+dup_idxVs[i]);
			}

		}
	}
	
	//Match edge points with verticies
	for(size_t i=0;i<edgesVE.size();i++){
		for(size_t i2=verticiesVVp.size();i2-->0;){
			//The reason for the inverse iteration is that if the vertex type is cell and steiner point simultaneously, it is considered to be a cell type. 
			if(edgesVE[i].endPointsP22[0].equal(verticiesVVp[i2]->centerP2)){
				edgesVE[i].endVerticies[0]=verticiesVVp[i2];
			}
			if(edgesVE[i].endPointsP22[1].equal(verticiesVVp[i2]->centerP2)){
				edgesVE[i].endVerticies[1]=verticiesVVp[i2];
			}
		}
	}

	//Match verticies with its original pin
	for(size_t i=0;i<verticiesVVp.size();i++){
		point2F &centerP2=verticiesVVp[i]->centerP2;
		net* cur_netp=def.netsMSNp[netnameS];
		unordered_map<string,string> cellPinMSS=cur_netp->cellPinMSS;
		for(auto& it_cp:cellPinMSS){
			string comp_name=it_cp.first;
			string pin_name=it_cp.second;
			
			if(comp_name=="PIN"){
				IOpin* cur_IOpin_p=def.IOpinsMSpp[pin_name];
				point2F center_pin_p2=cur_IOpin_p->centerP2;

				if(centerP2.equal_p2(center_pin_p2)&& verticiesVVp[i]->idxI<STEINER_VERTEX){
					verticiesVVp[i]->IOpinp=cur_IOpin_p;
					verticiesVVp[i]->pinTypeI=IOPIN;
					break;
				}
			}
			else{
				comp* cur_comp_p=def.compsMSCp[comp_name];
				pin* pinp=cur_comp_p->cellCp->pinsMSPp[pin_name];

				point2F center_pin_P2={cur_comp_p->cellCp->pinsMSPp[pinp->nameS]->centerP2.xF,cur_comp_p->locP2.yF+cur_comp_p->cellCp->sizeP2.yF/2};
				if(centerP2.equal_p2(center_pin_P2) && verticiesVVp[i]->idxI<STEINER_VERTEX){
					verticiesVVp[i]->pinp=pinp;
					verticiesVVp[i]->pinTypeI=PIN;
					break;
				}
			}

		}
		//Test the pin if the vertex is assigned to it.
		if(verticiesVVp[i]->pinTypeI==STEINER_PIN){
			if(verticiesVVp[i]->typeI==STEINER){
				//verticiesVVp[i]->pinTypeI==STEINER_PIN;
				continue;
			}
			else{
				//This error message is printed when the coordinates from flute file do not match with the coordinate of net(in c++ code).
				cout<<"ERROR"<<endl; 
			}
		}
	}
	
	//Build graph
	graph *netGraphp= new graph(verticiesVVp,edgesVE);
	return netGraphp;

}

graph *direct_connection(Def &def, net *netp){
	//Return graph consists of only two elements.
	//map<string,string> &cellPinMSS=netp->cellPinMSS;
	unordered_map<string,string> &cellPinMSS=netp->cellPinMSS;
	
	comp* first_compp=NULL;
	pin* first_pinp=NULL;
	IOpin* first_IOpinp=NULL;
	comp* second_compp=NULL;
	pin* second_pinp=NULL;
	IOpin* second_IOpinp=NULL;
	
	bool is_first=true;
	int IOPIN_idx=0;//0 : no IO pin has found. 1: first element is IO pin, 2: second element is IO pin, 3: all elements are IO pin.
	//for(auto it=cellPinMSS.begin();it!=cellPinMSS.end();it++){
	for(auto& it:cellPinMSS){
		//if(it->first!="PIN")
		if(it.first!="PIN"){
			if(is_first){
				//first_compp=def.compsMSCp[it->first];
				//first_pinp=first_compp->cellCp->pinsMSPp[it->second]; 
				first_compp=def.compsMSCp[it.first];
				first_pinp=first_compp->cellCp->pinsMSPp[it.second]; 
				is_first=false;
			}
			else{
				//second_compp=def.compsMSCp[it->first];
				//second_pinp=second_compp->cellCp->pinsMSPp[it->second]; 
				second_compp=def.compsMSCp[it.first];
				second_pinp=second_compp->cellCp->pinsMSPp[it.second]; 
			}
		}
		else{
			if(is_first){
				//first_IOpinp=def.IOpinsMSpp[it->second];
				first_IOpinp=def.IOpinsMSpp[it.second];
				IOPIN_idx=1;
				is_first=false;
			}
			else{
				//second_IOpinp=def.IOpinsMSpp[it->second];
				second_IOpinp=def.IOpinsMSpp[it.second];
				IOPIN_idx=(IOPIN_idx==0)?2:3;
			}
		}
	}
	
	if(IOPIN_idx==3){
		graph * graphp=NULL;
		return graphp;
	}
	
	vertex *firstVp=new vertex;
	vertex *secondVp=new vertex;
	bool first_pin=true; //true if first element is pin, false if IOpin
	bool second_pin=true;
	if(IOPIN_idx==0){
		first_pin=true;
		second_pin=true;
	}
	else if (IOPIN_idx==1){
		first_pin=false;
		second_pin=true;
	}
	else if (IOPIN_idx==2){
		first_pin=true;
		second_pin=false;
	}
	
	firstVp->idxI=1;
	firstVp->typeI=CELL;
	secondVp->idxI=2;
	secondVp->typeI=CELL;
	
	if(first_pin){
		firstVp->pinTypeI=PIN;
		firstVp->pinp=first_pinp;
		firstVp->centerP2=first_pinp->centerP2;
	}
	else{
		firstVp->pinTypeI=IOPIN;
		firstVp->IOpinp=first_IOpinp;
		firstVp->centerP2=first_IOpinp->centerP2;
	}
	if(second_pin){
		secondVp->pinTypeI=PIN;
		secondVp->pinp=second_pinp;
		secondVp->centerP2=second_pinp->centerP2;
	}
	else{
		secondVp->pinTypeI=IOPIN;
		secondVp->IOpinp=second_IOpinp;
		secondVp->centerP2=second_IOpinp->centerP2;
	}
	edge edgeE(firstVp,secondVp);
	vector<vertex*> verticiesVVp={firstVp,secondVp};
	vector<edge> edgeVE={edgeE};
	graph *netGraphp=new graph(verticiesVVp,edgeVE);
	return netGraphp;
}

//v5
gcell_grid* edge_shifting(Def &def){
	gcell_grid* gcell_gridp=new gcell_grid(def);

	for(auto it_net=def.netsMSNp.begin();it_net!=def.netsMSNp.end();it_net++){
		net* netp=it_net->second;
		graph* netgp=netp->steinerG;

		if(netgp==NULL) continue;
		vertex *** &verticiesppp=netgp->verticiesppp;
		map<int,int> &verticiesMII=netgp->verticiesMII; //vertex idx to number of edges
		map<int,size_t> &verIdxIdxMII=netgp->verIdxIdxMII; // vertex idx to verticiesppp idx
		vector<edge*> edgesVEp=netgp->edgesVEp;
	
		//Collect all Steiner points from the net
		vector<steiner_container> steiner_containerV;
		for(auto it_ver=verticiesMII.begin();it_ver!=verticiesMII.end();it_ver++){
			vertex* vertexp=verticiesppp[verIdxIdxMII[it_ver->first]][0];
			int vertex_degree=it_ver->second;
			if(vertex_degree>=3 && vertexp->typeI==STEINER){
				steiner_container sc(vertexp,netgp,it_ver->first);
				steiner_containerV.push_back(sc);
			}
		}
		
		//Connect steiner container to each stainer container
		for(size_t i1=0;i1<steiner_containerV.size();i1++){
			steiner_container& sc1=steiner_containerV[i1];
			vector<vertex*> &adjacent_steiner=sc1.adjacent_steiner;
			vector<steiner_container> &adjacent_steiner_container=sc1.adjacent_steiner_container;

			for(size_t i=0;i<adjacent_steiner.size();i++){
				vertex* steinerp=adjacent_steiner[i];
				for(size_t i2=0;i2<steiner_containerV.size();i2++){
					if(i1!=i2){
						steiner_container steiner_container=steiner_containerV[i2];
						if(steinerp==steiner_container.steiner_vertexp){
							adjacent_steiner_container.push_back(steiner_container);
						}
					}
				}
			}
			assert(adjacent_steiner.size()==adjacent_steiner_container.size());
		}

		
		//Find overlapped sliding edge
		for(size_t i1=0;i1<steiner_containerV.size();i1++){
			steiner_container& sc1=steiner_containerV[i1];
			edge & sliding_edge1_0=sc1.sliding_edge[0];
			edge & sliding_edge1_1=sc1.sliding_edge[1];
			vector<steiner_container> sc1_adjacent_steiner_containerV=sc1.adjacent_steiner_container;
			for(size_t i2=0;i2<sc1_adjacent_steiner_containerV.size();i2++){
				steiner_container& sc2=sc1_adjacent_steiner_containerV[i2];
				//edge & sliding_edge2_0=sc2.sliding_edge[0];
				//edge & sliding_edge2_1=sc2.sliding_edge[1];
				for(int i_xy=0;i_xy<=1;i_xy++){
					if(i_xy==0){
						bool is_same_coord1_00=sliding_edge1_0.endPointsP22[0].equal(sc2.steiner_vertexp->centerP2);
						bool is_same_coord1_01=sliding_edge1_0.endPointsP22[1].equal(sc2.steiner_vertexp->centerP2);
						if(is_same_coord1_00 || is_same_coord1_01){
							sc1.extend_sliding_edge(sc2,i_xy);
							sc2.extend_sliding_edge(sc1,i_xy);
						}
					}
					else{
						bool is_same_coord2_10=sliding_edge1_1.endPointsP22[0].equal(sc2.steiner_vertexp->centerP2);
						bool is_same_coord2_11=sliding_edge1_1.endPointsP22[1].equal(sc2.steiner_vertexp->centerP2);
						if(is_same_coord2_10 || is_same_coord2_11){
							sc1.extend_sliding_edge(sc2,i_xy);
							sc2.extend_sliding_edge(sc1,i_xy);
						}				
					}
				}
			}
		}

		//Delete graph from gcell_grid
		gcell_gridp->insert_graph_info_into_gcell_grid(netgp,true);

		//Shift edge (edge shifting)
		for(int i_xy=0;i_xy<=1;i_xy++){
			//i_xy==0: fix y, sweep x -> sliding range is horizontal, Steiner points should have the same x-coordinate.
			//i_xy==1:                -> sliding range is vertical, Steiner points should have the same y-coordinate.
			for(size_t i1=0;i1<steiner_containerV.size();i1++){
				steiner_container& sc1=steiner_containerV[i1];
				edge & sliding_edge1_0=sc1.sliding_edge[0];
				edge & sliding_edge1_1=sc1.sliding_edge[1];
				float& down_left_sc1=(i_xy==0)?sliding_edge1_0.endPointsP22[0].xF:sliding_edge1_1.endPointsP22[0].yF; // hor edge
				float& up_right_sc1=(i_xy==0)?sliding_edge1_0.endPointsP22[1].xF:sliding_edge1_1.endPointsP22[1].yF; // ver edge
				float sc1_line_org[2]={down_left_sc1,up_right_sc1};
				bool is_steiner_not_endpoint1=!((sc1.steiner_vertexp->centerP2.equal(sliding_edge1_0.endPointsP22[0])) ||
												(sc1.steiner_vertexp->centerP2.equal(sliding_edge1_0.endPointsP22[1])) ||
												(sc1.steiner_vertexp->centerP2.equal(sliding_edge1_1.endPointsP22[0])) ||
												(sc1.steiner_vertexp->centerP2.equal(sliding_edge1_1.endPointsP22[1])));
				

				vector<steiner_container> sc1_adjacent_steiner_containerV=sc1.adjacent_steiner_container;
				for(size_t i2=0;i2<sc1_adjacent_steiner_containerV.size();i2++){
					steiner_container& sc2=sc1_adjacent_steiner_containerV[i2];
					bool is_same_xy1=(i_xy==0 && sc2.steiner_vertexp->centerP2.equalX(sc1.steiner_vertexp->centerP2));				
					bool is_same_xy2=(i_xy==1 && sc2.steiner_vertexp->centerP2.equalY(sc1.steiner_vertexp->centerP2));
					if(!(is_same_xy1 || is_same_xy2)) continue; // Although sliding edges are overlapped with each other, we skip the case that the steiner point does not align to each other. This is because all the adjacent Steiner points have same x or y coordinate. 
					//If this condition is not conserved, it occurs a problem.
					//Ex. Steiner point A and B are adjacent in y direction (same x coordinate), then, their sliding edges must be overlapped since they share the same edge. However, this is not the case that we want to move the steiner points.

					edge & sliding_edge2_0=sc2.sliding_edge[0];
					edge & sliding_edge2_1=sc2.sliding_edge[1];
					float& down_left_sc2=(i_xy==0)?sliding_edge2_0.endPointsP22[0].xF:sliding_edge2_1.endPointsP22[0].yF;
					float& up_right_sc2=(i_xy==0)?sliding_edge2_1.endPointsP22[1].xF:sliding_edge2_1.endPointsP22[1].yF;
					float sc2_line_org[2]={down_left_sc2,up_right_sc2};
					bool is_steiner_not_endpoint2=!((sc2.steiner_vertexp->centerP2.equal(sliding_edge2_0.endPointsP22[0])) ||
												(sc2.steiner_vertexp->centerP2.equal(sliding_edge2_0.endPointsP22[1])) ||
												(sc2.steiner_vertexp->centerP2.equal(sliding_edge2_1.endPointsP22[0])) ||
												(sc2.steiner_vertexp->centerP2.equal(sliding_edge2_1.endPointsP22[1])));
					
					bool xy_condition=(i_xy==0 && sc1.steiner_vertexp->centerP2.xF==sc2.steiner_vertexp->centerP2.xF) ||
									(i_xy==1 && sc1.steiner_vertexp->centerP2.yF==sc2.steiner_vertexp->centerP2.yF);		

					bool is_steiner_not_endpoint=is_steiner_not_endpoint1 && is_steiner_not_endpoint2; //We do not consider the Steiner points laid on the end of the line.
					
					if(is_lines_intersect(sc1_line_org,sc2_line_org) && is_steiner_not_endpoint && xy_condition){
						//Find common range including extended sliding ranges
						float common_min_value=(i_xy==0)?max(sc1.extended_sliding_range[0].endPointsP22[0].xF,sc2.extended_sliding_range[0].endPointsP22[0].xF):max(sc1.extended_sliding_range[1].endPointsP22[0].yF,sc2.extended_sliding_range[1].endPointsP22[0].yF);
						float common_max_value=(i_xy==0)?min(sc1.extended_sliding_range[0].endPointsP22[1].xF,sc2.extended_sliding_range[0].endPointsP22[1].xF):min(sc1.extended_sliding_range[1].endPointsP22[1].yF,sc2.extended_sliding_range[1].endPointsP22[1].yF);
	
						point2F tmp_p1_1={common_min_value,min(sc1.steiner_vertexp->centerP2.yF,sc2.steiner_vertexp->centerP2.yF)};
						point2F tmp_p1_2={min(sc1.steiner_vertexp->centerP2.xF,sc2.steiner_vertexp->centerP2.xF),common_min_value};
						point2F tmp_p2_1={common_max_value,max(sc1.steiner_vertexp->centerP2.yF,sc2.steiner_vertexp->centerP2.yF)};
						point2F tmp_p2_2={max(sc1.steiner_vertexp->centerP2.xF,sc2.steiner_vertexp->centerP2.xF),common_max_value};
						
						point2F start_pointPF=(i_xy==0)?tmp_p1_1:tmp_p1_2;
						point2F end_pointPF=(i_xy==0)?tmp_p2_1:tmp_p2_2;
						
						edge min_edge=gcell_gridp->search_min_demand_edge(start_pointPF,end_pointPF,i_xy);

						point2F min_start_pointPF=min_edge.endPointsP22[0];
						point2F min_end_pointPF=min_edge.endPointsP22[1];
						
						assert((i_xy==0 && min_start_pointPF.equalX(min_end_pointPF)) || (i_xy==1 && min_start_pointPF.equalY(min_end_pointPF)));

						bool is_start_org_range1=(i_xy==0)?min_start_pointPF.xF>=sc1_line_org[0]:min_start_pointPF.yF>=sc1_line_org[0];
						bool is_start_org_range2=(i_xy==0)?min_start_pointPF.xF>=sc2_line_org[0]:min_start_pointPF.yF>=sc2_line_org[0];
						bool is_end_org_range1=(i_xy==0)?min_end_pointPF.xF<=sc1_line_org[1]:min_end_pointPF.yF<=sc1_line_org[1];
						bool is_end_org_range2=(i_xy==0)?min_end_pointPF.xF<=sc2_line_org[1]:min_end_pointPF.yF<=sc2_line_org[1];

						assert(is_start_org_range1 || is_end_org_range1);
						assert(is_start_org_range2 || is_end_org_range2);
						
						int changed_sp_idx1=sc1.vertex_idx;
						int changed_sp_idx2=sc2.vertex_idx;
						
						if(!is_start_org_range1 || !is_start_org_range2 || !is_end_org_range1 || !is_end_org_range2){
							//map<int,edge*> extended_sliding_edgeM=(!is_start_org_range1 || !is_end_org_range1)?sc1.extended_sliding_edgeM[0]:sc2.extended_sliding_edgeM[0];
							map<int,edge> extended_sliding_edgeM=(!is_start_org_range1 || !is_end_org_range1)?sc1.extended_sliding_edgeM[0]:sc2.extended_sliding_edgeM[0];
							for(auto it=extended_sliding_edgeM.begin();it!=extended_sliding_edgeM.end();it++){
								//float &edge_start=(i_xy==0)?it->second->endPointsP22[0].xF:it->second->endPointsP22[0].yF;
								//float &edge_end=(i_xy==0)?it->second->endPointsP22[1].xF:it->second->endPointsP22[1].yF;
								float &edge_start=(i_xy==0)?it->second.endPointsP22[0].xF:it->second.endPointsP22[0].yF;
								float &edge_end=(i_xy==0)?it->second.endPointsP22[1].xF:it->second.endPointsP22[1].yF;
								//if(!(is_start_org_range1 && is_end_org_range1)){
								if(!is_start_org_range1 || !is_end_org_range1){
									//if the determined sliding edge is beyond the small range, change the connection of the steiner points.
									bool x_cond=(i_xy==0 && edge_start<=min_start_pointPF.xF && edge_end>=min_start_pointPF.xF);
									bool y_cond=(i_xy==1 && edge_start<=min_start_pointPF.yF && edge_end>=min_start_pointPF.yF);
									if(x_cond || y_cond){
										netgp->change_edge(sc1.vertex_idx,it->first,sc2.vertex_idx);
									}
									changed_sp_idx1=it->first;

								}
								//if(!(is_start_org_range2 && is_end_org_range2)){
								if(!is_start_org_range2 || !is_end_org_range2){
									bool x_cond=(i_xy==0 && edge_start<=min_start_pointPF.xF && edge_end>=min_start_pointPF.xF);
									bool y_cond=(i_xy==1 && edge_start<=min_start_pointPF.yF && edge_end>=min_start_pointPF.yF);
									if(x_cond || y_cond){
										netgp->change_edge(sc2.vertex_idx,it->first,sc1.vertex_idx);
									}
									changed_sp_idx2=it->first;

								}


							}
						}
						//Move Steiner Points
						point2F sp1_new_coord=sc1.steiner_vertexp->centerP2;
						point2F sp2_new_coord=sc2.steiner_vertexp->centerP2;
						if(i_xy==0){
							//if i_xy==0 , both xFs of min_start_pointPF and max_start_pointPF are the same.
							sp1_new_coord.xF=min_start_pointPF.xF;
							sp2_new_coord.xF=min_start_pointPF.xF;
						}
						else{
							sp1_new_coord.yF=min_start_pointPF.yF;
							sp2_new_coord.yF=min_start_pointPF.yF;
						}
						netgp->move_vertex(changed_sp_idx1,sp1_new_coord);
						netgp->move_vertex(changed_sp_idx2,sp2_new_coord);

					}
				}
			}
		}

		//Set the shape of diagonal edges
		int dummy_idx=1;

		for(size_t i=0;i<edgesVEp.size();i++){
			edge* edgep=edgesVEp[i];
			point2F start_point=edgep->endPointsP22[0];
			point2F end_point=edgep->endPointsP22[1];
			if(!start_point.equalX(end_point) && !start_point.equalY(end_point)){
				//find suitable location
				point2F mid_point1={start_point.xF,end_point.yF};
				point2F mid_point2={end_point.xF,start_point.yF};
				
				float demand_case1=0;
				float demand_case2=0;
				
				//case1
				demand_case1+=gcell_gridp->search_demand_for_line(start_point,mid_point1);
				demand_case1+=gcell_gridp->search_demand_for_line(mid_point1,end_point);

				//case2
				demand_case2+=gcell_gridp->search_demand_for_line(start_point,mid_point2);
				demand_case2+=gcell_gridp->search_demand_for_line(mid_point2,end_point);

				//add dummy vertex
				vertex *newVp=new vertex;
				vector<vertex*> connected_verticies;
				newVp->idxI=DUMMY_VERTEX+dummy_idx++;
				newVp->typeI=DUMMY;
				newVp->pinTypeI=DUMMY_PIN;
				if(demand_case1<=demand_case2){
					newVp->centerP2=mid_point1;
				}
				else{
					newVp->centerP2=mid_point2;
				}
				connected_verticies.push_back(edgep->endVerticies[0]);
				connected_verticies.push_back(edgep->endVerticies[1]);

				netgp->delete_edge_and_add_vertex(edgep,newVp,connected_verticies);
			}
		}
		//Reassign graph from gcell_grid
		gcell_gridp->insert_graph_info_into_gcell_grid(netgp);
	}
	return gcell_gridp;
}

void split_steiner_tree(Def &def){
	//Split Steiner trees and set a target for each pin
	for(auto it_net=def.netsMSNp.begin();it_net!=def.netsMSNp.end();it_net++){
		net* netp=it_net->second;
		graph* netgp=netp->steinerG;
		
		if(netgp==NULL) continue;
		vector<edge*> &edgesVEp=netgp->edgesVEp;
		for(size_t i=0; i<edgesVEp.size();i++){
			vertex* start_vertex=edgesVEp[i]->endVerticies[0];
			int &start_idxI=start_vertex->idxI;
			int &start_typeI=start_vertex->typeI;
			int &start_vertypeI=start_vertex->pinTypeI;
			pin *&start_pinp=start_vertex->pinp;
			IOpin *&start_IOpinp=start_vertex->IOpinp;
			point2F &start_center=start_vertex->centerP2;

			vertex* end_vertex=edgesVEp[i]->endVerticies[1];
			int &end_idxI=end_vertex->idxI;
			int &end_typeI=end_vertex->typeI;
			int &end_vertypeI=end_vertex->pinTypeI;
			pin *&end_pinp=end_vertex->pinp;
			IOpin *&end_IOpinp=end_vertex->IOpinp;
			point2F &end_center=end_vertex->centerP2;
			
			//Set target
			if(start_vertypeI==PIN && end_vertypeI!=PIN){
				if(!(start_pinp->targetTypeI==PIN)){//PIN-to-PIN is major priority. If the target of the pin has been assigned to pin, then though there is a segment closer than the pin, we ignore that.
					bool set_start_target=false;
					
					if(start_pinp->targetP2==NULL) set_start_target=true;
					else if(distance_sq(start_center,end_center)<distance_sq(start_center,*start_pinp->targetP2)) set_start_target=true;
					
					if(set_start_target){
						start_pinp->targetTypeI=end_vertypeI;
						start_pinp->target=end_vertex;
						start_pinp->targetP2=&end_center;
					}
				}
			}
			else if(start_vertypeI!=PIN && end_vertypeI==PIN){
				if(!(end_pinp->targetTypeI==PIN)){
					bool set_end_target=false;
					if(end_pinp->targetP2==NULL) set_end_target=true;
					else if(distance_sq(end_center,start_center)<distance_sq(end_center,*end_pinp->targetP2)) set_end_target=true;
					
					if(set_end_target){
						end_pinp->targetTypeI=start_vertypeI;
						end_pinp->target=start_vertex;
						end_pinp->targetP2=&start_center;
					}
				}
			}
			else if(start_vertypeI==PIN && end_vertypeI==PIN){
				bool set_start_target=false;
				bool set_end_target=false;

				if(start_pinp->targetP2==NULL) set_start_target=true;
				else if(distance_sq(start_center,end_center)<distance_sq(start_center,*start_pinp->targetP2)) set_start_target=true;

				if(end_pinp->targetP2==NULL) set_end_target=true;
				else if(distance_sq(end_center,start_center)<distance_sq(end_center,*end_pinp->targetP2)) set_end_target=true;

				if(set_start_target){
					start_pinp->targetTypeI=end_vertypeI;
					start_pinp->target=end_pinp;
					start_pinp->targetP2=&end_pinp->centerP2;
				}
				if(set_end_target){
					end_pinp->targetTypeI=start_vertypeI;
					end_pinp->target=start_pinp;
					end_pinp->targetP2=&start_pinp->centerP2;
				}
			}
		}
	}
}