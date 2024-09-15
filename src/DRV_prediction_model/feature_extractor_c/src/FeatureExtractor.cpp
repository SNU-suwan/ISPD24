#include "FeatureExtractor.h"
#include "global.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <string>
#include <chrono>

#define NAN_M2_PITCH 64

void FeatureExtractor::accumulate_unfriendly(std::unordered_map<std::string, int>& used_type, std::unordered_map<std::string, int>& drv_type){
    auto start = std::chrono::high_resolution_clock::now();
    init();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Init : " << duration.count() << " ms" << std::endl;

    // calculate the number of cells of each type
    for (const auto& comp : def.compsMSCp) {
        std::string comp_type = comp.second->cellCp->nameS;
        if (used_type.count(comp_type) == 0) {
            used_type[comp_type] = 1;
        }
        else used_type[comp_type]++;
    }
    std::vector<std::string> intsec_cell_list;

    // calculate the number of cells of each type on drv
    for (auto& adrv : drv.drvV) {
        if(adrv.typeS == "Polygon not rectangle") continue;
        BBox<int> drv_bbox = convert_int3(adrv.bbox) - coreI[0];
        BBox<int> drv_grid = convert_to_grid(drv_bbox, gridXI, gridYI);
        for (int x = std::max(drv_grid.lb.x, 0); x <= std::min(drv_grid.rt.x, gridNum.x - 1); x++) {
            for (int y = std::max(drv_grid.lb.y, 0); y <= std::min(drv_grid.rt.y, gridNum.y - 1); y++) {
                assert(x >= 0 && x < gridNum.x);
                assert(y >= 0 && y < gridNum.y);
                for (const auto& cell_name : cellMap[y][x]) {
                    if (std::find(intsec_cell_list.begin(), intsec_cell_list.end(), cell_name) == intsec_cell_list.end()) {
                        intsec_cell_list.push_back(cell_name);
                    }
                }
            }
        }
    }
    for (const auto& comp_name : intsec_cell_list) {
        std::string comp_type = def.compsMSCp[comp_name]->cellCp->nameS;
        if (drv_type.count(comp_type) == 0) {
            drv_type[comp_type] = 1;
        }
        else drv_type[comp_type]++;
    }
}

void FeatureExtractor::run() {

    auto start = std::chrono::high_resolution_clock::now();
    init();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Init : " << duration.count() << " ms" << std::endl;


    if (!drv.drvDirS.empty()) {
        start = end;
        check_drv(false);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Label : " << duration.count() << " ms" << std::endl;
    }
    
	start = end;
	std::cout << "Parsing congestion map";
	congestion_map();
		
	end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << " : " << duration.count() << " ms" << std::endl;

/*
    start = end;
    m2_short(false);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "M2 short Label : " << duration.count() << " ms" << std::endl;
    */
    
    
    start = end;
    pin_rudy();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "PIN_RUDY : " << duration.count() << " ms" << std::endl;
    
    start = end;
    rudy_map();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "RUDY : " << duration.count() << " ms" << std::endl;
    
   
    start = end;
    hv_net_density();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "HV net density : " << duration.count() << " ms" << std::endl;
    
    
    start = end;
    avg_pin_access();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "pin access : " << duration.count() << " ms" << std::endl;
    
    start = end;
    pin_density();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "pin density : " << duration.count() << " ms" << std::endl;
    
    if (patSize != gcellSize) {
        start = end;
        pin_pattern();
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "pin pattern : " << duration.count() << " ms" << std::endl;
    }
    

    start = end;
    minimum_proximity();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "min proximity : " << duration.count() << " ms" << std::endl;
    
  
    start = end;
    weighted_unfriendly();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "unfriendly : " << duration.count() << " ms" << std::endl;
    
   
    start = end;
    local_global_self_crossing_net();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Local/global/self-crossing : " << duration.count() << " ms" << std::endl;
    
	start = end;
    macro_density();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Macro density : " << duration.count() << " ms" << std::endl;
    
	start = end;
    routing_capacity();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Routing capacity : " << duration.count() << " ms" << std::endl;

    if(is_ugnet){
        start = end;
        generate_graph();
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Graph : " << duration.count() << " ms" << std::endl;
    }

}

void FeatureExtractor::init() {
    
    // set die
    die.lb.x = (std::round(static_cast<double>(def.die_sizeP22[0].xF) * 1000)) / 1000.0;
    die.lb.y = (std::round(static_cast<double>(def.die_sizeP22[0].yF) * 1000)) / 1000.0;
    die.rt.x = (std::round(static_cast<double>(def.die_sizeP22[1].xF) * 1000)) / 1000.0;
    die.rt.y = (std::round(static_cast<double>(def.die_sizeP22[1].yF) * 1000)) / 1000.0;
 
    // set core
    core.lb = die.lb + def.rowsVp[0]->rowX;
    core.rt = die.rt - def.rowsVp[0]->rowX;
    // manage integer die, core
    dieI = convert_int3(die);
    coreI = convert_int3(core);

    // cal grid area
    patSizeI = std::round(patSize * 1000);
    gcellSizeI = std::round(gcellSize * 1000);
    patXI = std::round(patX * 1000);
    patYI = std::round(patY * 1000);
    gridXI = std::round(gridX * 1000);
    gridYI = std::round(gridY * 1000);

    patArea = patX * patY;
    gridArea = gridX * gridY;
    patAreaI = patXI * patYI;
    gridAreaI = gridXI * gridYI;

    // calculate patNum and gridNum
    auto calNum = [&coreI = coreI](int grid_size_x, int grid_size_y) {
        int xNum = coreI.xWidth() / grid_size_x;
        if (coreI.xWidth() % grid_size_x != 0) xNum++;
        int yNum = coreI.yWidth() / grid_size_y;
        if (coreI.yWidth() % grid_size_y != 0) yNum++;
        return Point<int>(xNum, yNum);
    };
    patNum = calNum(patXI, patYI);
    gridNum = calNum(gridXI, gridYI);
    // create cell map
    cellMap.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) cellMap[i].resize(gridNum.x);

    //std::cout << gridNum << std::endl;
    for (const auto& compM : def.compsMSCp) {
        auto comp_p = compM.second;
        auto cell_p = comp_p->cellCp;
        BBox<double> comp_bbox(comp_p->locP2.xF, comp_p->locP2.yF, comp_p->locP2.xF + cell_p->sizeP2.xF, comp_p->locP2.yF + cell_p->sizeP2.yF);
        BBox<int> comp_bboxI = convert_int3(comp_bbox);
        comp_bboxI = comp_bboxI - coreI[0];
        BBox<int> bbox_grid = convert_to_grid(comp_bboxI, gridXI, gridYI);
        for (int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++) {
            for (int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++) {
                //std::cout << "(" << x << ", " << y << ")" << std::endl;
                cellMap[y][x].push_back(compM.first);
            }
        }        
    }

	row_height=def.rowsVp[1]->rowHeightF-def.rowsVp[0]->rowHeightF;
}

void FeatureExtractor::printInfo() {
    std::cout << "patSize : " << patSize << ", gcellSize : " << gcellSize << std::endl;
    std::cout << "die : " << die << ", core : " << core << std::endl;
    std::cout << "patNum : " << patNum << ", gridNum : " << gridNum << std::endl;
}

void FeatureExtractor::print_def() {

    int count = 1;

    for (auto it = def.compsMSCp.cbegin(); it != def.compsMSCp.cend(); ++it) {
        auto comp_pointer = it->second;
        std::cout << "name : " << comp_pointer->nameS << ", dir : " << comp_pointer->dirS << std::endl;
        std::cout << "loc : (" << comp_pointer->locP2.xF << ", " << comp_pointer->locP2.yF << ")" << std::endl;

        auto cell_pointer = comp_pointer->cellCp;
        std::cout << "name : " << cell_pointer->nameS << ", dir : " << cell_pointer->dirS << std::endl;
        std::cout << "origin : (" << cell_pointer->originP2p->xF << ", " << cell_pointer->originP2p->yF << ")" << std::endl;
        std::cout << "size : (" << cell_pointer->sizeP2.xF << ", " << cell_pointer->sizeP2.yF << ")" << std::endl;
    

        for (auto pt = cell_pointer->pinsMSPp.cbegin(); pt != cell_pointer->pinsMSPp.cend(); ++pt) {
            std::cout << "pin name : " << pt->first;
            for (auto seg : pt->second->segmentsVp) {
                for (auto p : seg->pointsVp) {
                    std::cout << " (" << p->xF << ", " << p->yF << "), ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;

        }

        if (count++ > 10) break;
    }
}

void FeatureExtractor::pin_density() {
    pinDensity = cal_pin_density(gridNum, gridXI, gridYI);
    //convert_to_33(pinDensity[0], pinDensity33MinM1, pinDensity33MaxM1);
    //convert_to_33(pinDensity[1], pinDensity33MinM2, pinDensity33MaxM2);
    write_csv("pin_density_M1", pinDensity[0]);
    write_csv("pin_density_M2", pinDensity[1]);
    //write_csv("pin_density_M1_min", pinDensity33MinM1);
    //write_csv("pin_density_M1_max", pinDensity33MaxM1);
    //write_csv("pin_density_M2_min", pinDensity33MinM2);
    //write_csv("pin_density_M2_max", pinDensity33MaxM2);
}

void FeatureExtractor::pin_pattern() {
    pinPattern = cal_pin_density(patNum, patXI, patYI);
    write_csv("pin_pattern_M1", pinPattern[0]);
    write_csv("pin_pattern_M2", pinPattern[1]);
}
FeatureExtractor::Grid3D FeatureExtractor::cal_pin_density(Point<int> grid_num, int grid_x, int grid_y) {
    
    Grid3D pin_density;
    
    // initialize pin density matrix
    pin_density.clear();
    
    pin_density.resize(2);
    for (int k = 0; k < 2; k++) {
        pin_density[k].resize(grid_num.y);
        for (int i = 0; i < grid_num.y; i++) {
            pin_density[k][i].resize(grid_num.x, 0);
        }
    }

    // Calculate features
    for (const auto& compM : def.compsMSCp) {
        auto comp_p = compM.second;
        auto cell_p = comp_p->cellCp;
        //Point<double> cellLoc_d(comp_p->locP2.xF, comp_p->locP2.yF);
        //Point<int> cellLoc = convert_int3(cellLoc_d);

        for (const auto& pinM : cell_p->pinsMSPp) {
            if (pinM.first == "VDD" || pinM.first == "VSS") continue;

            std::vector<BBox<int>> rect_list_M1;
            std::vector<BBox<int>> rect_list_M2;
            for (auto seg : pinM.second->segmentsVp){
                BBox<double> seg_bbox_d((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
                BBox<int> seg_bbox = convert_int3(seg_bbox_d);
                seg_bbox = seg_bbox - coreI[0];
                if (seg->layerS == "M1") rect_list_M1.push_back(seg_bbox);
                else rect_list_M2.push_back(seg_bbox);
                
            } 

            auto make_grid_list = [&](vector<BBox<int>>& rect_list) {
                // reduce rect_list
                std::vector<int> hanan_x, hanan_y;
                for (const auto& bbox : rect_list) {
                    if (std::find(hanan_x.begin(), hanan_x.end(), bbox.lb.x) == hanan_x.end())
                        hanan_x.push_back(bbox.lb.x);
                    if (std::find(hanan_x.begin(), hanan_x.end(), bbox.rt.x) == hanan_x.end())
                        hanan_x.push_back(bbox.rt.x);
                    if (std::find(hanan_y.begin(), hanan_y.end(), bbox.lb.y) == hanan_y.end())
                        hanan_y.push_back(bbox.lb.y);
                    if (std::find(hanan_y.begin(), hanan_y.end(), bbox.rt.y) == hanan_y.end())
                        hanan_y.push_back(bbox.rt.y);
                }
                std::sort(hanan_x.begin(), hanan_x.end());
                std::sort(hanan_y.begin(), hanan_y.end());            

                std::vector<BBox<int>> grid_list;
                int x_size = hanan_x.size();
                int y_size = hanan_y.size();
                for (int i = 0; i < x_size - 1; i++) {
                    for (int j = 0; j < y_size - 1; j++) {
                        BBox<int> grid(hanan_x[i], hanan_y[j], hanan_x[i + 1], hanan_y[j + 1]);
                        for (const auto& rect_bbox : rect_list) {
                            if ((grid & rect_bbox) == grid) {
                                grid_list.push_back(grid);
                                break;
                            }
                        }
                    }
                }
                return grid_list;
            };
            
            std::vector<BBox<int>> grid_list_M1 = make_grid_list(rect_list_M1);
            std::vector<BBox<int>> grid_list_M2 = make_grid_list(rect_list_M2);

            int grid_area = grid_x * grid_y;

            // add pin density to grid
            for (auto& bbox : grid_list_M1) {
                BBox<int> bbox_grid = convert_to_grid(bbox, grid_x, grid_y);
                //std::cout << bbox << ", " << bbox_grid << std::endl;
                for (int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++) {
                    for (int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++) {
                        BBox<int> grid(x * grid_x, y * grid_y, (x + 1) * grid_x, (y + 1) * grid_y);
                        int intersect = (bbox & grid).getArea();
                        pin_density[0][grid_num.y - y - 1][x] += static_cast<double>(intersect) / (grid_area);
                    }
                }
            }
            for (auto& bbox : grid_list_M2) {
                BBox<int> bbox_grid = convert_to_grid(bbox, grid_x, grid_y);
                //std::cout << bbox << ", " << bbox_grid << std::endl;
                for (int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++) {
                    for (int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++) {
                        BBox<int> grid(x * grid_x, y * grid_y, (x + 1) * grid_x, (y + 1) * grid_y);
                        int intersect = (bbox & grid).getArea();
                        pin_density[1][grid_num.y - y - 1][x] += static_cast<double>(intersect) / (grid_area);
                    }
                }
            }
        }                
    }

    return pin_density;

}

void FeatureExtractor::avg_pin_access() {
    Grid2D pinCount, pinAccessCount;
    double pitch = 0.064f;
    int pitchI = 64;

    pinAccessInfo.clear();
    pinCount.clear();
    pinAccessCount.clear();
    pinAccessInfo.resize(gridNum.y);
    pinCount.resize(gridNum.y);
    pinAccessCount.resize(gridNum.y);

    for(int i = 0; i < gridNum.y; i++){
        pinAccessInfo[i].resize(gridNum.x, 0);
        pinCount[i].resize(gridNum.x, 0);
        pinAccessCount[i].resize(gridNum.x, 0);
    }

    for(auto it = def.compsMSCp.cbegin(); it != def.compsMSCp.cend(); ++it){
        auto comp_p = it->second;
        auto stdcell_p = comp_p->cellCp;
        Point<double> cellLoc(comp_p->locP2.xF, comp_p->locP2.yF);
        for(auto pt = stdcell_p->pinsMSPp.cbegin(); pt != stdcell_p->pinsMSPp.cend(); ++pt){
            std::string pin_name = pt->first;
            if(pin_name == "VDD" || pin_name == "VSS") continue;

            auto segments = pt->second;
            auto sep_bbox = seperate_bbox(segments);
            auto sep_bbox_M1 = sep_bbox.first;
            auto sep_bbox_M2 = sep_bbox.second;
            vector<pair<int, int>> pin_exist;
            pin_exist.clear();
            

            for(auto bboxes : sep_bbox_M1){
                for(auto single_bbox : bboxes){
                    BBox<int> single_bbox_coreI = convert_int3(single_bbox) - coreI[0];
                    BBox<int> bbox_grid = convert_to_grid(single_bbox_coreI, gridXI, gridYI);

                    int num_access;

                    num_access = std::floor((single_bbox_coreI.rt.y - bbox_grid.lb.y * gridYI) / pitchI) - std::floor((single_bbox_coreI.lb.y - bbox_grid.lb.y * gridYI) / pitchI);

                    if(num_access > 0){
                        for(int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
                            pinAccessCount[gridNum.y - bbox_grid.lb.y - 1][x] += num_access;
                            if(find(pin_exist.begin(), pin_exist.end(), pair<int, int>(x, bbox_grid.lb.y)) == pin_exist.end()){
                                pin_exist.push_back(pair<int, int>(x, bbox_grid.lb.y));
                            }
                            
                        }
                    }
                }                    
            }

            for(auto bboxes : sep_bbox_M2){
                for(auto single_bbox : bboxes){
                    BBox<int> single_bbox_coreI = convert_int3(single_bbox) - coreI[0];
                    BBox<int> bbox_grid = convert_to_grid(single_bbox_coreI, gridXI, gridYI);
                    int num_access;
                    for(int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
                        num_access = std::floor((std::min(single_bbox_coreI.rt.x, (x + 1) * gridXI) - x * gridXI) / pitchI) - std::floor((std::max(single_bbox_coreI.lb.x, x * gridXI) - x * gridXI) / pitchI);
                        if(num_access != 0){
                            pinAccessCount[gridNum.y - bbox_grid.lb.y - 1][x] += num_access;
                            if(find(pin_exist.begin(), pin_exist.end(), pair<int, int>(x, bbox_grid.lb.y)) == pin_exist.end()){
                                pin_exist.push_back(pair<int, int>(x, bbox_grid.lb.y));
                            }
                        }
                    }
                }
            }
            for(auto grid_coordinate : pin_exist){
                pinCount[gridNum.y - grid_coordinate.second - 1][grid_coordinate.first] ++;
            }

        }
    }
    for(int x = 0; x < gridNum.x; x++){
        for(int y = 0; y < gridNum.y; y++){
            pinAccessCount[gridNum.y - y - 1][x] -= calculate_overlap(x, y);
            if(pinCount[gridNum.y - y - 1][x] != 0)
                pinAccessInfo[gridNum.y - y - 1][x] = pinAccessCount[gridNum.y - y - 1][x] / pinCount[gridNum.y - y - 1][x];
        }
    }
   //convert_to_33(pinAccessInfo, pinAccessInfo33Min, pinAccessInfo33Max);
    write_csv("pin_access_info", pinAccessInfo);
    //write_csv("pin_access_info_min", pinAccessInfo33Min);
    //write_csv("pin_access_info_max", pinAccessInfo33Max);
}

void FeatureExtractor::rudy_map() {
    double metalWidth = 0.028f;
    shortRudy.clear();
    longRudy.clear();
    shortRudy.resize(gridNum.y);
    longRudy.resize(gridNum.y);
    for(int i = 0; i < gridNum.y; i++){
        shortRudy[i].resize(gridNum.x, 0);
        longRudy[i].resize(gridNum.x, 0);
    }
    for(auto it = def.netsMSNp.cbegin(); it != def.netsMSNp.cend(); ++it){
        auto net_name = it->first;
        auto net_p = it->second;
        auto comp_list = net_p->cellPinMSS;
        BBox<double> net_bbox(0, 0, 0, 0);
        BBox<double> IOpin_bbox(0, 0, 0, 0);
        double bbox_size;
        double rudy_value;

        for(auto it2 = comp_list.cbegin(); it2 != comp_list.cend(); ++it2){
            string comp_name = it2->first;
            string pin_name = it2->second;
            if(comp_name == "PIN"){
                double x = def.IOpinsMSpp[pin_name]->centerP2.xF;
                double y = def.IOpinsMSpp[pin_name]->centerP2.yF;

                /*if(x < core[0][0]) x = core[0][0];
                else if(x > core[1][0]) x = core[1][0];
                if(y < core[0][1]) y = core[0][1];
                else if(y > core[1][1]) y = core[1][1];*/
                IOpin_bbox = BBox<double>(x, y, x, y);
                if(net_bbox.lb.x == 0) net_bbox = IOpin_bbox;
                else net_bbox = net_bbox + IOpin_bbox;
                continue;
            }
			bool isMacro=false;
			for(auto it_macro=def.comp_macrosMSCp.begin();it_macro!=def.comp_macrosMSCp.end();it_macro++){
				string macro_name=it_macro->first;
				if(macro_name==comp_name){
					isMacro=true;
					break;
				}
			}
			pin* single_pin;
			if(!isMacro){
				single_pin = def.compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name];
			}
			else{
				single_pin = def.comp_macrosMSCp[comp_name]->macroMp->pinsMSPp[pin_name];
			}
            for(auto seg : single_pin->segmentsVp){
                BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
                if(net_bbox.lb.x == 0) net_bbox = seg_bbox;
                else net_bbox = net_bbox + seg_bbox;
            }
        }
        
        double w = (net_bbox.rt.x - net_bbox.lb.x);
        double h = (net_bbox.rt.y - net_bbox.lb.y);
        bbox_size = w + h;
        rudy_value = (w + h - metalWidth) * metalWidth / (w * h);

        net_bbox = net_bbox - core[0];
        BBox<int> bbox_grid = convert_to_grid(net_bbox, gridX, gridY);


        for(int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
            if(x < 0 || x >= gridNum.x) continue;
            for(int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++){
                if(y < 0 || y >= gridNum.y) continue;
                BBox<double> bbox(x * gridX, y * gridY, (x+1) * gridX, (y+1) * gridY);
                BBox<double> bbox_intersect = bbox & net_bbox;

                
                double intersect_ratio = bbox_intersect.getArea() / bbox.getArea();
                if(bbox_size < rudy_threshold * gcellSize) shortRudy[gridNum.y - y - 1][x] += (rudy_value * intersect_ratio);
                else longRudy[gridNum.y - y - 1][x] += (rudy_value * intersect_ratio);
            }
        }
    }
    //convert_to_33(shortRudy, shortRudy33Min, shortRudy33Max);
    //convert_to_33(longRudy, longRudy33Min, longRudy33Max);
    write_csv("short_rudy", shortRudy);
    write_csv("long_rudy", longRudy);
    //write_csv("short_rudy_min", shortRudy33Min);
    //write_csv("short_rudy_max", shortRudy33Max);
    //write_csv("long_rudy_min", longRudy33Min);
    //write_csv("long_rudy_max", longRudy33Max);
}

void FeatureExtractor::check_drv(bool check_num){
    int sum = 0;
    numDrv.clear();
    pinDRV.clear();
    congDRV.clear();
    numDrv.resize(gridNum.y);
    pinDRV.resize(gridNum.y);
    congDRV.resize(gridNum.y);
    for(int i = 0; i < gridNum.y; i++){
        numDrv[i].resize(gridNum.x);
        pinDRV[i].resize(gridNum.x);
        congDRV[i].resize(gridNum.x);
    }
    for(auto& single_drv : drv.drvV){
        string drv_name = single_drv.typeS;
        if(drv_name == "Polygon not rectangle") continue;
        
        BBox<double> bbox_core = single_drv.bbox - core[0];
        BBox<int> bbox_coreI = convert_int3(bbox_core);
        BBox<int> bbox_grid = convert_to_grid(bbox_coreI, gridXI, gridYI);
        for(int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
            if(x < 0 || x >= gridNum.x) continue;
            for(int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++){
                if(y < 0 || y >= gridNum.y) continue;
                numDrv[gridNum.y - y - 1][x]++;
                if (single_drv.layerS == "M1" || single_drv.layerS == "MINT1") pinDRV[gridNum.y - y - 1][x]++;
                else congDRV[gridNum.y - y - 1][x]++;
            }
        }
    }
    write_csv("drv_reg", numDrv);
    if(!check_num){
        for(int x = 0; x < gridNum.x; x++){
            for(int y = 0; y < gridNum.y; y++){
                if(numDrv[y][x] > 1) numDrv[y][x] = 1;
                if(pinDRV[y][x] > 1) pinDRV[y][x] = 1;
                if(congDRV[y][x] > 1) congDRV[y][x] = 1;
            }
        }
        write_csv("drv_check", numDrv);
    }
    //write_csv("drv_check", numDrv);
    //write_csv("pin_drv", pinDRV);
    //write_csv("cong_drv", congDRV);
}


void FeatureExtractor::m2_short(bool check_num){
    int sum = 0;
    M2Short.clear();
    M2Short.resize(gridNum.y);
    for(int i = 0; i < gridNum.y; i++){
        M2Short[i].resize(gridNum.x);
    }
    for(auto& single_drv : drv.drvV){
        
        //std::cout << single_drv.typeS << " " << single_drv.layerS << std::endl;
        if (!(single_drv.typeS == "Short" && single_drv.layerS == "MINT1")) continue;
        
        BBox<double> bbox_core = single_drv.bbox - core[0];
        BBox<int> bbox_coreI = convert_int3(bbox_core);
        BBox<int> bbox_grid = convert_to_grid(bbox_coreI, gridXI, gridYI);
        for(int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
            if(x < 0 || x >= gridNum.x) continue;
            for(int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++){
                if(y < 0 || y >= gridNum.y) continue;
                M2Short[gridNum.y - y - 1][x] ++;
            }
        }
    }
    if(!check_num){
        for(int x = 0; x < gridNum.x; x++){
            for(int y = 0; y < gridNum.y; y++){
                if(M2Short[y][x] > 1) M2Short[y][x] = 1;
            }
        }
    }
    write_csv("M2_short", M2Short);
}


void FeatureExtractor::minimum_proximity() {
    // initialize pin density matrix
    minProximity.clear();
    minProximity.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        minProximity[i].resize(gridNum.x, 0);
    }    

    // calculate the feature
    for (int x = 0; x < gridNum.x; x++) {
        for (int y = 0; y < gridNum.y; y++) {
            //std::cout << "(" << gridNum.y - y - 1 << ", " << x << ")" << std::endl;
            BBox<int> grid_bbox(x * gridXI, y * gridYI, (x + 1) * gridXI, (y + 1) * gridYI);

            std::vector<BBox<int>> pin_bbox_list;

            if (cellMap[y][x].size() == 0) {
                minProximity[gridNum.y - y - 1][x] = static_cast<double>(gridXI + gridYI) / 1000;
                continue;
            }
            for (auto comp_name : cellMap[y][x]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // pin 별로 bounding box를 만들어봐서 grid와 겹치면 넣는다
                for (const auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;

                    BBox<int> pin_bbox;
                    bool isFirst = true;
                    for (auto seg : pinM.second->segmentsVp) {
                        BBox<double> seg_bbox_d((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
                        BBox<int> seg_bbox = convert_int3(seg_bbox_d);
                        seg_bbox = seg_bbox - coreI[0];
                        if (isFirst) {
                            pin_bbox = seg_bbox;
                            isFirst = false;
                        } 
                        else pin_bbox = pin_bbox + seg_bbox;
                    }

                    BBox<int> inter_bbox = grid_bbox & pin_bbox;
                    if (inter_bbox.lb.x != -1 && inter_bbox.lb.y != -1 && inter_bbox.rt.x != -1 && inter_bbox.rt.y != -1) {
                        pin_bbox_list.push_back(pin_bbox);
                    }
                }
            }

            //int min_hpwl = 2 * gridSizeI;
            int max_hpwl = gridXI + gridYI;
            int min_hpwl = gridXI + gridYI;
            int num_pin = pin_bbox_list.size();
            for (int i = 0; i < num_pin - 1; i++) {
                for (int j = i + 1; j < num_pin; j++) {
                    BBox<int> two_pin_bbox = (pin_bbox_list[i] + pin_bbox_list[j]) & grid_bbox;
                    int hpwl = two_pin_bbox.xWidth() + two_pin_bbox.yWidth();
                    if (hpwl < min_hpwl) min_hpwl = hpwl;
                }
            }
            minProximity[gridNum.y - y - 1][x] = static_cast<double>(max_hpwl - min_hpwl) / 1000;

        }
    }
    //convert_to_33(minProximity, minProximity33Min, minProximity33Max);
    write_csv("min_proximity", minProximity);
    //write_csv("min_proximity_min", minProximity33Min);
    //write_csv("min_proximity_max", minProximity33Max);
}

// Calculate RUDY, but apply only the gcells with pins
void FeatureExtractor::pin_rudy() {

    int metalWidthI = 28;

    pinRudy.clear();
    pinRudy.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) pinRudy[i].resize(gridNum.x, 0);

    for (auto& net : def.netsMSNp) {
        auto net_name = net.first;
        auto net_p = net.second;
        auto& pin_list = net_p->cellPinMSS;
        BBox<int> net_bbox{0, 0, 0, 0};
        std::vector<Point<int>> pin_grid_list;

        for (auto& pin_info : pin_list) {
            auto comp_name = pin_info.first;
            auto pin_name = pin_info.second;

            if (comp_name == "PIN") {
                int x = std::round(def.IOpinsMSpp[pin_name]->centerP2.xF * 1000);
                int y = std::round(def.IOpinsMSpp[pin_name]->centerP2.yF * 1000);

                BBox<int> IOpin_bbox{x, y, x, y};
                if (net_bbox.lb.x == 0) net_bbox = IOpin_bbox;
                else net_bbox = net_bbox + IOpin_bbox;
                continue;
            }

            bool isMacro=false;
            for (auto& macro : def.comp_macrosMSCp) {
                auto macro_name = macro.first;
                if (macro_name == comp_name) {
                    isMacro = true;
                    break;
                }
            }

            pin* pin_p;
            if (!isMacro) pin_p = def.compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name];
            else pin_p = def.comp_macrosMSCp[comp_name]->macroMp->pinsMSPp[pin_name];
            
            for (auto& seg : pin_p->segmentsVp) {
                BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
                BBox<int> seg_bboxI = convert_int3(seg_bbox) - coreI[0];
                if (net_bbox.lb.x == 0) net_bbox = seg_bboxI;
                else net_bbox = net_bbox + seg_bboxI;

                BBox<int> grid_list = convert_to_grid(seg_bboxI, gridXI, gridYI);
                for (int x = grid_list.lb.x; x <= grid_list.rt.x; x++) {
                    for (int y = grid_list.lb.y; y <= grid_list.rt.y; y++) {
                        bool isExist = false;
                        for (auto& grid_loc : pin_grid_list) {
                            if (grid_loc.x == x && grid_loc.y == y) {
                                isExist = true;
                                break;
                            }
                        }
                        if (!isExist) pin_grid_list.emplace_back(x, y);
                    }
                }
            }
        }

        int w = net_bbox.rt.x - net_bbox.lb.x;
        int h = net_bbox.rt.y - net_bbox.lb.y;
        int bbox_size = w + h;
        double rudy_value = static_cast<double>((w + h - metalWidthI) * metalWidthI) / w / h;
        
        //net_bbox = net_bbox - coreI[0];

        for (auto& grid_loc : pin_grid_list) {
            if (grid_loc.x < 0 || grid_loc.x >= gridNum.x) {
                std::cout << "Out of range (x)" << std::endl;
                continue;
            }
            if (grid_loc.y < 0 || grid_loc.y >= gridNum.y) {
                std::cout << "Out of range (y)" << std::endl;
                continue;
            }
            BBox<int> grid_bbox(grid_loc.x * gridXI, grid_loc.y * gridYI, (grid_loc.x + 1) * gridXI, (grid_loc.y + 1) * gridYI);
            BBox<int> bbox_intersect = grid_bbox & net_bbox;
            double intersect_ratio = static_cast<double>(bbox_intersect.getArea()) / grid_bbox.getArea();
            // bbox threshold?
            if (bbox_size > rudy_threshold * gcellSizeI) {
                pinRudy[gridNum.y - grid_loc.y - 1][grid_loc.x] += (rudy_value * intersect_ratio);
            }
        }

    }
    write_csv("pin_rudy", pinRudy);

}


vector<vector<vector<pin*>>> FeatureExtractor::assign_pin_to_grid(){
    vector<vector<vector<pin*>>> pinMap;

	unordered_map<string,bool> is_done_map;
	int scale=1000;

    pinMap.clear();
    pinMap.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        pinMap[i].resize(gridNum.x);
    }

	for(size_t i_y=0;i_y<cellMap.size();i_y++){
		for(size_t i_x=0;i_x<cellMap[i_y].size();i_x++){
			BBox<double> grid_bbox_scaled(i_x*gridX*scale,i_y*gridY*scale,(i_x+1)*gridX*scale,(i_y+1)*gridY*scale);

            for(size_t i_c=0;i_c<cellMap[i_y][i_x].size();i_c++){
                string comp_name=cellMap[i_y][i_x][i_c];
                comp* compp=def.compsMSCp[comp_name];
                unordered_map<string,pin*> pinsMSPp=compp->cellCp->pinsMSPp;
                for(auto &it_p:pinsMSPp){
                    string pin_name=it_p.first;
                    pin* pinp=it_p.second;
                    BBox<double> pin_bbox=pin_to_bbox(pinp);
                    
                    pin_bbox.lb.x-=core[0][0];
                    pin_bbox.lb.y-=core[0][1];
                    pin_bbox.rt.x-=core[0][0];
                    pin_bbox.rt.y-=core[0][1];
                    
                    pin_bbox*=1000;
                    
                    //If partially covered,
                    if (grid_bbox_scaled.is_intersect(pin_bbox)){
                        pinMap[i_y][i_x].push_back(pinp);
                    }   
                }
            }
		}
	}

    for (const auto& comp_macroM : def.comp_macrosMSCp) {
        comp_macro* comp_macrop = comp_macroM.second;
        macro* macrop = comp_macrop->macroMp;
		
		unordered_map<string,pin*> pinsMSPp;
		for(auto &it:pinsMSPp){
			pin* pinp= it.second;
			vector<segment*> &segmentsVp=pinp->segmentsVp;
			
			double left=MAXF;
			double bottom=MAXF;
			double right=MINF;
			double top=MINF;

			for(size_t it2=0;it2<segmentsVp.size();it2++){
				vector<point2F*>& pointsVp=segmentsVp[it2]->pointsVp;
				for(size_t it3=0;it3<pointsVp.size();it3++){
					if(pointsVp[it3]->xF<left) left=pointsVp[it3]->xF;
					if(pointsVp[it3]->yF<bottom) bottom=pointsVp[it3]->yF;
					if(pointsVp[it3]->xF>right) right=pointsVp[it3]->xF;
					if(pointsVp[it3]->yF>top) top=pointsVp[it3]->yF;
				}
			}

			BBox<double> pin_bbox(left,bottom,right,top);
			BBox<int> pin_bboxI=convert_int3(pin_bbox);
			pin_bboxI=pin_bboxI-coreI[0];
			BBox<int> bbox_grid = convert_to_grid(pin_bboxI, gridXI, gridYI);
			
			for (int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++) {
				for (int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++) {
                    pinMap[y][x].push_back(pinp);
				}
			}

		}
        
    }
    return pinMap;
}

void FeatureExtractor::weighted_unfriendly() {

    // initialize pin density matrix
    unfriendly.clear();
    unfriendly.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        unfriendly[i].resize(gridNum.x, 0);
    }
        
    // calculate weighted-unfriendly
    for (const auto &compM : def.compsMSCp){
        auto comp_p = compM.second;
        auto cell_p = comp_p->cellCp;
        BBox<double> comp_bbox(comp_p->locP2.xF, comp_p->locP2.yF, comp_p->locP2.xF + cell_p->sizeP2.xF, comp_p->locP2.yF + cell_p->sizeP2.yF);
        BBox<int> comp_bboxI = convert_int3(comp_bbox);
        comp_bboxI = comp_bboxI - coreI[0];
        BBox<int> bbox_grid = convert_to_grid(comp_bboxI, gridXI, gridYI);
        for (int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
            for (int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++){
                //std::cout << "(" << gridNum.y - y - 1 << ", " << x << ")" << std::endl;
                BBox<int> grid_bbox(x * gridXI, y * gridYI, (x + 1) * gridXI, (y + 1) * gridYI);
                BBox<int> intersect = comp_bboxI & grid_bbox;
                if (intersect.lb.x == -1 || intersect.lb.y == -1 || intersect.rt.x == -1 || intersect.rt.y == -1)
                    continue;
                assert(stored_type_weight.count(cell_p->nameS) > 0);
                unfriendly[gridNum.y - y - 1][x] += stored_type_weight[cell_p->nameS] * intersect.getArea() / gridAreaI;
            }
        }
    }
    //convert_to_33(unfriendly, unfriendly33Min, unfriendly33Max);
    write_csv("unfriendly", unfriendly);
    //write_csv("unfriendly_min", unfriendly33Min);
    //write_csv("unfriendly_max", unfriendly33Max);
}


void FeatureExtractor::local_global_self_crossing_net() {
    vector<vector<vector<pin*>>> pinMap=assign_pin_to_grid();

    localNet.clear();
    localNet.resize(gridNum.y);
    
    globalNet.clear();
    globalNet.resize(gridNum.y);
    
    selfCrossingNet.clear();
    selfCrossingNet.resize(gridNum.y);
    
    for (int i = 0; i < gridNum.y; i++) {
        localNet[i].resize(gridNum.x);
        globalNet[i].resize(gridNum.x);
        selfCrossingNet[i].resize(gridNum.x);
    }
    

    for(size_t i_y=0;i_y<pinMap.size();i_y++){
		for(size_t i_x=0;i_x<pinMap[i_y].size();i_x++){
            //Local and global nets
            vector<string> included_net_namesV;
            vector<string> local_net_namesV;
            vector<pin*> local_net_pinV;

            included_net_namesV.clear();
            local_net_namesV.clear();
            local_net_pinV.clear();

            for(size_t i_p=0;i_p<pinMap[i_y][i_x].size();i_p++){
                pin* pinp=pinMap[i_y][i_x][i_p];
				if(pinp->netp==NULL) continue;
                string netNameS=pinp->netp->nameS;

                if (std::count(included_net_namesV.begin(),included_net_namesV.end(),netNameS)==1){
                    local_net_namesV.push_back(netNameS);
                }
                included_net_namesV.push_back(netNameS);
            }
            for(size_t i_p=0;i_p<pinMap[i_y][i_x].size();i_p++){
                pin* pinp=pinMap[i_y][i_x][i_p];
				if(pinp->netp==NULL) continue;
                string netNameS=pinp->netp->nameS;
				if(std::count(included_net_namesV.begin(),included_net_namesV.end(),netNameS)>1){
					local_net_pinV.push_back(pinp);
				}
			}

            int num_total_pins=pinMap[i_y][i_x].size();
            int num_local_nets=local_net_namesV.size();
            int num_local_pins=local_net_pinV.size();
            int num_global_nets=num_total_pins-num_local_pins;
			
            localNet[gridNum.y-i_y-1][i_x]=num_local_nets;
            globalNet[gridNum.y-i_y-1][i_x]=num_global_nets;

            //Self crossing net
            if(num_local_nets>=2){

                //Initialize self_crossing_net_bound
                unordered_map<int,unordered_map<string,vector<double>>> self_crossing_net_bound; //row: net_name: min,max
                for(pin* &pinp: local_net_pinV){
                    vector<double> bound={MAXF,MINF}; // min value, max value
                    string net_name=pinp->netp->nameS;

					double row=-1;
					if(pinp->cellp!=NULL){
						row=pinp->cellp->compp->locP2.yF;
					}
					else if(pinp->macrop!=NULL){
						//row=pinp->macrop->compp->locP2.yF;
						row=row_height* (int)((pinp->segmentsVp[0]->pointsVp[0]->yF-core.lb.y)/row_height)+core.lb.y;
					}
					row*=def.scaleI;

                    if(self_crossing_net_bound.find(row)==self_crossing_net_bound.end()){
						//Key not in map
                        unordered_map<string,vector<double>> net_bound;
                        net_bound.insert(make_pair(net_name,bound));
                        self_crossing_net_bound.insert(make_pair(row,net_bound));    
                    }
                    else {
                        self_crossing_net_bound[row].insert(make_pair(net_name,bound));
                    }
                }


                for(pin* &pinp: local_net_pinV){
                    string net_name=pinp->netp->nameS;
					
					double row=-1;
					if(pinp->cellp!=NULL){
						row=pinp->cellp->compp->locP2.yF;
					}
					else if(pinp->macrop!=NULL){
						//row=pinp->macrop->compp->locP2.yF;
						//row=pinp->segmentsVp[0]->pointsVp[0]->yF;
						row=row_height* (int)((pinp->segmentsVp[0]->pointsVp[0]->yF-core.lb.y)/row_height)+core.lb.y;
					}
					row*=def.scaleI;

                    //Lower bound
                    vector<double> lb={MAXF,MAXF};
                    for(segment* &segmentp: pinp->segmentsVp){
                        for(point2F* &coordsp: segmentp->pointsVp){
                            if(coordsp->xF<lb[0] || (coordsp->xF==lb[0] && coordsp->yF<lb[1]) ){
                                lb[0]=coordsp->xF;
                                lb[1]=coordsp->yF;
                            }
                        }
                    }

                    if(lb[0]<self_crossing_net_bound[row][net_name][0]){
                        self_crossing_net_bound[row][net_name][0]=lb[0];
                    }
                    if(lb[0]>self_crossing_net_bound[row][net_name][1]){
                        self_crossing_net_bound[row][net_name][1]=lb[0];
                    }
                }


                vector<string> self_crossing_netV;
                self_crossing_netV.clear();
                for(auto &it:self_crossing_net_bound){
                    double row=it.first;
                    int i_count=-1;
                    int j_count=-1;
                    for(auto &it1:self_crossing_net_bound[row]){
                        i_count++;
                        string net_name1=it1.first;
                        double lb_1=it1.second[0];
                        double ub_1=it1.second[1];
                        for(auto &it2:self_crossing_net_bound[row]){
                            j_count++;
                            string net_name2=it2.first;
                            if(i_count>=j_count || net_name1==net_name2) continue;
                            double lb_2=it2.second[0];
                            double ub_2=it2.second[1];

                            //if line intersect
                            if(!((min(lb_1,ub_1)<min(lb_2,ub_2) && max(lb_1,ub_1)<min(lb_2,ub_2)) || \
                            min(lb_1,ub_1)>min(lb_2,ub_2) && min(lb_1,ub_1)>max(lb_2,ub_2))){
                                if(std::count(self_crossing_netV.begin(),self_crossing_netV.end(),net_name1)==0){
                                    self_crossing_netV.push_back(net_name1);
                                }
                                if(std::count(self_crossing_netV.begin(),self_crossing_netV.end(),net_name2)==0){
                                    self_crossing_netV.push_back(net_name2);
                                }
                            }
                        }
                    }
                    selfCrossingNet[gridNum.y-i_y-1][i_x]=self_crossing_netV.size();
                }
            }
        }
    }
    write_csv("local_net", localNet);
    write_csv("global_net", globalNet);
    write_csv("self_crossing_net", selfCrossingNet);
	/*
    convert_to_33(localNet, localNet33Min, localNet33Max);
    convert_to_33(globalNet, globalNet33Min, globalNet33Max);
    convert_to_33(selfCrossingNet, selfCrossingNet33Min, selfCrossingNet33Max);
	
    write_csv("local_net_min", localNet33Min);
    write_csv("local_net_max", localNet33Max);
    write_csv("global_net_min", globalNet33Min);
    write_csv("global_net_max", globalNet33Max);
    write_csv("self_crossing_net_min", selfCrossingNet33Min);
    write_csv("self_crossing_net_max", selfCrossingNet33Max);
	*/
}

void FeatureExtractor::macro_density() {
    // initialize pin density matrix
    macroDensity.clear();
    macroDensity.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        macroDensity[i].resize(gridNum.x, 0);
    }

    // Calculate features
    for (const auto& comp_macroM : def.comp_macrosMSCp) {
        comp_macro* comp_macrop = comp_macroM.second;
        macro* macrop = comp_macrop->macroMp;
		point2F& macro_sizeP2=macrop->sizeP2;
		point2F& macro_locP2=comp_macrop->locP2;
		BBox<double> macro_bbox(macro_locP2.xF,macro_locP2.yF,macro_locP2.xF+macro_sizeP2.xF,macro_locP2.yF+macro_sizeP2.yF);

		BBox<int> macro_bboxI=convert_int3(macro_bbox);
		macro_bboxI=macro_bboxI-coreI[0];

		BBox<int> bbox_grid = convert_to_grid(macro_bboxI, gridXI, gridYI);
		//std::cout << bbox << ", " << bbox_grid << std::endl;
		for (int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++) {
			for (int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++) {
				BBox<int> grid(x * gridXI, y * gridYI, (x + 1) * gridXI, (y + 1) * gridYI);
				bool is_overlapped=macro_bboxI.is_intersect(grid);
				if(is_overlapped){
					//macroDensity[gridNum.y - y - 1][x] += 1;
                    int intersect = (macro_bboxI & grid).getArea();
					macroDensity[gridNum.y - y - 1][x] += static_cast<double>(intersect)/gridAreaI;
					//assert(macroDensity[gridNum.y-y-1][x]==1);
				}
			}
		}               
    }
    //convert_to_33(macroDensity, macroDensity33Min, macroDensity33Max);
    write_csv("macro_density", macroDensity);
    //write_csv("macro_density_min", macroDensity33Min);
    //write_csv("macro_density_max", macroDensity33Max);
}

void FeatureExtractor::routing_capacity() {
    routing_capacity_hor.clear();
    routing_capacity_ver.clear();
    
   routing_capacity_hor.resize(gridNum.y);
   routing_capacity_ver.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        routing_capacity_hor[i].resize(gridNum.x, 0);
        routing_capacity_ver[i].resize(gridNum.x, 0);
    }
   
   // Calculate capacity wrt a size of a side of a gcell.
   if(bottom_most_routing_layer == "NONE" || bottom_most_routing_layer == ""){
      bottom_most_routing_layer="MINT1";
   }
   if(upper_most_routing_layer=="NONE" || upper_most_routing_layer==""){
      if(name.find("MINT2",0) != std::string::npos){
         upper_most_routing_layer="MINT2";
      }
      else if(name.find("MINT4",0) != std::string::npos){
         upper_most_routing_layer="MINT4";
      }
      else{
         cout<<"Invalid input"<<endl;
         exit(1);
      }
   }


   int bottommost_layer_idx=METAL_LAYER_IDX[bottom_most_routing_layer];
   int topmost_layer_idx=METAL_LAYER_IDX[upper_most_routing_layer];
   /*
   string bottommmost_layer;
   string topmost_layer;
   for(auto &it=METAL_LAYER_IDX.begin();it!=METAL_LAYER_IDX.end();it++){
      if (it->second==bottommost_layer_idx){
         bottommost_layer=it->first;
      }
      if (it->second==topmost_layer_idx){
         topmost_layer=it->first;
      }
   }   
   */
   int horizontal_tracks=0;
   int vertical_tracks=0;
   for(int i=bottommost_layer_idx;i<=topmost_layer_idx;i++){
      if(METAL_ROUTING_DIRECTION[i]=="HOR"){
         horizontal_tracks+=METAL_ROUTING_TRACK[i];
      }
      if(METAL_ROUTING_DIRECTION[i]=="VER"){
         vertical_tracks+=METAL_ROUTING_TRACK[i];
      }
      if(METAL_ROUTING_DIRECTION[i]=="2D"){
         horizontal_tracks+=METAL_ROUTING_TRACK[i];
         vertical_tracks+=METAL_ROUTING_TRACK[i];
      }
   }

   int horizontal_capacity=std::round(gcellSize*1000)/std::round(SITE_Y*1000)*horizontal_tracks;
   int vertical_capacity=std::round(gcellSize*1000)/std::round(SITE_Y*1000)*vertical_tracks;
   for (int x=0;x<gridNum.x;x++){
      for(int y=0;y<gridNum.y;y++){
         routing_capacity_hor[gridNum.y - y - 1][x]=horizontal_capacity;
         routing_capacity_ver[gridNum.y - y - 1][x]=vertical_capacity;
      }
   }
   //cout<<"hor "<<horizontal_tracks<<endl;
   //cout<<"ver "<<vertical_tracks<<endl;
    write_csv("routing_capacity_horizontal", routing_capacity_hor);
    write_csv("routing_capacity_vertical", routing_capacity_ver);
}
void FeatureExtractor::congestion_map() {
    congDemHor.clear();
    congCapHor.clear();
    congDemHor.resize(gridNum.y);
    congCapHor.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        congDemHor[i].resize(gridNum.x + 1, 0);
        congCapHor[i].resize(gridNum.x + 1, 0);
    }
    congDemVer.clear();
    congCapVer.clear();
    congDemVer.resize(gridNum.y + 1);
    congCapVer.resize(gridNum.y + 1);
    for (int i = 0; i < gridNum.y + 1; i++) {
        congDemVer[i].resize(gridNum.x, 0);
        congCapVer[i].resize(gridNum.x, 0);
    }
    
    ifstream file(congPath);
    string line;

    bool start_parse = false;
    int maxLayer;
    int layer_count = 0;
    while (getline(file, line)){
        vector<string> tokens = split_by_space(line);
        if (tokens.size() <= 0) continue;
        if (tokens[0] == "Layer") {
            layer_count++;
            continue;
        }
        if (layer_count == 1 && tokens[0] == "Max" && tokens[1] == "Routing") {
            maxLayer = std::stoi(tokens[3].substr(tokens[3].length() - 1));
        }
        if (layer_count < 2) continue;

        if (tokens.size() != 6) {
            if (tokens.size() == 1 && tokens[0] == "1") {
                //std::cout << "End parsing." << std::endl;
            }
            else {
                for (auto& token : tokens) {
                    std::cout << token << " " << token.length() << std::endl;
                }
                std::cout << "Unexpected sentences! (" << tokens.size() << ")" << std::endl;
            }
            break;
        }

        std::string layer = tokens[0];
        int curLayer = std::stoi(layer.substr(layer.length() - 1));
        if (curLayer > maxLayer) continue;
        
        double lb_x = std::stod(tokens[1].substr(2));
        double lb_y = std::stod(tokens[2].substr(0, tokens[2].size() - 1));

        Point<int> lb;
        lb.x = std::round(lb_x * SCALE);
        lb.y = std::round(lb_y * SCALE);

        lb = lb - coreI.lb;

        if (!(lb.x >= 0 && lb.x <= gridNum.x * gridXI)) continue;
        if (!(lb.y >= 0 && lb.y <= gridNum.y * gridYI)) continue;

        //lb = (lb - coreI.lb.x) / gridSizeI;
        int gx = lb.x / gridXI;
        int gy = lb.y / gridYI;

        auto pos = tokens[5].find_first_of("/");
        int demand = std::stoi(tokens[5].substr(0, pos));
        int capacity = std::stoi(tokens[5].substr(pos + 1));
        
        // vertical layer
        if (curLayer % 2 == 0) {
            if (!(gridNum.y - gy >= 0 && gy >= 0)) continue;
            if (!(gx >= 0 && gx < gridNum.x)) continue;
            congCapVer[gridNum.y - gy][gx] += capacity;
            congDemVer[gridNum.y - gy][gx] += demand;
        }
        // Horizontal layer
        else {
            if (!(gridNum.y - gy > 0 && gy >= 0)) continue;
            if (!(gx >= 0 && gx <= gridNum.x)) continue;
            congCapHor[gridNum.y - gy - 1][gx] += capacity;
            congDemHor[gridNum.y - gy - 1][gx] += demand;
        }
    }

    ovflHor.clear();
    ovflVer.clear();
    ovflHor.resize(gridNum.y);
    ovflVer.resize(gridNum.y);
    for (int i = 0; i < gridNum.y; i++) {
        ovflHor[i].resize(gridNum.x, 0);
        ovflVer[i].resize(gridNum.x, 0);
    }

    int nan_track = 12;

    int hor_track, ver_track;
    hor_track = (maxLayer + 1) / 2 * nan_track;
    ver_track = maxLayer / 2 * nan_track;

    for (int i = 0; i < gridNum.y; i++) {
        for (int j = 0; j < gridNum.x; j++) {
            ovflVer[i][j] = static_cast<double>((congDemVer[i][j] + congDemVer[i + 1][j]) - (congCapVer[i][j] + congCapVer[i + 1][j])) / (2 * hor_track) + 1;
            ovflHor[i][j] = static_cast<double>((congDemHor[i][j] + congDemHor[i][j + 1]) - (congCapHor[i][j] + congCapHor[i][j + 1])) / (2 * ver_track) + 1;
        }
    }

    write_csv("h_congestion", ovflHor);
    write_csv("v_congestion", ovflVer);

}


void FeatureExtractor::hv_net_density() {
    HnetDensity.clear();
    VnetDensity.clear();
    HnetDensity.resize(gridNum.y);
    VnetDensity.resize(gridNum.y);
    for(int i = 0; i < gridNum.y; i++){
        HnetDensity[i].resize(gridNum.x, 0);
        VnetDensity[i].resize(gridNum.x, 0);
    }
    for(auto it = def.netsMSNp.cbegin(); it != def.netsMSNp.cend(); ++it){
        auto net_name = it->first;
        auto net_p = it->second;
        auto comp_list = net_p->cellPinMSS;
        BBox<double> net_bbox(0, 0, 0, 0);
        BBox<double> IOpin_bbox(0, 0, 0, 0);
        double bbox_size;
        double h_value, v_value;

        for(auto it2 = comp_list.cbegin(); it2 != comp_list.cend(); ++it2){
            string comp_name = it2->first;
            string pin_name = it2->second;
            if(comp_name == "PIN"){
                double x = def.IOpinsMSpp[pin_name]->centerP2.xF;
                double y = def.IOpinsMSpp[pin_name]->centerP2.yF;

                /*if(x < core[0][0]) x = core[0][0];
                else if(x > core[1][0]) x = core[1][0];
                if(y < core[0][1]) y = core[0][1];
                else if(y > core[1][1]) y = core[1][1];*/
                IOpin_bbox = BBox<double>(x, y, x, y);
                if(net_bbox.lb.x == 0) net_bbox = IOpin_bbox;
                else net_bbox = net_bbox + IOpin_bbox;
                continue;
            }
			bool isMacro=false;
			for(auto it_macro=def.comp_macrosMSCp.begin();it_macro!=def.comp_macrosMSCp.end();it_macro++){
				string macro_name=it_macro->first;
				if(macro_name==comp_name){
					isMacro=true;
					break;
				}
			}
			pin* single_pin;
			if(!isMacro){
				single_pin = def.compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name];
			}
			else{
				single_pin = def.comp_macrosMSCp[comp_name]->macroMp->pinsMSPp[pin_name];
			}
            for(auto seg : single_pin->segmentsVp){
                BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
                if(net_bbox.lb.x == 0) net_bbox = seg_bbox;
                else net_bbox = net_bbox + seg_bbox;
            }
        }
        
        double hor_length = (net_bbox.rt.x - net_bbox.lb.x) / gridX;
        double ver_length = (net_bbox.rt.y - net_bbox.lb.y) / gridY;
        h_value = 1 / hor_length;
        v_value = 1 / ver_length;

        net_bbox = net_bbox - core[0];
        BBox<int> bbox_grid = convert_to_grid(net_bbox, gridX, gridY);


        for(int x = bbox_grid.lb.x; x <= bbox_grid.rt.x; x++){
            if(x < 0 || x >= gridNum.x) continue;
            for(int y = bbox_grid.lb.y; y <= bbox_grid.rt.y; y++){
                if(y < 0 || y >= gridNum.y) continue;
                BBox<double> bbox(x * gridX, y * gridY, (x+1) * gridX, (y+1) * gridY);
                BBox<double> bbox_intersect = bbox & net_bbox;
                
                double intersect_ratio = bbox_intersect.getArea() / bbox.getArea();
                HnetDensity[gridNum.y - y - 1][x] += h_value * intersect_ratio;
                VnetDensity[gridNum.y - y - 1][x] += v_value * intersect_ratio;
            }
        }
    }
    
    write_csv("h_net_density", HnetDensity);
    write_csv("v_net_density", VnetDensity);
}

void FeatureExtractor::write_csv(std::string feature_name, Grid2D& feature) {
    int h = feature.size();
    int w = feature[0].size();

    fs::path csv_path = outPath / static_cast<fs::path>(name);
    if (!fs::exists(csv_path)) {
        fs::create_directories(csv_path);
    }

    fs::path feature_path = csv_path / fs::path(feature_name).replace_extension("csv");
    if (fs::exists(feature_path)) {
        fs::remove(feature_path);
    }

    std::ofstream writeFile(feature_path.string());

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (x != w - 1) writeFile << feature[y][x] << ",";
            else writeFile << feature[y][x] << "\n";
        }
    }
    writeFile.close();
}

Point<int> FeatureExtractor::convert_int3 (Point<double>& p) {
    return Point<int>(std::round(p[0] * 1000), std::round(p[1] * 1000));
}

BBox<int> FeatureExtractor::convert_int3 (BBox<double>& bbox) {
    return BBox<int>(std::round(bbox[0][0] * 1000), std::round(bbox[0][1] * 1000), std::round(bbox[1][0] * 1000), std::round(bbox[1][1] * 1000));
}

BBox<int> FeatureExtractor::convert_to_grid(BBox<double>& bbox, double grid_size_x, double grid_size_y) {
    auto bbox_int = convert_int3(bbox);

    int grid_size_x_scale = std::round(grid_size_x * 1000);
    int grid_size_y_scale = std::round(grid_size_y * 1000);

    int lb_gx = std::floor(bbox_int.lb.x / grid_size_x_scale);
    int lb_gy = std::floor(bbox_int.lb.y / grid_size_y_scale);
    int rt_gx = std::floor(bbox_int.rt.x / grid_size_x_scale);
    int rt_gy = std::floor(bbox_int.rt.y / grid_size_y_scale);
    
    if (bbox_int.rt.x % grid_size_x_scale == 0 && rt_gx > 0) rt_gx--;
    if (bbox_int.rt.y % grid_size_y_scale == 0 && rt_gy > 0) rt_gy--;

    return BBox<int>(lb_gx, lb_gy, rt_gx, rt_gy);
}

pair<FeatureExtractor::BBox2D,FeatureExtractor::BBox2D> FeatureExtractor::seperate_bbox(pin* segments){
    pair<BBox2D, BBox2D> sep_bbox;
    BBox2D sep_bbox_M1;
    BBox2D sep_bbox_M2;
    bool find_rects;
    
    for(auto seg : segments->segmentsVp){
        BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
        //seg_bbox = seg_bbox - core[0];
        if(seg->layerS == "MINT1"){
            vector<BBox<double>> newVec;
            newVec.push_back(seg_bbox);
            sep_bbox_M2.push_back(newVec);
            continue;
        }
        if(seg == (segments->segmentsVp)[0]){
            vector<BBox<double>> newVec;
            newVec.push_back(seg_bbox);
            sep_bbox_M1.push_back(newVec);
            
            continue;
        }
        find_rects = false;
        
        for(auto& bboxVec : sep_bbox_M1){
            for(auto& single_bbox : bboxVec){
                if((single_bbox & seg_bbox).lb.x != -1){
                    bboxVec.push_back(seg_bbox);
                    find_rects = true;
                    break;
                }
            }
            if(find_rects) break;
        }
        if(!find_rects){
            vector<BBox<double>> newVec;
            newVec.push_back(seg_bbox);
            sep_bbox_M1.push_back(newVec);
            continue;
        }
    }
    //auto it = sep_bbox_M1.insert(sep_bbox.end(), sep_bbox_M2.begin(), sep_bbox_M2.end());
    sep_bbox.first = sep_bbox_M1;
    sep_bbox.second = sep_bbox_M2;
    return sep_bbox;
}

int FeatureExtractor::calculate_overlap(int x, int y){
    vector<BBox<double>> M1, M2;
    BBox<double> grid_bbox(x*gridX, y*gridY, (x+1)*gridX, (y+1)*gridY);

    int num_overlap = 0;
    auto cellname_list = cellMap[y][x];
    for(auto cellname : cellname_list){

        auto& pin_list = def.compsMSCp[cellname]->cellCp->pinsMSPp;
        for(auto& single_pin : pin_list){
            std::string pin_name = single_pin.first;\
            if(pin_name == "VDD" || pin_name == "VSS") continue;

            for(auto& seg : single_pin.second->segmentsVp){
                BBox<double> seg_bbox((seg->pointsVp)[0]->xF, (seg->pointsVp)[0]->yF, (seg->pointsVp)[2]->xF, (seg->pointsVp)[2]->yF);
                if((seg_bbox - core[0]).overlap(grid_bbox)){
                    if(seg->layerS == "M1") M1.push_back(seg_bbox);
                    else if(seg->layerS == "MINT1") M2.push_back(seg_bbox);
                }
            }
        }
    }
    for(auto& M1_bbox : M1){
        for(auto& M2_bbox : M2){
            if(M1_bbox.overlap(M2_bbox)){
                num_overlap++; 
            }
        }
    }
    return num_overlap;
}

BBox<int> FeatureExtractor::convert_to_grid(BBox<int>& bbox, int grid_size_x, int grid_size_y) {
    
    int lb_gx = std::floor(bbox.lb.x / grid_size_x);
    int lb_gy = std::floor(bbox.lb.y / grid_size_y);
    int rt_gx = std::floor(bbox.rt.x / grid_size_x);
    int rt_gy = std::floor(bbox.rt.y / grid_size_y);
    
    if (bbox.rt.x % grid_size_x == 0 && rt_gx > 0) rt_gx--;
    if (bbox.rt.y % grid_size_y == 0 && rt_gy > 0) rt_gy--;

    return BBox<int>(lb_gx, lb_gy, rt_gx, rt_gy);
}

void FeatureExtractor::convert_to_33(Grid2D &feature, Grid2D &feature33Min, Grid2D &feature33Max){
    feature33Min.clear();
    feature33Max.clear();
    feature33Min.resize(gridNum.y);
    feature33Max.resize(gridNum.y);
    for(int i = 0; i < gridNum.y; i++){
        feature33Min[i].resize(gridNum.x);
        feature33Max[i].resize(gridNum.x);
    }
    for(int i = 0; i < gridNum.x; i++){
        for(int j = 0; j < gridNum.y; j++){
            double maxval = -99999999;
            double minval = 99999999;
            for(int m = -1; m <= 1; m++){
                for(int n = -1; n<= 1; n++){
                    if( i+m < 0 || i+m >= gridNum.x || j+n < 0 || j+n >= gridNum.y) continue;
                    if(feature[j+n][i+m] > maxval) maxval = feature[j+n][i+m];
                    if(feature[j+n][i+m] < minval) minval = feature[j+n][i+m];
                }
            }
            feature33Min[j][i] = minval;
            feature33Max[j][i] = maxval;
        }
    }
}


vector<string> FeatureExtractor::split_by_space(string line) {
    vector<string> tokens;
    size_t prev = line.find_first_not_of(" ");
    size_t pos;

    while ((pos = line.find_first_of(" ", prev)) != string::npos) {
        if (pos > prev) {
            tokens.push_back(line.substr(prev, pos - prev));
        }
        prev = pos + 1;
    }
    if (prev < line.length()) {
        tokens.push_back(line.substr(prev, string::npos));
    }

    return tokens;
}
/*
void FeatureExtractor::generate_graph() {

    fs::path feature_path = outPath / fs::path(name);

    fs::path node_feature_path = feature_path / fs::path("graph_x.txt");
    fs::path edge_index_path = feature_path / fs::path("graph_a.txt");
    fs::path edge_weight_path = feature_path / fs::path("graph_aw.txt");
    fs::path graph_index_path = feature_path / fs::path("graph_i.txt");

    // open output files
    std::ofstream nf_out, ei_out, ew_out, gi_out;    
    //std::ofstream dbg_out;

    nf_out.open(node_feature_path.string());
    ei_out.open(edge_index_path.string());
    ew_out.open(edge_weight_path.string());
    gi_out.open(graph_index_path.string());

    int node_index = 0;
    int graph_index = 0;
    

    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {

            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                for (int i = 0; i < 11; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue; 
            }

            std::vector<pin*> pin_list;
            std::unordered_map<std::string, std::vector<pin*>> net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;  // avoid pins that cannot estimate FLUTE direction
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;

                            if (net_list.count(net_name) == 1) {
                                net_list[net_name].push_back(pinM.second);
                            }
                            else {
                                std::vector<pin*> net_pin_list;
                                net_pin_list.push_back(pinM.second);
                                net_list.insert({net_name, net_pin_list});
                            }
                            break;
                        }
                    }
                }
            }
            

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 11; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue;                
            }

            // node feature representation
            std::vector<std::tuple<pin*, double, double, std::vector<int>>> pin_info;
            for (auto pin_p : pin_list) {
                if (pin_p->targetP2 == nullptr) continue;
                std::unordered_map<int, int> ap_x_dic;
                std::vector<int> digit_pin(9);
                int num_ap = 0;
                double sum_x_ap = 0;

                for (auto pin_ap_p : pin_p->apsVp) {
                    int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                    if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                        if (ap_x_dic.count(ap_x) > 0) ap_x_dic[ap_x]++;
                        else ap_x_dic[ap_x] = 1;           
                        
                        num_ap++;
                        sum_x_ap += static_cast<double>(ap_x - grid_bbox.lb.x) / gridXI;
                        int track = (std::round(pin_ap_p->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                        digit_pin[track - 2] = 1;        
                    }
                }
                // Calculate average relative x-coord of APs in grid 
                pin_p->avg_ap_x_ingrid = (num_ap == 0) ? 0 : sum_x_ap / num_ap;

                // Calculate dominated x/y-coord of APs in grid
                int pin_dx;
                int max_ap = -1;
                for (auto& ap_x_dic_iter : ap_x_dic) {
                    if (ap_x_dic_iter.second > max_ap) {
                        max_ap = ap_x_dic_iter.second;
                        pin_dx = ap_x_dic_iter.first;
                    }
                }

                int pin_dy = 0;
                for (auto pin_ap_p : pin_p->apsVp) {
                    int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                    if (ap_x == pin_dx) pin_dy += std::round(pin_ap_p->pointP2.yF * 1000);                            
                }
                pin_dy /= max_ap;

                pin_p->center_ingridP2.xF = static_cast<double>(pin_dx) / 1000;
                pin_p->center_ingridP2.yF = static_cast<double>(pin_dy) / 1000;

                double pin_dx_d = static_cast<double>(pin_dx - grid_bbox.lb.x) / gridXI;

                nf_out << pin_p->avg_ap_x_ingrid << " " << pin_dx_d << " ";

                for (int i = 0; i < 9; i++) {
                    nf_out << digit_pin[i] << " ";
                }
                nf_out << std::endl;

                pin_info.push_back(std::make_tuple(pin_p, pin_p->avg_ap_x_ingrid, pin_dx_d, digit_pin));

            }

            // edge feature representation
            int num_eff_pin = pin_info.size();

            for (int i = 0; i < num_eff_pin; i++) {
                auto pin_i_p = std::get<0>(pin_info[i]);
                auto& pin_i_digit = std::get<3>(pin_info[i]);
                
                for (int j = 0; j < num_eff_pin; j++) {
                    if (i == j) continue;

                    auto pin_j_p = std::get<0>(pin_info[j]);
                    auto& pin_j_digit = std::get<3>(pin_info[j]);
                    
                    // doesn't draw edge between local net pins
                    if (pin_i_p->netp->nameS == pin_j_p->netp->nameS) continue;

                    // calculate primary/secondary/last direction 
                    bool j_is_left = (std::get<2>(pin_info[i]) >= std::get<2>(pin_info[j])) ? true : false;
                    DIR dir;

                    int pin_cx = std::round(pin_i_p->center_ingridP2.xF * 1000);
                    int pin_cy = std::round(pin_i_p->center_ingridP2.yF * 1000);

                    int pin_tx = std::round(pin_i_p->targetP2->xF * 1000);
                    int pin_ty = std::round(pin_i_p->targetP2->yF * 1000);

                    int dx = pin_tx - pin_cx;
                    int dy = pin_ty - pin_cy;

                    if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                        dir = DIR::SEC;
                    }
                    else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                        dir = DIR::SEC;
                    }
                    else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                        if (j_is_left) dir = DIR::PRIM;
                        else dir = DIR::LAST;
                    }
                    else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                        if (j_is_left) dir = DIR::LAST;
                        else dir = DIR::PRIM;
                    }

                    // calculate distance between two pins
                    double dist = std::abs(std::get<1>(pin_info[i]) - std::get<1>(pin_info[j]));

                    // calculate pin ap overlap ratio
                    int num_grid_track = pin_i_digit.size();
                    int num_track_ap_i = 0;
                    int num_track_ap_both = 0;
                    
                    for (int k = 0; k < num_grid_track; k++) {
                        if (pin_i_digit[k] == 1) num_track_ap_i++;
                        if (pin_i_digit[k] == 1 && pin_j_digit[k] == 1) num_track_ap_both++;
                    }
                    double ap_overlap_ratio = static_cast<double>(num_track_ap_both) / num_track_ap_i;

                    // h-congestion
                    //double h_cong = ovflHor[y][x];
                    double h_cong = HnetDensity[y][x];

                    // print edge information (index & weight)
                    ei_out << node_index + j << " " << node_index + i << std::endl;
                    
                    if (dir == DIR::PRIM) ew_out << 1 << " " << 0 << " " << 0 << " ";
                    else if (dir == DIR::SEC) ew_out << 0 << " " << 1 << " " << 0 << " ";
                    else ew_out << 0 << " " << 0 << " " << 1 << " ";
                    ew_out << dist << " " << ap_overlap_ratio << " " << h_cong << std::endl;

                }
            }

            
            // graph index representation
            for (int i = 0; i < num_eff_pin; i++) {
                gi_out << graph_index << std::endl;
            }

            node_index += num_eff_pin;
            graph_index++;
        }
    }

    // close output files
    
    nf_out.close();
    ei_out.close();
    ew_out.close();
    gi_out.close();
    //dbg_out.close();

}
*/
void FeatureExtractor::generate_graph() {

    fs::path feature_path = outPath / fs::path(name);

    fs::path node_feature_path = feature_path / fs::path("graph_x.txt");
    fs::path edge_index_path = feature_path / fs::path("graph_a.txt");
    fs::path edge_weight_path = feature_path / fs::path("graph_aw.txt");
    fs::path graph_index_path = feature_path / fs::path("graph_i.txt");

    // open output files
    std::ofstream nf_out, ei_out, ew_out, gi_out;    
    //std::ofstream dbg_out;

    nf_out.open(node_feature_path.string());
    ei_out.open(edge_index_path.string());
    ew_out.open(edge_weight_path.string());
    gi_out.open(graph_index_path.string());

    int node_index = 0;
    int graph_index = 0;
    
    bool gi_printed = false;

    std::cout << "Size: " << gridNum.y << ", " << gridNum.x << std::endl;

    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {

            //if (graph_index != y * gridNum.y + x) {
            //   std::cout << "Wrong index: " << y << ", " << x << std::endl;
            //}

            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                for (int i = 0; i < 11; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;
                gi_printed = true;
                continue; 
            }

            std::vector<pin*> pin_list;
            std::unordered_map<std::string, std::vector<pin*>> net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;  // avoid pins that cannot estimate FLUTE direction
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;

                            if (net_list.count(net_name) == 1) {
                                net_list[net_name].push_back(pinM.second);
                            }
                            else {
                                std::vector<pin*> net_pin_list;
                                net_pin_list.push_back(pinM.second);
                                net_list.insert({net_name, net_pin_list});
                            }
                            break;
                        }
                    }
                }
            }
            

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 11; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;
                gi_printed = true;
                continue;                
            }

            // node feature representation
            std::vector<std::tuple<pin*, double, double, std::vector<int>>> pin_info;
            for (auto pin_p : pin_list) {
                if (pin_p->targetP2 == nullptr) continue;
                std::unordered_map<int, int> ap_x_dic;
                std::vector<int> digit_pin(9);
                int num_ap = 0;
                double sum_x_ap = 0;

                for (auto pin_ap_p : pin_p->apsVp) {
                    int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                    if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                        if (ap_x_dic.count(ap_x) > 0) ap_x_dic[ap_x]++;
                        else ap_x_dic[ap_x] = 1;           
                        
                        num_ap++;
                        sum_x_ap += static_cast<double>(ap_x - grid_bbox.lb.x) / gridXI;
                        int track = (std::round(pin_ap_p->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                        digit_pin[track - 2] = 1;        
                    }
                }
                // Calculate average relative x-coord of APs in grid 
                pin_p->avg_ap_x_ingrid = (num_ap == 0) ? 0 : sum_x_ap / num_ap;

                // Calculate dominated x/y-coord of APs in grid
                int pin_dx;
                int max_ap = -1;
                for (auto& ap_x_dic_iter : ap_x_dic) {
                    if (ap_x_dic_iter.second > max_ap) {
                        max_ap = ap_x_dic_iter.second;
                        pin_dx = ap_x_dic_iter.first;
                    }
                }

                int pin_dy = 0;
                for (auto pin_ap_p : pin_p->apsVp) {
                    int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                    if (ap_x == pin_dx) pin_dy += std::round(pin_ap_p->pointP2.yF * 1000);                            
                }
                pin_dy /= max_ap;

                pin_p->center_ingridP2.xF = static_cast<double>(pin_dx) / 1000;
                pin_p->center_ingridP2.yF = static_cast<double>(pin_dy) / 1000;

                double pin_dx_d = static_cast<double>(pin_dx - grid_bbox.lb.x) / gridXI;

                nf_out << pin_p->avg_ap_x_ingrid << " " << pin_dx_d << " ";


                // y_avg, pin length
                /*
                int y_min = 200;
                int y_max = -100;
                int pin_l = 0;

                for (int i = 0; i < 9; i++) {
                    if (digit_pin[i]) {
                        if (i < y_min) y_min = i;
                        if (i > y_max) y_max = i;
                        pin_l++;
                    }
                }
                if (y_min == 200) y_min = 0;
                if (y_max == -100) y_max = 0;
                nf_out << static_cast<double>(y_max + y_min + 4) / 24 << " " << static_cast<double>(pin_l) / 12 << std::endl;
                */
                for (int i = 0; i < 9; i++) nf_out << digit_pin[i] << " ";
                nf_out << endl;

                pin_info.push_back(std::make_tuple(pin_p, pin_p->avg_ap_x_ingrid, pin_dx_d, digit_pin));

            }

            // edge feature representation
            int num_eff_pin = pin_info.size();
            if (num_eff_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 11; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;
                gi_printed = true;
                continue;                
            }

            for (int i = 0; i < num_eff_pin; i++) {
                auto pin_i_p = std::get<0>(pin_info[i]);
                auto& pin_i_digit = std::get<3>(pin_info[i]);
                
                for (int j = 0; j < num_eff_pin; j++) {
                    if (i == j) continue;

                    auto pin_j_p = std::get<0>(pin_info[j]);
                    auto& pin_j_digit = std::get<3>(pin_info[j]);
                    
                    // doesn't draw edge between local net pins
                    if (pin_i_p->netp->nameS == pin_j_p->netp->nameS) continue;

                    // calculate primary/secondary/last direction 
                    bool j_is_left = (std::get<2>(pin_info[i]) >= std::get<2>(pin_info[j])) ? true : false;
                    DIR dir;

                    int pin_cx = std::round(pin_i_p->center_ingridP2.xF * 1000);
                    int pin_cy = std::round(pin_i_p->center_ingridP2.yF * 1000);

                    int pin_tx = std::round(pin_i_p->targetP2->xF * 1000);
                    int pin_ty = std::round(pin_i_p->targetP2->yF * 1000);

                    int dx = pin_tx - pin_cx;
                    int dy = pin_ty - pin_cy;

                    if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                        dir = DIR::SEC;
                    }
                    else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                        dir = DIR::SEC;
                    }
                    else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                        if (j_is_left) dir = DIR::PRIM;
                        else dir = DIR::LAST;
                    }
                    else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                        if (j_is_left) dir = DIR::LAST;
                        else dir = DIR::PRIM;
                    }

                    // calculate distance between two pins
                    double dist = std::abs(std::get<1>(pin_info[i]) - std::get<1>(pin_info[j]));

                    // calculate pin ap overlap ratio
                    int num_grid_track = pin_i_digit.size();
                    int num_track_ap_i = 0;
                    int num_track_ap_both = 0;
                    
                    for (int k = 0; k < num_grid_track; k++) {
                        if (pin_i_digit[k] == 1) num_track_ap_i++;
                        if (pin_i_digit[k] == 1 && pin_j_digit[k] == 1) num_track_ap_both++;
                    }
                    double ap_overlap_ratio = static_cast<double>(num_track_ap_both) / num_track_ap_i;

                    // h-congestion
                    //double h_cong = ovflHor[y][x];
                    double h_cong = HnetDensity[y][x];

                    // print edge information (index & weight)
                    ei_out << node_index + j << " " << node_index + i << std::endl;
                    
                    if (dir == DIR::PRIM) ew_out << 1 << " " << 0 << " " << 0 << " ";
                    else if (dir == DIR::SEC) ew_out << 0 << " " << 1 << " " << 0 << " ";
                    else ew_out << 0 << " " << 0 << " " << 1 << " ";
                    ew_out << dist << " " << ap_overlap_ratio << " " << h_cong << std::endl;

                }
            }

            
            // graph index representation
            for (int i = 0; i < num_eff_pin; i++) {
                gi_out << graph_index << std::endl;
                gi_printed = true;
            }

            node_index += num_eff_pin;
            graph_index++;

            if (!gi_printed) std::cout << "missed: " << y << ", " << x << std::endl;
        }
    }

    // close output files
    
    nf_out.close();
    ei_out.close();
    ew_out.close();
    gi_out.close();
    //dbg_out.close();

}


// for ASAP7
/*
void FeatureExtractor::generate_graph() {

    fs::path feature_path = outPath / fs::path(name);

    fs::path node_feature_path = feature_path / fs::path("node_feature.txt");
    fs::path edge_index_path = feature_path / fs::path("edge_index.txt");
    fs::path edge_weight_path = feature_path / fs::path("edge_weight_hcong.txt");
    fs::path graph_index_path = feature_path / fs::path("graph_index.txt");

    // open output files
    std::ofstream nf_out, ei_out, ew_out, gi_out;    
    //std::ofstream dbg_out;

    nf_out.open(node_feature_path.string());
    ei_out.open(edge_index_path.string());
    ew_out.open(edge_weight_path.string());
    gi_out.open(graph_index_path.string());

    int node_index = 0;
    int graph_index = 0;
    

    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {

            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                for (int i = 0; i < 8; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue; 
            }

            std::vector<pin*> pin_list;
            std::unordered_map<std::string, std::vector<pin*>> net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;  // avoid pins that cannot estimate FLUTE direction
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;

                            if (net_list.count(net_name) == 1) {
                                net_list[net_name].push_back(pinM.second);
                            }
                            else {
                                std::vector<pin*> net_pin_list;
                                net_pin_list.push_back(pinM.second);
                                net_list.insert({net_name, net_pin_list});
                            }
                            break;
                        }
                    }
                }
            }
            

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 8; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue;                
            }

            // node feature representation
            std::vector<std::tuple<pin*, double, double, std::vector<int>>> pin_info;
            for (auto pin_p : pin_list) {
                if (pin_p->targetP2 == nullptr) continue;
                std::unordered_map<int, int> ap_x_dic;
                std::vector<int> digit_pin(6);
                int num_ap = 0;
                double sum_x_ap = 0;

                for (auto pin_ap_p : pin_p->apsVp) {
                    int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                    if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                        if (ap_x_dic.count(ap_x) > 0) ap_x_dic[ap_x]++;
                        else ap_x_dic[ap_x] = 1;           
                        
                        num_ap++;
                        sum_x_ap += static_cast<double>(ap_x - grid_bbox.lb.x) / gridXI;
                        //int track = (std::round(pin_ap_p->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                        int track = (std::round(pin_ap_p->pointP2.yF * 1000) - grid_bbox.lb.y) / 27; // ASAP7 pitch = 27
                        digit_pin[track - 1] = 1;        
                    }
                }
                // Calculate average relative x-coord of APs in grid 
                pin_p->avg_ap_x_ingrid = (num_ap == 0) ? 0 : sum_x_ap / num_ap;

                // Calculate dominated x/y-coord of APs in grid
                int pin_dx;
                int max_ap = -1;
                for (auto& ap_x_dic_iter : ap_x_dic) {
                    if (ap_x_dic_iter.second > max_ap) {
                        max_ap = ap_x_dic_iter.second;
                        pin_dx = ap_x_dic_iter.first;
                    }
                }

                int pin_dy = 0;
                for (auto pin_ap_p : pin_p->apsVp) {
                    int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                    if (ap_x == pin_dx) pin_dy += std::round(pin_ap_p->pointP2.yF * 1000);                            
                }
                pin_dy /= max_ap;

                pin_p->center_ingridP2.xF = static_cast<double>(pin_dx) / 1000;
                pin_p->center_ingridP2.yF = static_cast<double>(pin_dy) / 1000;

                double pin_dx_d = static_cast<double>(pin_dx - grid_bbox.lb.x) / gridXI;

                nf_out << pin_p->avg_ap_x_ingrid << " " << pin_dx_d << " ";

                for (int i = 0; i < 6; i++) {
                    nf_out << digit_pin[i] << " ";
                }
                nf_out << std::endl;

                pin_info.push_back(std::make_tuple(pin_p, pin_p->avg_ap_x_ingrid, pin_dx_d, digit_pin));

            }

            // edge feature representation
            int num_eff_pin = pin_info.size();

            for (int i = 0; i < num_eff_pin; i++) {
                auto pin_i_p = std::get<0>(pin_info[i]);
                auto& pin_i_digit = std::get<3>(pin_info[i]);
                
                for (int j = 0; j < num_eff_pin; j++) {
                    if (i == j) continue;

                    auto pin_j_p = std::get<0>(pin_info[j]);
                    auto& pin_j_digit = std::get<3>(pin_info[j]);
                    
                    // doesn't draw edge between local net pins
                    if (pin_i_p->netp->nameS == pin_j_p->netp->nameS) continue;

                    // calculate primary/secondary/last direction 
                    bool j_is_left = (std::get<2>(pin_info[i]) >= std::get<2>(pin_info[j])) ? true : false;
                    DIR dir;

                    int pin_cx = std::round(pin_i_p->center_ingridP2.xF * 1000);
                    int pin_cy = std::round(pin_i_p->center_ingridP2.yF * 1000);

                    int pin_tx = std::round(pin_i_p->targetP2->xF * 1000);
                    int pin_ty = std::round(pin_i_p->targetP2->yF * 1000);

                    int dx = pin_tx - pin_cx;
                    int dy = pin_ty - pin_cy;

                    if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                        dir = DIR::SEC;
                    }
                    else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                        dir = DIR::SEC;
                    }
                    else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                        if (j_is_left) dir = DIR::PRIM;
                        else dir = DIR::LAST;
                    }
                    else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                        if (j_is_left) dir = DIR::LAST;
                        else dir = DIR::PRIM;
                    }

                    // calculate distance between two pins
                    double dist = std::abs(std::get<1>(pin_info[i]) - std::get<1>(pin_info[j]));

                    // calculate pin ap overlap ratio
                    int num_grid_track = pin_i_digit.size();
                    int num_track_ap_i = 0;
                    int num_track_ap_both = 0;
                    
                    for (int k = 0; k < num_grid_track; k++) {
                        if (pin_i_digit[k] == 1) num_track_ap_i++;
                        if (pin_i_digit[k] == 1 && pin_j_digit[k] == 1) num_track_ap_both++;
                    }
                    double ap_overlap_ratio = static_cast<double>(num_track_ap_both) / num_track_ap_i;

                    // h-congestion
                    double h_cong = ovflHor[y][x];
                    //double h_cong = HnetDensity[y][x];

                    // print edge information (index & weight)
                    ei_out << node_index + j << " " << node_index + i << std::endl;
                    
                    if (dir == DIR::PRIM) ew_out << 1 << " " << 0 << " " << 0 << " ";
                    else if (dir == DIR::SEC) ew_out << 0 << " " << 1 << " " << 0 << " ";
                    else ew_out << 0 << " " << 0 << " " << 1 << " ";
                    ew_out << dist << " " << ap_overlap_ratio << " " << h_cong << std::endl;

                }
            }

            
            // graph index representation
            for (int i = 0; i < num_eff_pin; i++) {
                gi_out << graph_index << std::endl;
            }

            node_index += num_eff_pin;
            graph_index++;
        }
    }

    // close output files
    
    nf_out.close();
    ei_out.close();
    ew_out.close();
    gi_out.close();
    //dbg_out.close();

}
*/


//////////////////////////////////////////////// OLD /////////////////////////////////////////////////////////


/*
void FeatureExtractor::generate_graph() {

    std::vector<std::tuple<double, double, double, double>> params;

    params.emplace_back(0.5, 0.2, 0.1, 384);
    params.emplace_back(0.5, 0.2, 0.1, 192);
    params.emplace_back(0.5, 0.25, 0, 384);
    params.emplace_back(0.5, 0.25, 0, 192);
    params.emplace_back(1, 0, 0, 384);
    params.emplace_back(1, 0, 0, 192);

    int num_params = params.size();

    fs::path feature_path = outPath / fs::path(name);

    fs::path node_feature_path = feature_path / fs::path("node_feature_pin_digit.txt");
    fs::path edge_index_path = feature_path / fs::path("edge_index.txt");
    fs::path edge_weight_path = feature_path / fs::path("edge_weight.txt");
    fs::path graph_index_path = feature_path / fs::path("graph_index.txt");

    // open output files
    std::ofstream gi_out;    
    //std::ofstream dbg_out;

    //nf_out.open(node_feature_path.string());
    //ei_out.open(edge_index_path.string());
    //ew_out.open(edge_weight_path.string());
    gi_out.open(graph_index_path.string());

    int node_index = 0;
    int graph_index = 0;
    

    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {

            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                fs::path node_feature_pd_path = feature_path / fs::path("node_feature_pd.txt");
                std::ofstream npd_out(node_feature_pd_path, std::ios::app);
                for (int i = 0; i < 9; i++) npd_out << 0 << " ";
                npd_out << std::endl;
                npd_out.close();

                for (int k = 1; k <= num_params; k++) {
                    std::string file_name = "node_feature_ew" + to_string(k) + ".txt";
                    fs::path node_feature_path = feature_path / fs::path(file_name);
                    std::ofstream nf_out(node_feature_path, std::ios::app);
                    for (int i = 0; i < 4; i++) nf_out << 0 << " ";
                    nf_out << std::endl;
                    nf_out.close();
                } 

                node_index++;
                gi_out << graph_index++ << std::endl;

                continue; 
            }

            std::vector<pin*> pin_list;
            std::unordered_map<std::string, std::vector<pin*>> net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;  // avoid pins that cannot estimate FLUTE direction
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;

                            if (net_list.count(net_name) == 1) {
                                net_list[net_name].push_back(pinM.second);
                            }
                            else {
                                std::vector<pin*> net_pin_list;
                                net_pin_list.push_back(pinM.second);
                                net_list.insert({net_name, net_pin_list});
                            }
                            break;
                        }
                    }
                }
            }
            

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                fs::path node_feature_pd_path = feature_path / fs::path("node_feature_pd.txt");
                std::ofstream npd_out(node_feature_pd_path, std::ios::app);
                for (int i = 0; i < 9; i++) npd_out << 0 << " ";
                npd_out << std::endl;
                npd_out.close();

                for (int k = 1; k <= num_params; k++) {
                    std::string file_name = "node_feature_ew" + to_string(k) + ".txt";
                    fs::path node_feature_path = feature_path / fs::path(file_name);
                    std::ofstream nf_out(node_feature_path, std::ios::app);
                    for (int i = 0; i < 4; i++) nf_out << 0 << " ";
                    nf_out << std::endl;
                    nf_out.close();
                } 
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue;                
            }

            // net node representation using FLUTE
            std::vector<std::tuple<std::string, int, std::vector<pin*>, std::vector<int>>> net_info;
            int net_index = 0;
            for (auto& net : net_list) {
                auto net_name = net.first;
                auto& net_pin_v = net.second;
                std::vector<int> digit_pin(9);

                for (auto pin_p : net_pin_v) {
                        
                    if (pin_p->targetP2 == nullptr) continue;

                    //dbg_out << "Pin : " << "(" << pin_p->cellp->compp->nameS << " " << pin_p->nameS << ")" << std::endl;

                    std::unordered_map<int, int> ap_x_dic;

                    for (auto pin_ap_p : pin_p->apsVp) {
                        int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                        if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                            if (ap_x_dic.count(ap_x) > 0) ap_x_dic[ap_x]++;
                            else ap_x_dic[ap_x] = 1;           

                            int track = (std::round(pin_ap_p->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                            digit_pin[track - 2] = 1;        
                 
                        }
                    }

                    int pin_dx;
                    int max_ap = -1;
                    for (auto& ap_x_dic_iter : ap_x_dic) {
                        if (ap_x_dic_iter.second > max_ap) {
                            max_ap = ap_x_dic_iter.second;
                            pin_dx = ap_x_dic_iter.first;
                        }
                    }

                    int pin_dy = 0;
                    for (auto pin_ap_p : pin_p->apsVp) {
                        int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                        if (ap_x == pin_dx) pin_dy += std::round(pin_ap_p->pointP2.yF * 1000);                            
                    }
                    pin_dy /= max_ap;

                    pin_p->center_ingridP2.xF = static_cast<double>(pin_dx) / 1000;
                    pin_p->center_ingridP2.yF = static_cast<double>(pin_dy) / 1000;

                }        
                
                net_info.push_back(std::make_tuple(net_name, net_index++, net_pin_v, digit_pin));

                fs::path node_feature_pd_path = feature_path / fs::path("node_feature_pd.txt");
                std::ofstream npd_out(node_feature_pd_path, std::ios::app);
                for (int i = 0; i < 9; i++) {
                    npd_out << digit_pin[i] << " ";
                }
                npd_out << std::endl;
                npd_out.close();

            }
            
            std::vector<std::vector<std::vector<double>>> adj_matrix;

            adj_matrix.clear();
            adj_matrix.resize(num_params);
            for (int i = 0; i < num_params; i++) {
                adj_matrix[i].resize(num_net);
                for (int j = 0; j < num_net; j++) adj_matrix[i][j].resize(num_net);
            }

            for (int k = 0; k < num_params; k++)
                for (int i = 0; i < num_net; i++) 
                    for (int j = 0; j < num_net; j++) adj_matrix[k][i][j] = 0;
            

            // edge weight calculation
            for (auto& net_i : net_info) {
                double weight = 0;
                for (auto& net_j : net_info) {
                    if (std::get<1>(net_i) == std::get<1>(net_j)) continue;

                    for (auto& pin_i : std::get<2>(net_i)) {
                        if (pin_i->targetP2 == nullptr) continue;
                        for (auto& pin_j : std::get<2>(net_j)) {
                            bool j_is_left = (pin_i->centerP2.xF >= pin_j->centerP2.xF) ? true : false;

                            // calculate direction value
                            DIR dir;

                            int pin_cx = std::round(pin_i->center_ingridP2.xF * 1000);
                            int pin_cy = std::round(pin_i->center_ingridP2.yF * 1000);

                            int pin_tx = std::round(pin_i->targetP2->xF * 1000);
                            int pin_ty = std::round(pin_i->targetP2->yF * 1000);

                            int dx = pin_tx - pin_cx;
                            int dy = pin_ty - pin_cy;

                            if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                                dir = DIR::SEC;
                            }
                            else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                                dir = DIR::SEC;
                            }
                            else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                                if (j_is_left) dir = DIR::PRIM;
                                else dir = DIR::LAST;
                            }
                            else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                                if (j_is_left) dir = DIR::LAST;
                                else dir = DIR::PRIM;
                            }


                            int c_ix = std::round(pin_i->center_ingridP2.xF * 1000);
                            int c_jx = std::round(pin_j->center_ingridP2.xF * 1000);

                            if (c_ix < grid_bbox.lb.x) c_ix = grid_bbox.lb.x;
                            if (c_ix > grid_bbox.rt.x) c_ix = grid_bbox.rt.x;
                            
                            if (c_jx < grid_bbox.lb.x) c_jx = grid_bbox.lb.x;
                            if (c_jx > grid_bbox.rt.x) c_jx = grid_bbox.rt.x;
                            
                            int n_pin_i_ap = 0;
                            int n_int_ap = 0;
                        
                            for (auto& it_ap : pin_i->apsMIVAp) {
                                bool in_grid = false;
                                for (auto ap_p : it_ap.second) {
                                    int ap_x = std::round(ap_p->pointP2.xF * 1000);
                                    if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                                        in_grid = true;
                                        break;
                                    }
                                }

                                if (in_grid && pin_j->apsMIVAp.count(it_ap.first) > 0) {                                
                                    n_pin_i_ap++;
                                    for (auto pin_pp : pin_j->apsMIVAp[it_ap.first]) {
                                        int ap_x = std::round(pin_pp->pointP2.xF * 1000);
                                        if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                                            n_int_ap++;
                                            break;
                                        }
                                    }
                                } 
                            }

                            for (int k = 0; k < num_params; k++) {
                                double dir_val, dist_factor;

                                if (dir == DIR::PRIM) dir_val = std::get<0>(params[k]);
                                else if (dir == DIR::SEC) dir_val = std::get<1>(params[k]);
                                else if (dir == DIR::LAST) dir_val = std::get<2>(params[k]);
                                else std::cout << "Direction error!" << std::endl;

                                double dist_const = std::get<3>(params[k]);

                                if (j_is_left) dist_factor = static_cast<double>(std::abs(grid_bbox.lb.x - c_ix) - std::abs(c_ix - c_jx) +  dist_const) / (std::abs(grid_bbox.lb.x - c_ix) + dist_const);
                                else dist_factor = static_cast<double>(std::abs(grid_bbox.rt.x - c_ix) - std::abs(c_ix - c_jx) +  dist_const) / (std::abs(grid_bbox.rt.x - c_ix) + dist_const); 

                                double weight = dist_factor * dir_val * n_int_ap / n_pin_i_ap;
       
                                weight *= (ovflHor[y][x] + 1);


                                adj_matrix[k][std::get<1>(net_i)][std::get<1>(net_j)] += weight;
                            }
                        }
                    }
                }
            }
            

            // net - net node edge representation
            for (int k = 0; k < num_params; k++) {
                std::string nf_file_name = "node_feature_ew" + to_string(k + 1) + ".txt";
                fs::path node_feature_ew_path = feature_path / fs::path(nf_file_name);
                std::ofstream new_out(node_feature_ew_path, std::ios::app);

                std::string ei_file_name = "edge_index_" + to_string(k + 1) + ".txt";
                fs::path edge_index_path = feature_path / fs::path(ei_file_name);
                std::ofstream ei_out(edge_index_path, std::ios::app);

                std::string ew_file_name = "edge_weight_" + to_string(k + 1) + ".txt";
                fs::path edge_weight_path = feature_path / fs::path(ew_file_name);
                std::ofstream ew_out(edge_weight_path, std::ios::app);

                
                for (int i = 0; i < num_net; i++) {
                    double weight_sum = 0;
                    int num_edge = 0;
                    double max_weight = -999;
                    for (int j = 0; j < num_net; j++) {
                        if (i == j) continue;
                        if (adj_matrix[k][i][j] > 0) {
                            ei_out << node_index + j << " " << node_index + i << std::endl;
                            ew_out << adj_matrix[k][i][j] << std::endl;

                            num_edge++;
                            weight_sum += adj_matrix[k][i][j];
                            if (adj_matrix[k][i][j] > max_weight) max_weight = adj_matrix[k][i][j];
                        }
                    }
                    double avg_weight = (num_edge != 0) ? (weight_sum / num_edge) : 0;

                    if (num_edge == 0) new_out << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    else new_out << weight_sum << " " << avg_weight << " " << max_weight << " " << num_edge << std::endl;
                }

                new_out.close();
                ei_out.close();
                ew_out.close();
            }
            
            // graph index representation
            for (int i = 0; i < num_net; i++) {
                gi_out << graph_index << std::endl;
            }

            node_index += num_net;
            graph_index++;
        }
    }

    // close output files
    
    //nf_out.close();
    //ei_out.close();
    //ew_out.close();
    gi_out.close();
    //dbg_out.close();

}
*/



/*
void FeatureExtractor::generate_graph() {

    fs::path feature_path = outPath / fs::path(name);

    fs::path node_feature_path = feature_path / fs::path("node_feature_pin_digit.txt");
    fs::path edge_index_path = feature_path / fs::path("edge_index.txt");
    fs::path edge_weight_path = feature_path / fs::path("edge_weight.txt");
    fs::path graph_index_path = feature_path / fs::path("graph_index.txt");
    //fs::path debug_path = feature_path / fs::path("debug.txt");

    // open output files
    std::ofstream nf_out, ei_out, ew_out, gi_out;    
    //std::ofstream dbg_out;

    nf_out.open(node_feature_path.string());
    ei_out.open(edge_index_path.string());
    ew_out.open(edge_weight_path.string());
    gi_out.open(graph_index_path.string());
    //dbg_out.open(debug_path.string());

    int node_index = 0;
    int graph_index = 0;
    
    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {

            //dbg_out << "==================================" << std::endl;
            //dbg_out << "x: " << x << ", y: " << y << std::endl;
            //dbg_out << "==================================" << std::endl;

            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            //dbg_out << "Grid location : " << grid_bbox << std::endl; 

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                for (int i = 0; i < 2; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue; 
            }

            std::vector<pin*> pin_list;
            std::unordered_map<std::string, std::vector<pin*>> net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;  // avoid pins that cannot estimate FLUTE direction
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;

                            if (net_list.count(net_name) == 1) {
                                net_list[net_name].push_back(pinM.second);
                            }
                            else {
                                std::vector<pin*> net_pin_list;
                                net_pin_list.push_back(pinM.second);
                                net_list.insert({net_name, net_pin_list});
                            }
                            break;
                        }
                    }
                }
            }
            
            for (auto& net : net_list) {
                dbg_out << "Net : " << net.first << " / Pin : ";
                for (auto pin_p : net.second) {
                    dbg_out << "(" << pin_p->cellp->compp->nameS << " " << pin_p->nameS << ") ";
                }
                dbg_out << std::endl;
            }

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 2; i++) nf_out << 0 << " ";
                nf_out << std::endl;
                node_index++;
                gi_out << graph_index++ << std::endl;

                continue;                
            }

            // net node representation using FLUTE
            std::vector<std::tuple<std::string, int, std::vector<double>, std::vector<pin*>>> net_info;
            int net_index = 0;
            for (auto& net : net_list) {
                auto net_name = net.first;
                auto& net_pin_v = net.second;

                //dbg_out << "DIR calculation!" << std::endl;
                //dbg_out << "Net : " << net_name << std::endl;

                std::vector<double> app_dir(2); // LR
                for (int j = 0; j < 2; j++) app_dir[j] = 0; 

                for (auto pin_p : net_pin_v) {
                        
                    if (pin_p->targetP2 == nullptr) continue;

                    //dbg_out << "Pin : " << "(" << pin_p->cellp->compp->nameS << " " << pin_p->nameS << ")" << std::endl;

                    std::unordered_map<int, int> ap_x_dic;

                    for (auto pin_ap_p : pin_p->apsVp) {
                        int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                        if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                            if (ap_x_dic.count(ap_x) > 0) ap_x_dic[ap_x]++;
                            else ap_x_dic[ap_x] = 1;                            
                        }
                    }

                    int pin_dx;
                    int max_ap = -1;
                    for (auto& ap_x_dic_iter : ap_x_dic) {
                        if (ap_x_dic_iter.second > max_ap) {
                            max_ap = ap_x_dic_iter.second;
                            pin_dx = ap_x_dic_iter.first;
                        }
                    }

                    int pin_dy = 0;
                    for (auto pin_ap_p : pin_p->apsVp) {
                        int ap_x = std::round(pin_ap_p->pointP2.xF * 1000);
                        if (ap_x == pin_dx) pin_dy += std::round(pin_ap_p->pointP2.yF * 1000);                            
                    }
                    pin_dy /= max_ap;

                    pin_p->center_ingridP2.xF = static_cast<double>(pin_dx) / 1000;
                    pin_p->center_ingridP2.yF = static_cast<double>(pin_dy) / 1000;

                    //dbg_out << "Dominant : (" << pin_dx << ", " << pin_dy << ")" << std::endl;

                    int pin_cx = std::round(pin_p->centerP2.xF * 1000);
                    int pin_cy = std::round(pin_p->centerP2.yF * 1000);

                    int pin_tx = std::round(pin_p->targetP2->xF * 1000);
                    int pin_ty = std::round(pin_p->targetP2->yF * 1000);

                    //dbg_out << "Center : (" << pin_cx << ", " << pin_cy << "), Target : (" << pin_tx << ", " << pin_ty << ")" << std::endl;

                    int dx = pin_tx - pin_dx;
                    int dy = pin_ty - pin_dy;

                    if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                        //if (pin_ty >= grid_bbox.rt.y) {
                            //dbg_out << "UP DIRECTION" << std::endl;
                            app_dir[0] += SEC_DIR; 
                            app_dir[1] += SEC_DIR;
                        //}
                        //else dbg_out << "LOCAL" << std::endl;
                    } 
                    else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                        //if (pin_ty <= grid_bbox.lb.y) {
                            //dbg_out << "DOWN DIRECTION" << std::endl;
                            app_dir[0] += SEC_DIR; 
                            app_dir[1] += SEC_DIR;
                        //} 
                        //else dbg_out << "LOCAL" << std::endl;
                    }
                    else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                        //if (pin_tx <= grid_bbox.lb.x) {
                            //dbg_out << "LEFT DIRECTION" << std::endl;
                            app_dir[0] += PRIM_DIR;
                            app_dir[1] += LAST_DIR;
                        //}
                        //else dbg_out << "LOCAL" << std::endl;
                    }
                    else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                        //if (pin_tx >= grid_bbox.rt.x) {
                            //dbg_out << "RIGHT DIRECTION" << std::endl;
                            app_dir[0] += LAST_DIR;
                            app_dir[1] += PRIM_DIR;
                        //}
                        //else dbg_out << "LOCAL" << std::endl;
                    }

                }
                // Normalize
                double ap_sum = app_dir[0] + app_dir[1];
                if (ap_sum == 0) {
                    app_dir[0] = 0.5;
                    app_dir[1] = 0.5;
                }
                else {
                    app_dir[0] = app_dir[0] / ap_sum;
                    app_dir[1] = app_dir[1] / ap_sum;
                }
                net_info.push_back(std::make_tuple(net_name, net_index++, app_dir, net_pin_v));
                
                //for (int j = 0; j < 2; j++) nf_out << app_dir[j] << " ";
                //nf_out << std::endl;
                std::vector<int> digit_pin(9);
                for (auto pin_ : net_pin_v) {
                    for (auto pin_ap : pin_->apsVp) {
                        int x_ap = std::round(pin_ap->pointP2.xF * 1000);
                        if (grid_bbox.lb.x <= x_ap && x_ap <= grid_bbox.rt.x) {
                            int track = (std::round(pin_ap->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                            digit_pin[track - 2] = 1;        
                        }
                    }
                }
                
                for (int i = 0; i < 9; i++) {
                    nf_out << digit_pin[i] << " ";
                }
                nf_out << std::endl;

                //dbg_out << "Node Feature : " << app_dir[0] << ", " << app_dir[1] << std::endl;
            }
            


            std::vector<std::vector<double>> adj_matrix(num_net, std::vector<double>(num_net));

            for (int i = 0; i < num_net; i++) {
                for (int j = 0; j < num_net; j++) adj_matrix[i][j] = 0;
            }

            // edge weight calculation
            //dbg_out << "Edge Weight Calculation!!" << std::endl;
            for (auto& net_i : net_info) {
                double weight = 0;
                for (auto& net_j : net_info) {
                    if (std::get<1>(net_i) == std::get<1>(net_j)) continue;
                    //dbg_out << std::get<0>(net_i) << " <-- " << std::get<0>(net_j) << std::endl;

                    for (auto& pin_i : std::get<3>(net_i)) {
                        if (pin_i->targetP2 == nullptr) continue;
                        for (auto& pin_j : std::get<3>(net_j)) {
                            bool j_is_left = (pin_i->centerP2.xF >= pin_j->centerP2.xF) ? true : false;

                            // calculate direction value
                            double dir_val = (j_is_left) ? std::get<2>(net_i)[0] : std::get<2>(net_i)[1];

                            int pin_cx = std::round(pin_i->center_ingridP2.xF * 1000);
                            int pin_cy = std::round(pin_i->center_ingridP2.yF * 1000);

                            int pin_tx = std::round(pin_i->targetP2->xF * 1000);
                            int pin_ty = std::round(pin_i->targetP2->yF * 1000);

                            int dx = pin_tx - pin_cx;
                            int dy = pin_ty - pin_cy;

                            if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                                dir_val = SEC_DIR;
                            }
                            else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                                dir_val = SEC_DIR;
                            }
                            else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                                if (j_is_left) dir_val = PRIM_DIR;
                                else dir_val = LAST_DIR;
                            }
                            else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                                if (j_is_left) dir_val = LAST_DIR;
                                else dir_val = PRIM_DIR;
                            }


                            int c_ix = std::round(pin_i->center_ingridP2.xF * 1000);
                            int c_jx = std::round(pin_j->center_ingridP2.xF * 1000);

                            if (c_ix < grid_bbox.lb.x) c_ix = grid_bbox.lb.x;
                            if (c_ix > grid_bbox.rt.x) c_ix = grid_bbox.rt.x;
                            
                            if (c_jx < grid_bbox.lb.x) c_jx = grid_bbox.lb.x;
                            if (c_jx > grid_bbox.rt.x) c_jx = grid_bbox.rt.x;
                            
                            int n_pin_i_ap = 0;
                            int n_int_ap = 0;
                        
                            for (auto& it_ap : pin_i->apsMIVAp) {
                                bool in_grid = false;
                                for (auto ap_p : it_ap.second) {
                                    int ap_x = std::round(ap_p->pointP2.xF * 1000);
                                    if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                                        in_grid = true;
                                        break;
                                    }
                                }

                                if (in_grid && pin_j->apsMIVAp.count(it_ap.first) > 0) {                                
                                    n_pin_i_ap++;
                                    for (auto pin_pp : pin_j->apsMIVAp[it_ap.first]) {
                                        int ap_x = std::round(pin_pp->pointP2.xF * 1000);
                                        if (ap_x >= grid_bbox.lb.x && ap_x <= grid_bbox.rt.x) {
                                            n_int_ap++;
                                            break;
                                        }
                                    }
                                } 
                            }

                            double dist_factor;

                            if (j_is_left) {
                                dist_factor = static_cast<double>(std::abs(grid_bbox.lb.x - c_ix) - std::abs(c_ix - c_jx) +  DIST_CONST) / (std::abs(grid_bbox.lb.x - c_ix) + DIST_CONST); 
                            }
                            else dist_factor = static_cast<double>(std::abs(grid_bbox.rt.x - c_ix) - std::abs(c_ix - c_jx) +  DIST_CONST) / (std::abs(grid_bbox.rt.x - c_ix) + DIST_CONST); 

                            //dbg_out << "Center x of pin " << pin_i->nameS << " : " << c_ix << ", Center x of pin " << pin_j->nameS << " : " << c_jx << std::endl;
                            //dbg_out << "Total pin track of pin " << pin_i->nameS << " : " << n_pin_i_ap << ", intersected pin : " << n_int_ap << std::endl;
                            //dbg_out << "Distance factor : " << dist_factor << std::endl;
                            //dbg_out << "Direction value : " << dir_val << std::endl;

                            double weight = dist_factor * dir_val * n_int_ap / n_pin_i_ap;
                            adj_matrix[std::get<1>(net_i)][std::get<1>(net_j)] += weight;

                            //dbg_out << "Calculated weight : " << weight << std::endl;
                        }
                    }
                }
            }

            // net - net node edge representation
            for (int i = 0; i < num_net; i++) {
                for (int j = 0; j < num_net; j++) {
                    if (i == j) continue;
                    if (adj_matrix[i][j] > 0) {
                        ei_out << node_index + j << " " << node_index + i << std::endl;
                        ew_out << adj_matrix[i][j] << std::endl;

                        //dbg_out << j << "-->" << i << " " << adj_matrix[i][j] << std::endl;
                    }
                }
            }
            
            // graph index representation
            for (int i = 0; i < num_net; i++) {
                gi_out << graph_index << std::endl;
            }

            node_index += num_net;
            graph_index++;
        }
    }

    // close output files
    nf_out.close();
    ei_out.close();
    ew_out.close();
    gi_out.close();
    //dbg_out.close();

}
*/


/*
void FeatureExtractor::generate_graph_new() {

    fs::path feature_path = outPath / fs::path(name);

    fs::path pin_net_edge_path = feature_path / fs::path("pin_net_edge.txt");
    fs::path net_net_edge_path = feature_path / fs::path("net_net_edge.txt");
    fs::path net_pin_edge_path = feature_path / fs::path("net_pin_edge.txt");
    fs::path pin_node_feature_path = feature_path / fs::path("pin_node_feature.txt");
    fs::path net_node_feature_path = feature_path / fs::path("net_node_feature.txt");
    fs::path net_graph_index_path = feature_path / fs::path("net_graph_index.txt");

    // open output files
    std::ofstream pne_out, nne_out, npe_out, pnf_out, nnf_out, ngi_out;
    pne_out.open(pin_net_edge_path.string());
    nne_out.open(net_net_edge_path.string());
    npe_out.open(net_pin_edge_path.string());
    pnf_out.open(pin_node_feature_path.string());
    nnf_out.open(net_node_feature_path.string());
    ngi_out.open(net_graph_index_path.string());

    int pin_node_index = 0;
    int net_node_index = 0;
    int graph_index = 0;
    
    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {
            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                for (int i = 0; i < 4; i++) nnf_out << 0 << " ";
                nnf_out << std::endl;
                net_node_index++;
                ngi_out << graph_index++ << std::endl;

                continue; 
            }

            std::vector<pin*> pin_list;
            std::vector<std::string> net_list;
            std::vector<std::string> local_net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;
                            if (std::find(net_list.begin(), net_list.end(), net_name) != net_list.end()){
                                if (std::find(local_net_list.begin(), local_net_list.end(), net_name) == local_net_list.end()) {
                                    local_net_list.push_back(net_name);
                                }
                            }
                            else net_list.push_back(net_name);
                            break;
                        }
                    }
                }
            }

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 4; i++) nnf_out << 0 << " ";
                nnf_out << std::endl;
                net_node_index++;
                ngi_out << graph_index++ << std::endl;

                continue;                
            }

            // pin node feature representation
            // Assume Nangate 15nm (12 track) but access point can be made in center 9 tracks
            std::vector<std::tuple<pin*, std::vector<int>, double>> pin_info;
            for (auto pin_ : pin_list) {
                //  dir (INPUT/OUTPUT) + pin access point vector + average x-coord
                int num_ap = 0;
                double sum_x_ap = 0;
                std::vector<int> digit_pin(9);
                for (auto pin_ap : pin_->apsVp) {
                    int x_ap = std::round(pin_ap->pointP2.xF * 1000);
                    if (grid_bbox.lb.x <= x_ap && x_ap <= grid_bbox.rt.x) {
                        num_ap++;
                        sum_x_ap += static_cast<double>(x_ap - grid_bbox.lb.x) / gridXI;

                        int track = (std::round(pin_ap->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                        digit_pin[track - 2] = 1;        
                    }
                }
                pin_info.push_back(std::make_tuple(pin_, digit_pin, sum_x_ap / num_ap));


            }

            std::sort(pin_info.begin(), pin_info.end(), 
                      [](std::tuple<pin*, std::vector<int>, double> const &t1, std::tuple<pin*, std::vector<int>, double> const &t2) {
                          return std::get<2>(t1)< std::get<2>(t2);
                      }
            );

            for (auto& pin_info_i : pin_info) {
                pnf_out << std::get<2>(pin_info_i) << " ";
                for (int i = 0; i < 9; i++) pnf_out << std::get<1>(pin_info_i)[i] << " ";
                pnf_out << std::endl;
            }


            // net node representation using FLUTE
            std::vector<std::tuple<std::string, std::vector<int>, std::vector<int>>> net_info;
            for (int i = 0; i < num_net; i++) {
                std::string net_name = net_list[i];

                std::vector<int> app_dir, con_pin_index;
                app_dir.resize(4);
                for (int j = 0; j < 4; j++) app_dir[j] = 0;

                for (int j = 0; j < num_pin; j++) {
                    auto pin_p = std::get<0>(pin_info[j]);
                    if (pin_p->netp->nameS == net_name) {
                        
                        if (pin_p->targetP2 == nullptr) continue;

                        con_pin_index.push_back(j);
                        int pin_cx = std::round(pin_p->centerP2.xF * 1000);
                        int pin_cy = std::round(pin_p->centerP2.yF * 1000);

                        int pin_tx = std::round(pin_p->targetP2->xF * 1000);
                        int pin_ty = std::round(pin_p->targetP2->yF * 1000);

                        int dx = pin_tx - pin_cx;
                        int dy = pin_ty - pin_cy;

                        if (dy > 0 && std::abs(dx) < std::abs(dy)) { // up
                            if (pin_ty >= grid_bbox.rt.y) app_dir[0] = 1;
                        } 
                        else if (dy < 0 && std::abs(dx) < std::abs(dy)) { // down
                            if (pin_ty <= grid_bbox.lb.y) app_dir[1] = 1;
                        }
                        else if (dx < 0 && std::abs(dx) > std::abs(dy)) { // left
                            if (pin_tx <= grid_bbox.lb.x) app_dir[2] = 1;
                        }
                        else if (dx > 0 && std::abs(dx) > std::abs(dy)) { // right
                            if (pin_tx >= grid_bbox.rt.x) app_dir[3] = 1;
                        }
                    }
                }
                net_info.push_back(std::make_tuple(net_name, app_dir, con_pin_index));
                
                for (int j = 0; j < 4; j++) nnf_out << app_dir[j] << " ";
                nnf_out << std::endl;
            }

            // pin - net / net - pin node edge representation
            std::vector<int> pin_net_index;
            pin_net_index.resize(num_pin);

            for (int i = 0; i < num_net; i++) {
                auto &net_info_i = net_info[i];
                int con_pin_num = std::get<2>(net_info_i).size();

                for (int j = 0; j < con_pin_num; j++) {
                    int pin_index = std::get<2>(net_info_i)[j];
                    pin_net_index[j] = i;
                    assert(std::get<0>(pin_info[pin_index])->netp->nameS == std::get<0>(net_info_i));

                    pne_out << pin_index + pin_node_index << " " << i + net_node_index << std::endl;
                    npe_out << i + net_node_index << " " << pin_index + pin_node_index << std::endl;
                }
            }



            // net - net node edge representation
            for (int i = 0; i < num_net - 1; i++) {

                BBox<int> net_i{0, 0, 0, 0};
                for (int k = 0; k < std::get<2>(net_info[i]).size(); k++) {
                    int pin_index = std::get<2>(net_info[i])[k];
                    auto pin_p = std::get<0>(pin_info[pin_index]);

                    int pin_cx = std::round(pin_p->centerP2.xF * 1000);
                    int pin_cy = std::round(pin_p->centerP2.yF * 1000);

                    int pin_tx = std::round(pin_p->targetP2->xF * 1000);
                    int pin_ty = std::round(pin_p->targetP2->yF * 1000);

                    BBox<int> pin_k{pin_cx, pin_cy, pin_tx, pin_ty};
                        
                    if (net_i.lb.x == 0) net_i = pin_k;
                    else net_i = net_i + pin_k;
                
                }
                for (int j = i + 1; j < num_net; j++) {
                    
                    BBox<int> net_j{0, 0, 0, 0};

                    for (int k = 0; k < std::get<2>(net_info[j]).size(); k++) {
                        int pin_index = std::get<2>(net_info[j])[k];
                        auto pin_p = std::get<0>(pin_info[pin_index]);

                        int pin_cx = std::round(pin_p->centerP2.xF * 1000);
                        int pin_cy = std::round(pin_p->centerP2.yF * 1000);

                        int pin_tx = std::round(pin_p->targetP2->xF * 1000);
                        int pin_ty = std::round(pin_p->targetP2->yF * 1000);

                        BBox<int> pin_k{pin_cx, pin_cy, pin_tx, pin_ty};          

                        if (net_j.lb.x == 0) net_j = pin_k;
                        else net_j = net_j + pin_k;
                    }

                    if ((net_i.lb.x <= net_j.lb.x && net_j.lb.x <= net_i.rt.x) || (net_j.lb.x <= net_i.lb.x && net_i.lb.x <= net_j.rt.x)) {
                        nne_out << i + net_node_index << " " << j + net_node_index << std::endl;
                        nne_out << j + net_node_index << " " << i + net_node_index << std::endl;
                    }


                }
            }
            
            // graph index representation
            for (int i = 0; i < num_net; i++) {
                ngi_out << graph_index << std::endl;
            }

            pin_node_index += num_pin;
            net_node_index += num_net;
            graph_index++;
        }
    }

    // close output files
    pne_out.close();
    nne_out.close();
    npe_out.close();
    pnf_out.close();
    nnf_out.close();
    ngi_out.close();

}

*/

/*
void FeatureExtractor::generate_graph() {

    fs::path feature_path = outPath / fs::path(name);

    fs::path pin_net_edge_path = feature_path / fs::path("pin_net_edge.txt");
    fs::path net_net_edge_path = feature_path / fs::path("net_net_edge.txt");
    fs::path net_pin_edge_path = feature_path / fs::path("net_pin_edge.txt");
    fs::path pin_node_feature_path = feature_path / fs::path("pin_node_feature.txt");
    fs::path net_node_feature_path = feature_path / fs::path("net_node_feature.txt");
    fs::path net_graph_index_path = feature_path / fs::path("net_graph_index.txt");

    // open output files
    std::ofstream pne_out, nne_out, npe_out, pnf_out, nnf_out, ngi_out;
    pne_out.open(pin_net_edge_path.string());
    nne_out.open(net_net_edge_path.string());
    npe_out.open(net_pin_edge_path.string());
    pnf_out.open(pin_node_feature_path.string());
    nnf_out.open(net_node_feature_path.string());
    ngi_out.open(net_graph_index_path.string());

    int pin_node_index = 0;
    int net_node_index = 0;
    int graph_index = 0;
    
    for (int y = 0; y < gridNum.y; y++) {
        for (int x = 0; x < gridNum.x; x++) {
            int lx = x;
            int ly = gridNum.y - y - 1;

            BBox<int> grid_bbox{lx * gridXI, ly * gridYI, (lx + 1) * gridXI, (ly + 1) * gridYI};
            grid_bbox = grid_bbox + coreI[0];

            if (cellMap[ly][lx].size() == 0) {
                // if no cells are in the grid
                // create empty node
                for (int i = 0; i < 4; i++) nnf_out << 0 << " ";
                nnf_out << std::endl;
                net_node_index++;
                ngi_out << graph_index++ << std::endl;

                continue; 
            }

            bool for_debug = false;
            for (auto cell_name : cellMap[ly][lx]) {
                if (cell_name == "P1/P1/U2142") {
                    for_debug = true;
                    break;
                }
            }

            std::vector<pin*> pin_list;
            std::vector<std::string> net_list;
            std::vector<std::string> local_net_list;

            // iterate pin list and save all pins overlapped with the grid
            for (auto comp_name : cellMap[ly][lx]) {
                auto comp_p = def.compsMSCp[comp_name];
                auto cell_p = comp_p -> cellCp;

                // find pin and nets accessing to the grid
                for (auto& pinM : cell_p->pinsMSPp) {
                    if (pinM.first == "VDD" || pinM.first == "VSS") continue;
                    for (auto aps : pinM.second->apsVp){
                        if (aps->layerS != "M1") continue; // just M1 pin
                        int ap_x = std::round(aps->pointP2.xF * 1000);
                        int ap_y = std::round(aps->pointP2.yF * 1000);
                        if (grid_bbox.lb.x <= ap_x && ap_x <= grid_bbox.rt.x && grid_bbox.lb.y <= ap_y && ap_y <= grid_bbox.rt.y) {
                            
                            if (pinM.second->netp == nullptr) continue;
                            pin_list.push_back(pinM.second);
                            std::string net_name = pinM.second->netp->nameS;
                            if (std::find(net_list.begin(), net_list.end(), net_name) != net_list.end()){
                                if (std::find(local_net_list.begin(), local_net_list.end(), net_name) == local_net_list.end()) {
                                    local_net_list.push_back(net_name);
                                }
                            }
                            else net_list.push_back(net_name);
                            break;
                        }
                    }
                }
            }

            int num_pin = pin_list.size();
            int num_net = net_list.size();
            
            if (num_pin == 0) {
                // if no pin is in the grid
                // create empty node
                for (int i = 0; i < 4; i++) nnf_out << 0 << " ";
                nnf_out << std::endl;
                net_node_index++;
                ngi_out << graph_index++ << std::endl;

                continue;                
            }

            // pin node feature representation
            // Assume Nangate 15nm (12 track) but access point can be made in center 9 tracks
            std::vector<std::tuple<bool, std::vector<int>, double>> pin_info;
            for (auto pin_ : pin_list) {
                //  dir (INPUT/OUTPUT) + pin access point vector + average x-coord
                int num_ap = 0;
                double sum_x_ap = 0;
                std::vector<int> digit_pin(9);
                for (auto pin_ap : pin_->apsVp) {
                    int x_ap = std::round(pin_ap->pointP2.xF * 1000);
                    if (grid_bbox.lb.x <= x_ap && x_ap <= grid_bbox.rt.x) {
                        num_ap++;
                        sum_x_ap += static_cast<double>(x_ap - grid_bbox.lb.x) / gridXI;

                        int track = (std::round(pin_ap->pointP2.yF * 1000) - grid_bbox.lb.y) / NAN_M2_PITCH;
                        digit_pin[track - 2] = 1;        
                    }
                }
                pnf_out << sum_x_ap / num_ap << " ";
                for (int i = 0; i < 9; i++) pnf_out << digit_pin[i] << " ";
                pnf_out << std::endl;

                if (for_debug) {
                    std::cout << sum_x_ap / num_ap << " ";
                    for (int i = 0; i < 9; i++) std::cout << digit_pin[i] << " ";
                    std::cout << std::endl;
                }

            }

            // net node representation
            int grid_cx = (grid_bbox.lb.x + grid_bbox.rt.x) / 2;
            int grid_cy = (grid_bbox.lb.y + grid_bbox.rt.y) / 2;

            for (int i = 0; i < num_net; i++) {
                std::string net_name = net_list[i];

                // Calculate approaching direction
                int con_pin_num = 0;
                double pin_dist[4] = {0};    // up, down, left, right

                for (auto& pin_info : def.netsMSNp[net_name]->cellPinMSS) {
                    auto comp_name = pin_info.first;
                    auto pin_name = pin_info.second;

                    if (std::find(cellMap[ly][lx].begin(), cellMap[ly][lx].end(), comp_name) != cellMap[ly][lx].end()) {
                        continue;
                    }

                    int pin_cx, pin_cy;
                    
                    if (comp_name == "PIN") {
                        pin_cx = std::round(def.IOpinsMSpp[pin_name]->centerP2.xF * 1000);
                        pin_cy = std::round(def.IOpinsMSpp[pin_name]->centerP2.yF * 1000);
                    } 
                    else {
                        pin_cx = std::round(def.compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name]->centerP2.xF * 1000);
                        pin_cy = std::round(def.compsMSCp[comp_name]->cellCp->pinsMSPp[pin_name]->centerP2.yF * 1000);
                    }

                    int dx = pin_cx - grid_cx;
                    int dy = pin_cy - grid_cy;

                    if (dy > 0 && std::abs(dx) < std::abs(dy)) pin_dist[0]++; // up
                    else if (dy < 0 && std::abs(dx) < std::abs(dy)) pin_dist[1]++; // down
                    else if (dx < 0 && std::abs(dx) > std::abs(dy)) pin_dist[2]++; // left
                    else if (dx > 0 && std::abs(dx) > std::abs(dy)) pin_dist[3]++; // right
                    con_pin_num++;
                }

                if (con_pin_num > 0) {
                    pin_dist[0] /= con_pin_num; 
                    pin_dist[1] /= con_pin_num;
                    pin_dist[2] /= con_pin_num;
                    pin_dist[3] /= con_pin_num;
                }

                for (int i = 0; i < 4; i++) nnf_out << pin_dist[i] << " ";
                nnf_out << std::endl;

                if (for_debug) {
                    for (int i = 0; i < 4; i++) std::cout << pin_dist[i] << " ";
                    std::cout << std::endl;
                }
            }

            // pin - net / net - pin node edge representation
            for (int i = 0; i < num_pin; i++) {
                auto pin_net_name = pin_list[i]->netp->nameS;
                for (int j = 0; j < num_net; j++) {
                    auto net_name = net_list[j];

                    if (pin_net_name == net_name) {
                        pne_out << i + pin_node_index << " " << j + net_node_index << std::endl;
                        npe_out << j + net_node_index << " " << i + pin_node_index << std::endl;


                        if (for_debug) {
                            std::cout << i << "(pin) " << j << "(net)" << std::endl;

                        }

                        break;
                    }
                }
            }

            // net - net node edge representation
            for (int i = 0; i < num_net - 1; i++) {
                for (int j = i + 1; j < num_net; j++) {
                    nne_out << i + net_node_index << " " << j + net_node_index << std::endl;
                    nne_out << j + net_node_index << " " << i + net_node_index << std::endl;
                }
            }
            
            // graph index representation
            for (int i = 0; i < num_net; i++) {
                ngi_out << graph_index << std::endl;
            }

            pin_node_index += num_pin;
            net_node_index += num_net;
            graph_index++;
        }
    }

    // close output files
    pne_out.close();
    nne_out.close();
    npe_out.close();
    pnf_out.close();
    nnf_out.close();
    ngi_out.close();

}

*/
