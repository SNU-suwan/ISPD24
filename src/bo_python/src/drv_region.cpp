#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <cmath>
#include <tuple>
#include <numeric>
#include <map>
#include <random>
#include <ctime>
#include "utils.h"
#include "parser.h"
#include "drv_region.h"

int MAX_RANGE_FOR_OVERFLOW = 100000;

RowsOfDRVRegion::RowsOfDRVRegion(){}
RowsOfDRVRegion::~RowsOfDRVRegion(){
    for(size_t i=0; i<rows_of_drv_region.size(); i++){
        delete rows_of_drv_region[i];
    }
}

RowsOfDRVRegion* ParseDrvRegionLog(
    const std::string& script_directory,
    const ParserWrapper& lef_def_parser,
    const int& gcell_scale
    ){
    
    vector<RowInDRVRegion*> drv_regions;
    int grid_side = lef_def_parser.def.row_height * gcell_scale;
    char line[1024];

    int counter = 0; // 0: label, y, 1: x_left, 2: x_right
   
    std::ifstream ifs;
    ifs.open(script_directory);
    if(ifs.is_open()){
		while(!ifs.eof()){
			ifs.getline(line,1024);
			char *token = strtok(line," ");
           
            counter = 0;

            int label = -1;
            int grid_y = -1;
            int grid_x_left = -1;
            int grid_x_right = -1;

			while(token !=NULL){
                std::string tmp_str(token);
                if(counter == 0){
                    label = atoi(token);
                    token = strtok(NULL," ");
                    grid_y = atoi (token);
                    counter = 1;
                }
                else if(counter == 1){
                    grid_x_left = atoi(token);
                    counter = 2;
                }
                else if (counter == 2){
                    grid_x_right = atoi(token);
                    counter = 1;

                    if(REGION_IDX == -1 || (REGION_IDX != -1 && label == REGION_IDX)){
                        RowInDRVRegion* new_row_in_drv_region = new RowInDRVRegion;
                        new_row_in_drv_region->label = label;
                        assert(label != -1);

                        new_row_in_drv_region->y = ConvertGridToChip(grid_y, gcell_scale, false, lef_def_parser);
                        new_row_in_drv_region->x_left = ConvertGridToChip(grid_x_left, gcell_scale, true, lef_def_parser);
                        new_row_in_drv_region->x_right = ConvertGridToChip(grid_x_right, gcell_scale, true, lef_def_parser) + grid_side;
                        new_row_in_drv_region->x_right = min(new_row_in_drv_region->x_right, lef_def_parser.def.chip_x + lef_def_parser.def.offset.x);
                        
                        new_row_in_drv_region->grid_y = grid_y;
                        new_row_in_drv_region->grid_x_left = grid_x_left;
                        new_row_in_drv_region->grid_x_right = grid_x_right + 1;

                        drv_regions.push_back(new_row_in_drv_region);
                    }
                    
                }
                token = strtok(NULL," ");
            }
            
        }
    }
    ifs.close();

    RowsOfDRVRegion *rows_in_drv_region = new RowsOfDRVRegion;
    rows_in_drv_region->rows_of_drv_region = drv_regions;
    return rows_in_drv_region;

}

void RowsOfDRVRegion::Preprocess(const ParserWrapper& lef_def_parser, const int& num_candidates, const bool& is_cts_design){
    for(size_t i = 0; i < rows_of_drv_region.size(); i++){
        rows_of_drv_region[i]->Initialize(lef_def_parser, num_candidates, is_cts_design);
    }
}



void RowsOfDRVRegion::MatchCellsWithDrvRegions(const ParserWrapper& lef_def_parser){
    vector<size_t> index_to_be_deleted;
    for(size_t i = 0; i < rows_of_drv_region.size(); i++){
        RowInDRVRegion*& cur_row_in_drv_region = rows_of_drv_region[i];
        
        const int& y = cur_row_in_drv_region->y;
        const int& x_left = cur_row_in_drv_region->x_left;
        const int& x_right = cur_row_in_drv_region->x_right;

        const std::vector<placed_cell*>& cells_in_a_row = lef_def_parser.def.placed_cells_row_map.at(y);
        for(size_t i_c = 0; i_c < cells_in_a_row.size(); i_c++){
            placed_cell* cur_cell = cells_in_a_row.at(i_c);
            int displacement_left = cur_cell->x - x_left;
            int displacement_right = x_right - (cur_cell->x + cur_cell->cell->cell_width);
            if(displacement_left >= 0 && displacement_right >= 0){ 
                cur_cell->is_fixed = 0;
                cur_row_in_drv_region->cells.push_back(cur_cell);
            }
            else {
                if(cur_cell->is_fixed == -1){
                    cur_cell->is_fixed = 1;
                }
            }

            //If a cell overlaps with the boundary of gcell, set the new boundary as the side of the cell.
            if(displacement_left < 0 && (cur_cell->x + cur_cell->cell->cell_width) - x_left > 0){    
                cur_row_in_drv_region->x_left += (cur_cell->x + cur_cell->cell->cell_width) - x_left;
            }
            if(displacement_right < 0 && x_right - cur_cell->x > 0){
                cur_row_in_drv_region->x_right -= x_right - cur_cell->x;
            }
        }

        if(cur_row_in_drv_region->cells.size() == 0){
            index_to_be_deleted.push_back(i);
        }
    }
    sort(index_to_be_deleted.begin(), index_to_be_deleted.end(), greater<size_t>());
    for(int i = 0 ; i < index_to_be_deleted.size(); i++){
        delete rows_of_drv_region[index_to_be_deleted[i]];
        rows_of_drv_region.erase(rows_of_drv_region.begin()+index_to_be_deleted[i]);
    }
}

void RowInDRVRegion::SeparateCTSCells(){
    vector<size_t> index_to_be_deleted;
    for(size_t i_c = 0; i_c < cells.size(); i_c++){
        placed_cell * cell = cells[i_c];
        if(cell->is_cts_cell){
            cts_cells.push_back(cell);
            index_to_be_deleted.push_back(i_c);
        }
    }

    sort(index_to_be_deleted.begin(), index_to_be_deleted.end(), greater<int>());
    for(int i_node_to_be_deleted = 0 ; i_node_to_be_deleted < index_to_be_deleted.size(); i_node_to_be_deleted++){
        cells.erase(cells.begin()+index_to_be_deleted[i_node_to_be_deleted]);
    }
}



void RowsOfDRVRegion::SeparateCTSCells(){
    for(size_t i = 0; i < rows_of_drv_region.size(); i++){
        rows_of_drv_region[i]->SeparateCTSCells();
    }
}

void RowInDRVRegion::AssignInvalidDistanceToCells(const ParserWrapper& lef_def_parser, const int& delta){
    
    vector<size_t> cts_cell_idx;
    for(size_t i = 0; i < cells.size(); i++){
        if(cells[i]->is_cts_cell) cts_cell_idx.push_back(i);
    }

    if(cts_cell_idx.size() != 0){
        for(size_t i = 0 ; i < cells.size(); i++){
            if(cells[i]->is_cts_cell) continue;
            placed_cell* cell = cells[i];

            for(size_t i_c = 0 ; i_c < cts_cell_idx.size(); i_c++){
                placed_cell* cts_cell = cells[cts_cell_idx[i_c]];
                int offset_whitespace = cts_cell->x - cell->x;
                for(int i_c_before = 0 ; i_c_before < cts_cell_idx.size(); i_c_before ++){
                    if(cells[cts_cell_idx[i_c_before]]->x < cell->x){
                        //CTS cells become whitespace
                        offset_whitespace += cells[cts_cell_idx[i_c_before]]->cell->cell_width;
                    }
                }
                offset_whitespace /= lef_def_parser.def.cpp;
                int original_cumulative_whitespace = 0;
                for(int i_w = 0; i_w <= i ; i_w ++){
                    original_cumulative_whitespace += cells_original_whitespace[i_w];
                }
                offset_whitespace += original_cumulative_whitespace;
                
                int left_offset = 1 - (cell->cell->cell_width / lef_def_parser.def.cpp);
                int right_offset = (cts_cell->cell->cell_width / lef_def_parser.def.cpp) - 1;

                for(int i_w = left_offset; i_w <= right_offset; i_w ++){
                    if (-1 * delta <= offset_whitespace + i_w - original_cumulative_whitespace 
                        && offset_whitespace + i_w - original_cumulative_whitespace <= delta){
                        cell->invalid_cumulative_whitespace.push_back(offset_whitespace + i_w); 
                    }
                }
            }
        }
    }
}

void RowsOfDRVRegion::SaveOriginalWhitespace(const ParserWrapper& lef_def_parser){
    for(size_t i = 0; i<rows_of_drv_region.size(); i++){
        rows_of_drv_region[i]->SaveOriginalWhitespace(lef_def_parser);
    }
}

void RowsOfDRVRegion::AssignInvalidDistanceToCells(const ParserWrapper& lef_def_parser, const int& delta){
    for(size_t i = 0; i<rows_of_drv_region.size(); i++){
        rows_of_drv_region[i]->AssignInvalidDistanceToCells(lef_def_parser, delta);
    }
}

RowsOfDRVRegion* DecodeDrvRegion(
    const std::string& script_directory,
    const ParserWrapper& lef_def_parser,
    const int& delta,
    const bool& is_cts_design,
    const int& gcell_scale
    ){

    RowsOfDRVRegion* drv_regions = ParseDrvRegionLog(script_directory, lef_def_parser, gcell_scale);
    drv_regions->MatchCellsWithDrvRegions(lef_def_parser);
    
    return drv_regions;
} 

void ClusteredDRVRegions1Dimension::Destroy(){
    for(size_t i = 0; i < rows_of_drv_regions.size(); i++){
        delete rows_of_drv_regions[i];
    }
}

void ClusteredDRVRegions1Dimension::Clear(){
    total_range = 1;
    num_cells = 0;
    ranges.clear();
    rows_of_drv_regions.clear();
    total_predicted_value = 0;
}


void ClusteredDRVRegions1Dimension::Preprocess(){
    num_cells = 0;
    ranges.clear();
    int num_cases;
    int debug_total_range = 1;
    for(size_t i = 0; i < rows_of_drv_regions.size(); i++){
        RowInDRVRegion* row_of_drv_region = rows_of_drv_regions[i];
        
        if(is_shift){
            num_cases = row_of_drv_region->cases_shift;
        }
        else{
            num_cases = row_of_drv_region->cases_flip;
        }
        debug_total_range *= num_cases;
        ranges.push_back(num_cases);

        num_cells += row_of_drv_region->cells.size();
    }
    assert(debug_total_range == total_range);
}

void ClusteredDRVRegions1Dimension::AddDRVRegion(RowInDRVRegion *row_in_drv_region){
    num_cells += row_in_drv_region->cells.size();
    rows_of_drv_regions.push_back(row_in_drv_region);
    if(is_shift){
        total_range *= row_in_drv_region->cases_shift;
    }
    else{
        total_range *= row_in_drv_region->cases_flip;
    }
}

RowInDRVRegion::RowInDRVRegion(RowInDRVRegion& _row_in_drv_region){
    y = _row_in_drv_region.y;
    x_left = _row_in_drv_region.x_left;
    x_right = _row_in_drv_region.x_right;

    grid_y = _row_in_drv_region.grid_y;
    grid_x_left = _row_in_drv_region.grid_x_left;
    grid_x_right = _row_in_drv_region.grid_x_right;

    cases_shift = _row_in_drv_region.cases_shift;
    cells = _row_in_drv_region.cells;
    
    predicted_value = _row_in_drv_region.predicted_value;
    delta_min = _row_in_drv_region.delta_min;
}


tuple<vector<int>, int, int> GetWhitespaceVector(const RowInDRVRegion* row_in_drv_region, const ParserWrapper& lef_def_parser){
    vector<int> whitespace_distribution_vector;
    int total_cell_width = 0;
    int total_whitespace = 0;

    const vector<placed_cell*> &cells = row_in_drv_region->cells;

    int previous_right = row_in_drv_region->x_left;
    for(size_t i = 0; i < cells.size(); i++){
        int cell_position = cells[i]->x;
        int whitespace = (cell_position - previous_right) / lef_def_parser.def.cpp;
        assert(whitespace >= 0);
        whitespace_distribution_vector.push_back(whitespace);
        previous_right = cells[i]->x + cells[i]->cell->cell_width;
        total_cell_width += cells[i]->cell->cell_width;
        total_whitespace += whitespace;
    }
    for(size_t i = 0; i < row_in_drv_region->cts_cells.size(); i++){
        total_cell_width += row_in_drv_region->cts_cells[i]->cell->cell_width;
        total_whitespace -= row_in_drv_region->cts_cells[i]->cell->cell_width;
        previous_right = max(previous_right, row_in_drv_region->cts_cells[i]->x + row_in_drv_region->cts_cells[i]->cell->cell_width);
    }
    int last_whitespace = (row_in_drv_region->x_right - previous_right) / lef_def_parser.def.cpp;
    assert(last_whitespace >= 0);
    whitespace_distribution_vector.push_back(last_whitespace);
    total_whitespace += last_whitespace;

    tuple<vector<int>, int, int> whitespace_information = make_tuple(whitespace_distribution_vector, total_cell_width, total_whitespace);
    return whitespace_information;
}

void RowInDRVRegion::SaveOriginalWhitespace(const ParserWrapper& lef_def_parser){
    cells_original_whitespace.clear();
    shift_graph.layers.clear();
    encode_decode_LUT.shift_decode_map.clear();
    encode_decode_LUT.shift_encode_map.clear();

    int total_space = x_right - x_left;
    total_cell_width = 0;
    total_whitespace = 0;
    
    tuple<vector<int>, int, int> whitespace_information = GetWhitespaceVector(this, lef_def_parser);
    cells_original_whitespace = get<0>(whitespace_information);
    total_cell_width = get<1>(whitespace_information);
    total_whitespace = get<2>(whitespace_information);

    //Get whitespace information
    int num_slots = 1;
    int layer_idx = 0;
    int cumulative_whitespace = cells_original_whitespace[0];
    for(int i = 0; i <= int(cells_original_whitespace.size())-2; i++){
        int &current_whitespace = cells_original_whitespace[i];
        int &next_whitespace = cells_original_whitespace[i+1];
        
        int cumulative_slots = i + 1;
        Layer* new_layer = new Layer(cumulative_whitespace, num_slots, cumulative_slots, layer_idx++);
        shift_graph.layers.push_back(new_layer);
        num_slots = 1;

        cumulative_whitespace += next_whitespace;
    }
    assert(cumulative_whitespace == total_whitespace);
    Layer* new_layer = new Layer(cumulative_whitespace, num_slots, cells_original_whitespace.size(), layer_idx++);
    shift_graph.layers.push_back(new_layer);

    encode_decode_LUT.shift_decode_map.insert(make_pair(0,cells_original_whitespace));
    encode_decode_LUT.shift_encode_map.insert(make_pair(cells_original_whitespace,0));

    assert(total_whitespace == (total_space - total_cell_width) / lef_def_parser.def.cpp);
}

void RowInDRVRegion::SaveOriginalFlip(){
    encode_decode_LUT.flip_decode_map.clear();
    encode_decode_LUT.flip_encode_map.clear();

    vector<bool> original_flip(cells.size(),false);

    encode_decode_LUT.flip_decode_map.insert(make_pair(0,original_flip));
    encode_decode_LUT.flip_encode_map.insert(make_pair(original_flip,0));
}

vector<vector<int>> BuildCombinationDistribution(int N, int K, int base)
{
    //Does not work when K = 0
    vector<vector<int>> combinations;
    vector<int> combination;
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
    
    // print integers and permute bitmask
    do {
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]){
                combination.push_back(i);
            }
        }
        combinations.push_back(combination);
        combination.clear();
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    vector<vector<int>> combination_vectors = ConvertToVector(N, combinations, base);
    return combination_vectors;
}

vector<vector<int>> ConvertToVector(int N, vector<vector<int>> combinations, int base){
    vector<vector<int>> combination_vectors;
    for(int i = 0; i < combinations.size(); i++){
        vector<int> new_vector (N,base);
        vector<int> combination = combinations[i];
        for(int i_elem = 0 ; i_elem < combination.size(); i_elem ++){
            new_vector[combination[i_elem]]++;
        }
        combination_vectors.push_back(new_vector);
    }
    return combination_vectors;
}

void Graph::InitializeGraph(const int& delta, const bool& is_cts_design, const RowInDRVRegion& drv_region_row){
    if(!IS_CHOI){
        int total_whitespace = drv_region_row.total_whitespace;
        int slots_whitespace = drv_region_row.cells.size() + 1;
        for(int i = 0; i < layers.size(); i++){
            Layer* whitespace_info = layers[i];
            int& whitespace = whitespace_info->cumulative_whitespace;
            int& num_slots = whitespace_info->num_slots;
            int& cumulative_slots = whitespace_info->cumulative_slots;
            if(i < layers.size()-1){
                for(int i_delta = -1 * delta; i_delta <= delta; i_delta++){
                    int node_value = whitespace + i_delta;
                    if(node_value < 0 || node_value > total_whitespace) continue;

                    Node *new_node = new Node(node_value, whitespace_info);
                    whitespace_info->nodes.push_back(new_node);
                }
            }
            else{
                int node_value = total_whitespace;
                Node *new_node = new Node(node_value, whitespace_info);
                whitespace_info->nodes.push_back(new_node);
            }
        }
    }
    else{
        //Choi
        int total_whitespace = drv_region_row.total_whitespace;
        int slots_whitespace = drv_region_row.cells.size() + 1;
        int base_whitespace = total_whitespace / slots_whitespace;
        for(int i = 0; i < layers.size(); i++){
            Layer* whitespace_info = layers[i];
            int& whitespace = whitespace_info->cumulative_whitespace;
            int& num_slots = whitespace_info->num_slots;
            int& cumulative_slots = whitespace_info->cumulative_slots;
            if(i < layers.size()-1){
                for(int i_delta = -1 * delta; i_delta <= delta; i_delta++){
                    int node_value = whitespace + i_delta;
                    if(node_value < 0 
                        || node_value > cumulative_slots * (base_whitespace + 1)
                        || ((base_whitespace > 0) && (node_value < base_whitespace * cumulative_slots))
                        )
                        continue;

                    Node *new_node = new Node(node_value, whitespace_info);
                    whitespace_info->nodes.push_back(new_node);
                }
            }
            else{
                int node_value = total_whitespace;
                Node *new_node = new Node(node_value, whitespace_info);
                whitespace_info->nodes.push_back(new_node);
            }
        }
    
    }
}

void Graph::BuildGraphAndRefine(const RowInDRVRegion& drv_region_row){
    int total_whitespace = drv_region_row.total_whitespace;
    int slots_whitespace = drv_region_row.cells.size() + 1;
    int base_whitespace = total_whitespace / slots_whitespace;

    //Sort layer by the value of nodes
    for(int i_layer = 0; i_layer < layers.size(); i_layer++){
        vector<Node*> node_layer = layers[i_layer]->nodes;
        sort(node_layer.begin(), node_layer.end(),
                [](Node* lhs, Node* rhs) ->bool{
                    return lhs->value < rhs->value;
                });
    }

    for(int i = layers.size()-2; i >= 0; i--){
        Layer* next_layer = layers[i + 1];
        Layer* cur_layer = layers[i];
        int next_cumulative_whitespace = next_layer->cumulative_whitespace;
        int next_slot = next_layer->num_slots;

        for(int i_next_node = 0; i_next_node < next_layer->nodes.size(); i_next_node++){
            Node* next_node = next_layer->nodes[i_next_node];
            int &next_whitespace = next_node->value;
            for(int i_node = 0; i_node < cur_layer->nodes.size(); i_node++){
                Node *cur_node = cur_layer->nodes[i_node];
                int &cur_whitespace = cur_node->value;
                if(IS_CHOI){
                    if (cur_whitespace + next_slot * base_whitespace <= next_whitespace
                        && next_whitespace <= cur_whitespace + next_slot * (1 + base_whitespace)){
                        next_node->parent_nodes.push_back(cur_node);
                        cur_node->child_nodes.push_back(next_node);
                    }
                }
                else{
                    if (cur_whitespace <= next_whitespace){
                        next_node->parent_nodes.push_back(cur_node);
                        cur_node->child_nodes.push_back(next_node);
                    }
                }
            }
        }

        //Search for the nodes to delete without connection.
        vector<int> node_idx_to_be_deleted;
        for(int i_node = 0; i_node < cur_layer->nodes.size(); i_node++){
            if(cur_layer->nodes[i_node]->child_nodes.size() == 0){
                node_idx_to_be_deleted.push_back(i_node);
            }
        }
        sort(node_idx_to_be_deleted.begin(), node_idx_to_be_deleted.end(), greater<int>());
        for(int i_node_to_be_deleted = 0 ; i_node_to_be_deleted < node_idx_to_be_deleted.size(); i_node_to_be_deleted++){
            delete cur_layer->nodes[node_idx_to_be_deleted[i_node_to_be_deleted]];
            cur_layer->nodes.erase(cur_layer->nodes.begin()+node_idx_to_be_deleted[i_node_to_be_deleted]);
        }
    }

    //Further delete nodes without having any connection with parent
    bool is_all_nodes_having_parents = false;
    while(!is_all_nodes_having_parents){
        for(int i = layers.size() - 1; i >= 1; i--){
            Layer* cur_layer = layers[i];
            vector<int> node_idx_to_be_deleted;
            for(int i_node = 0; i_node < cur_layer->nodes.size(); i_node++){
                Node* cur_node = cur_layer->nodes[i_node];
                vector<Node*> parent_nodes = cur_node->parent_nodes;

                if(parent_nodes.size()==0){
                    node_idx_to_be_deleted.push_back(i_node);
                    //Erase from parent list of child nodes
                    vector<Node*>& child_nodes = cur_node->child_nodes;
                    for(int i_child = 0; i_child < child_nodes.size(); i_child++){
                        Node* child_node = child_nodes[i_child];
                        vector<Node*> &parent_of_child_nodes = child_node->parent_nodes;
                        assert(parent_of_child_nodes.size() != 0);
                        
                        int idx_to_erase = -1;
                        for(int i_parent = 0; i_parent < parent_of_child_nodes.size(); i_parent++){
                            if(parent_of_child_nodes[i_parent] == cur_node){
                                idx_to_erase = i_parent;
                                break;
                            }
                        }
                        assert(idx_to_erase != -1);
                        parent_of_child_nodes.erase(parent_of_child_nodes.begin()+idx_to_erase);
                        
                    }

                }
            }
            sort(node_idx_to_be_deleted.begin(), node_idx_to_be_deleted.end(), greater<int>());
            for(int i_node_to_be_deleted = 0 ; i_node_to_be_deleted < node_idx_to_be_deleted.size(); i_node_to_be_deleted++){
                delete cur_layer->nodes[node_idx_to_be_deleted[i_node_to_be_deleted]];
                cur_layer->nodes.erase(cur_layer->nodes.begin()+node_idx_to_be_deleted[i_node_to_be_deleted]);
            }
        }

        bool check_parents = true;
        for(int i = layers.size() - 1; i >= 1; i--){
            Layer* cur_layer = layers[i];
            for(int i_node = 0; i_node < cur_layer->nodes.size(); i_node++){
                Node* cur_node = cur_layer->nodes[i_node];
                vector<Node*> parent_nodes = cur_node->parent_nodes;
                if(parent_nodes.size() == 0){
                    check_parents = false;
                    break;
                }
            }
            if(!check_parents) break;
        }
        is_all_nodes_having_parents = check_parents;
    }
}

void Graph::GetAllPaths(){
    //The input graph should not have verticies located on the mid-layers(not the start or the end) without parents. 
    if(layers[layers.size()-1]->nodes.size()==0){
        return;
    }
    Node* last_node = layers[layers.size()-1]->nodes[0];
    vector<Node*> node_stream;
    vector<Node*> stack = {last_node};

    int i_d = 0;
    while (stack.size() != 0){
        Node *current_node = stack[0];
        pop_front(stack);

        while(node_stream.size() != 0 
            && node_stream.size() > layers.size() - 1 - current_node->contained_layer->idx){
            pop_front(node_stream);
        }
        node_stream.insert(node_stream.begin(), current_node);
        
        vector<Node*> parent_nodes = current_node->parent_nodes;
        if(parent_nodes.size() == 0){
            if(node_stream.size() == layers.size()){
                NodeStream new_node_stream(node_stream);
                node_streams.insert(node_streams.begin(),new_node_stream);
                i_d ++;
            }
            pop_front(node_stream);
        }
        else{
            for(int i = 0; i < parent_nodes.size(); i++){
                stack.insert(stack.begin(),parent_nodes[i]);
            }
        }
    }
}

void Graph::CheckNumberOfPaths(){
    //Checking whether the full combinations have been correctly considered.
    num_cases = 0; // For the situation where the only one case exists.
    if(!(layers[layers.size()-1]->nodes.size() == 0)){
        layers[layers.size()-1]->nodes[0]->num_cases = 1;
        for(int i = layers.size() - 2; i >= 0; i--){
            vector<Node*> cur_nodes = layers[i]->nodes;
            for(int i_node = 0; i_node < cur_nodes.size(); i_node++){
                Node* cur_node = cur_nodes[i_node];
                vector<Node*> child_nodes = cur_node -> child_nodes;
                int cumulative_num_cases = 0;
                for(int i_child = 0; i_child < child_nodes.size(); i_child++){
                    if(child_nodes[i_child]->num_cases == -1){
                        cumulative_num_cases = -1;
                        break;
                    }
                    else{
                        cumulative_num_cases += child_nodes[i_child]->num_cases;
                    }
                }
                cumulative_num_cases = (cumulative_num_cases < 0)?-1:cumulative_num_cases; // if overflow, num_cases = -1
                cur_node -> num_cases = cumulative_num_cases;

                if(i == 0){
                    if(num_cases != -1){
                        if(cur_node->num_cases == -1){ 
                            num_cases = -1;
                        }
                        else{
                            num_cases += cur_node -> num_cases;
                        }
                    }
                }
            }
        }
    }
}



void Graph::CountNumberOfCombinations(const RowInDRVRegion& drv_region_row){
    int base_whitespace = drv_region_row.total_whitespace / (drv_region_row.cells.size() + drv_region_row.cts_cells.size() + 1);
    num_all_combinations = 0;
    for(int i = 0; i < node_streams.size(); i++){
        //vector<Node*>& node_stream = node_streams[i];
        NodeStream& node_stream = node_streams[i];
        int num_combinations = 1;
        int previous_whitespace = 0;
        for(int i_node = 0; i_node < node_stream.stream.size(); i_node++){
            Node* node = node_stream.stream[i_node];
            int n = node->contained_layer->num_slots;
            int cumulative_whitespace = node->value;
            int whitespace = cumulative_whitespace - previous_whitespace;
            int k = (n - whitespace) >=0 ? whitespace : whitespace - base_whitespace * n;

            int num_combination = nChoosek(n,k);
    
            node_stream.num_combination_for_each_layer.push_back(num_combination);
            num_combinations *= num_combination;

            previous_whitespace = cumulative_whitespace;
        }
        node_stream.num_combinations = num_combinations;
        num_all_combinations += num_combinations;

        assert(node_stream.stream.size() == node_stream.num_combination_for_each_layer.size());
    }
}



vector<int> BuildACombinationOfIdx(int idx, int num_slots, int whitespace, int base_whitespace){
    vector<int> combination(num_slots, base_whitespace);
    
    int tmp_idx = idx;
    int tmp_slots = num_slots;
    int tmp_whitespace = whitespace;
    for(int i = num_slots - 1; i >= 0; i--){
        if(tmp_slots == tmp_whitespace){
            combination[i]++;
            tmp_slots -= 1;
            tmp_whitespace -= 1;
        }
        else{
            int num_candidates_prefix_0 = nChoosek(tmp_slots - 1, tmp_whitespace);
            if(tmp_idx < num_candidates_prefix_0){
                tmp_slots -= 1;
            }
            else{
                tmp_idx -= num_candidates_prefix_0;
                combination[i] ++;
                tmp_slots -= 1;
                tmp_whitespace -= 1;
            }
        }
    }
    return combination;
}


vector<int> UniformNumberGenerator(int N, int num_candidates){
    vector<int> numbers_with_uniform_interval;
    if(N >= num_candidates || N < 0){
        int interval = 1;
        for(int i = 0; i < num_candidates; i++){
            numbers_with_uniform_interval.push_back(interval * i);
        }
    }
    else{
        for(int i = 0; i < N; i++){
            numbers_with_uniform_interval.push_back(i);
        }
    }
    return numbers_with_uniform_interval;
}

int SampleFromDifferentProbability(const vector<float>& probability){
    int random_number = rand()%100;

    float ref = 0;
    int return_idx = probability.size() - 1;
    for(int i = 0; i < probability.size(); i++){
        ref += probability[i] * 100;
        if (random_number < ref){
            return_idx = i;
            break;
        }
    }
    return return_idx;
}

vector<float> GenerateProabilityDistribution(vector<Node*> nodes){
    assert(nodes.size() != 0);

    vector<float> probability;
    int total_cases = 0;
    for(int i_c = 0; i_c < nodes.size(); i_c++){
        if(nodes[i_c]->num_cases != -1 && total_cases >= 0){
            total_cases += nodes[i_c]->num_cases;
        }
        else{
            total_cases = -1;
        }
    }
    
    if(total_cases > 0){
        for(int i_c = 0; i_c < nodes.size(); i_c++){
            probability.push_back(float(nodes[i_c]->num_cases) / total_cases);
        }
    }
    else{
        int n = nodes.size();
        total_cases = (2*n*n*n + 6*n*n + 4*n)/12; // Sigma Sigma n = Sigma n(n+1)/2  
        for(int i_c = 0; i_c < nodes.size(); i_c++){
            float value = (n-i_c)*((n-i_c)+1)/2; // Sigma n 
            probability.push_back(value / total_cases);
        }
    }
    
    return probability;
}


vector<int> Graph::GetPathWithIdx(const int& idx, const int& idx_idx, const int& num_candidates, const map<vector<int>,int>& shift_encode_map, const RowInDRVRegion& drv_region_row){
    assert((idx < num_cases || num_cases < 0) && "The input must not exceed the maximum input range.");
    bool is_random_sampling = (num_cases < num_candidates && num_cases != -1)?false:true; 
    bool is_lottery_sampling = (num_cases * 0.8 < num_candidates && num_candidates < num_cases && num_cases != -1)?true:false; // The constant 0.8 is added to prevent from sampling a few number of samples from all samples.

    vector<int> whitespace_distribution;
    int tmp_idx = idx;
    int before_whitespace = 0;
    Node* selected_node = NULL;

    if(!is_random_sampling){
        for(int i = 0; i < layers.size(); i++){
            vector<Node*> nodes;
            if(i == 0){
                nodes = layers[i]->nodes;
            }
            else{
                assert(selected_node != NULL);
                nodes = selected_node->child_nodes;
                if(i != layers.size()-1){
                    assert(nodes.size() > 0);
                }
            }
            for(int i_n = 0; i_n < nodes.size(); i_n++){
                if(tmp_idx - nodes[i_n]->num_cases < 0){
                    assert(nodes[i_n]->value - before_whitespace >= 0);
                    whitespace_distribution.push_back(nodes[i_n]->value - before_whitespace);
                    before_whitespace = nodes[i_n]->value;
                    selected_node = nodes[i_n];
                    break;
                }
                else tmp_idx -= nodes[i_n]->num_cases;
                assert(tmp_idx < idx || idx == 0);
            }
        }
    }
    else{
        //Iteratively generate random node combinations that is unique.
        while(true){
            whitespace_distribution.clear();
            before_whitespace = 0;
            for(int i = 0; i < layers.size(); i++){
                vector<Node*> nodes;
                if(i == 0){
                    nodes = layers[i]->nodes;
                }
                else{
                    assert(selected_node != NULL);
                    nodes = selected_node->child_nodes;
                    if(i != layers.size()-1){
                        assert(nodes.size() > 0);
                    }
                }
                vector<float> probability = GenerateProabilityDistribution(nodes);
                int i_n = SampleFromDifferentProbability(probability);
                assert(nodes[i_n]->value - before_whitespace >= 0);
                whitespace_distribution.push_back(nodes[i_n]->value - before_whitespace);
                before_whitespace = nodes[i_n]->value;
                selected_node = nodes[i_n];
            }
            if(shift_encode_map.find(whitespace_distribution) == shift_encode_map.end() || is_lottery_sampling){
                break;
            }
        }
    }
    
    assert(whitespace_distribution.size() == layers.size());
    return whitespace_distribution;
}

bool RowInDRVRegion::CheckOverlapWithCTSCell(const vector<int>& whitespace_vector, const ParserWrapper& lef_def_parser){
    if(cts_cells.size() == 0){
        return false;
    }
    bool is_overlap = false;
    int cumulative_whitespace = 0;
    for(int i_c = 0; i_c < cells.size(); i_c++){
        cumulative_whitespace += whitespace_vector[i_c];
        int new_cell_left = cells[i_c]->x + cumulative_whitespace * lef_def_parser.def.cpp;
        int new_cell_right = new_cell_left + cells[i_c]->cell->cell_width;
        for(int i_cts = 0; i_cts < cts_cells.size(); i_cts++){
            int cts_cell_left = cts_cells[i_cts]->x;
            int cts_cell_right = cts_cells[i_cts]->x + cts_cells[i_cts]->cell->cell_width;

            if(!
                ((cts_cell_left >= new_cell_right && cts_cell_left > new_cell_left) 
                || (cts_cell_right <= new_cell_left && cts_cell_right < new_cell_right))
            ){
                is_overlap = true;
                break;
            }
        }
        if(is_overlap) break;
    }
    return is_overlap;
}

void RowInDRVRegion::BuildShiftLUT(const ParserWrapper& lef_def_parser, const int& delta, const int& num_candidates, const bool& is_cts_design){
    
    cases_shift = 1;
    if (total_whitespace == 0){
        return;
    }

    shift_graph.InitializeGraph(delta, is_cts_design, *this);
    shift_graph.BuildGraphAndRefine(*this);
    shift_graph.CheckNumberOfPaths();
    
    if(shift_graph.num_cases == 0) return;

    int encoded_value = 0;
    vector<int> random_number_idx = UniformNumberGenerator(shift_graph.num_cases, num_candidates - 1); // -1 is due to remove the original whitespace distribution from counting.

    map<vector<int>, int> tmp_shift_encode_map;
    vector<vector<int>> tmp_shift_encode_vector;
    tmp_shift_encode_vector.push_back(cells_original_whitespace);
    int tmp_encoded_value = 1;

    for(int i = 0; i < random_number_idx.size(); i++){
        const int& idx = random_number_idx[i];
        vector<int> new_whitespace_vector = shift_graph.GetPathWithIdx(idx, i, num_candidates, tmp_shift_encode_map, *this);
      
        assert(new_whitespace_vector.size() == cells.size() + 1);
        if(new_whitespace_vector != cells_original_whitespace){
            tmp_shift_encode_map.insert(make_pair(new_whitespace_vector,tmp_encoded_value));
            tmp_shift_encode_vector.push_back(new_whitespace_vector);
        }
    }

    std::sort(tmp_shift_encode_vector.begin(), tmp_shift_encode_vector.end(), 
        [](const vector<int> &lhs, const vector<int> &rhs){
            bool result = false;
            for(size_t i = 0; i < lhs.size(); i++){
                if (lhs[i] < rhs[i]){
                    result = true;
                    break;
                }
                else if (lhs[i] > rhs[i]){
                    result = false;
                    break;
                }
            }
            return result;
        }
    );
    
    encode_decode_LUT.shift_encode_map.clear();
    encode_decode_LUT.shift_decode_map.clear();
    for(int i = 0 ; i < tmp_shift_encode_vector.size(); i++){
        encode_decode_LUT.shift_encode_map.insert(make_pair(tmp_shift_encode_vector[i], encoded_value));
        encode_decode_LUT.shift_decode_map.insert(make_pair(encoded_value, tmp_shift_encode_vector[i]));
        encoded_value++;
    }
    
    cases_shift = encoded_value;
}

vector<bool> ConvertDecimalToBinary(int decimal, int size){
    vector<bool> binary(size,0);
    for(int i = 0; decimal > 0; i++){
        binary[i] = decimal % 2;
        decimal /= 2;
    }
    return binary;
}

void RowInDRVRegion::BuildFlipLUT(const ParserWrapper& lef_def_parser, const int& delta, const int& num_candidates, const bool& is_cts_design){
    
    cases_flip = 1;
    int encoded_value = 1;
    
    size_t num_cells = cells.size();
    int num_flip_cases = pow(2,num_cells);

    vector<bool> zero_vector(num_cells,0);
    vector<int> random_number_idx = UniformNumberGenerator(num_flip_cases, num_candidates - 1); // -1 is due to remove the original whitespace distribution from counting.
    
    map<vector<bool>, int> tmp_flip_encode_map;
    vector<vector<bool>> tmp_flip_encode_vector;
    int tmp_encoded_value = 1;
    for(int i = 0; i < random_number_idx.size(); i++){
        const int& idx = random_number_idx[i];
        vector<bool> new_flip_vector;
        if(num_flip_cases <= num_candidates && num_flip_cases > 0){
            new_flip_vector = ConvertDecimalToBinary(idx, cells.size());
        }
        else{
            //Iteratively generate random node combinations that is unique.
            while(true){
                new_flip_vector.clear();
                for(int i = 0; i < cells.size(); i++){
                    bool i_n = rand()%2;
                    new_flip_vector.push_back(i_n);
                }
                if(tmp_flip_encode_map.find(new_flip_vector) == tmp_flip_encode_map.end()){
                    break;
                }
            }
        }
        if(new_flip_vector != zero_vector){
            tmp_flip_encode_map.insert(make_pair(new_flip_vector,tmp_encoded_value));
            tmp_flip_encode_vector.push_back(new_flip_vector);
        }
    }
    std::sort(tmp_flip_encode_vector.begin(), tmp_flip_encode_vector.end(), 
        [](const vector<bool> &lhs, const vector<bool> &rhs){
            bool result = false;
            for(size_t i = 0; i < lhs.size(); i++){
                if (lhs[i] < rhs[i]){
                    result = true;
                    break;
                }
                else if (lhs[i] > rhs[i]){
                    result = false;
                    break;
                }
            }
            return result;
        }
    );
    for(int i = 0; i < tmp_flip_encode_vector.size(); i++){
        encode_decode_LUT.flip_encode_map.insert(make_pair(tmp_flip_encode_vector[i], encoded_value));
        encode_decode_LUT.flip_decode_map.insert(make_pair(encoded_value, tmp_flip_encode_vector[i]));
        encoded_value ++;
    }
    cases_flip = encoded_value;
}

void RowInDRVRegion::BuildLUT(const ParserWrapper& lef_def_parser, const int& delta, const int& num_candidates, const bool& is_cts_design){
    BuildShiftLUT(lef_def_parser, delta, num_candidates, is_cts_design);
    BuildFlipLUT(lef_def_parser, delta, num_candidates, is_cts_design);

}
void RowInDRVRegion::Initialize(const ParserWrapper& lef_def_parser, const int& num_candidates, const bool& is_cts_design, const int& delta){
    SaveOriginalWhitespace(lef_def_parser);
    SaveOriginalFlip();
    BuildLUT(lef_def_parser, delta, num_candidates, is_cts_design);
}


void InputForBO::Clear(){
    input.clear();
    output_value = 0;
}


vector<RowInDRVRegion*> ClipDrvRegions(RowInDRVRegion* drv_region, const int& gcell_scale ,const ParserWrapper& lef_def_parser, const bool& is_cts_design, const bool& is_performing_flip){
    vector<RowInDRVRegion*> clipped_drv_regions;

    if (is_cts_design){
        int previous_right_coordinate = drv_region->x_left;
        int previous_right_idx = -1;
        int previous_right_grid = -1;    
        for(size_t i = 0; i < drv_region->cells.size(); i++){
            vector<placed_cell*> cells = drv_region->cells;

            if(cells[i]->is_cts_cell || i == drv_region->cells.size() - 1){
                int cumulative_whitespace = 0;
                int previous_right_cell = previous_right_coordinate;
                for(int i_c = previous_right_idx + 1; i_c <= i; i_c++){
                    cumulative_whitespace += cells[i_c]->x - previous_right_cell;
                    previous_right_cell = cells[i_c]->x + cells[i_c]->cell->cell_width;
                }
                if(i - previous_right_idx > 1
                    && ((cumulative_whitespace != 0 || is_performing_flip) || i == drv_region->cells.size() -1)                
                    ){
                    RowInDRVRegion* new_drv_region = new RowInDRVRegion();
                    new_drv_region->label = drv_region->label;
                    new_drv_region->y = drv_region->y;
                    new_drv_region->x_left = previous_right_coordinate;
                    if(cells[i]->is_cts_cell){
                        new_drv_region->x_right = cells[i]->x;
                    } 
                    else {
                        new_drv_region->x_right = drv_region->x_right;
                    }

                    new_drv_region->grid_y = drv_region->grid_y;
                    new_drv_region->grid_x_left = ConvertChipToGrid(new_drv_region->x_left, gcell_scale, true, lef_def_parser);  
                    if(new_drv_region->grid_x_left == previous_right_grid){
                        new_drv_region->grid_x_left++;
                    }
                    if(i == drv_region->cells.size()-1){
                        new_drv_region->grid_x_right = drv_region->grid_x_right;                  
                    }
                    else{
                        new_drv_region->grid_x_right = ConvertChipToGrid(new_drv_region->x_right,gcell_scale,true,lef_def_parser);
                    }
                    assert(new_drv_region->grid_x_left <= new_drv_region->grid_x_right);
                    
                    int right_i = (cells[i]->is_cts_cell)?i:i+1;
                    new_drv_region->total_cell_width = 0;
                    for(int i_c = previous_right_idx + 1; i_c < right_i; i_c++){
                        new_drv_region->cells.push_back(drv_region->cells[i_c]);
                        new_drv_region->total_cell_width += drv_region->cells[i_c]->cell->cell_width;
                    }    

                    previous_right_grid = new_drv_region->grid_x_right;
                    clipped_drv_regions.push_back(new_drv_region);
                }
                previous_right_coordinate = cells[i]->x + cells[i]->cell->cell_width;
                previous_right_idx = i;
            }
        }
    }
    if(clipped_drv_regions.size() == 0){
        clipped_drv_regions.push_back(drv_region);
    }

    return clipped_drv_regions;
}

void ClusteredDRVRegionsMultipleDimensions::Destroy(){
    for(size_t i = 0; i < inputs_for_BO.size(); i++){
        InputForBO & input_for_BO = inputs_for_BO[i];
        vector<ClusteredDRVRegions1Dimension>& tmp_vector = input_for_BO.input;
        for(size_t i2 = 0; i2< tmp_vector.size(); i2++){
            tmp_vector[i2].Destroy();
        }
    }
}

void RowInDRVRegion::CalculateMinDelta(){
    int whitespace_base = total_whitespace / cells_original_whitespace.size();

    //Find max num slots to set delta_offset
    int max_num_slots = -1;
    for(int i = 0; i < shift_graph.layers.size(); i++){
        int num_slots = shift_graph.layers[i]->num_slots;
        if(num_slots > max_num_slots){
            max_num_slots = num_slots;
        }
    }
    
    int delta = max_num_slots * whitespace_base;

    Layer* current_whitespace_info = shift_graph.layers[shift_graph.layers.size()-1];
    int current_whitespace = current_whitespace_info->cumulative_whitespace;
    for(int i = int(shift_graph.layers.size())-2; i >= 0; i--){
        Layer* whitespace_info = shift_graph.layers[i];
        int current_num_slots = current_whitespace_info->num_slots;
        int previous_whitespace = whitespace_info->cumulative_whitespace;

        if(current_whitespace - (previous_whitespace + delta) > current_num_slots){
            int new_delta = current_whitespace - previous_whitespace - current_num_slots;
            assert(new_delta > delta);
            delta = new_delta;

            current_whitespace = previous_whitespace + delta;
        }
        else{
            int min_available_whitespace = 100000;
            for(int i_delta = -1 * delta; i_delta <= delta; i_delta++){
                int new_previous_whitespace = previous_whitespace + i_delta;
                if(current_whitespace - new_previous_whitespace <= current_num_slots){
                    if(new_previous_whitespace < min_available_whitespace){
                        min_available_whitespace = new_previous_whitespace;
                    }
                }
            }
            current_whitespace = min_available_whitespace;
            assert(current_whitespace != 100000);
        }
        current_whitespace_info = whitespace_info;
    }
    delta_min = delta;
}

ClusteredDRVRegionsMultipleDimensions ClusterDRVRegions(    
        RowsOfDRVRegion*& drv_regions, 
        const int& max_dimension, 
        const ParserWrapper& lef_def_parser,
        const int& max_num_candidates,
        const int& max_num_cases,
        const bool& is_cts_design,
        const bool& is_performing_shift,
        const bool& is_performing_flip,
        const int& gcell_scale
        ){
    
    vector<RowInDRVRegion*>& drv_region_row = drv_regions->rows_of_drv_region;
    
    
    //Clipping DRV Regions
    //Clip DRV regions to have maxmimum max_num_cells cells in the drv region.
    vector<RowInDRVRegion*> clipped_drv_regions;
    for(size_t i = 0; i < drv_region_row.size(); i++){
        vector<RowInDRVRegion*> clipped_drv_regions_tmp = ClipDrvRegions(drv_region_row[i], gcell_scale, lef_def_parser, is_cts_design, is_performing_flip);
        for(size_t i_drv_region = 0; i_drv_region < clipped_drv_regions_tmp.size(); i_drv_region++){
            clipped_drv_regions_tmp[i_drv_region]->Initialize(lef_def_parser, max_num_candidates, is_cts_design);
            if(clipped_drv_regions_tmp[i_drv_region]->cases_shift == 1 && !is_performing_flip) continue;
            clipped_drv_regions.push_back(clipped_drv_regions_tmp[i_drv_region]);
        }
    }
    
    

    //Regroup DRV regions
    InputForBO input_for_BO;
    ClusteredDRVRegionsMultipleDimensions container_of_clustered_drv_regions;

    ClusteredDRVRegions1Dimension portion_of_clustered_drv_regions_shift;
    portion_of_clustered_drv_regions_shift.is_shift = true;
    portion_of_clustered_drv_regions_shift.AddDRVRegion(clipped_drv_regions[0]);

    ClusteredDRVRegions1Dimension portion_of_clustered_drv_regions_flip;
    portion_of_clustered_drv_regions_flip.is_shift = false;
    portion_of_clustered_drv_regions_flip.AddDRVRegion(clipped_drv_regions[0]);

    int remaining_dimension = max_dimension - int(is_performing_shift) - int(is_performing_flip);
    for(int i = 0; i < clipped_drv_regions.size()-1; i++){
        bool is_skip_shift = false;
        bool is_skip_flip = false;
        bool is_shift_new_dim = false;
        bool is_flip_new_dim = false;

        bool is_different_cluster = (clipped_drv_regions[i]->label == clipped_drv_regions[i+1]->label)?false:true;
        if(clipped_drv_regions[i+1]->cases_shift == 1) is_skip_shift = true;
        if(portion_of_clustered_drv_regions_shift.total_range * clipped_drv_regions[i+1]->cases_shift > max_num_cases && !is_skip_shift) is_shift_new_dim = true;
        if(portion_of_clustered_drv_regions_flip.total_range * clipped_drv_regions[i+1]->cases_flip > max_num_cases) is_flip_new_dim = true;

        if(!is_performing_shift){
            is_skip_shift = true;
            is_shift_new_dim = false;
        }
        if(!is_performing_flip){
            is_skip_flip = true;
            is_flip_new_dim = false;
        }

        int requiring_dimension = int(is_shift_new_dim) + int(is_flip_new_dim);

        if(remaining_dimension < requiring_dimension || is_different_cluster){
            if(is_performing_shift && portion_of_clustered_drv_regions_shift.num_cells != 0 && portion_of_clustered_drv_regions_shift.total_range != 1) input_for_BO.input.push_back(portion_of_clustered_drv_regions_shift);
            if(is_performing_flip && portion_of_clustered_drv_regions_flip.num_cells != 0) input_for_BO.input.push_back(portion_of_clustered_drv_regions_flip);
            container_of_clustered_drv_regions.inputs_for_BO.push_back(input_for_BO);
            input_for_BO.Clear();

            portion_of_clustered_drv_regions_shift.Clear();
            portion_of_clustered_drv_regions_flip.Clear();

            remaining_dimension = max_dimension;
        }

        remaining_dimension -= requiring_dimension;

        //Shift
        if(!is_skip_shift){
            if(is_shift_new_dim && portion_of_clustered_drv_regions_shift.num_cells != 0 && portion_of_clustered_drv_regions_shift.total_range != 1){ // Right after adding to container_of_clustered_drv_regions, skip this stage
                portion_of_clustered_drv_regions_shift.Preprocess();
                input_for_BO.input.push_back(portion_of_clustered_drv_regions_shift);
                portion_of_clustered_drv_regions_shift.Clear();
            }
            portion_of_clustered_drv_regions_shift.AddDRVRegion(clipped_drv_regions[i+1]);
        }
        //Flip
        if(!is_skip_flip){
            if(is_flip_new_dim && portion_of_clustered_drv_regions_flip.num_cells != 0){
                portion_of_clustered_drv_regions_flip.Preprocess();
                input_for_BO.input.push_back(portion_of_clustered_drv_regions_flip);
                portion_of_clustered_drv_regions_flip.Clear();
            }
            portion_of_clustered_drv_regions_flip.AddDRVRegion(clipped_drv_regions[i+1]);
        }

        if(i == clipped_drv_regions.size()-2){
            if(is_performing_shift && !is_skip_shift) input_for_BO.input.push_back(portion_of_clustered_drv_regions_shift);
            if(is_performing_flip) input_for_BO.input.push_back(portion_of_clustered_drv_regions_flip);
            container_of_clustered_drv_regions.inputs_for_BO.push_back(input_for_BO);

            input_for_BO.Clear();
            portion_of_clustered_drv_regions_shift.Clear();
            portion_of_clustered_drv_regions_flip.Clear();
        }
    }
    return container_of_clustered_drv_regions;
}


bool GetDirection(const int &num_flipped){
    return (num_flipped%2==0)?true:false; // false: MSB<-LSB , true: LSB->MSB
}

int EncodeShift(const RowInDRVRegion* row_in_drv_region, const ParserWrapper& lef_def_parser){
    tuple<vector<int>, int, int> whitespace_information = GetWhitespaceVector(row_in_drv_region, lef_def_parser);
    return row_in_drv_region->encode_decode_LUT.shift_encode_map.at(get<0>(whitespace_information));
}
int EncodeFlip(const RowInDRVRegion* row_in_drv_region){
    vector<bool> flip_information;
    for(int i = 0; i < row_in_drv_region->cells.size(); i++){
        placed_cell* cell = row_in_drv_region->cells[i];
        flip_information.push_back(cell->is_flipped);
    }
    return row_in_drv_region->encode_decode_LUT.flip_encode_map.at(flip_information);
}

void DecodeShift(int encoded_input, RowInDRVRegion* row_in_drv_region, const ParserWrapper& lef_def_parser){
    
    assert(row_in_drv_region->encode_decode_LUT.shift_decode_map.find(encoded_input) != row_in_drv_region->encode_decode_LUT.shift_decode_map.end()); // The encoded key must exist in LUT.
    vector<int> whitespace_distribution = row_in_drv_region->encode_decode_LUT.shift_decode_map.at(encoded_input);
   
    vector<placed_cell*> &cells = row_in_drv_region->cells;
    int previous_right = row_in_drv_region->x_left;
    for(int i = 0; i < cells.size(); i++){
        cells[i]->x = previous_right + whitespace_distribution[i] * lef_def_parser.def.cpp;
        previous_right = cells[i]->x + cells[i]->cell->cell_width;
    }
    assert(previous_right + whitespace_distribution[whitespace_distribution.size()-1]* lef_def_parser.def.cpp == row_in_drv_region->x_right);
}
void DecodeFlip(int encoded_input, RowInDRVRegion* row_in_drv_region){
    assert(row_in_drv_region->encode_decode_LUT.flip_decode_map.find(encoded_input) != row_in_drv_region->encode_decode_LUT.flip_decode_map.end()); // The encoded key must exist in LUT.
    vector<bool> flip_distribution = row_in_drv_region->encode_decode_LUT.flip_decode_map.at(encoded_input);

    vector<placed_cell*> &cells = row_in_drv_region->cells;
    for(int i = 0; i < cells.size(); i++){
        cells[i]->is_flipped = flip_distribution[i];
    }
}


int Encode(const ClusteredDRVRegions1Dimension& clustered_drv_regions, const ParserWrapper& lef_def_parser){
    vector<int> encoded_outputs;
    vector<int> total_num_cases;
    for(size_t i = 0; i < clustered_drv_regions.rows_of_drv_regions.size(); i++){
        RowInDRVRegion* current_row_of_drv_regions = clustered_drv_regions.rows_of_drv_regions[i];
        int encoded_output_row = -1;
        if(clustered_drv_regions.is_shift){
            encoded_output_row = EncodeShift(current_row_of_drv_regions, lef_def_parser);
            total_num_cases.push_back(current_row_of_drv_regions->cases_shift);
        }
        else{
            encoded_output_row = EncodeFlip(current_row_of_drv_regions);
            total_num_cases.push_back(current_row_of_drv_regions->cases_flip);
        }
        encoded_outputs.push_back(encoded_output_row);
    }

    int encoded_output = 0;
    for(int i1 = total_num_cases.size()-1; i1 >= 0 ; i1--){
        int weight = 1;
        for(int i2 = 0; i2 < i1; i2++){
            weight *= total_num_cases[i2];
        }
        weight *= encoded_outputs[i1];
        encoded_output += weight;
    }
    return encoded_output;
};

void Decode(int encoded_input, ClusteredDRVRegions1Dimension& clustered_drv_regions, const ParserWrapper& lef_def_parser){
    vector<int> total_num_cases;
    int max_num_decode = 1;
    for(size_t i = 0; i < clustered_drv_regions.rows_of_drv_regions.size(); i++){
        RowInDRVRegion* current_row_of_drv_regions = clustered_drv_regions.rows_of_drv_regions[i];
        if(clustered_drv_regions.is_shift){
            total_num_cases.push_back(current_row_of_drv_regions->cases_shift);
            max_num_decode *= current_row_of_drv_regions->cases_shift;
        }
        else{
            total_num_cases.push_back(current_row_of_drv_regions->cases_flip);
            max_num_decode *= current_row_of_drv_regions->cases_flip;
        }
    }

    if(encoded_input >= max_num_decode){
        cout<<encoded_input<<" / "<<max_num_decode<<endl;
        assert(encoded_input < max_num_decode && "The input to the Decode must not exceed the maximum input range.");
    }

    for(int i1 = total_num_cases.size()-1; i1 >= 0 ; i1--){
        RowInDRVRegion * current_row_of_drv_regions = clustered_drv_regions.rows_of_drv_regions[i1];
        int weight = 1;
        for(int i2 = 0; i2 < i1; i2++){
            weight *= total_num_cases[i2];
        }
        
        int quotient = encoded_input / weight;
        if(i1 != 0){
            if(clustered_drv_regions.is_shift){
                DecodeShift(quotient, current_row_of_drv_regions, lef_def_parser);
            }
            else{
                DecodeFlip(quotient, current_row_of_drv_regions);
            }
        }
        else{
            assert(encoded_input < total_num_cases[0]);
            if(clustered_drv_regions.is_shift){
                DecodeShift(encoded_input, current_row_of_drv_regions, lef_def_parser);
            }
            else{
                DecodeFlip(encoded_input, current_row_of_drv_regions);
            }
        }
        encoded_input -= quotient * weight;
    }
}

void RunReadLogit(const string& logit_directory, const string& read_logit_python_directory, const string& script_directory, const string& script_out_directory, const bool& is_updating_value_by_sum){
    std::string script;
    if(!is_updating_value_by_sum){
        script = "python3 " 
                + read_logit_python_directory + " "
                + logit_directory + " "
                + script_directory + " "
                + script_out_directory + " "
                + to_string(false);
    }
    else{
        script = "python3 " 
                + read_logit_python_directory + " "
                + logit_directory + " "
                + script_directory + " "
                + script_out_directory + " "
                + to_string(true);
    }
    system(script.c_str());
}

void UpdatePredictionValue(ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const string& logit_directory, const string& read_logit_python_directory, const string& script_directory, const string& script_out_directory, const bool& is_updating_value_by_sum){
    if(!is_updating_value_by_sum){
        std::ofstream ofs;
        ofs.open(script_directory);

        map<tuple<string,string,string>, bool> for_checking_grid_overlap;
        for(size_t i = 0; i < clustered_drv_regions.inputs_for_BO.size(); i++){
            InputForBO & input_for_BO = clustered_drv_regions.inputs_for_BO[i];
            vector<ClusteredDRVRegions1Dimension> & input_cluster = input_for_BO.input;

            for(size_t i_input = 0; i_input < input_cluster.size(); i_input++){
                vector<RowInDRVRegion*>& rows_in_drv_region = input_cluster[i_input].rows_of_drv_regions;
                for(size_t i_row = 0 ; i_row < rows_in_drv_region.size(); i_row++){
                    RowInDRVRegion* row_in_drv_region = rows_in_drv_region[i_row];
                    
                    string grid_y = to_string(row_in_drv_region->grid_y);
                    string grid_x_left = to_string(row_in_drv_region->grid_x_left);
                    string grid_x_right = to_string(row_in_drv_region->grid_x_right);

                    tuple<string,string,string> new_grid_tuple{grid_y, grid_x_left, grid_x_right};
                    if(for_checking_grid_overlap.find(new_grid_tuple) == for_checking_grid_overlap.end()){
                        for_checking_grid_overlap.insert(make_pair(new_grid_tuple,true));
                        string script = grid_y + " " + grid_x_left + " " + grid_x_right;
                        ofs<<script<<endl;
                    }
                }
            }
        }
        ofs.close();

        RunReadLogit(logit_directory, read_logit_python_directory, script_directory, script_out_directory, is_updating_value_by_sum);

        std::ifstream ifs;
        ifs.open(script_out_directory);
        assert(ifs.is_open() && "[ERROR] Prediction file does not exist!");
        
        for_checking_grid_overlap.clear();
        for(size_t i = 0; i < clustered_drv_regions.inputs_for_BO.size(); i++){
            vector<ClusteredDRVRegions1Dimension>& input_cluster = clustered_drv_regions.inputs_for_BO[i].input;

            double total_prediction_value_cluster = 0;
            for(size_t i_input = 0; i_input < input_cluster.size(); i_input++){
                vector<RowInDRVRegion*>& rows_in_drv_region = input_cluster[i_input].rows_of_drv_regions;
                double total_predicted_value = 0;
                for(size_t i_row = 0 ; i_row < rows_in_drv_region.size(); i_row++){
                    RowInDRVRegion* row_in_drv_region = rows_in_drv_region[i_row];
                    string grid_y = to_string(row_in_drv_region->grid_y);
                    string grid_x_left = to_string(row_in_drv_region->grid_x_left);
                    string grid_x_right = to_string(row_in_drv_region->grid_x_right);

                    tuple<string,string,string> new_grid_tuple{grid_y, grid_x_left, grid_x_right};
                    if(for_checking_grid_overlap.find(new_grid_tuple) == for_checking_grid_overlap.end()){
                        for_checking_grid_overlap.insert(make_pair(new_grid_tuple,true));
                    }
                    else{
                        continue;
                    }

                    char line[1024];

                    if(ifs.is_open()){
                        ifs.getline(line,1024);
                        char *token=strtok(line," ");
                
                        string tmp_grid_y;
                        string tmp_grid_x_left;
                        string tmp_grid_x_right;

                        string tmp_str(token);
                        
                        tmp_grid_y = tmp_str;
                        token=strtok(NULL," ");
                        tmp_str.assign(token);
                        tmp_grid_x_left = tmp_str;
                        token=strtok(NULL," ");
                        tmp_str.assign(token);
                        tmp_grid_x_right = tmp_str;

                        assert(grid_y == tmp_grid_y && grid_x_left == tmp_grid_x_left && grid_x_right == tmp_grid_x_right);
                        
                        token=strtok(NULL," ");
                        tmp_str.assign(token);
                        double predicted_value = stod(tmp_str);
                        total_predicted_value += predicted_value;
                        row_in_drv_region->predicted_value = predicted_value;
                        
                    }
                }
                input_cluster[i_input].total_predicted_value = total_predicted_value;
                total_prediction_value_cluster += total_predicted_value;
            }
            clustered_drv_regions.inputs_for_BO[i].output_value = total_prediction_value_cluster;
        }
        ifs.close();
    }
    else{
        RunReadLogit(logit_directory, read_logit_python_directory, script_directory, script_out_directory, is_updating_value_by_sum);
        
        double sum_total;
        std::ifstream ifs;
        ifs.open(script_out_directory);
        char line[1024];
        if(ifs.is_open()){
            ifs.getline(line,1024);
            char *token=strtok(line," ");
            string tmp_str(token);
            sum_total = stod(tmp_str);
        }
        ifs.close();

        for(size_t i = 0; i < clustered_drv_regions.inputs_for_BO.size(); i++){
            vector<ClusteredDRVRegions1Dimension> & input_cluster = clustered_drv_regions.inputs_for_BO[i].input;
            clustered_drv_regions.inputs_for_BO[i].output_value = sum_total;

        }
    }
}


double calculate_output(const vector<InputForBO>& input_vectors, const int& i_input, const bool& is_updating_value_by_sum){
    double output = 0;
    if(!is_updating_value_by_sum){
        const vector<ClusteredDRVRegions1Dimension>& input_vector = input_vectors[i_input].input;
        for(int i = 0; i < input_vector.size(); i++){
            output += input_vector[i].total_predicted_value;
        }        
    }
    else{
        output = input_vectors[i_input].output_value;
    }
    if(output<=0){
        cout<<output<<endl;
        assert(output <= 0);
    }
    return output;
}


vector<vector<int>> GetInputRange(const ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser){
    vector<vector<int>> max_input_range;
    const vector<InputForBO>& input_vectors = clustered_drv_regions.inputs_for_BO;

    for(int i_input = 0; i_input < input_vectors.size(); i_input++){
        const vector<ClusteredDRVRegions1Dimension>& input_vector = input_vectors[i_input].input;
        vector<int> input_range;

        for(int i = 0; i < input_vector.size(); i++){
            
            if(input_vector[i].total_range < 0){
                cout<<"[ERROR] The number of total range is negetive."<<endl;
                cout<<input_vector[i].num_cells<<endl;
                cout<<i_input<<" "<<i<<endl;
                cout<<input_vector[i].total_range<<endl;
            }
            input_range.push_back(input_vector[i].total_range);
        }
        max_input_range.push_back(input_range);
    }

    return max_input_range;
}

vector<vector<int>> GetEncodedInputs(const ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser){
    const vector<InputForBO>& input_vectors = clustered_drv_regions.inputs_for_BO;
    vector<vector<int>> encoded_inputs;
    
    for(int i_input = 0; i_input < input_vectors.size(); i_input++){
        const vector<ClusteredDRVRegions1Dimension>& input_vector = input_vectors[i_input].input;
        int input_dimension = input_vector.size();
        vector<int> encoded_input;
        
        for(int i = 0; i < input_vector.size(); i++){
            encoded_input.push_back(Encode(input_vector[i], lef_def_parser)); //Shift
        }
        encoded_inputs.push_back(encoded_input);
    }
    return encoded_inputs;
}

vector<double> GetPredictionValues(
        ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, 
        const ParserWrapper& lef_def_parser, 
        const bool& is_updating_value_by_sum,
        const string& drv_prediction_directory,
        const string& drv_prediction_model_directory,
        const string& lef_directory,
        const string& def_output_directory,
        const string& output_directory,
        const string& iterative_drv_script_directory,

        const string& new_logit_directory,
        const string& read_logit_directory,
        const string& iterative_coord_log_directory,
        const string& iterative_read_logit_script_directory,

        const bool& is_performing_classification_also
        ){    

    if(is_performing_classification_also){
        int inference_type;
        if (INFERENCE_TYPE == UGNET_UNSEEN){
            inference_type = UGNET_UNSEEN;
            IS_STOP_UGNET_UNSEEN = false;
        }
        else inference_type = UGNET;
        RunDrvPrediction(
            drv_prediction_directory, 
            lef_directory,
            def_output_directory, 
            output_directory, 
            iterative_drv_script_directory,
            inference_type);
    }
    RunDrvPrediction(
            drv_prediction_directory, 
            lef_directory,
            def_output_directory, 
            output_directory, 
            iterative_drv_script_directory,
            INFERENCE_TYPE
            );
    UpdatePredictionValue(clustered_drv_regions, new_logit_directory, read_logit_directory, iterative_coord_log_directory, iterative_read_logit_script_directory, is_updating_value_by_sum);

    vector<double> prediction_values;
    
    const vector<InputForBO>& input_vectors = clustered_drv_regions.inputs_for_BO;
    for(int i_input = 0; i_input < input_vectors.size(); i_input++){
        const vector<ClusteredDRVRegions1Dimension>& input_vector = input_vectors[i_input].input;
        double output = calculate_output(input_vectors, i_input, is_updating_value_by_sum);
        prediction_values.push_back(output);
    }
    return prediction_values;
}


void RefinePlacement(vector<vector<int>> inputs, ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser){
    vector<InputForBO>& input_vectors = clustered_drv_regions.inputs_for_BO;

    for(int i_input = 0; i_input < input_vectors.size(); i_input++){
        for(int i_sample = 0; i_sample < input_vectors[i_input].input.size(); i_sample++){
            ClusteredDRVRegions1Dimension& drv_region = input_vectors[i_input].input[i_sample];
            Decode(inputs[i_input][i_sample], drv_region, lef_def_parser);
        }
    }
}

void PrintWSDistribution(RowsOfDRVRegion* drv_regions, const ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser, const int& num_candidates){
    ofstream ofs;
    ofs.open("whitespace.log");
    const vector<InputForBO>& input_vectors = clustered_drv_regions.inputs_for_BO;
    cout<<"y    #Cells  #WS    #shift/flip   #total_shift/flip   is_shift  y   x_left   x_right   BO_idx  Dim  label"<<endl;
    ofs<<"y     #Cells  #WS    #shift/flip   #total_shift/flip   is_shift   y   x_left   x_right   BO_idx  Dim  label"<<endl;
    map<int,float> avg_ws;
    int label_ws_sum = 0;
    int label_num_cells = 0;
    int prev_label = input_vectors[0].input[0].rows_of_drv_regions[0]->label;
    for(int i_input = 0; i_input < input_vectors.size(); i_input++){
        const vector<ClusteredDRVRegions1Dimension>& input_vector = input_vectors[i_input].input;
        for(int i_rows = 0; i_rows < input_vector.size(); i_rows++){
            const std::vector<RowInDRVRegion*>& rows_of_drv_regions = input_vector[i_rows].rows_of_drv_regions;
            for(int i_row = 0; i_row < rows_of_drv_regions.size(); i_row++){
                const RowInDRVRegion* cur_row = rows_of_drv_regions[i_row];
                const vector<placed_cell*>& cells = cur_row->cells;
                vector<int> ws_vec;
                int ws_sum = 0;
                int prev_right = cur_row->x_left;
                for(int i_c=0;i_c<cells.size();i_c++){
                    int ws= (cells[i_c]->x - prev_right)/lef_def_parser.def.cpp;
                    ws_vec.push_back(ws);
                    prev_right = cells[i_c]->x + cells[i_c]->cell->cell_width;
                    ws_sum += ws;
                }
                
                int ws = (cur_row->x_right - prev_right)/lef_def_parser.def.cpp;
                ws_vec.push_back((cur_row->x_right - prev_right)/lef_def_parser.def.cpp);
                ws_sum+=ws;

                string script = to_string(cur_row->y)
                            + " " + to_string(cur_row->cells.size())
                            + " "+ to_string(ws_sum);
                if(input_vector[i_rows].is_shift){ 
                    script  +=  " " + to_string(cur_row->cases_shift)
                                + " " + to_string(cur_row->shift_graph.num_cases + 1);
                }
                else{ 
                    script  +=  " " + to_string(cur_row->cases_flip)
                                + " " + to_string(int(pow(2, cur_row->cells.size())));
                }
                script +=   " "+ to_string(input_vector[i_rows].is_shift)
                            + " "+ to_string(cur_row->grid_y)
                            + " "+ to_string(cur_row->grid_x_left)
                            + " "+ to_string(cur_row->grid_x_right)
                            + " "+ to_string(i_input)
                            + " "+ to_string(i_rows)
                            + " "+ to_string(cur_row->label);
                cout<<script<<endl;
                ofs<<script;
                ofs<<endl;

                if(input_vector[i_rows].is_shift){
                    if(prev_label != cur_row->label){
                        avg_ws.insert(make_pair(prev_label,float(label_ws_sum)/float(label_num_cells)));
                        prev_label = cur_row->label;
                        label_ws_sum = 0;
                        label_num_cells = 0;
                    }
                    label_ws_sum += ws_sum;
                    label_num_cells += cur_row->cells.size();
                }
            }

        }
    }

    ofs.close();
    cout<<"Label Avg_ws"<<endl;
    for(auto &it: avg_ws){
        cout<<it.first<<" "<<it.second<<endl;
    }
}