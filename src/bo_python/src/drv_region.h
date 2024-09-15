#ifndef PLACE_DRV_REGION_H_
#define PLACE_DRV_REGION_H_

#include <vector>
#include <tuple>
#include "cell.h"
#include "parser.h"
#include "map"

enum {
    UNET,
    UGNET,
    CNN_HALF,
    CNN_FULL,
    UGNET_UNSEEN,
    GCNN,
    UGNET_REG
};

extern int GCELL_SCALE;
extern int INFERENCE_TYPE;
extern std::string INFERENCE_NAME;
extern int IS_UPDATING_VALUE_BY_SUM;
extern int DELTA;
extern bool IS_CHOI;
extern bool IS_CADENCE;
extern bool IS_STOP_UGNET_UNSEEN;
extern int REGION_IDX;


struct Node;
struct Layer;
struct Graph;
struct RowInDRVRegion;

struct Node{
  Node(int value_, Layer* contained_layer_) : value(value_), contained_layer(contained_layer_) {
    parent_nodes.clear();
    child_nodes.clear();
    num_cases = 0;
  }
 
  int value;
  int num_cases;
  Layer* contained_layer;
  vector<Node*> parent_nodes;
  vector<Node*> child_nodes;
};

struct Layer{
  Layer(int cumulative_whitespace_, int num_slots_, int cumulative_slots_, int idx_)
      : cumulative_whitespace(cumulative_whitespace_), num_slots(num_slots_), cumulative_slots(cumulative_slots_), idx(idx_){
      nodes.clear();  
    }
  int cumulative_whitespace;
  int num_slots;
  int cumulative_slots;
  int idx;
  vector<Node*> nodes;
};


struct NodeStream{
    NodeStream(){
      num_combination_for_each_layer.clear();
    };
    NodeStream(vector<Node*> stream_) : stream(stream_){
      num_combination_for_each_layer.clear();
    };
    vector<Node*> stream;
    vector<int> num_combination_for_each_layer;
    int num_combinations;
};


class Graph{
  public:
    Graph(){
      layers.clear();
      num_all_combinations = -1;
      num_cases = -1;
    }
    ~Graph(){
      for(int i_layer = 0; i_layer < layers.size(); i_layer++){
        Layer* layer = layers[i_layer];
        for(int i_node = 0; i_node < layer->nodes.size(); i_node++){
          delete layer->nodes[i_node];
        }
        delete layer;
      }
    }
    vector<Layer*> layers;
    vector<NodeStream> node_streams;
    int num_cases;
    int num_all_combinations;
    
  void InitializeGraph(const int& delta, const bool& is_cts_design, const RowInDRVRegion& drv_region_row);
  void BuildGraphAndRefine(const RowInDRVRegion& drv_region_row);
  void GetAllPaths();
  void CheckNumberOfPaths();
  void CountNumberOfCombinations(const RowInDRVRegion& drv_region_row);
  vector<int> GetPathWithIdx(const int& idx, const int& idx_idx, const int& num_candidates, const map<vector<int>,int>& shift_encode_map, const RowInDRVRegion& drv_region_row);
};

struct LUT{
  public:
    map<vector<int>, int> shift_encode_map;
    map<int, vector<int>> shift_decode_map;
    
    map<vector<bool>, int> flip_encode_map;
    map<int, vector<bool>> flip_decode_map;
};


class RowInDRVRegion{
 public:
  RowInDRVRegion(){
    label = -1;
    y = -1;
    x_left = -1;
    x_right = -1;
    grid_y = -1;
    grid_x_left = -1;
    grid_x_right = -1;
    predicted_value = -1;
    cases_shift = 1;
    cases_flip = -1;
    delta_min = -1;
    cells.clear();
    cells_original_whitespace.clear();
    encode_decode_LUT.shift_encode_map.clear();
    encode_decode_LUT.shift_decode_map.clear();
    encode_decode_LUT.flip_encode_map.clear();
    encode_decode_LUT.flip_decode_map.clear();
  };
  RowInDRVRegion(RowInDRVRegion& _row_in_drv_region);

  int label;
  int y;
  int x_left;
  int x_right;

  int total_whitespace; // unit of cpp
  int total_cell_width;

  int grid_y, grid_x_left, grid_x_right;
  double predicted_value;

  int delta_min;
  int cases_shift;
  int cases_flip;
  
  std::vector<placed_cell*> cells; // If the input from the placement stage, cts cells can be moved. If the input is from the CTS stage, the CTS cells are removed from the vector and they are included in the cts_cells vector.
  std::vector<placed_cell*> cts_cells; // For after CTS designs only. 
  Graph shift_graph;
  std::vector<int> cells_original_whitespace;
  LUT encode_decode_LUT;
  

  void SaveOriginalCellPosition();
  void SaveOriginalWhitespace(const ParserWrapper& lef_def_parser);
  void SaveOriginalFlip();
  void Initialize(const ParserWrapper& lef_def_parser, const int& num_candidates, const bool& is_cts_design, const int& delta = DELTA);
  void BuildLUT(const ParserWrapper& lef_def_parser, const int& delta, const int& num_candidates, const bool& is_cts_design);
  void BuildShiftLUT(const ParserWrapper& lef_def_parser, const int& delta, const int& num_candidates, const bool& is_cts_design);
  void BuildFlipLUT(const ParserWrapper& lef_def_parser, const int& delta, const int& num_candidates, const bool& is_cts_design);

  bool CheckOverlapWithCTSCell(const vector<int>& whitespace_vector, const ParserWrapper& lef_def_parser);

  void CalculateMinDelta();
  void SeparateCTSCells();
  void AssignInvalidDistanceToCells(const ParserWrapper& lef_def_parser, const int& delta);
 private:
  
};


//Wrapper for row_of_drv_region (to destroy)
class RowsOfDRVRegion{
 public:
  RowsOfDRVRegion();
  ~RowsOfDRVRegion();
  std::vector<RowInDRVRegion*> rows_of_drv_region;
  
  void Preprocess(const ParserWrapper& lef_def_parser, const int& num_candidates, const bool& is_cts_design);
  void MatchCellsWithDrvRegions(const ParserWrapper& lef_def_parser);
  void SeparateCTSCells();
  void AssignInvalidDistanceToCells(const ParserWrapper& lef_def_parser, const int& delta);
  void SaveOriginalWhitespace(const ParserWrapper& lef_def_parser);
};


class ClusteredDRVRegions1Dimension{
 public:
  ClusteredDRVRegions1Dimension(){
    Clear();
  };
  ~ClusteredDRVRegions1Dimension(){};
  
  bool is_shift; //True: shift , False: flip

  int total_range;
  int num_cells;
  std::vector<int> ranges;
  std::vector<RowInDRVRegion*> rows_of_drv_regions;

  double total_predicted_value;

  void Clear();
  void Destroy();
  void Preprocess();
  void AddDRVRegion(RowInDRVRegion *row_in_drv_region);
};

class InputForBO{
  public:
    InputForBO(){
      input.clear();
      output_value = -1;
    }
    std::vector<ClusteredDRVRegions1Dimension> input;
    double output_value;

    void Clear();
};

// This class contains "max_dimension" clustered_drv_regions.
class ClusteredDRVRegionsMultipleDimensions{
 public:
  std::vector<InputForBO> inputs_for_BO;

  void Destroy();
};



RowsOfDRVRegion* ParseDrvRegionLog(const std::string& script_directory, const ParserWrapper& lef_def_parser, const int& gcell_scale = GCELL_SCALE);
RowsOfDRVRegion* DecodeDrvRegion(const std::string& script_directory, const ParserWrapper& lef_def_parser, const int& delta, const bool& is_cts_design, const int& gcell_scale = GCELL_SCALE);

ClusteredDRVRegionsMultipleDimensions ClusterDRVRegions(RowsOfDRVRegion*& drv_regions, const int& max_dimension, const ParserWrapper& lef_def_parser, const int& max_num_candidates, const int& max_num_inputs, const bool& is_cts_design, const bool& is_performing_shift, const bool& is_performing_flip, const int& gcell_scale = GCELL_SCALE);

void UpdatePredictionValue(ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const string& logit_directory, const string& read_logit_python_directory, const string& script_directory, const string& script_out_directory, const bool& is_updating_value_by_sum);
double calculate_output(const vector<InputForBO>& input_vectors, const int& i_input, const bool& is_updating_value_by_sum);

vector<vector<int>> GetInputRange(const ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser);

vector<vector<int>> GetEncodedInputs(const ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser);
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
        );

      
void RefinePlacement(vector<vector<int>> inputs, ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser);
void PrintWSDistribution(RowsOfDRVRegion* drv_regions, const ClusteredDRVRegionsMultipleDimensions& clustered_drv_regions, const ParserWrapper& lef_def_parser, const int& num_candidates);

tuple<vector<int>, int, int> GetWhitespaceVector(const RowInDRVRegion* row_in_drv_region, const ParserWrapper& lef_def_parser);
int EncodeShift(const RowInDRVRegion* row_in_drv_region, const ParserWrapper& lef_def_parser);
void DecodeShift(int encoded_input, RowInDRVRegion* row_in_drv_region, const ParserWrapper& lef_def_parser);
int Encode(const ClusteredDRVRegions1Dimension& clustered_drv_regions, const ParserWrapper& lef_def_parser);
void Decode(int encoded_input, ClusteredDRVRegions1Dimension& clustered_drv_regions, const ParserWrapper& lef_def_parser);

int SampleFromDifferentProbability(const vector<float>& probability);
vector<float> GenerateProabilityDistribution(vector<Node*> nodes);

vector<vector<int>> ConvertToVector(int N, vector<vector<int>> combinations, int base);



#endif
