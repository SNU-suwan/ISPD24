#ifndef PLACE_UTILS_H_
#define PLACE_UTILS_H_

#include "parser.h"

void RunDrvPrediction(
        const std::string& drv_prediction_directory, 
        const std::string& lef_directory,
        const std::string& def_directory,
        const std::string& prediction_directory, 
        const std::string& script_directory,
        const int &inference_type
        );

void RunDrvRegion(
    const std::string& drv_region_directory,
    const std::string& logit_directory,
    const std::string& script_directory,
    const std::string& blockages_directory,
    double & threshold,
    const std::string& lef_directory,
    const std::string& def_directory,
    const std::string& track_pitch,
    const std::string& blockage_layers
    );

int ConvertGridToChip(const int& grid_value, const int& gcell_scale, const bool& is_xy, const ParserWrapper& lef_def_parser);
int ConvertChipToGrid(const int& chip_coordinate, const int& gcell_scale, const bool& is_xy, const ParserWrapper& lef_def_parser);
void RemoveScript(const std::string& output_directory, const std::string& design_name_with_postfix, const bool& is_updating_value_by_sum);
void RemoveDefCsvNpz(const std::string& def_directory, const std::string& output_directory, const std::string& design_name_with_postfix, const bool& is_updating_value_by_sum);
vector<string> Parse(string text);
void PrintVector(const vector<vector<double>> vector_in_vector);

template <typename T>
void pop_front(std::vector<T>& vec){
    assert(!vec.empty());
    vec.erase(vec.begin());
}

void PrintVector(const vector<vector<double>> vector_input);
void PrintVector(const vector<vector<int>> vector_input);
void PrintVector(const vector<double> vector_input);

void RemoveLog(const string& log_dir);

unsigned nChoosek( unsigned n, unsigned k );
bool StartsWith(const string& s, const string& reference);

#endif
