#include <iostream>
#include <fstream>
#include <cstring>
#include "parser.h"
#include "drv_region.h"
#include "utils.h"
using namespace std;

bool IS_STOP_UGNET_UNSEEN = false;
bool IS_INITIALIZATION = true;
std::vector<std::string> benchmark_lists{"aes","b18","b19","ecg","eth","jpeg","ldpc","nova","tate","vga"};

void RunDrvPrediction(
        const std::string& drv_prediction_directory, 
        const std::string& lef_directory,
        const std::string& def_directory,
        const std::string& prediction_directory, 
        const std::string& script_directory,
        const int &inference_type
        ){
    std::string drv_prediction_model_directory;
    std::string script;
    cout<<" INF TYPE: "<<inference_type<<endl;
    if (IS_INITIALIZATION && !(inference_type == UGNET_UNSEEN)){
        if (!IS_CADENCE){
           
           drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
            script = " -i " + drv_prediction_directory + "/v20_clas_info.pkl"
                + " -c " + drv_prediction_directory + "/v20_clas.pth ";
        
        }
        else{
            drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
            script = " -i " + drv_prediction_directory + "/PGNN_asap_clas.pkl"
                + " -c " + drv_prediction_directory + "/PGNN_asap_clas.pth"
                + " -is_cadence 1";
        }
        IS_INITIALIZATION = false;
        
    }
    else if (inference_type == UGNET_UNSEEN && !IS_STOP_UGNET_UNSEEN){
        drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
        string bench_name;
        for(size_t i = 0 ; i<benchmark_lists.size(); i++){
            if(def_directory.find(benchmark_lists[i]) != std::string::npos){
                bench_name = benchmark_lists[i];
                break;
            }
        }
        script = " -i " + drv_prediction_directory + "/unseen/" + bench_name + ".pkl"
            + " -c " + drv_prediction_directory + "/unseen/"+ bench_name + ".pth ";
        IS_STOP_UGNET_UNSEEN = true;
    }
    else if (inference_type == UNET || (inference_type == UGNET_UNSEEN && IS_STOP_UGNET_UNSEEN)){
        drv_prediction_model_directory = drv_prediction_directory + "/inference_unet.py";
        script = " -i " + drv_prediction_directory + "/unet_info_reg.pkl"
            + " -c " + drv_prediction_directory + "/unet_reg.pth ";
    }
    else if (inference_type == UGNET){
        if(!IS_CADENCE){
            drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
            script = " -i " + drv_prediction_directory + "/v20_clas_info.pkl"
                + " -c " + drv_prediction_directory + "/v20_clas.pth ";
        }
        else{
            drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
            script = " -i " + drv_prediction_directory + "/PGNN_asap_clas.pkl"
                + " -c " + drv_prediction_directory + "/PGNN_asap_clas.pth"
                + " -is_cadence 1";
        }
    }
    else if (inference_type == CNN_HALF){
        drv_prediction_model_directory = drv_prediction_directory + "/inference_cnn.py";
        script = " -i " + drv_prediction_directory + "/PinUGnet_info_clas.pkl"
            + " -c " + drv_prediction_directory + "/cnn_training_half.pth ";
    }
    else if (inference_type == CNN_FULL){
        drv_prediction_model_directory = drv_prediction_directory + "/inference_cnn.py";
        script = " -i " + drv_prediction_directory + "/PinUGnet_info_clas.pkl"
            + " -c " + drv_prediction_directory + "/cnn_all.pth ";
    }
    else if (inference_type == GCNN){
        drv_prediction_model_directory = drv_prediction_directory + "/inference_gcnn.py";
        script = " -i " + drv_prediction_directory + "/gcnn.pkl"
            + " -c " + drv_prediction_directory + "/gcnn.pth ";
    }
    else if (inference_type == UGNET_REG){
        if(!IS_CADENCE){
            drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
            script = " -i " + drv_prediction_directory + "/v20_reg_info.pkl"
                + " -c " + drv_prediction_directory + "/v20_reg.pth ";
        }
        else{
            drv_prediction_model_directory = drv_prediction_directory + "/inference_graph.py";
            script = " -i " + drv_prediction_directory + "/PGNN_asap_reg.pkl"
                + " -c " + drv_prediction_directory + "/PGNN_asap_reg.pth"
                + " -is_cadence 1";
        }
    }
    script = "python3 -u " + drv_prediction_model_directory + " -d " + def_directory + " -o " + prediction_directory + " -l " + lef_directory + script;
    
    std::ofstream ofs;
    ofs.open(script_directory);
    cout<<script<<endl;
    ofs<<script<<endl;
    ofs.close();

    std::string tmp_script = "source " + script_directory;
    system(tmp_script.c_str());
}

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
    ){

    std::string tmp_script = "python3 " 
        + drv_region_directory + " " 
        + logit_directory + " " 
        + to_string(threshold) + " " 
        + script_directory + " "
        + blockages_directory + " "
        + lef_directory + " "
        + def_directory + " "
        + track_pitch + " "
        + blockage_layers;
    
    system(tmp_script.c_str());

}

void RemoveScript(const std::string& output_directory, const std::string& design_name_with_postfix, const bool& is_updating_value_by_sum){
    std::string script_directory = output_directory + "/scripts";

    std::string script_1 = script_directory + "/run_drv_prediction_" + design_name_with_postfix;
    std::string script_2 = script_directory + "/coords_" + design_name_with_postfix;
    std::string script_3 = script_directory + "/" + design_name_with_postfix;
  
    script_1 += ".bashrc";
    script_2 += ".log";
    script_3 += ".log";

    std::string remove_script_1 = "rm " + script_1;
    std::string remove_script_2 = "rm " + script_2;
    std::string remove_script_3 = "rm " + script_3;
    cout<<remove_script_1<<endl;
    cout<<remove_script_2<<endl;
    cout<<remove_script_3<<endl;
    system(remove_script_1.c_str());
    if(!is_updating_value_by_sum){
        system(remove_script_2.c_str());
    }
    system(remove_script_3.c_str());
}

void RemoveDefCsvNpz(const std::string& def_directory, const std::string& output_directory, const std::string& design_name_with_postfix, const bool& is_updating_value_by_sum){
    std::string csv_directory = output_directory + "/csv/" + design_name_with_postfix;
    std::string npz_directory;
    std::string prediction_directory;
    
    if (INFERENCE_TYPE == UNET || INFERENCE_TYPE == UGNET_UNSEEN){
        npz_directory = output_directory + "/npz/" + design_name_with_postfix + ".npz";
        prediction_directory = output_directory + "/predictions/" + design_name_with_postfix + "_reg_logit.npy";
    }
    else if (INFERENCE_TYPE == UGNET){
        npz_directory = output_directory + "/npz/" + design_name_with_postfix + "_class.npz";
        prediction_directory = output_directory + "/predictions/" + design_name_with_postfix + "_class_logit.npy";
    }
    else if (INFERENCE_TYPE == CNN_HALF){
        npz_directory = output_directory + "/npz/" + design_name_with_postfix + "_cnn_half.npz";
        prediction_directory = output_directory + "/predictions/" + design_name_with_postfix + "_cnn_half_logit.npy";
    }
    else if (INFERENCE_TYPE == CNN_FULL){
        npz_directory = output_directory + "/npz/" + design_name_with_postfix + "_cnn_full.npz";
        prediction_directory = output_directory + "/predictions/" + design_name_with_postfix + "_cnn_full_logit.npy";
    }
    else if (INFERENCE_TYPE == GCNN){
        npz_directory = output_directory + "/npz/" + design_name_with_postfix + "_gcnn.npz";
        prediction_directory = output_directory + "/predictions/" + design_name_with_postfix + "_cnn_full_logit.npy";
    }
    else if (INFERENCE_TYPE == UGNET_REG){
        npz_directory = output_directory + "/npz/" + design_name_with_postfix + "_reg_graph.npz";
        prediction_directory = output_directory + "/predictions/" + design_name_with_postfix + "_reg_graph_logit.npy";
    }

    std::string prediction_directory2 = output_directory + "/predictions/" + design_name_with_postfix + "_class_pred.npy";

    std::string remove_def_script = "rm " + def_directory;
    std::string remove_csv_script = "rm -r " + csv_directory;
    std::string remove_npz_script = "rm " + npz_directory;
    std::string remove_prediction_script = "rm " + prediction_directory;
    std::string remove_prediction_script2 = "rm " + prediction_directory2;

    std::string script_directory = output_directory + "/scripts";
    
    cout<<"Remove .def, .csv, .npz, for " + design_name_with_postfix <<endl;
    system(remove_def_script.c_str());
    system(remove_csv_script.c_str());
    system(remove_npz_script.c_str());
    system(remove_prediction_script.c_str());
    system(remove_prediction_script2.c_str());
}

int ConvertGridToChip(const int& grid_value, const int& gcell_scale, const bool& is_xy, const ParserWrapper& lef_def_parser){
    //gcell_scale : the unit of gcell in terms of a row height
    //is_xy : true for x, false for y
    
    const int& row_height = lef_def_parser.def.row_height;
    const Point<int>& offset = lef_def_parser.def.offset;

    int chip_coordinate = (is_xy) ? 
        grid_value * (row_height * gcell_scale) + offset.x 
        : grid_value * (row_height * gcell_scale) + offset.y;
    
    return chip_coordinate;
}

int ConvertChipToGrid(const int& chip_coordinate, const int& gcell_scale, const bool& is_xy, const ParserWrapper& lef_def_parser){
    //gcell_scale : the unit of gcell in terms of a row height
    //is_xy : true for x, false for y
    
    const int& row_height = lef_def_parser.def.row_height;
    const Point<int>& offset = lef_def_parser.def.offset;

    int grid_value = (is_xy) ? 
        (chip_coordinate - offset.x ) / (row_height * gcell_scale) 
        : (chip_coordinate - offset.y) / (row_height * gcell_scale);
    
    return grid_value;
}

vector<string> Parse(string text){
    size_t n = text.length();
    char char_array[n+1];
    strcpy(char_array, text.c_str());    
    char *token=strtok(char_array," ");

    vector<string> words;
    while (token != NULL){
		string tmp_str(token);
        words.push_back(tmp_str);
        token=strtok(NULL," ");
    }
    return words;
}


void PrintVector(const vector<vector<double>> vector_input){
    for(int i1 = 0; i1 < vector_input.size(); i1++){
        cout<<"idx"<<i1<<" ";
        for(int i2 = 0; i2 < vector_input[i1].size(); i2++){
            cout<<vector_input[i1][i2]<<" ";
        }
    }
    cout<<endl;
}
void PrintVector(const vector<vector<int>> vector_input){
    for(int i1 = 0; i1 < vector_input.size(); i1++){
        cout<<"idx"<<i1<<" ";
        for(int i2 = 0; i2 < vector_input[i1].size(); i2++){
            cout<<vector_input[i1][i2]<<" ";
        }
    }
    cout<<endl;
}
void PrintVector(const vector<double> vector_input){
    for(int i1 = 0; i1 < vector_input.size(); i1++){
        cout<<"idx"<<i1<<" "<<vector_input[i1]<<" ";
    }
    cout<<endl;
}

unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

bool StartsWith(const string& s, const string& reference){
    if (s.rfind(reference,0) == 0) return true;
    return false;
}

void RemoveLog(const string& log_dir){
    string remove_log_script = "rm " + log_dir;
    system(remove_log_script.c_str());
}
