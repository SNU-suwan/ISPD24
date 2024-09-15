#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include "Place.h"

using namespace std;


int GCELL_SCALE;
int MAX_ITER;
int DELTA;
int INFERENCE_TYPE = 1; // 0: UNET, 1: UGNET, 2: CNN_HALF, 3: CNN_FULL, 4: UGNET + UNSEEN, 5: Graph + CNN, 6: UGNET (Regression)
std::string INFERENCE_NAME;
std::string TRACK_PITCH;
bool IS_CHOI = true;
bool IS_CADENCE = false;
int REGION_IDX = -1;

vector<string> CLOCK_CELLS = {
    //"BUF",
    //"CLK",
    "DFF",
    "SDFF",
    //"TBUF"
};
vector<string> CLOCK_CELLS_NAME = {
    "cts",
    "ccd",
    "cto"
};
string INVERTER_NAME = "INV";


int main(int argc, char *argv[]){
    
    ///// INITIALIZE VARIALBES /////
    std::chrono::steady_clock::time_point time_0 = std::chrono::steady_clock::now();
    GCELL_SCALE = 1;
    DELTA = 5;
    
    MAX_ITER = 1;
    int MAX_DIMENSION = 5; //->20
    int MAX_NUM_CANDIDATES = 100; // The maximum number for single input space domain.
    int MAX_NUM_CASES = 1000; // The number for concatenating adjacent rows if their multiplication of input space is smaller than "MAX_NUM_CASES".
    int GCELL_SIZE = 1;
    
    bool is_performing_shift = true;
    bool is_performing_flip = false;
    bool is_cts_from_place = false;
    bool is_updating_value_by_sum = false;

    std::string project_directory = ".";
    std::string drv_prediction_directory = project_directory + "/src/DRV_prediction_model";
    std::string drv_prediction_model_directory;
    std::string drv_region_directory = project_directory + "/src/DRV_region/get_drv_region.py";
    std::string read_logit_directory = project_directory + "/src/read_logit_result/read_logit.py";
    std::string script_directory = project_directory + "/outputs/scripts";
    std::string direct_input_def = "" ;
    std::string direct_input_lef = "" ;
    std::string blockage_layers = "" ;
    bool is_cts_design = false;

    is_cts_design = true; //True: CTS, False: placement
    double threshold = 0.4985;

    std::string design_name;
    
    // Parse arguments
    std::vector<std::string> arg;
    for(int i = 1; i < argc; i++){
        arg.push_back(std::string(argv[i]));
    }   
    for(int i = 0; i < arg.size(); i++){
        if (arg[i]=="-m"  || arg[i]=="-max_iter"){
            MAX_ITER = stoi(arg[i+1]);
        }
        else if (arg[i]=="-c" || arg[i]=="-choi"){
            if(arg[i+1] == "0"){
                IS_CHOI = false;
            }
            else if(arg[i+1] == "1"){
                IS_CHOI = true;
            }
        }
        else if (arg[i] == "-cts"){
            if (arg[i+1] == "1"){
                is_cts_design = true;
            }
            else if (arg[i+1] == "0"){
                is_cts_design = false;
            }
        }
        else if (arg[i]=="-d"  || arg[i]=="-dimension"){
            MAX_DIMENSION = stoi(arg[i+1]);
        } 
        else if (arg[i]=="-p"){
            //MAX_NUM_CELLS = stoi(arg[i+1]);
            MAX_NUM_CASES = stoi(arg[i+1]);
        }
        else if (arg[i] == "-name"){
            design_name = arg[i+1];
        }
        else if (arg[i] == "-delta"){
            DELTA = stoi(arg[i+1]);
        }
        else if (arg[i] == "-threshold"){
            threshold = stod(arg[i+1]);
        }
        else if (arg[i] == "-shift"){
            if (arg[i+1] == "1"){
                is_performing_shift = true;
            }
            else if (arg[i+1] == "0"){
                is_performing_shift = false;
            }
        }
        else if (arg[i] == "-flip"){
            if (arg[i+1] == "1"){
                is_performing_flip = true;
            }
            else if (arg[i+1] == "0"){
                is_performing_flip = false;
            }
        }
        else if(arg[i] == "-candidate" || arg[i] == "-cand"){
            MAX_NUM_CANDIDATES = stoi(arg[i+1]);
        }
        else if (arg[i] == "-r" || arg[i] == "-region"){
            REGION_IDX = stoi(arg[i+1]);
        }
        else if (arg[i] == "-inference_type" || arg[i] == "-inf" || arg[i] == "-INF"){
            INFERENCE_TYPE = stoi(arg[i+1]);
        }
        else if (arg[i] == "-cadence"){
            IS_CADENCE = true;
        }
        else if (arg[i] == "-def"){
            direct_input_def = arg[i+1];
        }
        else if (arg[i] == "-lef"){
            direct_input_lef = arg[i+1];
        }
        else if (arg[i] == "-track_pitch"){
            TRACK_PITCH = arg[i+1];
        }
        else if (arg[i] == "-blockage_layers"){
            blockage_layers = arg[i+1];
        }
        else if (arg[i] == "-cfp"){
            if (arg[i+1] == "0"){
                is_cts_from_place = false;
            }
            else if (arg[i+1] == "1"){
                is_cts_from_place = true;
                is_cts_design = true;
            }
            else{
                cout<<"Invalid input " << arg[i+1] <<endl;
                return -1;
            }
            
        }
    }
    
    if(!is_performing_shift && !is_performing_flip){
        cout<<"[ERROR] One of operations, shift and flip, must be turned on."<<endl;
        exit(-1);
    }

    if (INFERENCE_TYPE == UNET){
        INFERENCE_NAME = "reg";
    }
    else if (INFERENCE_TYPE == UGNET){
        INFERENCE_NAME = "class";
    }
    else if (INFERENCE_TYPE == CNN_HALF){
        INFERENCE_NAME = "cnn_half";
    }
    else if (INFERENCE_TYPE == CNN_FULL){
        INFERENCE_NAME = "cnn_full";
    }
    else if (INFERENCE_TYPE == UGNET_UNSEEN){
        INFERENCE_NAME = "unseen";
    }
    else if (INFERENCE_TYPE == GCNN){
        INFERENCE_NAME = "gcnn";
    }
    else if (INFERENCE_TYPE == UGNET_REG){
        INFERENCE_NAME = "reg_graph";
    }

    std::string postfix = "M" + to_string(MAX_ITER)  //If you want to change this, you must change the python code as well.
        ;
    std::string def_postfix = postfix;
    if (is_cts_from_place){ 
        def_postfix = postfix 
            + "_CTS0"
            + "_DELTA" + to_string(DELTA)
            + "_D" + to_string(MAX_DIMENSION)
            + "_P" + to_string(MAX_NUM_CASES)
            + "_CAND" + to_string(MAX_NUM_CANDIDATES)
            + "_TH" + to_string(int(threshold*10000));
        postfix = postfix 
            + "_CTS0" 
            + "_DELTA" + to_string(DELTA)
            + "_D" + to_string(MAX_DIMENSION)
            + "_P" + to_string(MAX_NUM_CASES)
            + "_CAND" + to_string(MAX_NUM_CANDIDATES)
            + "_TH" + to_string(int(threshold*10000));
    }
    else{
        postfix = postfix
        + "_CTS" + to_string(is_cts_design)
        + "_DELTA" + to_string(DELTA)
        + "_D" + to_string(MAX_DIMENSION)
        + "_P" + to_string(MAX_NUM_CASES)
        + "_CAND" + to_string(MAX_NUM_CANDIDATES)
        + "_TH" + to_string(threshold);
    }
    postfix = postfix \
        + "_SHIFT" + to_string(is_performing_shift)
        + "_FLIP" + to_string(is_performing_flip)
        + "_CHOI" + to_string(IS_CHOI)
        + "_INF" + to_string(INFERENCE_TYPE)
        ;
    def_postfix = def_postfix \
        + "_SHIFT" + to_string(is_performing_shift)
        + "_FLIP" + to_string(is_performing_flip)
        + "_CHOI" + to_string(IS_CHOI)
        + "_INF" + to_string(INFERENCE_TYPE)
        ;
    
    if (is_cts_from_place){
        postfix = "place_" + postfix;
        postfix += "_3stages_CTS";

        def_postfix = "place_" + def_postfix;
        def_postfix += "_3stages_CTS";
    } 
    if (IS_CADENCE){
        postfix += "_cadence";
    }

    if(REGION_IDX != -1){
        postfix += "_REGION" + to_string(REGION_IDX);
    }
    
    std::cout << design_name + "_" + postfix << endl;

    std::string lef_directory = direct_input_lef;
    std::string def_directory;
    std::string input_directory;
    std::string output_directory;
    if (!IS_CADENCE){
        input_directory = project_directory + "/inputs";
        output_directory = project_directory + "/outputs";
        if(lef_directory == ""){
            lef_directory = input_directory+"/lef/nangate15nm.lef";
        }
        
        if(!is_cts_design && !is_cts_from_place){
            def_directory = input_directory + "/def/placements/" + design_name + ".def";
        }
        else{
            if (is_cts_design && !is_cts_from_place){
                def_directory = input_directory + "/def/cts/" + design_name + ".def";
            }
            else if (is_cts_design && is_cts_from_place){
                def_directory = input_directory + "/def/cts/" + design_name + "_" + def_postfix + ".def";
            }
        }
    }
    else{
        input_directory = project_directory + "/inputs";
        output_directory = project_directory + "/outputs";
        if(lef_directory == ""){
            lef_directory = input_directory+"/lef/asap7.lef";
        }
        
        def_directory = input_directory + "/def/cts/" + design_name + ".def";
    }
    if(direct_input_def !=""){
        def_directory=direct_input_def;
    }
    
    std::string drv_region_log_directory = script_directory + "/drv_region_"  + design_name + "_" + postfix + ".log";
    std::string blockages_directory = script_directory + "/blockages_"  + design_name + "_" + postfix + ".tcl";
    string prediction_base = (is_cts_from_place)?design_name + "_" + postfix: design_name;
    std::string regression_logit_directory = output_directory + "/predictions/" + prediction_base + "_" + INFERENCE_NAME + "_logit.npy";
    ///// INITAILIZING VARIABLES DONE /////

    //Parsing Lef, Def
    ParserWrapper lef_def_parser(
        lef_directory, 
        def_directory,
        is_cts_design,
        IS_CADENCE
        );
    lef_def_parser.ParseLefDef(); 
    
    //Run DRV Prediction Model for "Classification"
    int drv_region_inference_type;
    if (INFERENCE_TYPE == UGNET_UNSEEN){
        drv_region_inference_type = UGNET_UNSEEN;
    }
    else{
        drv_region_inference_type = UGNET;
    }
    std::string drv_prediction_script_dir = script_directory + "/run_drv_prediction_class_initial_" + design_name + "_" + postfix + ".bashrc";
    cout<<"DRV pred: " <<drv_prediction_directory<<endl;
    RunDrvPrediction(
        drv_prediction_directory, 
        lef_directory,
        def_directory, 
        output_directory, 
        drv_prediction_script_dir,
        drv_region_inference_type);
    std::chrono::steady_clock::time_point time_0_1 = std::chrono::steady_clock::now();
    RemoveLog(drv_prediction_script_dir);
    

    //Get DRV Region
    RunDrvRegion(
        drv_region_directory,
        output_directory + "/predictions/" + prediction_base + "_class_logit.npy",
        drv_region_log_directory,
        blockages_directory,
        threshold,
        lef_directory,
        def_directory,
        TRACK_PITCH,
        blockage_layers
        );    
    std::chrono::steady_clock::time_point time_0_2 = std::chrono::steady_clock::now();
    
    //Run DRV Prediction Model for "Regression"
    drv_prediction_script_dir = script_directory + "/run_drv_prediction_reg_intial_"+ design_name + "_" + postfix + ".bashrc";
    RunDrvPrediction(
        drv_prediction_directory, 
        lef_directory,
        def_directory, 
        output_directory, 
        drv_prediction_script_dir,
        INFERENCE_TYPE);
    RemoveLog(drv_prediction_script_dir);
    
    std::chrono::steady_clock::time_point time_0_3 = std::chrono::steady_clock::now();
    
    RowsOfDRVRegion* drv_regions = DecodeDrvRegion(drv_region_log_directory, lef_def_parser, DELTA, is_cts_design);
    ClusteredDRVRegionsMultipleDimensions clustered_drv_regions = ClusterDRVRegions(drv_regions, MAX_DIMENSION, lef_def_parser, MAX_NUM_CANDIDATES, MAX_NUM_CASES, is_cts_design, is_performing_shift, is_performing_flip);
    std::chrono::steady_clock::time_point time_0_4 = std::chrono::steady_clock::now();
    RemoveLog(drv_region_log_directory);

    vector<InputForBO>& input_vectors = clustered_drv_regions.inputs_for_BO;

    //Get Initial DRV values
    std::string drv_coord_dir = script_directory + "/coord_initial_" + design_name + "_" + postfix +".log";
    std::string initial_drv_region_value_dir = script_directory + "/" + design_name + '_' + postfix + "_initial.log";
    UpdatePredictionValue(clustered_drv_regions, regression_logit_directory, read_logit_directory, drv_coord_dir, initial_drv_region_value_dir, is_updating_value_by_sum);   
    vector<double> initial_prediction_values;
    RemoveLog(drv_coord_dir);
    RemoveLog(initial_drv_region_value_dir);
    
    for(int i_input = 0; i_input < input_vectors.size(); i_input++){
        vector<ClusteredDRVRegions1Dimension>& input_vector = input_vectors[i_input].input;
        double output = calculate_output(input_vectors, i_input, is_updating_value_by_sum);
        initial_prediction_values.push_back(output);
    }
    
    std::cout<<"initial_drv ";
    PrintVector(initial_prediction_values);
    std::cout<<"done"<<endl;

    std::chrono::steady_clock::time_point time_1 = std::chrono::steady_clock::now();
    std::cout << "Initialize (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(time_1 - time_0).count()) /1000000.0  <<std::endl;
    std::cout << "  DRV_pred_class (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(time_0_1 - time_0).count()) /1000000.0  <<std::endl;
    std::cout << "  DRV_region (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(time_0_2 - time_0_1).count()) /1000000.0  <<std::endl;
    std::cout << "  DRV_pred_reg (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(time_0_3 - time_0_2).count()) /1000000.0  <<std::endl;
    std::cout << "  Build input space (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(time_0_4 - time_0_3).count()) /1000000.0  <<std::endl;
    
    std::cout << "Initialization done." << std::endl;
    bool is_running = true;
    while(is_running){
        string code;
        std::getline(std::cin, code);    
        vector<string> info = Parse(code);
        
        if(info[0] == "get_input_range"){
            vector<vector<int>> input_range = GetInputRange(clustered_drv_regions, lef_def_parser);
            std::cout<<"input_range ";
            PrintVector(input_range);
        }
        else if(info[0] == "get_encode"){
            vector<vector<int>> encoded_inputs = GetEncodedInputs(clustered_drv_regions, lef_def_parser);
            std::cout<<"encoded_value ";
            PrintVector(encoded_inputs);
        }
        else if(info[0] == "get_drv"){ 
            //get_drv ${iteration}
            std::string iterative_design_name = design_name + "_" + postfix + "_"+ info[1];
            std::string def_output_directory = output_directory + "/def/" + iterative_design_name + ".def";
            std::string new_logit_directory = output_directory + "/predictions/" + iterative_design_name + "_" + INFERENCE_NAME + "_logit.npy";
            std::string iterative_drv_script_directory = script_directory + "/run_drv_prediction_"+ iterative_design_name +".bashrc";
            std::string iterative_coord_log_directory = script_directory + "/coords_" + iterative_design_name + ".log";
            std::string iterative_read_logit_script_directory = script_directory + "/" + iterative_design_name + ".log";

            lef_def_parser.def.ValidationTest(); 
            lef_def_parser.def.WriteDef(def_output_directory);

            bool is_performing_classification_also = (info[1] == "final")?true:false;

            vector<double> prediction_values = GetPredictionValues(
                clustered_drv_regions,
                lef_def_parser,
                is_updating_value_by_sum,

                drv_prediction_directory,
                drv_prediction_model_directory,
                lef_directory,
                def_output_directory,
                output_directory,
                iterative_drv_script_directory,

                new_logit_directory,
                read_logit_directory,
                iterative_coord_log_directory,
                iterative_read_logit_script_directory,

                is_performing_classification_also
            );
            std::cout<<"drv ";
            PrintVector(prediction_values);

            if (!(info[1] == "0" || info[1] == "final")){
                RemoveDefCsvNpz(def_output_directory, output_directory, iterative_design_name, is_updating_value_by_sum);
            }
            RemoveScript(output_directory, iterative_design_name, is_updating_value_by_sum);
        }
        else if(info[0] == "decode"){ 
            //decode ${values to be decoded}
            vector<int> single_input;
            vector<vector<int>> decoded_inputs;
            int i_BO = 0;
            for(int i=1; i<info.size(); i++){
                single_input.push_back(stoi(info[i]));
                if(single_input.size() == clustered_drv_regions.inputs_for_BO[i_BO].input.size() || i == info.size()-1){
                    decoded_inputs.push_back(single_input);
                    single_input.clear();
                    i_BO ++;
                }
            }
            RefinePlacement(decoded_inputs, clustered_drv_regions, lef_def_parser);
            lef_def_parser.def.ValidationTest();
        }
        else if(info[0] == "exit" || info[0] == "quit"){
            is_running=false;
            std::cout<<"Terminatie the program."<<endl;
        }
        else{
            std::cout<<"[Warning] Invalid input. "<<endl;
        }
        std::cout<<"done"<<endl;
    }    
}
