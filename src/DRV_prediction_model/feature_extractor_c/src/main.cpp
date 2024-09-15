#include "global.h"
#include "FeatureExtractor.h"
#include "run_flute.h"

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
	namespace po = boost::program_options;

	std::string def_path_s, lef_path_s, drv_path_s, cong_path_s, weight_name_s, csv_path_s;//, design_name_s;
    bool is_ugnet;
	int multi_thread;
	
	int divide_unit = 104;
	int exit_threshold = 30;
    
    std::string project_directory =
        "./";

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("def,D", po::value<std::string>(&def_path_s)->default_value(project_directory+"/inputs/def"), "specify a def file.")
		("lef,L", po::value<std::string>(&lef_path_s)->default_value(project_directory+"/inputs/lef/nangate15nm.lef"),"specify a lef file")
		("drv,V", po::value<std::string>(&drv_path_s)->default_value(project_directory+"/inputs/drv"), "specify the path of the folder that contains drv log files")
		("cong,C", po::value<std::string>(&cong_path_s)->default_value(project_directory+"/inputs/cong"), "specify the path of the folder that contains cong files")
		("unfriendly-weight,W", po::value<std::string>(&weight_name_s)->default_value(project_directory+"/src/DRV_prediction_model/feature_extractor_c/type_weight_opc_expand.csv"), "specify the path of the unfriendly weight file")
		("output,O", po::value<std::string>(&csv_path_s)->default_value(project_directory+"/outputs/csv"), "specify the path of the output csv folder")
		("multi,M", po::value<int>(&multi_thread)->default_value(0), "multi-threading")
		("graph,G", po::value<bool>(&is_ugnet)->default_value(0), "graph")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if(!def_path_s.compare(def_path_s.length()-4,4,".def")==0){
        cout<<"A single def file should be inserted."<<endl;
        return 1;
	}
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}


	fs::path def_path = fs::absolute(fs::path(def_path_s));
	fs::path lef_path = fs::absolute(fs::path(lef_path_s));
	fs::path drv_path = fs::absolute(fs::path(drv_path_s));
	fs::path cong_path = fs::absolute(fs::path(cong_path_s));
	fs::path type_weight_path = fs::absolute(fs::path(weight_name_s));
	fs::path feature_path = fs::absolute(fs::path(csv_path_s));

	fs::create_directories(csv_path_s);

	if (!fs::exists(def_path)) {
		std::cout << "DEF file folder doesn't exist!" << std::endl;
		return 1;		
	}

	if (!fs::exists(lef_path)) {
		std::cout << "LEF file doesn't exist!" << std::endl;
		return 1;
	}
	    
    auto start = std::chrono::high_resolution_clock::now();

	std::unordered_map<std::string, int> used_type;
	std::unordered_map<std::string, int> drv_type;
	std::unordered_map<std::string, double> type_weight;



	if(!fs::exists(type_weight_path)){
		cout << "###  Type weight file doesn't exist! ###";
		cout << "###  Calculate weight per cell type  ###" << endl;

		if (!fs::exists(drv_path)) {
			std::cout << "DRV log file folder doesn't exist!" << std::endl;
			std::cout << "Terminated.." << std::endl;
			return 1;
		}

		int count_i=1;
        
        std::cout << "Parsing LEF";
    	Lef lef(lef_path);
    	lef.parse();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << " : " << duration.count() << " ms" << std::endl;

		for (auto& p : fs::directory_iterator(def_path)) {
			cout << "Num design : "<< count_i++ <<endl;
			auto design_name = p.path().stem();
			auto def_dir = def_path / fs::path(design_name.string() + ".def");

			auto drv_dir = drv_path / fs::path(design_name.string() + ".log");

			if (!fs::exists(drv_dir)) {
				std::cout << design_name << " DRV log file doesn't exist!, skipped.." << std::endl;
				continue;
			}

			//if (design_name.string().substr(0, 3) == "mgc") continue;

			std::cout << "Design : " << design_name.string() << std::endl;
			std::cout << "Parsing DEF" << std::endl;
			Def def(def_dir, lef, false);
			def.parse();
	
			std::cout << "Parsing DRV" << std::endl;
			Drv drv_log(drv_dir);
			drv_log.parse();
			auto cong_dir = cong_path / fs::path(design_name.string() + ".cong");
	    	
			FeatureExtractor fe(design_name.string(), feature_path, cong_dir, def, drv_log, 0.768, 0.768, 1, 1, is_ugnet);
			fe.accumulate_unfriendly(used_type, drv_type);
		}
    	for (const auto &type : used_type){
        	if (drv_type.count(type.first) == 0){
            	type_weight[type.first] = 0;
        	}
        	else
            	type_weight[type.first] = static_cast<double>(drv_type[type.first]) / used_type[type.first];
    	}
		std::ofstream writeFile(type_weight_path.string());
		for(auto it = type_weight.begin(); it != type_weight.end(); it++){
			writeFile << it->first << "," << it->second << "\n";
		}
		writeFile.close();
		cout << "Done!" <<endl;

		for (const auto& type : type_weight) {
			auto type_name = type.first;
			std::cout << "(" << type_name << ") : " << "used " << used_type[type_name] << ", drv " << drv_type[type_name] << ", weight " << type.second << std::endl;
        }
	} else {
		string str_buf, cellname;
		size_t pos;
		double cellweight;
		std::ifstream readFile(type_weight_path.string());
		while(!readFile.eof()){
			getline(readFile, str_buf);
			if(str_buf.size() == 0) continue;
			pos = str_buf.find(',');
			cellname = str_buf.substr(0, pos);
			cellweight = stod(str_buf.substr(pos+1, str_buf.size()));
			type_weight[cellname] = cellweight;
			
		}
		readFile.close();
	}
	
	int design_num = 1;

	int count_i=0;

	int real_count = 0;

	fs::directory_entry p(def_path);
    //cout << "Num design : "<< ++count_i <<endl;
    //if (multi_thread > 0 && !(count_i >= (multi_thread - 1) * divide_unit && count_i  <= multi_thread * divide_unit)) continue;

    auto design_name = p.path().stem();
    /*
    if(design_name_s == "NULL"){
        design_name_s = design_name;
    }
    */
    auto csv_path = feature_path / design_name.string();
    auto feature_csv_path = feature_path / fs::path(design_name);
    
    if(fs::exists(feature_csv_path)){
		std::cout <<"Feature exists. Skip design " << design_name.string() << std::endl;
		return 1;
	}
    
    
    std::cout << "Parsing LEF";
    Lef lef(lef_path);
    lef.parse();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << " : " << duration.count() << " ms" << std::endl;

    
    auto debug_path = csv_path / fs::path("edge_weight_hnetden.txt");

    /*
    if (fs::exists(debug_path)) {
        std::cout << design_name.string() << " already exists!" << std::endl;
        continue;
    }

    if (real_count++ >= exit_threshold) break;
    */
    
    //auto def_dir = def_path / fs::path(design_name.string() + ".def");
    auto def_dir = def_path;
    
    std::cout << "Design : " << design_name.string() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    
    std::cout << "Parsing DEF";
    Def def(def_dir, lef, false);
    def.parse();

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << " : " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();

    if(is_ugnet){
        std::cout << "Run FLUTE";

        // FLUTE
        initialize_net_v2(def);
        gcell_grid *gcell_gridp = edge_shifting(def);
        split_steiner_tree(def);

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << " : " << duration.count() << " ms" << std::endl;
    }

    auto drv_dir = drv_path / fs::path(design_name.string() + ".log");
    Drv drv_log;

    if (fs::exists(drv_dir)) {
        start = std::chrono::high_resolution_clock::now();

        std::cout << "Parsing DRV log" << std::endl;
        drv_log = Drv(drv_dir);
        drv_log.parse();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "DRV log parsing : " << duration.count() << " ms" << std::endl;
    }
    else {
        std::cout << "DRV file doesn't exist.. skipped.." << std::endl;
    }

    auto cong_dir = cong_path / fs::path(design_name.string() + ".cong");
    if (!fs::exists(cong_dir)) {
        std::cout << design_name.string() << " of congestion file doesn't exist!, skipped.." << std::endl;
    }

    std::cout << "Extracting Features" << std::endl;
    FeatureExtractor fe(design_name, feature_path, cong_dir, def, drv_log, 0.768, 0.768, 1, 1, is_ugnet);
    fe.stored_type_weight = type_weight;
    fe.run();

	return 0;
}
