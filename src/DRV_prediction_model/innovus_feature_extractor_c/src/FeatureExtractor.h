#pragma once

#include "BBox.h"
#include "parser/def.h"
#include "parser/drv.h"

#define SCALE 1000

#include <filesystem>
namespace fs = std::filesystem;

enum class DIR {
    PRIM, SEC, LAST
};


class FeatureExtractor {
public:

    using Grid2D = vector<vector<double>>;
    using Grid3D = vector<vector<vector<double>>>;
    using BBox2D = vector<vector<BBox<double>>>;

    FeatureExtractor(std::string name_, fs::path out_, fs::path cout_, Def& def_, Drv& drv_, double pat_, double gcell_, int Ysize_, int Xsize_, bool is_innovus_) 
    : name(name_), outPath(out_), congPath(cout_), def(def_), drv(drv_), patSize(pat_), gcellSize(gcell_), patY(pat_ * Ysize_), patX(pat_ * Xsize_), gridY(gcell_ * Ysize_), gridX(gcell_ * Xsize_), is_innovus(is_innovus_) {}
    
    void run();
    void accumulate_unfriendly(std::unordered_map<std::string, int>& used_type, std::unordered_map<std::string, int>& drv_type);

    // process functions
    void init();
    void printInfo();
    void print_def();

    // Write feature extraction functions here // 
    void pin_density();
    void pin_pattern();
    void avg_pin_access();
    void rudy_map();
    void m2_short(bool check_num);
    void check_drv(bool check_num);
    void minimum_proximity();
    void weighted_unfriendly();
    void local_global_self_crossing_net();
	void macro_density();
    void pin_rudy();
	void routing_capacity();
    void congestion_map();
    void hv_net_density();

    void generate_graph();
    void generate_graph_new();
    void generate_2Dgraph();
    void generate_2Dgraph_2Dedge();
    void net_direction();
    void flute_accuracy();


    // utils
    Grid3D cal_pin_density(Point<int> grid_num, int grid_x, int grid_y);
	vector<vector<vector<pin*>>> assign_pin_to_grid();
    void write_csv (std::string feature_name, Grid2D& feature);
    pair<BBox2D, BBox2D> seperate_bbox (pin* pin_name);
    int calculate_overlap (int x, int y);
    Point<int> convert_int3 (Point<double>& p);
    BBox<int> convert_int3 (BBox<double>& bbox);
    BBox<int> convert_to_grid (BBox<double>& bbox, double grid_size_x, double grid_size_y);
    BBox<int> convert_to_grid (BBox<int>& bbox, int grid_size_x, int grid_size_y);
    void convert_to_33 (Grid2D &feature, Grid2D &feature33Min, Grid2D &feature33Max);
    vector<string> split_by_space(string line);
    Grid2D direction_map(Grid2D &pin_image, string dir);


public:

    Def& def;
    Drv& drv;
    std::string name;
    fs::path outPath, congPath;

    double patSize, gcellSize;
    double patX, patY, gridX, gridY;
    double patArea, gridArea;
    double rudy_threshold = 15;
	double row_height;

    //graph 2D node feature
    int pinGridNum;

    int patSizeI, gcellSizeI;
    int patXI, patYI, gridXI, gridYI;
    int patAreaI, gridAreaI;
    bool is_innovus;

    BBox<double> die, core;
    BBox<int> dieI, coreI;

    Point<int> patNum, gridNum;

    vector<vector<vector<string>>> cellMap;

    std::unordered_map<std::string, double> stored_type_weight;

    Grid2D pinAccessInfo;
    Grid2D shortRudy, longRudy, pinRudy;
    Grid2D numDrv, M2Short, pinDRV, congDRV;
    Grid3D pinDensity, pinPattern;
    Grid2D minProximity, unfriendly;
	Grid2D localNet, globalNet, selfCrossingNet;
	Grid2D macroDensity;
    Grid2D routing_capacity_hor, routing_capacity_ver, HnetDensity, VnetDensity;
	/*
    Grid2D localNet33Min, globalNet33Min, selfCrossingNet33Min, shortRudy33Min, longRudy33Min, pinDensity33Min, minProximity33Min;
    Grid2D localNet33Max, globalNet33Max, selfCrossingNet33Max, shortRudy33Max, longRudy33Max, pinDensity33Max, minProximity33Max;
	*/

    Grid2D congDemHor, congCapHor, congDemVer, congCapVer;
    Grid2D congOvflHor, congOvflVer;
    Grid2D ovflHor, ovflVer;


};
