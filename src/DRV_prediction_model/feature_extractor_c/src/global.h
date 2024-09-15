#pragma once

#include "parser/global_variables.h"
#include "parser/default_classes.h"
#include "parser/structures.h"
#include "parser/lef.h"
#include "parser/def.h"
#include "parser/drv.h"
#include "parser/utils.h"

#include "BBox.h"

#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <math.h>
#include <map>
#include <memory>
#include <cassert>
#include <tuple>

#include <boost/program_options.hpp>

namespace fs = std::filesystem;


extern string bottom_most_routing_layer;
extern string upper_most_routing_layer;

extern map<string,int> METAL_LAYER_IDX;
//extern map<string,int> METAL_ROUTING_TRACK;
//extern map<string,string> METAL_ROUTING_DIRECTION;
extern map<int,int> METAL_ROUTING_TRACK;
extern map<int,string> METAL_ROUTING_DIRECTION;
extern vector<string> CLOCK_UNITS;

/*
	SITE information
*/
#define SITE_Y 0.768

/*
   Prefix of buffers or FFs to avoid when considering routing layers
*/
#define BUF1 "BUF"
#define BUF2 "CLKBUF"
#define INV1 "INV"
#define DFF1 "DFF"
#define NAND1 "NAND"
#define AOI1 "AOI" //"AOI21_X1"
#define OAI1 "OAI" //"OAI22_X1"

/*
   Names of routing metal layers
*/
#define METAL1 "M1"
#define METAL2 "MINT1"
#define METAL3 "MINT2"
#define METAL4 "MINT3"
#define METAL5 "MINT4"
#define METAL6 "MINT5"

/*
   Routing tracks with respect to the y-size of SITE for each layer
   (int)
*/
#define M1_TRACKS 0 
#define M2_TRACKS 12 
#define M3_TRACKS 12 
#define M4_TRACKS 12 
#define M5_TRACKS 12 
#define M6_TRACKS 12 

/*
	Routing direction for each layer 
	(HOR -> horizontal, VER -> vertical, 2D -> both horizontal and vertical direction)
*/
#define M1_ROUTING_DIRECTION "VER"
#define M2_ROUTING_DIRECTION "HOR"
#define M3_ROUTING_DIRECTION "VER"
#define M4_ROUTING_DIRECTION "HOR"
#define M5_ROUTING_DIRECTION "VER"
#define M6_ROUTING_DIRECTION "HOR"
