#ifndef __UTILS__
#define __UTILS__

class comp;
class cell;
class pin;
class ap;
#include <vector>
#include <string>
#include "default_classes.h"
#include "structures.h"


using namespace std;
vector<string> split(string);
point2I convert_coords_to_idx(point2F);
point2F convert_idx_to_coords(point2I);
string get_direction(const point2F& locP2, const point2F& targetP2);
bool file_exist_test(const string&);
string get_origin_name(const string &);
float distance_sq(const point2F &, const point2F &);
void debug_test(comp* & compp);
void debug_test(cell* & cellp);
void begin_timer(struct timeval* begin);
double end_timer(struct timeval* begin, struct timeval* end);
bool is_lines_intersect(const float* point1, const float* point2);
BBox<double> pin_to_bbox(pin* pinp);
map<int,string> build_metal_idx_layer();
#endif
