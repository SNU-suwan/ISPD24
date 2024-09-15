#pragma once

#include "global.h"
#include <algorithm>

extern "C" {
    #include "flute-3.1/flute.h"
}

#define NET_NUM_THRESHOLD 1000
#define STEINER_VERTEX 10000
#define DUMMY_VERTEX 1000


void initialize_net_v2(Def &def);

graph *run_flute_for_each_net(Def&, vector<point2F>, string);
graph *direct_connection(Def&, net *);

gcell_grid* edge_shifting(Def&);

void split_steiner_tree(Def &);
