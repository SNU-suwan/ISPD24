#ifndef PLACE_CELL_H_
#define PLACE_CELL_H_

#include <iostream>
#include <vector>
#include <tuple>

struct standard_cell{
	standard_cell(){}
	standard_cell(std::string type_, int cell_width_, int cell_height_)
		: type(type_), cell_width(cell_width_), cell_height(cell_height_){}
	std::string type="NULL";
	int cell_width;
	int cell_height;
};

struct placed_cell{
	placed_cell(){}
	placed_cell(int x_, int y_, int is_fixed_, bool is_flipped_, std::string name_, std::string direction_, standard_cell * std_cell_, bool is_cts_cell_)
		: x(x_), y(y_), is_fixed(is_fixed_), is_flipped(is_flipped_), name(name_), direction(direction_), original_direction(direction_), cell(std_cell_), is_cts_cell(is_cts_cell_){}

	int x = -1;
	int y = -1;
	int is_fixed = -1; 
	int idx = -1;
	bool is_flipped=false;
	bool is_cts_cell=false;
	
	std::string name="NULL";
	std::string direction;
	std::string original_direction;
	standard_cell* cell;

	std::vector<int> invalid_cumulative_whitespace; 

	void Flip();
};

struct net{
	net(std::string name_): name(name_){}
	std::string name;
	std::vector<std::tuple<std::string,std::string>> cell_pin_tuple_vec;
};
#endif
