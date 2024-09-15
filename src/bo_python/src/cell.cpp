#include "cell.h"

void placed_cell::Flip(){
    if(is_flipped){
        if(direction == "S") direction = "FS";
        else if(direction == "N") direction = "FN";
        else if(direction == "FS") direction = "S";
        else if(direction == "FN") direction = "N";
    }
}