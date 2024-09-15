import numpy as np
import os
import sys
import time
import decimal

import copy
from def_parser import *
from lef_parser import *

######################## DRV REGIONS ##########################

def initialize_drv_region(logit_np, threshold, low_bound):
    drv_region=np.zeros([logit_np.shape[0],logit_np.shape[1]]) # 0: white, 1: black
    for row_idx in range(logit_np.shape[0]):
        queue=[]
        color_flag=0
        for column_idx in range(logit_np.shape[1]):
            gcell=logit_np[row_idx][column_idx]
            if gcell>threshold:
                color_flag=1
                for coord in queue:
                    drv_region[coord[0]][coord[1]]=color_flag
                if len(queue)>0:
                    if queue[0][1]>0:
                        drv_region[row_idx][queue[0][1]-1]=color_flag # margin left
                queue=[]

            if gcell>=low_bound:
                drv_region[row_idx][column_idx]=color_flag
                if not color_flag==1:
                    queue.append([row_idx,column_idx])

            if gcell<low_bound:
                queue=[]
                if color_flag==1 and column_idx<logit_np.shape[1]:
                    drv_region[row_idx][column_idx]=color_flag # margin right
                color_flag=0
    return drv_region

def cluster_drv_region(drv_region, x_offset = 3, y_offset = 3):
    label_drv_region = copy.deepcopy(drv_region)
    drv_points_location = np.where(drv_region == 1)

    cluster_idx = 2 # starts with 2
    for i in range(len(drv_points_location[0])):
        cur_x = drv_points_location[0][i]
        cur_y = drv_points_location[1][i]

        if label_drv_region[cur_x][cur_y] > 1:
            continue

        if label_drv_region[cur_x][cur_y] == 1:
            queue = [[cur_x,cur_y]]

            while len(queue) != 0:
                [x, y] = queue.pop(0)            

                for x_i in range(x_offset * 2 + 1):
                    new_x = x + (x_i - x_offset)
                    if new_x < 0 or new_x > np.shape(drv_region)[0] - 1:
                        continue
                    for y_i in range(y_offset * 2 + 1):
                        new_y = y + (y_i - y_offset)
                        if new_y < 0 or new_y > np.shape(drv_region)[1] - 1:
                            continue
                        if label_drv_region[new_x][new_y] == 1:
                            label_drv_region[new_x][new_y] = cluster_idx
                            queue.append([new_x,new_y])
            cluster_idx += 1
            
    return label_drv_region.astype(int)

def build_label_row_dict(label_drv_region):
    label_row_dict = {}
    for label in range(int(np.max(label_drv_region))):
        if label == 0 or label == 1:
            continue
        clustered_drv_region = np.where(label_drv_region == label)
        label_row_dict[label] = {}

        for i, y in enumerate(clustered_drv_region[0]):
            x = clustered_drv_region[1][i]
            if y not in label_row_dict[label]:
                label_row_dict[label][y] = [x,x]
            label_row_dict[label][y] = [min(x,label_row_dict[label][y][0]), max(x,label_row_dict[label][y][1])]
            assert(label_drv_region[y][x]>1)
            
    return label_row_dict

def add_marginal_rows(drv_region, label_drv_region, label_row_dict, margin = 1):
    for label in label_row_dict:
        label_row_list = list(label_row_dict[label].keys())
        for idx, y in enumerate(label_row_dict[label]):
            if len(label_row_dict[label]) > 1 and idx == len(label_row_dict[label])-1:
                continue
            else:
                if len(label_row_dict[label]) == 1 or label_row_list[idx+1] - label_row_list[idx] != 1:
                    for adjacent_row in range(y-margin, y+margin +1):
                        if adjacent_row < 0 or np.shape(drv_region)[0]-1 < adjacent_row:
                            continue
                        for x in range(label_row_dict[label][y][0], label_row_dict[label][y][1]+1):
                            drv_region[adjacent_row][x] = 1
                            label_drv_region[adjacent_row][x] = label

                            
def write_row_dict(drv_region):
    row_dict = {}
    for y in range(np.shape(drv_region)[0]):
        is_writing = False
        for x in range(np.shape(drv_region)[1]):
            if drv_region[y][x] != 0 and not is_writing: #234 117 -> 234 118
                is_writing = True
                start_x = x
            elif (drv_region[y][x] == 0 or x == np.shape(drv_region)[0] - 1) and is_writing:
                is_writing = False
                end_x = x - 1

                if y not in row_dict:
                    row_dict[y] = []
                row_dict[y].append(np.array([[start_x,y],[end_x,y]]))
    return row_dict
                            
def add_sides_if_lacks_whitespace(drv_region, label_drv_region, def_parser):
    label_row_dict = build_label_row_dict(label_drv_region)
    for label, row_dict in label_row_dict.items():
        for row, segments in row_dict.items():
            left_coord = [segments[0], row]
            right_coord = [segments[1], row]

            new_left_coord = ConvertLogitCoordToGrid(left_coord, np.shape(drv_region)[0])
            new_right_coord = ConvertLogitCoordToGrid(right_coord, np.shape(drv_region)[0])
            amount_whitespace = cal_whitespace_bbox([new_left_coord, new_right_coord], def_parser)

            if amount_whitespace == 0:
                size = new_right_coord[0] - new_left_coord[0]
                nearest_whitespace_loc_left = SearchWhitespace(new_left_coord, def_parser, True, drv_region)
                nearest_whitespace_loc_right = SearchWhitespace(new_right_coord, def_parser, False, drv_region)

                dist_left = None
                dist_right = None
                if nearest_whitespace_loc_left is not None:
                    dist_left = new_left_coord[0] - nearest_whitespace_loc_left
                if nearest_whitespace_loc_right is not None:
                    dist_right = nearest_whitespace_loc_right -  new_right_coord[0]

                if dist_left is not None and dist_left < size/2 : #234 117 -> 234 118
                    for x in range(new_left_coord[0], nearest_whitespace_loc_left -1, -1):
                        drv_region[left_coord[1]][x] = 1
                        label_drv_region[left_coord[1]][x] = label

                if dist_right is not None and dist_right < size/2:
                    for x in range(new_right_coord[0], nearest_whitespace_loc_right+1, 1):
                        drv_region[left_coord[1]][x] = 1
                        label_drv_region[left_coord[1]][x] = label

############################## UTILS ##########################

def ConvertGridToCoord(grid, offset, gcell_height):
    grid = np.array(grid)
    return (offset + grid*gcell_height) 


def cal_whitespace_grid(grid, def_parser):
    gcell_size = def_parser.rows[1].y-def_parser.rows[0].y
    coord = ConvertGridToCoord(grid, def_parser.offset[0], gcell_size)
    
    x_left = coord[0]
    x_right = coord[0] + gcell_size
    whitespaces_in_row = def_parser.whitespace_dict[coord[1]]
    
    whitespace_width = 0
    for whitespace in whitespaces_in_row:
        if not ((whitespace.x_left < x_left and whitespace.x_right <= x_left) or (whitespace.x_right > x_right and whitespace.x_left >= x_right)):
            whitespace_width = min(whitespace.x_right, x_right) - max(whitespace.x_left, x_left)
    
    return whitespace_width
        
def cal_whitespace_bbox(bbox, def_parser):
    total_whitespace = 0
    for x in range(bbox[0][0], bbox[1][0]+1):
        grid = [x, bbox[0][1]]
        total_whitespace += cal_whitespace_grid(grid, def_parser)
    return total_whitespace


def ConvertLogitCoordToGrid(logit_coord, drv_region_y_shape):
    return np.array([logit_coord[0], drv_region_y_shape - logit_coord[1] - 1])

def SearchWhitespace(boundary_coord, def_parser, is_left, drv_region):
    if is_left:
        start_coord = boundary_coord[0]
        end_coord = 0
        incr = -1
    else:
        start_coord = boundary_coord[0]
        end_coord = np.shape(drv_region)[0]-1
        incr = 1
        
    for x in range(start_coord, end_coord, incr):
        whitespace = cal_whitespace_grid([x, boundary_coord[1]], def_parser)
        if whitespace > 0 :
            return x
    return None



def initialize_parsers(def_dir, lef_dir, is_cadence=False):
    def_parser = DefParser(def_dir, is_cadence)
    def_parser.get_scale()
    lef_parser = LefParser(lef_dir, def_parser.get_scale())
    def_parser.parse(lef_parser)
    def_parser.get_whitespace()    
    return def_parser, lef_parser


########################### BLOCKAGES ###############################


def write_blockage_tcl(is_cadence, blockage_layers, track_pitch, side_length, drv_coords, def_parser, output_dir, start_idx = 3, box_size = 2):
    track_pitch = decimal.Decimal(track_pitch)
    gcell_height = def_parser.rows[1].y - def_parser.rows[0].y
    offset_x, offset_y = def_parser.rows[0].x, def_parser.rows[0].y
    
    if is_cadence:
        pre_script = 'createRouteBlk -layer {'
    else:
        pre_script = 'create_routing_blockage -layers {'
    
    for layer in blockage_layers:
        pre_script += layer + " "
    pre_script += "}"

    if is_cadence:
        pre_script += " -box "
    else:
        pre_script += " -boundary "

    with open(output_dir,'w') as f:
        for i in range(len(drv_coords[0])):
            bbox_grid_x = drv_coords[0][i]
            bbox_grid_y = drv_coords[1][i]
            bbox_x = bbox_grid_x * gcell_height + offset_x
            bbox_y = bbox_grid_y * gcell_height + offset_y

            bbox_offset = int(start_idx * track_pitch * def_parser.scale)
            side = int(side_length * track_pitch * def_parser.scale * box_size)

            lb_x = bbox_x + bbox_offset
            lb_y = bbox_y + bbox_offset
            rt_x = lb_x + side
            rt_y = lb_y + side
            if is_cadence:
                script = pre_script + '{' + str(lb_x/def_parser.scale) + ' ' + str(lb_y/def_parser.scale) + ' ' + str(rt_x/def_parser.scale) + ' ' + str(rt_y/def_parser.scale) + '}\n' 
            else:
                script = pre_script + '{{{' + str(lb_x/def_parser.scale) + ' ' + str(lb_y/def_parser.scale) + '} {' + str(rt_x/def_parser.scale) + ' ' + str(rt_y/def_parser.scale) + '}}}\n'
            f.write(script)
