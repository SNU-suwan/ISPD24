import numpy as np
import os
import sys
import time
from def_parser import *
from lef_parser import *
from utils import *

#track_pitch = 0.064
side_length = 1
start_idx = 3
#output_dir = 'test_blockages.log'


def get_region(logit_dir, def_parser, threshold, is_stop_verbose=False, low_bound=None):
    if low_bound is None:
        low_bound=threshold/2
    logit_np=np.load(logit_dir)
    
    if not is_stop_verbose:
        print("Gcell size : ", logit_np.shape)
    
    time0 = time.time()
    
    #Initailize DRV Region
    time1 = time.time()
    drv_region = initialize_drv_region(logit_np, threshold, low_bound)
    if not is_stop_verbose:
        print("  Initializing DRV regions has done in ",time.time() - time1, "s.")
    
    #Cluster(labeling) DRV region
    time1 = time.time()
    label_drv_region = cluster_drv_region(drv_region)
    label_row_dict = build_label_row_dict(label_drv_region)
    if not is_stop_verbose:
        print("  Clustering(labeling) DRV regions has done in",time.time() - time1, "s.")
    
    #Add marginal rows to drv regions consist of a single row
    time1 = time.time()
    add_marginal_rows(drv_region, label_drv_region, label_row_dict, margin = 1)
    if not is_stop_verbose:
        print("  Adding marginal rows to small DRV regions has done in",time.time() - time1, "s.")
    
    #Add DRV region if the row lacks whitespace
    time1 = time.time()
    add_sides_if_lacks_whitespace(drv_region, label_drv_region, def_parser)
    if not is_stop_verbose:
        print("  Increasing areas of DRV regions which lack whitespaces has done in ",time.time() - time1, "s.")
        print("Extracting DRV regions has done in ", time.time() - time0, "s.")
    
        
    return [np.where(drv_region==1)[1], drv_region.shape[0]-1-np.where(drv_region==1)[0], label_drv_region[np.where(drv_region==1)]], drv_region, label_drv_region


def arrange_info(drv_coords):
    x_coord = drv_coords[0]
    y_coord = drv_coords[1]
    labels = drv_coords[2]

    assert len(x_coord) == len(y_coord) == len(labels)
    
    label_row_dict = {}
    for i in range(len(x_coord)):
        if labels[i] not in label_row_dict:
            label_row_dict[labels[i]] = {}
        if y_coord[i] not in label_row_dict[labels[i]]:
            label_row_dict[labels[i]][y_coord[i]] = []
        label_row_dict[labels[i]][y_coord[i]].append(x_coord[i])
    
    return label_row_dict

def write_drv_region_log(output_dir, label_row_dict):
    with open(output_dir,'w') as f:
        for label, row_dict in label_row_dict.items():
            for key,vals in row_dict.items():
                f.write(str(label))
                f.write(' ')
                f.write(str(key))
                f.write(' ')
                #print(key)
                prev_val = vals[0]
                for val in vals:
                    cur_val = val
                    if cur_val - prev_val == 0 or cur_val == vals[-1]:
                        #print(cur_val)
                        f.write(str(cur_val))
                        f.write(' ')
                    elif cur_val - prev_val != 1:
                        #print(prev_val)
                        #print(cur_val)
                        f.write(str(prev_val))
                        f.write(' ')
                        f.write(str(cur_val))
                        f.write(' ')
                    prev_val = cur_val    
                f.write('\n')
    
if __name__=="__main__":
    print("Extracting DRV region")
    if len(sys.argv)<8:
        print('Usage: python3 get_drv_region.py $logit_np_dir $threshold $output_dir $blockages_dir $lef_dir $def_dir $track_pitch $blockage_layers')
        sys.exit()
    else:
        
        is_cadence=False
        if 'cadence' in sys.argv[6]:
            is_cadence=True
        
        def_parser, lef_parser = initialize_parsers(sys.argv[6],sys.argv[5], is_cadence)
        
        logit_dir=sys.argv[1]
        threshold=float(sys.argv[2])
        output_dir=sys.argv[3]
        blockages_dir=sys.argv[4]
        track_pitch = float(sys.argv[7])
        #blockage_layers = sys.argv[8].split(',')
        blockage_layers = sys.argv[8:]
        
        print("Logit dir : ",logit_dir)
        print("Threshold : ",threshold)
        print("Lower bound : ",threshold/2)
        print("Output dir : ",output_dir)
        print("Blockages dir : ",blockages_dir)
        print("Def dir : ", sys.argv[6])
        print("Track Pitch : ", sys.argv[7])
        
        drv_coords, drv_region, label_drv_region=get_region(logit_dir, def_parser, threshold)
        
        label_row_dict = arrange_info(drv_coords)
        write_drv_region_log(output_dir, label_row_dict)
        write_blockage_tcl(is_cadence, blockage_layers, track_pitch, side_length, drv_coords, def_parser, blockages_dir, start_idx)
