import os
import numpy as np
import sys

if __name__=="__main__":
    #print(sys.argv)
    if len(sys.argv)!=5:
        print('Usage: python3 read_logit.py $logit_np_dir $script_dir $script_out_dir $is_updating_value_by_sum(0 or 1)')
        sys.exit()
    else:
        logit_np_dir = sys.argv[1]
        script_dir = sys.argv[2]
        script_out_dir = sys.argv[3]
        is_updating_value_by_sum = int(sys.argv[4])
        
        if not is_updating_value_by_sum:
            coords = []
            with open(script_dir, 'r') as f:
                line = f.readline()
                while line:
                    info = line.split()
                    assert(len(info)==3)
                    coords.append([int(x) for x in info])
                    line = f.readline()

            logit_np = np.load(logit_np_dir)
            for coord in coords:
                y = coord[0]
                logit_sum = 0
                for x in range(coord[1],coord[2]+1):
                    new_x = np.shape(logit_np)[0] - 1 - y #np.shape(logit_np)[0] - y #len(logit_np[0]) - 1 - y
                    new_y = x
                    if new_y < np.shape(logit_np)[1]:
                        logit_sum += logit_np[new_x][new_y] # We sum all the values in the grid. It can be changed by taking maximum value in the grid
                coord.append(logit_sum)


            with open(script_out_dir, 'w') as f:
                for coord in coords:
                    f.write(str(coord[0]) + " " + str(coord[1]) + " " + str(coord[2]) + " " + str(coord[3]) + "\n")
        else:
            logit_np = np.load(logit_np_dir)
            sum_value = np.sum(logit_np)
            with open(script_out_dir, 'w') as f:
                f.write(str(sum_value))
        
        
