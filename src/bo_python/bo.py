import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from skopt import BO_SW
from skopt.plots import plot_gaussian_process
import time
import sys
import csv
import os

PROJECT_DIR="."

def convert_input_to_UTF(input_string):
    return bytes(input_string + '\n', 'UTF-8')

def send_stdin(input_string, process, is_wait=False):
    process.stdin.write(convert_input_to_UTF(input_string))
    process.stdin.flush()
    print(input_string)
    if is_wait:
        result = ""
        while result != "done": #Wait until the process finishes its job.
            result = read_stdout(process)
    
def read_stdout(process):
    result = process.stdout.readline().strip()
    result = result.decode('ascii')
    print(result)
    return result

def read_until(process, message):
    output = [""]
    while output[0] != message:
        output = read_stdout(process).split()
        if len(output)==0:
            print("[ERROR] The program halts.")
            sys.exit(-1)
        elif output[0] == "[Error]":
            process.terminate()
            process.wait()
            sys.exit(-1)
        '''
        while len(output) == 0:
            output = read_stdout(process).split()
        '''
    result = ""
    while result != "done": #Wait until the process finishes its job.
        result = read_stdout(process)
    output.pop(0)
    return output

def parse_output(output, output_type):
    output_items = []
    output_item = []
    for idx, item in enumerate(output):
        if item.startswith('idx'):
            if idx != 0:
                output_items.append(output_item)
                output_item = []
            continue
        elif idx == len(output)-1:
            output_items.append(output_item)
        if output_type == "input_range":
            max_range = int(item)-1
            new_range = (0, max_range)
            output_item.append(new_range)
        elif output_type == "encoded_value":
            output_item.append(int(item))
        elif output_type == "drv":
            output_item.append(float(item))
        elif output_type == "initial_drv":
            output_item.append(float(item))
    return output_items

def read_input_range(process):
    send_stdin("get_input_range", process)
    output_type = "input_range"
    output = read_until(process, output_type)
    return parse_output(output, output_type)

def read_encoded_input(process):
    send_stdin("get_encode", process)
    output_type = "encoded_value"
    output = read_until(process, output_type)
    return parse_output(output, output_type)

def read_drv_prediction(process, i):
    send_stdin("get_drv " + str(i), process)
    output_type = "drv"
    output = read_until(process, output_type)
    return parse_output(output, output_type)

def cal_delta(initial_value, improved_value):
    return (float(initial_value) - float(improved_value)) / float(initial_value)

class bo_process:
    def __init__(self, binary_dir, result_dir, design_name, post_fix, max_iter, n_initial_points, max_dimension,\
            max_delta, is_cts_design, threshold, max_num_cases, max_num_candidates, is_shift, is_flip, is_choi,\
            inference_type, is_cts_from_place, def_dir, is_cadence, lef_dir, track_pitch, blockage_layers):
        self.binary_dir = binary_dir
        self.result_dir = result_dir
        
        self.design_name = design_name
        self.post_fix = post_fix
        self.max_iter = max_iter
        self.n_initial_points = n_initial_points
        self.max_dimension = max_dimension
        self.max_delta = max_delta
        self.is_cts_design = is_cts_design
        self.max_num_cases = max_num_cases
        self.max_num_candidates = max_num_candidates
        self.threshold = threshold
        self.is_shift = is_shift
        self.is_flip = is_flip
        self.is_choi = is_choi
        self.inference_type = inference_type
        self.is_cts_from_place = is_cts_from_place
        self.def_dir = def_dir
        self.is_cadence = is_cadence
        self.lef_dir = lef_dir
        self.track_pitch = track_pitch
        self.blockage_layers = ','.join(blockage_layers)

        self.num_models = None
        self.input_ranges = None
        self.BO_models = []
        self.initial_inputs = None
        self.initial_outputs = None
        self.final_outputs = None

    def run(self):
        main_dir = self.binary_dir
        binary_script = main_dir + ' -m ' + str(self.max_iter) + ' -d ' + str(self.max_dimension) + ' -p ' + \
            str(self.max_num_cases) + ' -name ' + str(self.design_name) + ' -threshold ' + str(self.threshold) + ' -cand ' + \
            str(self.max_num_candidates) + ' -lef ' + self.lef_dir + ' -track_pitch ' + self.track_pitch + ' -blockage_layers ' + self.blockage_layers
        binary_script += ' -delta ' + str(self.max_delta)
        binary_script += ' -shift ' + str(self.is_shift)
        binary_script += ' -flip ' + str(self.is_flip) 
        binary_script += ' -choi ' + str(self.is_choi)
        binary_script += ' -inf ' + str(self.inference_type)
        binary_script += ' -cts ' + str(self.is_cts_design)
        binary_script += ' -cfp ' + str(int(self.is_cts_from_place))
        if self.def_dir is not None:
            binary_script += ' -def ' + self.def_dir
        if self.is_cadence:
            binary_script += ' -cadence '
        
        #binary_script += ' | tee ' + self.log_dir
            
        print(binary_script)
        process = Popen(binary_script, shell=True, stdout=PIPE, stdin=PIPE)            

        #Initialization
        result = ""
        output = read_until(process, "initial_drv")
        self.initial_outputs = parse_output(output, "initial_drv")
        while result != "Initialization done.":
            result = read_stdout(process)

        self.input_ranges = read_input_range(process)
        self.initial_inputs = read_encoded_input(process)

        self.num_models = len(self.initial_outputs)

        print("Initialize model")
        t_initialization = time.time()
        for i in range(self.num_models):
            print("Initial sample :", self.initial_inputs[i] , "/", np.array(self.input_ranges[i])[:,1], " Output :", self.initial_outputs[i])
            self.BO_models.append(
                BO_SW(
                    self.input_ranges[i],      # the bounds on each dimension of x
                    acq_func = "EI",      # the acquisition function
                    n_initial_points = self.n_initial_points,
                    n_calls = self.max_iter,
                    random_state = 1234,
                    x0 = self.initial_inputs[i],
                    y0 = self.initial_outputs[i]))
        print("Initalization done, it takes ", time.time()-t_initialization, "s.")


        #Start iteration
        for i_iter in range(self.max_iter):
            print("Iteration ", i_iter)
            t_iter = time.time()
            next_samples = []
            decode_script = "decode "
            
            for i_model in range(self.num_models):
                next_sample = self.BO_models[i_model].query_next_sample()
                next_samples.append(next_sample)
                for item in next_sample:
                    decode_script += str(item) + " "
            
            t_query = time.time()
            send_stdin(decode_script, process, True)
            t_placement_refinement = time.time()
            prediction_results = read_drv_prediction(process, i_iter+1)
            t_pred = time.time()
            for i_model in range(self.num_models):
                print(i_model, "Next sample :", next_samples[i_model], "/", np.array(self.input_ranges[i_model])[:,1], " Output :", prediction_results[i_model])
                self.BO_models[i_model].optimize(next_samples[i_model], prediction_results[i_model][0])
            t_opt = time.time()
            print("Iteration done. It takes ", time.time()-t_iter)
            print("  1. Acquisition function : ", t_query - t_iter)
            print("  2. Placement refinement : ", t_placement_refinement - t_query)
            print("  3. Prediction : ", t_pred - t_placement_refinement)
            print("  4. Optimization : ", t_opt - t_pred)
            print("  BO (1+4) = ", t_opt - t_pred + t_query - t_iter)

        best_samples = []
        decode_script = "decode "    
        for i_model in range(self.num_models):
            best_sample = self.BO_models[i_model].result.x
            best_samples.append(best_sample)
            for item in best_sample:
                decode_script += str(item) + " "
                
        send_stdin(decode_script, process, True)
        self.final_outputs = read_drv_prediction(process, 'final')
        process.terminate()
        process.wait()

    def write_results(self):
        # Write results

        with open(self.result_dir + '/' +self.design_name + self.post_fix + '.csv', 'w') as f:    
            writer = csv.writer(f)
            writer.writerow(['Idx','Range','Best x','Initial','BO opt', 'Imp.', 'final', 'Imp.'])
            for i_model in range(self.num_models):
                script = []
                script.append(i_model)
                script.append(self.input_ranges[i_model][0][1])
                script.append(self.BO_models[i_model].result.x[0])
                script.append(self.initial_outputs[i_model][0])
                script.append(self.BO_models[i_model].result.fun)
                script.append(cal_delta(self.initial_outputs[i_model][0], self.BO_models[i_model].result.fun))
                script.append(self.final_outputs[i_model][0])
                script.append(cal_delta(self.initial_outputs[i_model][0], self.final_outputs[i_model][0]))
                writer.writerow(script)
                

if __name__ == '__main__':
    if len(sys.argv) < 16:
        print("Usage: python3 bo.py \
            $design_name or def\
            $is_cts_design \
            $max_iter \
            $max_dimension \
            $threshold \
            $max_num_cases \
            $max_num_candidates \
            $max_delta \
            $is_shift \
            $is_flip \
            $is_choi \
            $inference_type \
            $is_cadence \
            $lef_directory \
            $track_pitch \
            $output_dir \
            $blockage_layers")
        print(len(sys.argv))
        print(sys.argv)
        sys.exit()
    else:
        is_cts_from_place = False

        design_name = sys.argv[1]
        def_dir = None
        if sys.argv[1].endswith('.def'):
            def_dir = sys.argv[1]
            design_name= def_dir[def_dir.rfind('/')+1:def_dir.rfind('.def')]
        is_cts_design = int(sys.argv[2])
        max_iter = int(sys.argv[3])
        max_dimension = int(sys.argv[4])
        threshold = sys.argv[5]
        while len(threshold) != 6:
            threshold += '0'
        
        max_num_cases = int(sys.argv[6])
        max_num_candidates = int(sys.argv[7])
        max_delta = int(sys.argv[8])
        is_shift = int(sys.argv[9])
        is_flip = int(sys.argv[10])
        is_choi = int(sys.argv[11])
        inference_type = int(sys.argv[12])
        is_cadence = int(sys.argv[13])
        lef_dir = sys.argv[14]
        track_pitch = sys.argv[15]
        #output_dir = sys.argv[16]
        blockage_layers = sys.argv[16:]
        

        n_initial_points = 10

        binary_dir = PROJECT_DIR + '/src/bo_python/src/main'
        result_dir = PROJECT_DIR
        output_dir = PROJECT_DIR + '/outputs/def'

        
        if is_cts_from_place:
            post_fix = '_place_M' + str(max_iter) + '_CTS0' + '_DELTA' + str(max_delta) + '_D' + str(max_dimension) + '_P' + str(max_num_cases) + '_CAND' + str(max_num_candidates) +  '_TH' + str(int(float(threshold)*10000)) + '_SHIFT' + str(is_shift) + '_FLIP' + str(is_flip) + '_CHOI' + str(is_choi) + '_INF' + str(inference_type)
            post_fix += '_3stages_CTS'
            print(design_name + post_fix)
        else:
            post_fix = '_M' + str(max_iter) + '_CTS' + str(is_cts_design) + '_DELTA' + str(max_delta) + '_D' + str(max_dimension) + '_P' + str(max_num_cases) + '_CAND' + str(max_num_candidates) +  '_TH' + str(threshold) + '00_SHIFT' + str(is_shift) + '_FLIP' + str(is_flip) + '_CHOI' + str(is_choi) + '_INF' + str(inference_type)

        output_def_name = output_dir + '/' + design_name + post_fix + '.def'
        print(output_def_name)

        is_skip = False
        '''
        for item in os.listdir(output_dir):
            if (design_name + post_fix in item) and ('final' in item):
                is_skip = True
        ''' 

        if not is_skip:
            bo = bo_process(binary_dir, result_dir, design_name, post_fix, max_iter, n_initial_points, max_dimension,\
                    max_delta, is_cts_design, threshold, max_num_cases, max_num_candidates, is_shift, is_flip, is_choi,\
                    inference_type, is_cts_from_place, def_dir, is_cadence, lef_dir, track_pitch, blockage_layers)
            bo.run()
        else:
            print(output_def_name, ' exists! Skip!')
        
