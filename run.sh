#!/bin/bash

benchmarks=(b18) # Make sure that input .def files are located in ./inputs/def
output_dir="./outputs"

#########################################

log_dir=$output_dir/log
if [[ ! -d $log_dir ]]; then
    mkdir -p $log_dir
fi
if [[ ! -d $output_dir/def ]]; then
    mkdir -p $output_dir/def
fi
if [[ ! -d $output_dir/scripts ]]; then
    mkdir -p $output_dir/scripts
fi
if [[ ! -d $output_dir/csv ]]; then
    mkdir -p $output_dir/scripts
fi
if [[ ! -d $output_dir/npz ]]; then
    mkdir -p $output_dir/scripts
fi


#source /etc/profile.d/conda.sh
#conda activate pytorch_1.8.1 #Env that contains pytorch_1.8.1 package
#conda activate pytorch

IS_SYNOPSYS=true # true: Nangate 15nm + Synopsys / false: ASAP 7nm + Innovus
if [[ $IS_SYNOPSYS == "true" ]]; then
    lef_dir=./inputs/lef/nangate15nm.lef
    is_cadence=0
    track_pitch=0.064
    blockage_layers="M1 V1 MINT1 VINT1 MINT2 VINT2 MINT3 VINT3 MINT4 VINT4 MINT5 VINT5"
else
    lef_dir=./inputs/lef/asap7.lef
    is_cadence=1
    track_pitch=0.144
    blockage_layers="M1 V1 M2 V2 M3 V3 M4 V4 M5 V5 M6 V6 M7 V7 M8 V8 M9 V9 Pad"
fi

threshold=0.4984
max_iter=40
max_dim=5
max_range=100
max_delta=5
max_num_candidates=100
max_num_cases=$max_num_candidates
inference_type=1
prediction_postfix=M${max_iter}_CTS1_DELTA${max_delta}_D${max_dim}_P${max_num_cases}_CAND${max_num_candidates}_TH${threshold}00_SHIFT1_FLIP0_CHOI1_INF${inference_type}

benchmark=b18



for benchmark in ${benchmarks[*]}
do  
    start_time=`date +%s`
    def_dir=./inputs/def/${benchmark}.def
    
    python3 -u ./src/bo_python/bo.py $def_dir 1 $max_iter $max_dim $threshold $max_num_cases $max_num_candidates $max_delta 1 0 1 $inference_type 0 $lef_dir $track_pitch $output_dir $blockage_layers | tee $log_dir/${benchmark}_ISPD24.log

    echo "Routing blockage file location: ./outputs/scripts/blockages_${benchmark}_M${max_iter}_CTS1_DELTA${max_delta}_D${max_dim}_P${max_range}_CAND${max_num_candidates}_TH${threshold}00_SHIFT1_FLIP0_CHOI1_INF${inference_type}.tcl"

    end_time=`date +%s`
    echo Total execution time took `expr $end_time - $start_time` s. | tee -a ${log_dir}/log/runtime_${benchmark}.log

done

#conda deactivate
