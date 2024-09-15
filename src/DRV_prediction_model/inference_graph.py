import numpy as np
import os
import pickle
import argparse
from tqdm import tqdm
import torch
from torch_geometric.data import Data

from GNN_new import PinUGnet_new

#torch.cuda.set_device(1)
#print("GPU :", torch.cuda.current_device())
'''
import tensorflow as tf
conf = tf.compat.v1.ConfigProto()
conf.gpu_options.per_process_gpu_memory_fraction = 0.4
conf.gpu_options.allow_growth = True
os.environ['CUDA_VISIBLE_DEVICES']='1'
session = tf.compat.v1.Session(config=conf)
'''

PROJECT_DIR="."
npz_base = ""
csv_path = ""
csv_dir = ""
inference_type = ""
design = ""
is_unseen = None

def normalize(feature, normal, mode, type): #min/max, mean/std || type = 'pin_pattern' or 'others'
    epsilon = 1e-6
    if mode == 'minmax':
        if type == 'pin_pattern':
            feature_ = (feature - normal['min_p'].unsqueeze(1).unsqueeze(1)) / (normal['max_p'].unsqueeze(1).unsqueeze(1) - normal['min_p'].unsqueeze(1).unsqueeze(1) + epsilon)
        elif type == 'others':
            feature_ = (feature - normal['min_o'].unsqueeze(1).unsqueeze(1)) / (normal['max_o'].unsqueeze(1).unsqueeze(1) - normal['min_o'].unsqueeze(1).unsqueeze(1) + epsilon)
        else:
            print('Error in normalize!')
    elif mode == 'meanstd':
        if type == 'pin_pattern':
            feature_ = (feature - normal['mean_p'].unsqueeze(1).unsqueeze(1)) / (normal['std_p'].unsqueeze(1).unsqueeze(1) + epsilon)
        elif type == 'others':
            feature_ = (feature - normal['mean_o'].unsqueeze(1).unsqueeze(1)) / (normal['std_o'].unsqueeze(1).unsqueeze(1) + epsilon)
        else:
            print('Error in normalize!')

    else:
        print('Error in normalize!')

    return feature_

def inference(net, device, pickle_data, out_path):
    
    global inference_type
    global npz_path
    assert(inference_type != "")
    
    npz_data = design + '.npz'

    npz_name = npz_data.rstrip('.npz')
    if not is_unseen:
        npz_name_with_type = npz_name + '_' + inference_type
    else:
        npz_name_with_type = npz_name + '_unseen'
    npz_feature = np.load(npz_path)
    print('Inferencing starts. Read npz '+npz_path)
    
    feature_array = []
    for feature in pickle_data['gcell_features']:
        if not feature.startswith('drv_check') and not feature.startswith('graph'):
            feature_array.append(npz_feature[feature])
        
    input_feature = np.stack(feature_array)
    
    if input_feature.ndim != 4:
        input_feature = np.expand_dims(input_feature, axis=0)

    if input_feature.ndim != 4:
        input_feature = np.expand_dims(input_feature, axis=0)

    input_feature = torch.from_numpy(input_feature).type(torch.FloatTensor)

    input_feature = normalize(input_feature, pickle_data, 'meanstd', 'others')
    input_feature = input_feature.to(device=device, dtype=torch.float32)
    
    graph_feature = Data(x=torch.from_numpy(npz_feature['graph_x']).type(torch.FloatTensor), 
                 edge_index=torch.from_numpy(npz_feature["graph_a"]).type(torch.LongTensor).t().contiguous(), 
                 edge_attr=torch.from_numpy(npz_feature["graph_aw"]).type(torch.FloatTensor), 
                 graph_index=torch.from_numpy(npz_feature["graph_i"]).type(torch.LongTensor))


    graph_feature = graph_feature.to(device=device)

    if not os.path.isdir(out_path):
        os.makedirs(out_path)
        
    with torch.no_grad():
        pred = net(input_feature, graph_feature)
        if pickle_data['type'] == 'clas':
            logit = torch.sigmoid(pred)
        else:
            logit = pred

        logit_np = logit.cpu().detach().numpy()
        logit_np = np.squeeze(logit_np, axis=0)
        logit_np = np.squeeze(logit_np, axis=0)
        np.save(out_path + f'/{npz_name_with_type}_logit.npy', logit_np)

        if pickle_data['type'] == 'clas':
            pred_np = (logit_np >= pickle_data['threshold']) * 1
            np.save(out_path + f'/{npz_name_with_type}_pred.npy', pred_np)


    print(f'End inferencing!')

def csv_to_npz(gcell_features, graph_features):
    global npz_path
    assert(design != "")
    if os.path.isfile(npz_path):
        return
    np_data = {}
    for feature in gcell_features: 
        feature_path = csv_dir + '/' + design + '/' + feature + '.csv'
        if not os.path.isfile(feature_path):
            print(feature_path + " does not exists!")
            break
        np_data[feature] = np.loadtxt(feature_path, delimiter=',')

    for feature in graph_features:
        feature_path = csv_dir + '/' + design + '/' + feature + '.txt'
        if not os.path.isfile(feature_path):
            print(feature_path + " does not exists!")
            break
        np_data[feature] = np.loadtxt(feature_path)                

    if not os.path.isdir(npz_base):
        os.makedirs(npz_base)
    np.savez_compressed(npz_path, **np_data)
    print("npz saved to "+npz_path)

def get_args():
    parser = argparse.ArgumentParser(description='Inferencing: DEF to prediction',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-m', '--model', metavar='M', type=str, default='unet',
                        help='prediction model', dest='model')
    parser.add_argument('-c', '--checkpoint(.pth)', metavar='C', type=str,
                        help='model checkpoint file', dest='cp')
    parser.add_argument('-i', '--feature-info(.pkl)', metavar='I', type=str,
                        help='feature info & normalization file', dest='info')
    parser.add_argument('-f', '--feature-extractor-path', metavar='F', type=str, default=PROJECT_DIR + '/src/DRV_prediction_model/feature_extractor_c',
                        help='feature extractor path', dest='fe_path')
    parser.add_argument('-d', '--def-path', metavar='D', type=str, default='def',
                    help='def path', dest='def_path')
    parser.add_argument('-l', '--lef-path', metavar='L', type=str, default='./inputs/lef/nangate15nm.lef',
                    help='lef path', dest='lef_path')
    parser.add_argument('-o', '--out-path', metavar='O', type=str, default='out',
                    help='output path', dest='out_path')
    parser.add_argument('-is_cadence', '--is_cadence?', metavar='is_cadence', type=str, default=False,
                    help='is inference from cadence?', dest='is_cadence')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    info_path = args.info
    with open(info_path, 'rb') as f:
        pickle_data = pickle.load(f)
    
    if pickle_data['type'] == 'clas':
        inference_type = 'class'
    else:
        inference_type = 'reg_graph'
    
    if 'unseen' in args.cp:
        is_unseen = True
    else:
        is_unseen = False
    csv_path = args.out_path + "/csv"
    csv_dir = csv_path
    npz_base = args.out_path + "/npz"
    def_path = args.def_path
    out_path = args.out_path + "/predictions"
    drv_path = args.out_path + "/ground_truth"
    lef_path = args.lef_path
    is_cadence = bool(int(args.is_cadence))
    if is_cadence:
        print("Tool compatibility has set for Cadence.")

    if not is_cadence:
        extractor_path = os.path.join(args.fe_path, 'main')
    else:
        path_splitted = args.fe_path.split('/')
        cadence_extract_path = '/'.join(path_splitted[:-1])+'/innovus_'+path_splitted[-1]
        extractor_path = os.path.join(cadence_extract_path, 'main')

    design = def_path.split('/')[-1][:def_path.split('/')[-1].rfind('.')]
    npz_path = npz_base + '/' + design + '_'+ inference_type + '.npz'
    if not is_cadence:
        print(f'{extractor_path} -D {def_path} -O {csv_path} -V {drv_path} -G 1 -L {lef_path}')
        os.system(f'{extractor_path} -D {def_path} -O {csv_path} -V {drv_path} -G 1 -L {lef_path}')
    else:
        print(f'{extractor_path} -D {def_path} -O {csv_path} -V {drv_path} -L {lef_path}')
        os.system(f'{extractor_path} -D {def_path} -O {csv_path} -V {drv_path} -L {lef_path}')

    
    csv_to_npz(pickle_data['gcell_features'], pickle_data['graph_features'])

    device = torch.device('cpu')
    if not is_cadence:
        net = PinUGnet_new(11, 6, len(pickle_data['gcell_features']), 1)
    else:
        net = PinUGnet_new(8, 6, len(pickle_data['gcell_features']), 1)

    weight_path = os.path.join(os.getcwd(), args.cp)
    net.load_state_dict(torch.load(weight_path, map_location=device))
    net.to(device=device)
    net.eval()
    
    print (f'Complete net loading!')

    inference(net, device, pickle_data, out_path)


    
