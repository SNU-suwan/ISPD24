import numpy as np

class CongParser:
    #Warning: Only the congestions on the layers MINT1 and MINT2 are considered.
    
    
    def __init__(self, cong_dir, scale, gcell_height, die_lb, die_rt, drv_region_shape, capacity_map = None):
        die_rt -= gcell_height
        
        self.cong_dir = cong_dir
        self.scale = scale
        self.gcell_height = gcell_height
        self.die_lb = die_lb
        self.die_rt = die_rt
        self.drv_region_shape = drv_region_shape
        
        self.cong_map = None
        self.demand_map = None
        self.capacity_map = capacity_map
        self.parse()
        
    def parse(self):
        if self.capacity_map is None:
            is_update_capacity = True
        else:
            is_update_capacity = False
        
        self.cong_map = np.zeros(self.drv_region_shape)
        self.demand_map = np.zeros(self.drv_region_shape)
        if is_update_capacity:
            self.capacity_map = [[{} for x in range(self.drv_region_shape[1])] for y in range(self.drv_region_shape[0])]
        
        is_writing_cong = False
        with open(self.cong_dir, 'r') as f:
            line = f.readline()
            while line:
                info = line.split()
                if len(info) > 1:
                    if not is_writing_cong and info[0] == 'Layer' and info[1] == 'Boundary':
                        is_writing_cong = True
                    elif is_writing_cong:
                        layer = info[0]
                        
                        if self.cong_dir.count('MINT4') == 0:
                            if layer == 'MINT3' or layer == 'MINT4':
                                #print('Skip above MINT3')
                                break
                        
                        bbox = np.array([[float(info[1][2:]), float(info[2][:-1])], [float(info[3][1:]), float(info[4][:-2])]]) * self.scale
                        if (bbox[0][0] < self.die_lb[0] or bbox[0][1] < self.die_lb[1]) or (bbox[0][0] > self.die_rt[0] or bbox[0][1] > self.die_rt[1]):
                            line = f.readline()
                            continue
                        
                        demand = int(info[5][:info[5].find('/')])
                        capacity = int(info[5][info[5].find('/') + 1:])

                        grid_idx = [int((bbox[0][0] - self.die_lb[0])/self.gcell_height), int((bbox[0][1] - self.die_lb[1])/self.gcell_height)]
    
                        #self.cong_map[grid_idx[0]][grid_idx[1]] += demand - capacity
                    
                        #Only left & bottom edges are consitdered
                        
                        try:
                            if is_update_capacity:
                                self.cong_map[grid_idx[1]][grid_idx[0]] += demand - capacity
                            else:
                                self.cong_map[grid_idx[1]][grid_idx[0]] += demand - self.capacity_map[grid_idx[1]][grid_idx[0]][layer]
                            self.demand_map[grid_idx[1]][grid_idx[0]] += demand

                            if is_update_capacity:
                                self.capacity_map[grid_idx[1]][grid_idx[0]][layer] = capacity
                        except:
                            pass
                        
                        '''
                        test_idx = np.array([207,253])
                        
                        if grid_idx[0] == test_idx[0] and grid_idx[1] == test_idx[1]:
                            print(self.cong_map[grid_idx[1]][grid_idx[0]])
                            print(demand - capacity, layer)
                            print(bbox)
                            print()
                        ''' 
                        '''
                        if grid_idx[0] > 0 and (layer=='MINT1' or layer == 'MINT3'):
                            #self.cong_map[grid_idx[0]-1][grid_idx[1]] += demand - capacity
                            self.cong_map[grid_idx[1]][grid_idx[0]-1] += demand - capacity
                            self.demand_map[grid_idx[1]][grid_idx[0]-1] += demand
                            self.capacity_map[grid_idx[1]][grid_idx[0]-1] += capacity
                        
                            if grid_idx[0] - 1 == test_idx[0] and grid_idx[1] == test_idx[1]:
                                print(self.cong_map[grid_idx[1]][grid_idx[0]-1])
                                print(demand - capacity, layer)
                                print(bbox)
                                print()

                        if grid_idx[1] > 0 and (layer=='MINT2' or layer == 'MINT4'):
                            #self.cong_map[grid_idx[0]][grid_idx[1]-1] += demand - capacity
                            self.cong_map[grid_idx[1]-1][grid_idx[0]] += demand - capacity
                            self.demand_map[grid_idx[1]-1][grid_idx[0]-1] += demand
                            self.capacity_map[grid_idx[1]-1][grid_idx[0]-1] += capacity
                            if grid_idx[0] == test_idx[0] and grid_idx[1] - 1 == test_idx[1]:
                                print(self.cong_map[grid_idx[1]-1][grid_idx[0]])
                                print(demand - capacity, layer)
                                print(bbox)
                                print()
                        '''
                            

                line = f.readline()
        self.cong_map = self.cong_map / 2
        self.demand_map = self.demand_map / 2
        #self.capacity_map = self.capacity_map / 2

'''
class CongParser:
    #Warning: Only the congestions on the layers MINT1 and MINT2 are considered.
    
    
    def __init__(self, cong_dir, scale, gcell_height, die_lb, die_rt, drv_region_shape):
        die_rt -= gcell_height
        
        self.cong_dir = cong_dir
        self.scale = scale
        self.gcell_height = gcell_height
        self.die_lb = die_lb
        self.die_rt = die_rt
        self.drv_region_shape = drv_region_shape
        
        self.cong_map = None
        self.parse()
        
    def parse(self):
        self.cong_map = np.zeros(self.drv_region_shape)
        
        is_writing_cong = False
        with open(self.cong_dir, 'r') as f:
            line = f.readline()
            while line:
                info = line.split()
                if len(info) > 1:
                    if not is_writing_cong and info[0] == 'Layer' and info[1] == 'Boundary':
                        is_writing_cong = True
                    elif is_writing_cong:
                        layer = info[0]
                        if layer == 'MINT3' or layer == 'MINT4':
                            #line = f.readline()
                            #continue
                            break

                        bbox = np.array([[float(info[1][2:]), float(info[2][:-1])], [float(info[3][1:]), float(info[4][:-2])]]) * self.scale
                        if (bbox[0][0] < self.die_lb[0] or bbox[0][1] < self.die_lb[1]) or (bbox[0][0] >= self.die_rt[0] or bbox[0][1] >= self.die_rt[1]):
                            line = f.readline()
                            continue
                        
                        demand = int(info[5][:info[5].find('/')])
                        capacity = int(info[5][info[5].find('/') + 1:])

                        grid_idx = [int((bbox[0][0] - self.die_lb[0])/self.gcell_height), int((bbox[0][1] - self.die_lb[1])/self.gcell_height)]

                        self.cong_map[grid_idx[0]][grid_idx[1]] += demand - capacity
                        
                        if grid_idx[0] > 0 and layer=='MINT1':
                            self.cong_map[grid_idx[0]-1][grid_idx[1]] += demand - capacity
                        if grid_idx[1] > 0 and layer=='MINT2':
                            self.cong_map[grid_idx[0]][grid_idx[1]-1] += demand - capacity
                            

                line = f.readline()
        self.cong_map = self.cong_map / 4
'''