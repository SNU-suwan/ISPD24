import sys

class Row:
    def __init__(self, row_name, x, y):
        self.name = row_name
        self.x = x
        self.y = y
        
class Comp:
    def __init__(self, comp_name, cell, x, y, dir_):
        self.name = comp_name
        self.cell = cell
        self.x = x
        self.y = y
        self.dir = dir_

class Whitespace:
    def __init__(self, row, x_left, x_right):
        self.row = row
        self.x_left = x_left
        self.x_right = x_right
        
class DefParser:
    def __init__(self, def_dir, is_cadence=False):
        self.def_dir = def_dir
        self.scale = None
        self.cpp = None
        self.die_size = []
        self.offset = []
        self.whitespace_dict = {}
        self.comp_dict = {}
        self.row_comp_dict = {}
        self.rows = []
        self.is_cadence = is_cadence

    def parse(self, lef):
        with open(self.def_dir, 'r') as f:
            line = f.readline()
            is_component = False
            while line:
                info = line.split()
                if len(info) > 0:
                    if info[0].startswith('ROW'):
                        new_row = Row(info[1], int(info[3]), int(info[4]))
                        if not new_row.y in self.row_comp_dict:
                            self.row_comp_dict[new_row.y] = []
                            self.whitespace_dict[new_row.y] = []
                        self.rows.append(new_row)
                        if not self.is_cadence:
                            self.cpp = int(info[-3])
                        else:
                            self.cpp = int(info[-2])
                    elif info[0].startswith('DIEAREA'):
                        if not self.is_cadence:
                            self.die_size = [int(info[10]), int(info[11])]
                        else:
                            self.die_size = [int(info[6]), int(info[7])]
                    elif info[0].startswith('COMPONENTS'):
                        is_component = True
                    elif info[0] == '-':
                        if is_component:
                            #new_comp = Comp(info[1], lef.cell_dict[info[2]], int(info[6]), int(info[7]), info[9])
                            new_comp = Comp(info[1], lef.cell_dict[info[2]], int(info[info.index('(')+1]), int(info[info.index('(')+2]), info[info.index('(')+4])
                            self.comp_dict[new_comp.name] = new_comp
                            self.row_comp_dict[new_comp.y].append(new_comp)
                    elif info[0] == 'END':
                        if info[1] == 'COMPONENTS':
                            is_component = False
                line = f.readline()
                
        for row, comps in self.row_comp_dict.items():
            self.row_comp_dict[row] = sorted(comps, key = lambda x: x.x)
        self.rows = sorted(self.rows, key= lambda x:x.x)
        self.offset = [self.rows[0].x, self.rows[0].y]
        
    def get_whitespace(self):
        offset_x = self.offset[0]
        for row, comps in self.row_comp_dict.items():
            previous_x_right = offset_x
            for comp in comps:
                if comp.x - previous_x_right > 0:
                    w_left = previous_x_right
                    w_right = comp.x
                    new_whitespace = Whitespace(row, w_left, w_right)
                    self.whitespace_dict[row].append(new_whitespace)
                previous_x_right = comp.x + comp.cell.width
            
            if previous_x_right != self.die_size[1] - offset_x:
                w_left = previous_x_right
                w_right = self.die_size[1] - offset_x
                new_whitespace = Whitespace(row, w_left, w_right)
                self.whitespace_dict[row].append(new_whitespace)
        
        for row, wss in self.whitespace_dict.items():
            self.whitespace_dict[row] = sorted(wss, key=lambda x:x.x_left)
    
    def get_scale(self):
        with open(self.def_dir, 'r') as f:
            line = f.readline()
            while line:
                info = line.split()
                if info[0].startswith('UNITS'):
                    self.scale = int(info[-2])
                    break
                line = f.readline()
        return self.scale