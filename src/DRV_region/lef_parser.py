class Cell:
    def __init__(self, name, width, height):
        self.name = name
        self.width = width
        self.height = height
        
class LefParser:
    def __init__(self, lef_dir, scale):
        self.lef_dir = lef_dir
        self.scale = scale
        self.cell_dict = {}
        self.parse()
        
    def parse(self):
        with open(self.lef_dir, 'r') as f:
            line = f.readline()
            while line:
                info = line.split()
                if len(info)>0:
                    if info[0] == "MACRO":
                        name = info[1]
                    elif info[0] == "SIZE":
                        width = int(float(info[1]) * self.scale)
                        height = int(float(info[3]) * self.scale)

                        new_cell = Cell(name, width, height)
                        self.cell_dict[name] = new_cell

                line = f.readline()