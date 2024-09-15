import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.modules.activation import ReLU
from torch.nn.modules.dropout import Dropout
from torch_scatter import scatter_mean, scatter_add
from torch_geometric.nn import GCNConv, SAGEConv, GraphConv, MessagePassing, HeteroConv, GATConv, Linear, to_hetero, global_add_pool, global_mean_pool, global_max_pool

class Conv3x3GNReLU(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        num_groups = 32
        if out_channels % num_groups != 0 or out_channels == num_groups:
            num_groups = out_channels // 2

        self.blocks = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, (3, 3), stride=1, padding=1, bias=False),
            nn.GroupNorm(num_groups, out_channels),
            nn.ReLU(inplace=True),
        )
    
    def forward(self, x):
        return self.blocks(x)

class DecoderBlock(nn.Module):
    def __init__(self, in_channels, out_channels, up_sample=True):
        super().__init__()
        self.up_sample = up_sample
        self.block = nn.Sequential(
            Conv3x3GNReLU(in_channels, out_channels),
            Conv3x3GNReLU(out_channels, out_channels),
        )

    def forward(self, x):
        if isinstance(x, list) or isinstance(x, tuple):
            x, skip = x
        else:
            skip = None

        if self.up_sample:
            x = F.interpolate(x, scale_factor=2, mode='nearest')
            
        if skip is not None:
            diffY = skip.size()[2] - x.size()[2]
            diffX = skip.size()[3] - x.size()[3]
            x = F.pad(x, [diffX // 2, diffX - diffX // 2,
                          diffY // 2, diffY - diffY // 2])

            x = torch.cat([x, skip], dim=1)

        return self.block(x)

class ResBlock(nn.Module):
    def __init__(self, in_channels, out_channels, downsample=True):
        super().__init__()
        num_groups = 32
        if out_channels % num_groups != 0 or out_channels == num_groups:
            num_groups = out_channels // 2

        if downsample:
            self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=(3, 3), stride=(2, 2), padding=(1, 1), bias=False)
        else:
            self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=(3, 3), stride=(1, 1), padding=(1, 1), bias=False)

        self.gn1 = nn.GroupNorm(num_groups, out_channels)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=(3, 3), stride=(1, 1), padding=(1, 1), bias=False)
        self.gn2 = nn.GroupNorm(num_groups, out_channels)

        if downsample:
            self.down = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, kernel_size=(1, 1), stride=(2, 2), bias=False),
                nn.GroupNorm(num_groups, out_channels),
            )

        self.downsample = downsample

    def forward(self, x):
        identity = x
        if self.downsample:
            identity = self.down(x)

        out = self.conv1(x)
        out = self.gn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.gn2(out)

        out += identity
        out = self.relu(out)

        return out

        
class ResBlockX2(nn.Module):
    def __init__(self, in_channels, out_channels, downsample=True):
        super().__init__()    
        self.block1 = ResBlock(in_channels, out_channels, downsample=downsample)
        self.block2 = ResBlock(out_channels, out_channels, downsample=False)

    def forward(self, x):
        x = self.block1(x)
        x = self.block2(x)
        return x

# in_channels = [16, 32, 64, 128, 256, 512]
class EncoderX(nn.Module):
    def __init__(self, g_channels, in_channels=[32, 64, 128, 256, 512]):
        super().__init__()
        self.conv1 = nn.Conv2d(g_channels, 16, kernel_size=3, stride=1, padding=1, bias=False)
        self.gn1 = nn.GroupNorm(8, 16)
        self.relu = nn.ReLU(inplace=True)

        self.en0 = ResBlockX2(16, in_channels[0])
        
        self.en1 = ResBlockX2(in_channels[0], in_channels[1])

        self.en2 = ResBlockX2(in_channels[1], in_channels[2])
        self.en3 = ResBlockX2(in_channels[2], in_channels[3])
        self.en4 = ResBlockX2(in_channels[3], in_channels[4])

    def forward(self, x):
        x = self.conv1(x)
        x = self.gn1(x)
        x = self.relu(x)

        x0 = self.en0(x)
        x1 = self.en1(x0)
        
        x2 = self.en2(x1)
        x3 = self.en3(x2)
        x4 = self.en4(x3)

        return [x4, x3, x2, x1, x0]


class DecoderXNoF(nn.Module):
    def __init__(self, e_channels = [512, 256, 128, 64, 32], out_channels = [256, 128, 64, 32, 32]):
        super().__init__()
        self.center = DecoderBlock(in_channels=e_channels[0], out_channels=e_channels[0], up_sample=False)
        self.layer1 = DecoderBlock(e_channels[0] + e_channels[1], out_channels[0])
        self.layer2 = DecoderBlock(e_channels[2] + out_channels[0], out_channels[1])
        self.layer3 = DecoderBlock(e_channels[3] + out_channels[1], out_channels[2])
        self.layer4 = DecoderBlock(e_channels[4] + out_channels[2], out_channels[3])
        self.layer5 = DecoderBlock(out_channels[3], out_channels[4])
        '''
        self.final = nn.Sequential(
            nn.Dropout(dropout),
            nn.Conv2d(out_channels[-1], n_classes, kernel_size=1)
        )
        '''

    def forward(self, x, encodes):
        skips = encodes[1:]

        center = self.center(encodes[0])

        decodes = self.layer1([center, skips[0]])
        decodes = self.layer2([decodes, skips[1]])
        decodes = self.layer3([decodes, skips[2]])
        decodes = self.layer4([decodes, skips[3]])
        decodes = self.layer5([decodes, None])

        decodes4 = F.interpolate(decodes, x.size()[2:], mode='bilinear', align_corners=False)

        return decodes4

class PinGConv(MessagePassing):
    def __init__(self, in_channels, edge_channels, out_channels, neighbor_channels):
        super().__init__(aggr='add')

        self.in_channels = in_channels
        self.edge_channels = edge_channels
        self.out_channels = out_channels

        self.eff_mlp = Linear(2 * in_channels + edge_channels, neighbor_channels)
        self.out_mlp = Linear(in_channels + neighbor_channels, out_channels)

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index, x=x, edge_attr=edge_attr)
        x_out = torch.cat([x, out], dim=1)
        return self.out_mlp(x_out)

    def message(self, x_i, x_j, edge_attr):
        tmp = torch.cat([x_i, x_j, edge_attr], dim=1)
        return self.eff_mlp(tmp)


class PinGNN(nn.Module):
    def __init__(self, c_in, c_hidden, c_out, c_edge, c_neighbor, num_layers=2, dp_rate=0.1):
        super().__init__()

        layers = []
        in_channels, out_channels = c_in, c_hidden
        for l_idx in range(num_layers - 1):
            layers += [
                PinGConv(in_channels=in_channels, edge_channels=c_edge, out_channels=out_channels, neighbor_channels=c_neighbor),
                nn.ReLU(inplace=True)
            ]
            in_channels = c_hidden
        layers += [PinGConv(in_channels=in_channels, edge_channels=c_edge, out_channels=c_out, neighbor_channels=c_neighbor)]
        self.layers = nn.ModuleList(layers)


    def forward(self, graph_data):
        x, edge_index, edge_attr, graph_index = graph_data.x, graph_data.edge_index, graph_data.edge_attr, graph_data.graph_index

        for l in self.layers:
            if isinstance(l, MessagePassing):
                x = l(x, edge_index, edge_attr)
            else:
                x = l(x)

        x = global_mean_pool(x, graph_index)

        return x

class PinGNNWhole_new(nn.Module):
    def __init__(self, n_node_features, n_edge_features, n_classes, dropout=0.2):
        super().__init__()
        self.n_node_featuresraph_ = n_node_features
        self.n_edge_features = n_edge_features
        self.n_classes = n_classes
        self.graph_conv = PinGNN(n_node_features, 32, 32, n_edge_features, 16, num_layers=3)   

        self.out_conv = nn.Sequential(
            nn.Dropout(dropout),
            nn.Conv2d(32, 16, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Conv2d(16, 1, kernel_size=1)
        )

    def forward(self, x, graph_x):
        graph_feature = self.graph_conv(graph_x)
        graph_feature_re = graph_feature.T.view(-1, graph_feature.shape[1], x.shape[2], x.shape[3])

        out = self.out_conv(graph_feature_re)

        return out

class PinGUnet_new(nn.Module):
    def __init__(self, n_node_features, n_edge_features, g_channels, n_classes):
        super().__init__()
        self.n_node_features = n_node_features
        self.n_edge_features = n_edge_features
        self.n_gc_channels = g_channels
        self.n_classes = n_classes
        self.graph_conv = PinGNN(n_node_features, 32, 32, n_edge_features, 16, num_layers=3)

        self.g_scale_conv = nn.Conv2d(g_channels, 32, kernel_size=1, stride=1)

        self.scale_conv = nn.Conv2d(64, 32, kernel_size=1, stride=1)

        self.encoder = EncoderX(32)
        self.decoder = DecoderXNoF(n_classes)

    def forward(self, x, graph_x):

        graph_feature = self.graph_conv(graph_x)
        graph_feature_re = graph_feature.T.view(-1, graph_feature.shape[1], x.shape[2], x.shape[3])

        x = self.g_scale_conv(x)

        in_x = torch.cat([x, graph_feature_re], dim=1)
        in_x = self.scale_conv(in_x)

        encodes = self.encoder(in_x)
        return self.decoder(in_x, encodes)

class PinUGnet_new(nn.Module):
    def __init__(self, n_node_features, n_edge_features, g_channels, n_classes, dropout=0.2):
        super().__init__()
        self.n_node_features = n_node_features
        self.n_edge_features = n_edge_features        
        self.n_gc_channels = g_channels
        self.n_classes = n_classes
        self.graph_conv = PinGNN(n_node_features, 32, 32, n_edge_features, 16, num_layers=3)

        self.gconv = nn.Conv2d(32, 32, kernel_size=3, stride=1, padding=1) #added

        self.encoder = EncoderX(g_channels)
        self.decoder = DecoderXNoF()

        self.final = nn.Sequential(
            nn.Dropout(dropout),
            nn.Conv2d(64, 32, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Conv2d(32, 16, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Conv2d(16, n_classes, kernel_size=1)
        )

    def forward(self, x, graph_x):

        graph_feature = self.graph_conv(graph_x)
        graph_feature_re = graph_feature.T.view(-1, graph_feature.shape[1], x.shape[2], x.shape[3])
        graph_feature_re = self.gconv(graph_feature_re) # added

        encodes = self.encoder(x)
        u_out = self.decoder(x, encodes)

        fe_x = torch.cat([u_out, graph_feature_re], dim=1)
        out = self.final(fe_x)
        return out

class Unet_new(nn.Module):
    def __init__(self, g_channels, n_classes, dropout=0.2):
        super().__init__()
        self.n_gc_channels = g_channels
        self.n_classes = n_classes
        self.encoder = EncoderX(g_channels)
        self.decoder = DecoderXNoF()
        self.final = nn.Sequential(
            nn.Dropout(dropout),
            nn.Conv2d(32, 16, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Conv2d(16, 1, kernel_size=1)        
        )

    def forward(self, x):
        encodes = self.encoder(x)
        decode = self.decoder(x, encodes)
        return self.final(decode)

class PinGUnet_new(nn.Module):
    def __init__(self, n_node_features, n_edge_features, g_channels, n_classes):
        super().__init__()
        self.n_node_features = n_node_features
        self.n_edge_features = n_edge_features
        self.n_gc_channels = g_channels
        self.n_classes = n_classes
        self.graph_conv = PinGNN(n_node_features, 32, 32, n_edge_features, 16, num_layers=3)

        self.g_scale_conv = nn.Conv2d(g_channels, 32, kernel_size=1, stride=1)

        self.scale_conv = nn.Conv2d(64, 32, kernel_size=1, stride=1)

        self.encoder = EncoderX(32)
        self.decoder = DecoderXNoF(n_classes)

    def forward(self, x, graph_x):

        graph_feature = self.graph_conv(graph_x)
        graph_feature_re = graph_feature.T.view(-1, graph_feature.shape[1], x.shape[2], x.shape[3])

        x = self.g_scale_conv(x)

        in_x = torch.cat([x, graph_feature_re], dim=1)
        in_x = self.scale_conv(in_x)

        encodes = self.encoder(in_x)
        return self.decoder(in_x, encodes)

class PinUGnet_new_withnum(nn.Module):
    def __init__(self, n_node_features, n_edge_features, g_channels, n_classes, dropout=0.2):
        super().__init__()
        self.n_node_features = n_node_features
        self.n_edge_features = n_edge_features        
        self.n_gc_channels = g_channels
        self.n_classes = n_classes
        self.graph_conv = PinGNN(n_node_features, 32, 32, n_edge_features, 16, num_layers=3)

        self.gconv = nn.Conv2d(32, 32, kernel_size=3, stride=1, padding=1) #added

        self.encoder = EncoderX(g_channels)
        self.decoder = DecoderXNoF()

        self.final = nn.Sequential(
            nn.Dropout(dropout),
            nn.Conv2d(64, 32, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Conv2d(32, 16, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Conv2d(16, n_classes, kernel_size=1)
        )
        self.avgpool = nn.AdaptiveAvgPool2d((1, 1))
        self.fc = nn.Linear(512, 1)


    def forward(self, x, graph_x):

        graph_feature = self.graph_conv(graph_x)
        graph_feature_re = graph_feature.T.view(-1, graph_feature.shape[1], x.shape[2], x.shape[3])
        graph_feature_re = self.gconv(graph_feature_re) # added

        encodes = self.encoder(x)
        u_out = self.decoder(x, encodes)

        fe_x = torch.cat([u_out, graph_feature_re], dim=1)
        out = self.final(fe_x)

        num_drv = self.avgpool(encodes[0])
        num_drv = torch.flatten(num_drv, 1)
        num_drv = self.fc(num_drv)

        return out, num_drv

class PinGUnet_num(nn.Module):
    def __init__(self, n_node_features, n_edge_features, g_channels, n_classes):
        super().__init__()
        self.n_node_features = n_node_features
        self.n_edge_features = n_edge_features
        self.n_gc_channels = g_channels
        self.n_classes = n_classes
        self.graph_conv = PinGNN(n_node_features, 32, 32, n_edge_features, 16, num_layers=3)

        self.g_scale_conv = nn.Conv2d(g_channels, 32, kernel_size=1, stride=1)

        self.scale_conv = nn.Conv2d(64, 32, kernel_size=1, stride=1)

        self.encoder = EncoderX(32)

        self.avgpool = nn.AdaptiveAvgPool2d((1, 1))
        self.final = nn.Sequential(
            nn.Dropout(0.2),
            nn.Linear(512, 512),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(512, 1)
        )


    def forward(self, x, graph_x):
        graph_feature = self.graph_conv(graph_x)
        graph_feature_re = graph_feature.T.view(-1, graph_feature.shape[1], x.shape[2], x.shape[3])

        x = self.g_scale_conv(x)

        in_x = torch.cat([x, graph_feature_re], dim=1)
        in_x = self.scale_conv(in_x)

        encodes = self.encoder(in_x)
        avg_pooled = self.avgpool(encodes[0])
        flattened = torch.flatten(avg_pooled, 1)
        output = self.final(flattened)

        return output