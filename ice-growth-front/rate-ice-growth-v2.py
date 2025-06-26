import numpy as np
import mdtraj as md
from tqdm import tqdm  # 进度条库

'''
配合z-density-stat-partial.py使用
'''

# 设置常量和输入输出文件路径
trr_file = 'npt-200.trr'  # 轨迹文件路径
top_file = 'npt-200.gro'  # gro文件路径
output_file = '1-sigma.txt'  # 输出文件路径
num_layers = 200  # 分层数
box_height = 8.90420  # 模拟盒子总高度（单位nm）从gro文件中获取
std_dev = 34.0878   # 标准差,从z-density-stat-partial.py的运行结果中获取

# 读取轨迹文件
traj = md.load(trr_file, top=top_file)
time_steps = traj.time  # 轨迹时间
z_coords = traj.xyz[:, :, 2]  # 获取每个原子的z坐标
num_atoms = traj.n_atoms

# 计算每层的粒子数
def calculate_layer_density(z_coords, num_layers, box_height):
    layer_density = np.zeros(num_layers)
    layer_z = np.linspace(0, box_height, num_layers + 1)  # 计算层的z坐标范围
    for i in range(len(z_coords)):  # 修改为 len(z_coords)
        atom_z = z_coords[i]
        layer_index = int(atom_z / (box_height / num_layers))  # 根据z坐标分配层号
        # 确保 layer_index 不超过最大层数
        layer_index = min(layer_index, num_layers - 1)
        layer_density[layer_index] += 1  # 增加对应层的粒子数
    return layer_density, layer_z

# 判断冰层
def define_ice_layers(layer_density, num_layers):
    avg_density = np.mean(layer_density)  # 计算平均粒子数
    ice_layers = np.zeros(num_layers, dtype=bool)
    for i in range(num_layers):
        #if layer_density[i] > 1.5 * avg_density or layer_density[i] < 0.5 * avg_density:    # 0.5N~1.5N, 是暂时的推测性判据；
        if (layer_density[i] > (avg_density + np.sqrt(2) * std_dev)) or (layer_density[i] < (avg_density - np.sqrt(2) * std_dev)):
        #if (layer_density[i] > (avg_density + 1 * std_dev)) or (layer_density[i] < (avg_density - 1 * std_dev)):
            ice_layers[i] = True
    return ice_layers

# 寻找下冰层（从0到中值）
def find_lower_ice_layer(ice_layers, M):
    for i in range(M):
        if ice_layers[i]:
            return i
    return -1  # 如果没有找到，返回-1

# 寻找上冰层（从最大层减20到中值）
def find_upper_ice_layer(ice_layers, M, num_layers, layer_density):
    avg_density = np.mean(layer_density)  # 计算平均粒子数
    start_layer = max(num_layers - 20, 0)  # 起始层号为最大层-20，确保不小于0
    for i in range(start_layer, M - 1, -1):  # 从 start_layer 开始到 M 层逐层递减
        # 判断当前层是否满足冰层条件（粒子数 > 1.5 倍平均值或 < 0.5 倍平均值）
        if layer_density[i] > 1.5 * avg_density or layer_density[i] < 0.5 * avg_density:
            return i  # 返回满足条件的层号
    return M  # 如果没有找到，返回中值层号

# 处理每个时间步，输出符合要求的层序号和z坐标
def process_trajectory(traj, num_layers, box_height, output_file):
    with open(output_file, 'w') as f:
        # 写入文件头
        f.write("Time (ns)\tLower Ice Layer\tUpper Ice Layer\tLower Z Coordinate (nm)\tUpper Z Coordinate (nm)\n")

        # 使用 tqdm 添加进度条
        for step in tqdm(range(len(traj)), desc="Processing trajectory", unit="step"):
            # 获取当前时间步的z坐标
            z_coords_step = z_coords[step]
            
            # 计算每层的粒子数
            layer_density, layer_z = calculate_layer_density(z_coords_step, num_layers, box_height)
            
            # 判断哪些层是冰层
            ice_layers = define_ice_layers(layer_density, num_layers)
            
            # 计算几何中心的层号M
            M = num_layers // 2
            
            # 分别寻找下冰层和上冰层
            lower_ice_layer = find_lower_ice_layer(ice_layers, M)
            upper_ice_layer = find_upper_ice_layer(ice_layers, M, num_layers, layer_density)

            # 如果没有找到满足条件的冰层，输出0
            if lower_ice_layer == -1:
                lower_ice_layer = 0
            if upper_ice_layer == -1:
                upper_ice_layer = num_layers - 1
            
            # 计算z坐标并输出
            lower_z = layer_z[lower_ice_layer] if lower_ice_layer != 0 else 0
            upper_z = layer_z[upper_ice_layer]
            
            # 记录时间（单位ns）
            time_ns = traj.time[step] / 1000.0
            
            # 写入输出文件
            f.write(f"{time_ns:.3f}\t{lower_ice_layer + 1}\t{upper_ice_layer + 1}\t{lower_z:.3f}\t{upper_z:.3f}\n")

# 执行处理过程
process_trajectory(traj, num_layers, box_height, output_file)

print(f"结果已输出到 {output_file}")
