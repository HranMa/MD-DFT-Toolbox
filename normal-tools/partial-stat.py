'''
单次使用该脚本的时候，注意输出后检查。
如果不小心重复执行，需要到输出文件里删去不需要的内容。
目前脚本的逻辑是：创建txt，写upper-water，继续写ice-zone，继续写lower-water
所以如果重复执行，不会擦除原本已经存在的内容。
'''

import numpy as np

def read_particle_distribution(file_path):
    """读取粒子分布文件 (npt-400.txt)，返回每层的粒子数"""
    layers = []
    particle_counts = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # 跳过表头
        
        for line in lines:
            parts = line.split()
            if len(parts) >= 3:
                layers.append(int(parts[0]))  # 层号
                try:
                    # 尝试将粒子数转换为浮动数
                    particle_count = float(parts[2])
                    # 确保粒子数是有效的
                    if particle_count >= 0:
                        particle_counts.append(particle_count)
                    else:
                        print(f"Warning: Invalid particle count {particle_count} at layer {parts[0]}, skipping.")
                except ValueError:
                    print(f"Warning: Invalid data for particle count at layer {parts[0]}, skipping.")
    
    return np.array(layers), np.array(particle_counts)

def read_layer_info(stat_file):
    """读取统计结果文件 (stat-result-total.txt)，获取分层信息"""
    layer_info = {}
    
    with open(stat_file, 'r') as f:
        lines = f.readlines()
        
        for line in lines:
            # 检查并提取每个分区的范围
            if "Lower Water Layers" in line:
                lower_water_range = line.split(":")[-1].strip()  # 获取'1,66'部分
                lower_water_range = list(map(int, lower_water_range.split(',')))  # 提取并转换为整数
                layer_info['lower_water'] = lower_water_range
            elif "Ice Layers" in line:
                ice_range = line.split(":")[-1].strip()  # 获取'67,142'部分
                ice_range = list(map(int, ice_range.split(',')))  # 提取并转换为整数
                layer_info['ice'] = ice_range
            elif "Upper Water Layers" in line:
                upper_water_range = line.split(":")[-1].strip()  # 获取'143,200'部分
                upper_water_range = list(map(int, upper_water_range.split(',')))  # 提取并转换为整数
                layer_info['upper_water'] = upper_water_range
    
    return layer_info

def calculate_statistics(particle_counts):
    """计算均值、方差和标准差"""
    mean = np.mean(particle_counts)
    variance = np.var(particle_counts)
    std_dev = np.std(particle_counts)
    
    return mean, variance, std_dev

def output_statistics(mean, variance, std_dev, layer_name, output_file):
    """将统计数据输出到文件"""
    with open(output_file, 'a') as f:
        f.write(f"{layer_name} Mean: {mean:.4f}\n")
        f.write(f"{layer_name} Variance: {variance:.4f}\n")
        f.write(f"{layer_name} Standard Deviation: {std_dev:.4f}\n\n")

def main():
    """主函数"""
    
    # 输入文件
    npt_file = 'npt-400.txt'  # 从 z-density-stat-partial.py 输出的文件
    stat_file = 'stat-result-total.txt'  # 从 z-density-stat-partial.py 输出的文件
    
    # 输出文件
    output_file = 'partial-stat.txt'  # 输出的统计文件
    
    # 读取粒子分布数据
    layers, particle_counts = read_particle_distribution(npt_file)
    
    # 读取分层信息
    layer_info = read_layer_info(stat_file)
    
    # 获取上下水区和冰区的粒子数
    lower_water_layers = layer_info.get('lower_water', [])
    ice_layers = layer_info.get('ice', [])
    upper_water_layers = layer_info.get('upper_water', [])
    
    # 修正索引，确保不超出范围
    if lower_water_layers and lower_water_layers[1] > len(particle_counts):
        lower_water_layers[1] = len(particle_counts)
    if ice_layers and ice_layers[1] > len(particle_counts):
        ice_layers[1] = len(particle_counts)
    if upper_water_layers and upper_water_layers[1] > len(particle_counts):
        upper_water_layers[1] = len(particle_counts)
    
    # 分别提取粒子数
    lower_water_counts = particle_counts[lower_water_layers[0]-1 : lower_water_layers[1]] if lower_water_layers else []
    ice_counts = particle_counts[ice_layers[0]-1 : ice_layers[1]] if ice_layers else []
    upper_water_counts = particle_counts[upper_water_layers[0]-1 : upper_water_layers[1]] if upper_water_layers else []
    
    # 计算各区的统计数据
    if lower_water_counts.size > 0:
        lower_water_mean, lower_water_variance, lower_water_std_dev = calculate_statistics(lower_water_counts)
        output_statistics(lower_water_mean, lower_water_variance, lower_water_std_dev, 'Lower Water', output_file)
    
    if ice_counts.size > 0:
        ice_mean, ice_variance, ice_std_dev = calculate_statistics(ice_counts)
        output_statistics(ice_mean, ice_variance, ice_std_dev, 'Ice', output_file)
    
    if upper_water_counts.size > 0:
        upper_water_mean, upper_water_variance, upper_water_std_dev = calculate_statistics(upper_water_counts)
        output_statistics(upper_water_mean, upper_water_variance, upper_water_std_dev, 'Upper Water', output_file)
    
    # 最后输出简洁的消息
    print(f"Partial statistics saved to '{output_file}'")

if __name__ == "__main__":
    main()
