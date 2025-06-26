import numpy as np

'''
配合rate-ice-growth-v2.py使用
用来获取冰晶生长前沿
'''

def read_gro_file(gro_file):
    """读取GRO文件，返回原子坐标和盒子尺寸，以及fSOL残基的z坐标最小值和最大值"""
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    coordinates = []
    fSOL_z_coordinates = []  # 用于存储fSOL残基的z坐标
    
    # 读取盒子尺寸
    box_size = np.array([float(val) for val in lines[-1].split()])
    
    # 读取原子坐标
    for line in lines[2:]:  # 跳过前两行
        if line.strip():
            try:
                # 尝试提取坐标值 根据gro文件列宽定义
                res_name = line[5:10].strip()  # 提取残基名称
                x = float(line[20:28].strip())
                y = float(line[28:36].strip())
                z = float(line[36:44].strip())
                coordinates.append([x, y, z])
                
                # 如果是fSOL残基，记录其z坐标
                if res_name == 'fSOL':
                    fSOL_z_coordinates.append(z)
                    
            except ValueError:
                # 如果提取坐标失败，打印警告并跳过该行
                print(f"Warning: Skipping line due to invalid coordinate data: {line.strip()}")
    
    # 计算fSOL的z坐标的最小值和最大值
    z_fSOL_min = np.min(fSOL_z_coordinates) if fSOL_z_coordinates else None
    z_fSOL_max = np.max(fSOL_z_coordinates) if fSOL_z_coordinates else None
    
    return np.array(coordinates), box_size, z_fSOL_min, z_fSOL_max

def calculate_particle_distribution(coordinates, box_size, n_layers, z_fSOL_min, z_fSOL_max, output_file='particle_distribution.txt'):
    """计算z方向上粒子数分布并进行分层处理"""
    # 计算z轴的切分区间
    z_min, z_max = np.min(coordinates[:, 2]), np.max(coordinates[:, 2])
    z_step = (z_max - z_min) / n_layers
    
    # 初始化每一层的粒子数
    particle_distribution = np.zeros(n_layers)
    z_coordinates = np.zeros(n_layers)
    
    # 根据z坐标分配粒子到分层
    for coord in coordinates:
        z = coord[2]
        # 找到z坐标对应的层
        layer_index = int((z - z_min) // z_step)
        if layer_index >= n_layers:
            layer_index = n_layers - 1  # 处理边界情况
        particle_distribution[layer_index] += 1
    
    # 计算每一层的z坐标
    for i in range(n_layers):
        if i == 0:
            z_coordinates[i] = z_min
        elif i == n_layers - 1:
            z_coordinates[i] = z_max
        else:
            z_coordinates[i] = z_min + (i + 0.5) * z_step
    
    # 根据新的分层规范进行计算
    # 下水区：若z_fSOL_min所在层为n，则1~(n-1)层为下水区
    lower_water_layers = int((z_fSOL_min - z_min) // z_step)
    lower_water_end = lower_water_layers  # 下水区结束层
    
    # 冰区：z_fSOL_min到z_fSOL_max这段范围的层应被划为冰层
    ice_layers_start = lower_water_end + 1  # 冰区从下水区后面开始
    ice_layers_end = int(ice_layers_start + (z_fSOL_max - z_fSOL_min) // z_step)  # 冰区结束层
    
    # 上水区：若z_fSOL_max所在层为m，则(m+1)~z_max为上水区
    upper_water_layers = int((z_max - z_fSOL_max) // z_step)
    upper_water_start = ice_layers_end + 1  # 上水区起始层
    
    # 输出每层的粒子分布和统计数据
    with open(output_file, 'w') as txtfile:
        # 添加表头
        txtfile.write("Layer Number\tz Coordinate (nm)\tParticle Count\n")
        for i in range(n_layers):
            txtfile.write(f"{i + 1}\t{z_coordinates[i]:.4f}\t{particle_distribution[i]}\n")
    
    print(f"\nData saved to '{output_file}'")
    
    # 输出分层信息，调整为连续且不重叠
    print(f"\nLower Water Layers: 1 ~ {lower_water_end}")  # 输出下水区
    print(f"Ice Layers: {ice_layers_start} ~ {ice_layers_end}")  # 输出冰区
    print(f"Upper Water Layers: {upper_water_start} ~ {n_layers}")  # 输出上水区
    
    return lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers

def calculate_statistics(particle_distribution):
    """计算粒子数分布的统计数据：均值、方差、标准差"""
    mean = np.mean(particle_distribution)
    variance = np.var(particle_distribution)
    std_dev = np.std(particle_distribution)
    return mean, variance, std_dev

def output_statistics(mean, variance, std_dev, lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers, output_file="stat-result.txt"):
    """输出统计数据到文件"""
    with open(output_file, 'w') as f:
        f.write(f"Mean: {mean:.4f}\n")
        f.write(f"Variance: {variance:.4f}\n")
        f.write(f"Standard Deviation: {std_dev:.4f}\n")

        f.write(f"\nLower Water Layers: 1,{lower_water_end}")  # 输出下水区
        f.write(f"\nIce Layers: {ice_layers_start },{ice_layers_end}")  # 输出冰区
        f.write(f"\nUpper Water Layers: {upper_water_start},{n_layers}")  # 输出上水区
    print(f"Statistics saved to '{output_file}'\n")

def main():
    """主函数"""
    
    # 在这里集中的定义输入输出文件路径和其他参数
    gro_file = 'npt-10.gro'  # 输入你的GRO文件路径
    output_file = 'npt-10.txt'  # 输出文件名
    stat_file = 'stat-result-total.txt'  # 统计结果文件名
    n_layers = 200  # 设置分层数
    
    # 读取gro文件
    coordinates, box_size, z_fSOL_min, z_fSOL_max = read_gro_file(gro_file)
    
    if z_fSOL_min is None or z_fSOL_max is None:
        print("Error: No fSOL residues found in the GRO file.")
        return
    
    # 计算粒子分布并输出为txt
    lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers = calculate_particle_distribution(coordinates, box_size, n_layers=n_layers, z_fSOL_min=z_fSOL_min, z_fSOL_max=z_fSOL_max, output_file=output_file)
    
    # 计算统计数据
    particle_distribution = np.loadtxt(output_file, delimiter='\t', skiprows=1, usecols=2)  # 读取分布数据，跳过表头
    mean, variance, std_dev = calculate_statistics(particle_distribution)
    
    # 在输出到文件之前打印整体均值和标准差
    print(f"\nMean Particle Count: {mean:.4f}")
    print(f"Standard Deviation: {std_dev:.4f}")
    
    # 输出统计数据到文件
    output_statistics(mean, variance, std_dev, lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers, output_file=stat_file)

if __name__ == "__main__":
    main()
