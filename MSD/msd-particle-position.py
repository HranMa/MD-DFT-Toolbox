import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非图形化的Agg后端
import matplotlib.pyplot as plt

'''
配合particle-position.py使用
'''

# 设置周期性边界的盒子尺寸 (Lx, Ly, Lz)，单位：Ångström
Lx, Ly, Lz = 35, 35, 90  # 盒子的尺寸，以Å为单位

# 计算周期性边界修正的函数
def apply_periodic_boundary(x, L):
    return x - L * np.round(x / L)

# 读取文件，处理每个文件的粒子轨迹
def read_particle_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)  # 跳过第一行表头
    return data  # 返回 (n, 3) 的数组，n为时间步长，3为X, Y, Z坐标

# 计算均方位移（MSD）
def calculate_msd(trajectory):
    # 计算初始位置与每个时刻位置的位移
    initial_position = trajectory[0]  # 初始时刻位置
    displacements = trajectory - initial_position  # 所有时刻的位移
    # 处理周期性边界条件
    displacements[:, 0] = apply_periodic_boundary(displacements[:, 0], Lx)
    displacements[:, 1] = apply_periodic_boundary(displacements[:, 1], Ly)
    displacements[:, 2] = apply_periodic_boundary(displacements[:, 2], Lz)
    
    # 计算每个时刻的位移平方
    msd = np.mean(np.square(displacements), axis=1)
    return msd

# 获取所有txt文件路径
def get_txt_files():
    # 获取output目录下所有txt文件
    files = glob.glob("./output/*.txt")
    return files

# 主函数
def main():
    # 获取所有文件
    txt_files = get_txt_files()
    
    # 打印所有txt文件，询问是否继续
    print("读取以下文件：")
    for f in txt_files:
        print(f)
    
    proceed = input("是否继续处理这些文件? (y/n): ")
    if proceed.lower() != 'y':
        print("已取消处理")
        return
    
    msd_all = []  # 存储所有文件的MSD数据

    # 处理每一个txt文件
    for file_path in txt_files:
        print(f"正在处理文件: {file_path}")
        trajectory = read_particle_data(file_path)  # 读取粒子轨迹
        msd = calculate_msd(trajectory)  # 计算MSD
        msd_all.append(msd)
    
    # 计算所有文件的平均MSD
    msd_all = np.array(msd_all)
    mean_msd = np.mean(msd_all, axis=0)

    # 输出到msd.txt
    os.makedirs('./msd', exist_ok=True)  # 确保msd目录存在
    np.savetxt('./msd/msd.txt', mean_msd)

    # 绘制MSD曲线并保存为png图片
    plt.figure()
    plt.plot(mean_msd, label="MSD")
    plt.xlabel("Time Steps")
    plt.ylabel("Mean Squared Displacement (MSD) [Å²]")
    plt.title("Mean Squared Displacement (MSD)")
    plt.legend()
    plt.savefig('./msd/msd.png', dpi=300)  # 输出为png文件

    print("处理完毕，MSD已保存为msd/msd.txt 和 msd/msd.png")

# 调用主函数
if __name__ == "__main__":
    main()
