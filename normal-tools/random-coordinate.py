import random

'''
本脚本可以用来生成一定熟练的随机坐标
格式服从gmx的gro文件格式
'''

def generate_random_coordinates(x_range, y_range, z_range, num_points=10, output_file="random-coord.txt"):
    # 打开文件进行写入
    with open(output_file, "w") as file:
        for _ in range(num_points):
            # 生成三个坐标值，每个坐标值在指定的范围内
            x = round(random.uniform(x_range[0], x_range[1]), 3)
            y = round(random.uniform(y_range[0], y_range[1]), 3)
            z = round(random.uniform(z_range[0], z_range[1]), 3)
            # 强制保留3位小数，使用字符串格式化
            file.write(f"{x:.3f}   {y:.3f}   {z:.3f}\n")
# 硬编码的坐标范围
x_range = (0, 3.54623)  # x坐标范围 (1.0 到 5.0)
y_range = (0, 3.54623)  # y坐标范围 (-5.0 到 5.0)
z_range = (3.35, 8.5)  # z坐标范围 (2.0 到 8.0)

# 硬编码输出的随机坐标点数
num_points = 15  # 输出10个随机坐标点

# 生成随机坐标并写入文件
generate_random_coordinates(x_range, y_range, z_range, num_points, output_file="random-coord.txt")

print("随机坐标已生成并写入到 'random-coord.txt' 文件中。")
