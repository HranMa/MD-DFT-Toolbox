import re

"""
脚本工作逻辑：
处理 Gromacs 的 gro 文件，保留冰部分，删除位于 ice_zone 或超出盒子范围的水分子。
使用gmx solvate填充水分子的时候，可能会有一些水分子位于冰内部，这些水分子不应该被填充。
还有会一些超出模拟盒子范围的水分子，影响不大但可以删去。

ice_zone = 由genice生成的冰盒子
参数:
- input_file: 输入的 gro 文件路径/名称
- output_file: 输出的 gro 文件路径/名称
- box_dimensions: 手动设定的盒子大小 [x, y, z] (单位：nm)
- z_buffer: z 方向缓冲区大小，用于调整 ice_zone 范围

注意，当前脚本删除后水分子后，可能会导致原子编号不连续，但不影响后续模拟
建议额外执行一次EM能量最小化，并且使用输出文件进行后续模拟
"""

def process_gro_file(input_file, output_file, box_dimensions, z_buffer=0.0):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # Header lines
    header = lines[:2]
    box_line = lines[-1]
    atoms = lines[2:-1]

    # Parse atoms and separate into fSOL and SOL
    fsol_atoms = []
    sol_atoms = []
    for line in atoms:
        residue_name = line[5:10].strip()
        if residue_name == "fSOL":
            fsol_atoms.append(line)
        elif residue_name == "SOL":
            sol_atoms.append(line)
    
    # Get z-coordinate range for fSOL
    z_coords = [float(line[36:44].strip()) for line in fsol_atoms]
    z_min, z_max = min(z_coords), max(z_coords)

    # Define the ice zone with user-defined box dimensions
    ice_zone = {
        "x_min": 0.0,
        "x_max": box_dimensions[0],
        "y_min": 0.0,
        "y_max": box_dimensions[1],
        "z_min": max(0.0, z_min - z_buffer),  # z_min adjusted by buffer
        "z_max": min(box_dimensions[2], z_max + z_buffer),  # z_max adjusted by buffer
    }

    # Check and filter SOL atoms
    filtered_sol_atoms = []
    deleted_molecule_count = 0
    molecule_start = None
    for i, line in enumerate(sol_atoms):
        if i % 4 == 0:  # Start of a new water molecule
            molecule_start = i
        
        # Check if the atom is in ice_zone or out of box_dimensions
        x, y, z = map(float, (line[20:28].strip(), line[28:36].strip(), line[36:44].strip()))
        in_ice_zone = (
            ice_zone["x_min"] <= x <= ice_zone["x_max"] and
            ice_zone["y_min"] <= y <= ice_zone["y_max"] and
            ice_zone["z_min"] <= z <= ice_zone["z_max"]
        )
        out_of_box = (
            x < 0.0 or x > box_dimensions[0] or
            y < 0.0 or y > box_dimensions[1] or
            z < 0.0 or z > box_dimensions[2]
        )

        # If any atom of the molecule is in ice_zone or out of box, skip the whole molecule
        if in_ice_zone or out_of_box:
            if i % 4 == 3:  # At the end of a molecule
                deleted_molecule_count += 1
            continue

        # If molecule is valid, keep all its atoms
        if i % 4 == 3:  # End of a water molecule
            filtered_sol_atoms.extend(sol_atoms[molecule_start:molecule_start + 4])
    
    # Renumber the atoms
    combined_atoms = fsol_atoms + filtered_sol_atoms
    renumbered_atoms = []
    for atom_index, line in enumerate(combined_atoms, start=1):
        new_line = line[:15] + f"{atom_index:>5}" + line[20:]
        renumbered_atoms.append(new_line)

    # Update the total number of atoms
    total_atoms = len(renumbered_atoms)
    header[1] = f"{total_atoms}\n"

    # Write the output
    with open(output_file, 'w') as file:
        file.writelines(header)
        file.writelines(renumbered_atoms)
        file.write(box_line)
    
    # Print completion message
    print(f"处理完成！总共删除了 {deleted_molecule_count} 个水分子。输出文件为：{output_file}")


# ============= 配置区域 =============
# 输入和输出文件路径——建议把脚本和input放在一个路径下使用
input_gro = "ice-water.gro"         # 输入文件
output_gro = "ice-water-output.gro"  # 输出文件

# 手动设定盒子大小 (单位: nm)
box_dimensions = [4.0, 4.0, 10.0]  # gro文件最后的x, y, z 尺寸

# z 方向缓冲区大小 (单位: nm)
z_buffer = 0.1  # 根据需要调整，默认0.1即可

# 调用主函数
process_gro_file(input_gro, output_gro, box_dimensions, z_buffer)

