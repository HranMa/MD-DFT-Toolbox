# 打开原始 xvg 文件和新文件
input_file = 'water-ice.txt'  # 请替换为你的输入文件名
output_file = 'water-ice-output.txt'  # 新文件名

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for index, line in enumerate(infile):
        # 只写入奇数行（index 为偶数的行）
        if index % 2 == 0:
            outfile.write(line)

print(f'已将偶数行删除，并保存到 {output_file}')
