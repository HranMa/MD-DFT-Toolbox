'''
一个简单的三轴换四轴工具
'''

def convert_plane_to_four_axis(h, k, l):
    """三轴 (hkl) -> 四轴 (hkil)"""
    i = -(h + k)
    return h, k, i, l

def convert_plane_to_three_axis(h, k, i, l):
    """四轴 (hkil) -> 三轴 (hkl)"""
    return h, k, l  # 去掉 i

def convert_direction_to_four_axis(U, V, W):
    """三轴 [UVW] -> 四轴 [uvtw]"""
    u = (2 * U - V) / 3
    v = (2 * V - U) / 3
    t = -(u + v)
    w = W
    return u, v, t, w

def convert_direction_to_three_axis(u, v, t, w):
    """四轴 [uvtw] -> 三轴 [UVW]"""
    U = u - t
    V = v - t
    W = w
    return U, V, W

def main():
    print("欢迎使用三轴与四轴晶体学转换工具！")
    print("请选择转换类型：")
    print("1. 晶面转换 (hkl ↔ hkil)")
    print("2. 晶向转换 [UVW] ↔ [uvtw]")

    conversion_type = input("请输入选项（1或2）：")
    if conversion_type not in ["1", "2"]:
        print("输入无效，请输入 1 或 2。")
        return

    direction = input("请选择转换方向（输入 '3to4' 表示三轴到四轴，'4to3' 表示四轴到三轴）：")
    if direction not in ["3to4", "4to3"]:
        print("输入无效，请输入 '3to4' 或 '4to3'。")
        return

    if conversion_type == "1":  # 晶面转换
        if direction == "3to4":
            h = int(input("请输入三轴 h 值："))
            k = int(input("请输入三轴 k 值："))
            l = int(input("请输入三轴 l 值："))
            h, k, i, l = convert_plane_to_four_axis(h, k, l)
            print(f"四轴晶面指数为：({h} {k} {i} {l})")
        elif direction == "4to3":
            h = int(input("请输入四轴 h 值："))
            k = int(input("请输入四轴 k 值："))
            i = int(input("请输入四轴 i 值："))
            l = int(input("请输入四轴 l 值："))
            h, k, l = convert_plane_to_three_axis(h, k, i, l)
            print(f"三轴晶面指数为：({h} {k} {l})")

    elif conversion_type == "2":  # 晶向转换
        if direction == "3to4":
            U = int(input("请输入三轴 U 值："))
            V = int(input("请输入三轴 V 值："))
            W = int(input("请输入三轴 W 值："))
            u, v, t, w = convert_direction_to_four_axis(U, V, W)
            print(f"四轴晶向指数为：[ {u:.2f} {v:.2f} {t:.2f} {w} ]")
        elif direction == "4to3":
            u = float(input("请输入四轴 u 值："))
            v = float(input("请输入四轴 v 值："))
            t = float(input("请输入四轴 t 值："))
            w = int(input("请输入四轴 w 值："))
            U, V, W = convert_direction_to_three_axis(u, v, t, w)
            print(f"三轴晶向指数为：[ {U:.2f} {V:.2f} {W} ]")

if __name__ == "__main__":
    main()

