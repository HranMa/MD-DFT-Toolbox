"""
Hexagonal-Cubic Axis Conversion Tool

This script converts crystallographic indices between three-axis (hkl or [UVW]) 
and four-axis (hkil or [uvtw]) notations used in hexagonal and cubic crystal systems.
"""

def convert_plane_3_to_4_axis(h, k, l):
    """Convert three-axis plane indices (hkl) to four-axis (hkil)."""
    i = -(h + k)
    return h, k, i, l

def convert_plane_4_to_3_axis(h, k, i, l):
    """Convert four-axis plane indices (hkil) to three-axis (hkl)."""
    return h, k, l  # drop i

def convert_direction_3_to_4_axis(U, V, W):
    """Convert three-axis direction indices [UVW] to four-axis [uvtw]."""
    u = (2 * U - V) / 3
    v = (2 * V - U) / 3
    t = -(u + v)
    w = W
    return u, v, t, w

def convert_direction_4_to_3_axis(u, v, t, w):
    """Convert four-axis direction indices [uvtw] to three-axis [UVW]."""
    U = u - t
    V = v - t
    W = w
    return U, V, W

def main():
    print("Hexagonal-Cubic Axis Conversion Tool")
    print("Select conversion type:")
    print("1. Plane indices conversion (hkl ↔ hkil)")
    print("2. Direction indices conversion [UVW] ↔ [uvtw]")

    conversion_type = input("Enter option (1 or 2): ")
    if conversion_type not in ["1", "2"]:
        print("Invalid input. Please enter 1 or 2.")
        return

    direction = input("Choose conversion direction ('3to4' for 3-axis to 4-axis, '4to3' for 4-axis to 3-axis): ")
    if direction not in ["3to4", "4to3"]:
        print("Invalid input. Please enter '3to4' or '4to3'.")
        return

    if conversion_type == "1":  # Plane conversion
        if direction == "3to4":
            h = int(input("Enter three-axis h value: "))
            k = int(input("Enter three-axis k value: "))
            l = int(input("Enter three-axis l value: "))
            h, k, i, l = convert_plane_3_to_4_axis(h, k, l)
            print(f"Four-axis plane indices: ({h} {k} {i} {l})")
        elif direction == "4to3":
            h = int(input("Enter four-axis h value: "))
            k = int(input("Enter four-axis k value: "))
            i = int(input("Enter four-axis i value: "))
            l = int(input("Enter four-axis l value: "))
            h, k, l = convert_plane_4_to_3_axis(h, k, i, l)
            print(f"Three-axis plane indices: ({h} {k} {l})")

    elif conversion_type == "2":  # Direction conversion
        if direction == "3to4":
            U = int(input("Enter three-axis U value: "))
            V = int(input("Enter three-axis V value: "))
            W = int(input("Enter three-axis W value: "))
            u, v, t, w = convert_direction_3_to_4_axis(U, V, W)
            print(f"Four-axis direction indices: [ {u:.2f} {v:.2f} {t:.2f} {w} ]")
        elif direction == "4to3":
            u = float(input("Enter four-axis u value: "))
            v = float(input("Enter four-axis v value: "))
            t = float(input("Enter four-axis t value: "))
            w = int(input("Enter four-axis w value: "))
            U, V, W = convert_direction_4_to_3_axis(u, v, t, w)
            print(f"Three-axis direction indices: [ {U:.2f} {V:.2f} {W:.2f} ]")

if __name__ == "__main__":
    main()
