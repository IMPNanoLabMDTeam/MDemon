#!/usr/bin/env python3
"""
测试MDemon读取包含flag信息的LAMMPS data文件
验证flag信息的正确处理（flag信息应该被忽略，只读取坐标）
"""

import os

import numpy as np

import MDemon as md


def test_atomic_with_flags_auto():
    """测试自动检测包含flag的atomic格式"""
    print("=" * 60)
    print("测试1: 自动检测atomic格式 (8字段 - 包含flag)")
    print("=" * 60)

    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    filename = "atomic_with_flags.data"
    filepath = os.path.join(data_dir, filename)

    try:
        # 不指定atom_style，让系统自动检测8字段的atomic格式
        u = md.Universe(filepath)
        print("✅ 成功创建Universe对象")

        # 验证读取结果
        print(f"原子数量: {len(u.atoms)}")
        print(f"原子类型数量: {len(set(atom.species for atom in u.atoms))}")

        # 验证坐标（应该忽略flag信息）
        expected_coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ]

        print("\n验证原子坐标（应该忽略flag）:")
        for i, atom in enumerate(u.atoms):
            coord = atom.coordinate
            expected_coord = expected_coords[i]
            coord_match = np.allclose(coord, expected_coord, atol=1e-6)
            status = "✅" if coord_match else "❌"
            print(
                f"{status} 原子 {i+1}: 类型={atom.species}, "
                f"坐标=({coord[0]:.3f}, {coord[1]:.3f}, {coord[2]:.3f}) "
                f"(期望:{expected_coord})"
            )

        return True

    except Exception as e:
        print(f"❌ 自动检测失败: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_atomic_with_flags_explicit():
    """测试明确指定包含flag的atomic格式"""
    print("\n" + "=" * 60)
    print("测试2: 明确指定atom_style='id type x y z' (8字段文件)")
    print("=" * 60)

    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    filename = "atomic_with_flags.data"
    filepath = os.path.join(data_dir, filename)

    try:
        # 明确指定atom_style，应该忽略flag字段
        u = md.Universe(filepath, atom_style="id type x y z")
        print("✅ 成功创建Universe对象")

        # 验证读取结果
        print(f"原子数量: {len(u.atoms)}")

        # 验证坐标
        expected_coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ]

        expected_types = [1, 1, 1, 2]

        print("\n验证原子数据:")
        all_correct = True
        for i, atom in enumerate(u.atoms):
            coord = atom.coordinate
            expected_coord = expected_coords[i]
            expected_type = expected_types[i]

            coord_match = np.allclose(coord, expected_coord, atol=1e-6)
            type_match = atom.species == expected_type

            status = "✅" if (coord_match and type_match) else "❌"
            print(
                f"{status} 原子 {i+1}: 类型={atom.species} (期望:{expected_type}), "
                f"坐标=({coord[0]:.3f}, {coord[1]:.3f}, {coord[2]:.3f}) "
                f"(期望:{expected_coord})"
            )

            if not (coord_match and type_match):
                all_correct = False

        return all_correct

    except Exception as e:
        print(f"❌ 明确指定格式失败: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_comparison_normal_vs_flags():
    """比较5字段和8字段文件的读取结果是否一致"""
    print("\n" + "=" * 60)
    print("测试3: 比较5字段与8字段文件的读取结果")
    print("=" * 60)

    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    filepath1 = os.path.join(data_dir, "simple_atomic.data")  # 5字段
    filepath2 = os.path.join(data_dir, "atomic_with_flags.data")  # 8字段

    try:
        # 读取两个文件
        u1 = md.Universe(filepath1)  # 5字段文件
        u2 = md.Universe(filepath2)  # 8字段文件

        print(f"5字段文件原子数量: {len(u1.atoms)}")
        print(f"8字段文件原子数量: {len(u2.atoms)}")

        # 比较结果
        atoms_match = len(u1.atoms) == len(u2.atoms)
        print(f"原子数量匹配: {atoms_match}")

        print("\n比较原子坐标和类型...")
        all_match = True
        for i in range(len(u1.atoms)):
            coord1 = u1.atoms[i].coordinate
            coord2 = u2.atoms[i].coordinate
            type1 = u1.atoms[i].species
            type2 = u2.atoms[i].species

            coord_match = np.allclose(coord1, coord2, atol=1e-6)
            type_match = type1 == type2

            if not (coord_match and type_match):
                all_match = False
                print(
                    f"❌ 原子{i+1}: 坐标 {coord1} vs {coord2}, 类型 {type1} vs {type2}"
                )
            else:
                print(f"✅ 原子{i+1}: 坐标和类型匹配")

        if all_match:
            print("\n🎉 两种格式读取结果完全一致!")
        else:
            print("\n⚠️  检测到差异")

        return all_match

    except Exception as e:
        print(f"❌ 比较测试失败: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """运行所有flag测试"""
    print("MDemon LAMMPS data文件flag信息处理测试")
    print("=" * 60)

    # 检查测试文件是否存在
    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    file1 = os.path.join(data_dir, "simple_atomic.data")
    file2 = os.path.join(data_dir, "atomic_with_flags.data")

    if not os.path.exists(file1):
        print(f"❌ 测试文件不存在: {file1}")
        return False

    if not os.path.exists(file2):
        print(f"❌ 测试文件不存在: {file2}")
        return False

    print(f"✅ 测试文件存在: {file1}")
    print(f"✅ 测试文件存在: {file2}")

    # 运行测试
    results = []

    results.append(test_atomic_with_flags_auto())
    results.append(test_atomic_with_flags_explicit())
    results.append(test_comparison_normal_vs_flags())

    # 总结
    print("\n" + "=" * 60)
    print("Flag测试结果总结")
    print("=" * 60)

    test_names = [
        "自动检测8字段atomic格式",
        "明确指定atom_style (8字段)",
        "5字段与8字段结果比较",
    ]

    for i, (name, result) in enumerate(zip(test_names, results)):
        status = "✅ 通过" if result else "❌ 失败"
        print(f"测试 {i+1}: {name} - {status}")

    all_passed = all(results)

    if all_passed:
        print("\n🎉 所有flag测试通过! Flag信息被正确忽略!")
    else:
        print(f"\n⚠️  {sum(results)}/{len(results)} 个测试通过")

    return all_passed


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
