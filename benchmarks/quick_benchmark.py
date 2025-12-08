#!/usr/bin/env python3
"""
Quick Benchmark Script for CMTJ - Tests current version performance
No comparison, no matplotlib required - just fast timing results
"""

import time
import numpy as np
from typing import Callable

MAX_TIME = 100e-9
ITERATIONS_PER_TEST = 10


def timer(func: Callable, iterations: int = ITERATIONS_PER_TEST, warmup: int = 1) -> dict:
    """Time a function with warmup and multiple iterations"""
    # Warmup
    for _ in range(warmup):
        func()

    # Benchmark
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "mean": np.mean(times),
        "std": np.std(times),
        "min": np.min(times),
        "max": np.max(times),
    }


def main():
    print("="*70)
    print("CMTJ Quick Performance Benchmark")
    print("="*70)

    try:
        from cmtj import Layer, Junction, CVector, constantDriver, AxialDriver, SolverMode
        from cmtj.utils import FieldScan
    except ImportError as e:
        print(f"Error: Could not import cmtj: {e}")
        print("Please install: pip install -e .")
        return 1

    # Test 1: Single layer RK4
    print("\n[1/5] Single layer RK4 simulation...")

    def test1():
        layer = Layer(
            "free",
            mag=CVector(0.1, 0.1, 0.99),
            anis=CVector(0, 0, 1),
            Ms=1.0,
            thickness=1.5e-9,
            damping=0.01,
            cellSurface=100e-18,
            demagTensor=[CVector(0, 0, 0), CVector(
                0, 0, 0), CVector(0, 0, 1.0)],
        )
        layer.setAnisotropyDriver(constantDriver(400e3))
        j = Junction([layer])
        j.runSimulation(MAX_TIME, 1e-13, 1e-11)

    result1 = timer(test1, iterations=ITERATIONS_PER_TEST, warmup=1)
    print(f"  Time: {result1['mean']:.4f}s ± {result1['std']:.4f}s")

    # Test 2: Multi-layer
    print("\n[2/5] Multi-layer (3 layers) RK4 simulation...")

    def test2():
        demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
        l1 = Layer("l1", CVector(0, 0, 0.99), CVector(
            0, 0, 1), 1.03, 1e-9, 1e-9, demag, 0.01)
        l2 = Layer("l2", CVector(0, 0, 0.99), CVector(
            0, 0, 1), 1.03, 1e-9, 1e-9, demag, 0.01)
        l3 = Layer("l3", CVector(1, 0, 0), CVector(
            1, 0, 0), 0.12, 1e-9, 1e-9, demag, 0.015)
        l1.setAnisotropyDriver(constantDriver(489e3))
        l2.setAnisotropyDriver(constantDriver(514e3))
        l3.setAnisotropyDriver(constantDriver(17e3))
        j = Junction([l1, l2, l3])
        j.setIECDriver("l1", "l2", constantDriver(-1e-3))
        j.runSimulation(MAX_TIME, 5e-13, 1e-11)

    result2 = timer(test2, iterations=ITERATIONS_PER_TEST, warmup=1)
    print(f"  Time: {result2['mean']:.4f}s ± {result2['std']:.4f}s")

    # Test 3: Dormand-Prince
    print("\n[3/5] Adaptive Dormand-Prince solver...")

    def test3():
        layer = Layer(
            "free",
            mag=CVector(0.1, 0.1, 0.99),
            anis=CVector(0, 0, 1),
            Ms=1.0,
            thickness=1.5e-9,
            damping=0.01,
            cellSurface=100e-18,
            demagTensor=[CVector(0, 0, 0), CVector(
                0, 0, 0), CVector(0, 0, 1.0)],
        )
        layer.setAnisotropyDriver(constantDriver(400e3))
        j = Junction([layer])
        j.runSimulation(MAX_TIME, 1e-13, 1e-13,
                        solverMode=SolverMode.DormandPrince)

    result3 = timer(test3, iterations=ITERATIONS_PER_TEST, warmup=1)
    print(f"  Time: {result3['mean']:.4f}s ± {result3['std']:.4f}s")

    # Test 4: Field sweep
    print("\n[4/5] Field sweep (50 steps)...")

    def test4():
        layer = Layer(
            "free",
            mag=CVector(0, 0, 1),
            anis=CVector(0, 0, 1),
            Ms=1.0,
            thickness=1.5e-9,
            damping=0.01,
            cellSurface=100e-18,
            demagTensor=[CVector(0, 0, 0), CVector(
                0, 0, 0), CVector(0, 0, 1.0)],
        )
        layer.setAnisotropyDriver(constantDriver(400e3))
        j = Junction([layer])

        Hscan, Hvecs = FieldScan.amplitude_scan(-500e3, 500e3, 50, 0, 0)
        for Hv in Hvecs:
            j.clearLog()
            j.setLayerExternalFieldDriver("all", AxialDriver(*Hv))
            j.runSimulation(5e-9, 1e-13, 1e-11)

    result4 = timer(test4, iterations=ITERATIONS_PER_TEST, warmup=1)
    print(f"  Time: {result4['mean']:.4f}s ± {result4['std']:.4f}s")

    # Test 5: Tensor operations
    print("\n[5/5] Heavy tensor interaction calculations...")

    def test5():
        dipole = [CVector(0.1, 0.05, 0.02), CVector(
            0.05, 0.1, 0.02), CVector(0.02, 0.02, 0.05)]
        demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
        l1 = Layer("l1", CVector(0.1, 0.1, 0.99), CVector(
            0, 0, 1), 1.0, 1e-9, 1e-9, demag, 0.01)
        l2 = Layer("l2", CVector(0, 0, 1), CVector(
            0, 0, 1), 1.0, 1e-9, 1e-9, demag, 0.01)
        l1.setAnisotropyDriver(constantDriver(400e3))
        l2.setAnisotropyDriver(constantDriver(400e3))
        l1.setTopDipoleTensor(dipole)
        l2.setBottomDipoleTensor(dipole)
        j = Junction([l1, l2])
        j.runSimulation(MAX_TIME, 5e-14, 1e-11)

    result5 = timer(test5, iterations=ITERATIONS_PER_TEST, warmup=1)
    print(f"  Time: {result5['mean']:.4f}s ± {result5['std']:.4f}s")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    total_time = result1['mean'] + result2['mean'] + \
        result3['mean'] + result4['mean'] + result5['mean']
    print(f"Total benchmark time: {total_time:.4f}s")
    print("\nIndividual results:")
    print(
        f"  1. Single layer RK4:          {result1['mean']:>8.4f}s ± {result1['std']:.4f}s")
    print(
        f"  2. Multi-layer RK4:           {result2['mean']:>8.4f}s ± {result2['std']:.4f}s")
    print(
        f"  3. Dormand-Prince:            {result3['mean']:>8.4f}s ± {result3['std']:.4f}s")
    print(
        f"  4. Field sweep (50 steps):    {result4['mean']:>8.4f}s ± {result4['std']:.4f}s")
    print(
        f"  5. Tensor operations:         {result5['mean']:>8.4f}s ± {result5['std']:.4f}s")
    print("="*70)
    print("\nTo compare with a previous version, use:")
    print("  python benchmark_performance.py --compare HEAD~1")
    print("="*70)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
