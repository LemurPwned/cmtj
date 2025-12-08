#!/usr/bin/env python3
"""
Performance Benchmark Script for CMTJ
Compares optimized version performance against baseline (previous git commit)

Usage:
    python benchmark_performance.py                    # Benchmark current version
    python benchmark_performance.py --compare HEAD~1   # Compare with previous commit
    python benchmark_performance.py --compare <commit> # Compare with specific commit
"""

import argparse
import subprocess
import sys
import time
import json
from pathlib import Path
import numpy as np
from typing import Callable
import matplotlib.pyplot as plt

# Benchmark configuration
BENCHMARKS = {
    "single_layer_rk4": {
        "desc": "Single layer RK4 simulation",
        "iterations": 10,
        "warmup": 1,
    },
    "multi_layer_rk4": {
        "desc": "Multi-layer (3 layers) RK4 simulation", 
        "iterations": 10,
        "warmup": 1,
    },
    "dormand_prince": {
        "desc": "Single layer adaptive Dormand-Prince solver",
        "iterations": 100,
        "warmup": 1,
    },
    "field_sweep": {
        "desc": "Field sweep with 50 steps",
        "iterations": 10,
        "warmup": 1,
    },
    "tensor_operations": {
        "desc": "Heavy tensor interaction calculations",
        "iterations": 10,
        "warmup": 1,
    },
}


def run_benchmark(name: str, func: Callable, iterations: int = 5, warmup: int = 1) -> dict:
    """Run a single benchmark and return timing statistics"""
    times = []
    
    # Warmup runs
    for _ in range(warmup):
        func()
    
    # Actual benchmark runs
    for _ in range(iterations):
        start = time.perf_counter()
        func()
        end = time.perf_counter()
        times.append(end - start)
    
    return {
        "name": name,
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "min": np.min(times),
        "max": np.max(times),
        "median": np.median(times),
    }


def benchmark_single_layer_rk4():
    """Benchmark single layer RK4 simulation"""
    from cmtj import Layer, Junction, CVector, constantDriver
    
    Ms = 1.0
    Ku = 400e3
    damping = 0.01
    
    layer = Layer(
        "free",
        mag=CVector(0.1, 0.1, 0.99),
        anis=CVector(0, 0, 1),
        Ms=Ms,
        thickness=1.5e-9,
        damping=damping,
        cellSurface=100e-18,
        demagTensor=[CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)],
    )
    layer.setAnisotropyDriver(constantDriver(Ku))
    
    junction = Junction([layer])
    
    # Run simulation
    sim_time = 10e-9  # 10 nanoseconds
    dt = 1e-13
    write_freq = 1e-11
    
    junction.runSimulation(sim_time, dt, write_freq)


def benchmark_multi_layer_rk4():
    """Benchmark multi-layer RK4 simulation with coupling"""
    from cmtj import Layer, Junction, CVector, constantDriver, AxialDriver
    
    Ms1 = 1.03
    Ms2 = 1.03
    Ms3 = 0.12
    Ku1 = 489e3
    Ku2 = 514e3
    Ku3 = 17e3
    damping = 0.01
    
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
    
    l1 = Layer(
        "free",
        mag=CVector(0, 0, 0.99),
        anis=CVector(0, 0, 1),
        Ms=Ms1,
        thickness=1.0e-9,
        damping=damping,
        demagTensor=demag,
        cellSurface=1e-9,
    )
    l2 = Layer(
        "perpendicular",
        mag=CVector(0, 0, 0.99),
        anis=CVector(0, 0, 1),
        Ms=Ms2,
        thickness=1.0e-9,
        damping=damping,
        cellSurface=1e-9,
        demagTensor=demag,
    )
    l3 = Layer(
        "buffer",
        mag=CVector(1, 0, 0),
        anis=CVector(1, 0, 0),
        Ms=Ms3,
        thickness=1.0e-9,
        damping=0.015,
        cellSurface=1e-9,
        demagTensor=demag,
    )
    
    l1.setAnisotropyDriver(constantDriver(Ku1))
    l2.setAnisotropyDriver(constantDriver(Ku2))
    l3.setAnisotropyDriver(constantDriver(Ku3))
    
    junction = Junction([l1, l2, l3])
    junction.setIECDriver("free", "perpendicular", constantDriver(-1e-3))
    
    # External field
    Hext = 500e3
    junction.setLayerExternalFieldDriver("all", AxialDriver(constantDriver(0), constantDriver(0), constantDriver(Hext)))
    
    # Run simulation
    sim_time = 10e-9
    dt = 5e-13
    write_freq = 1e-11
    
    junction.runSimulation(sim_time, dt, write_freq)


def benchmark_dormand_prince():
    """Benchmark adaptive Dormand-Prince solver"""
    from cmtj import Layer, Junction, CVector, constantDriver, SolverMode
    
    Ms = 1.0
    Ku = 400e3
    damping = 0.01
    
    layer = Layer(
        "free",
        mag=CVector(0.1, 0.1, 0.99),
        anis=CVector(0, 0, 1),
        Ms=Ms,
        thickness=1.5e-9,
        damping=damping,
        cellSurface=100e-18,
        demagTensor=[CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)],
    )
    layer.setAnisotropyDriver(constantDriver(Ku))
    
    junction = Junction([layer])
    
    # Run with adaptive solver
    sim_time = 10e-9
    dt = 1e-13
    write_freq = 1e-11
    
    junction.runSimulation(sim_time, dt, write_freq, solverMode=SolverMode.DormandPrince)


def benchmark_field_sweep():
    """Benchmark field sweep simulation"""
    from cmtj import Layer, Junction, CVector, constantDriver, AxialDriver
    from cmtj.utils import FieldScan
    
    Ms = 1.0
    Ku = 400e3
    damping = 0.01
    
    layer = Layer(
        "free",
        mag=CVector(0, 0, 1),
        anis=CVector(0, 0, 1),
        Ms=Ms,
        thickness=1.5e-9,
        damping=damping,
        cellSurface=100e-18,
        demagTensor=[CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)],
    )
    layer.setAnisotropyDriver(constantDriver(Ku))
    
    junction = Junction([layer])
    
    # Field sweep
    Hmax = 500e3
    Hscan, Hvecs = FieldScan.amplitude_scan(
        start=-Hmax,
        stop=Hmax,
        steps=50,
        theta=0,
        phi=0,
    )
    
    sim_time = 5e-9
    dt = 1e-13
    write_freq = 1e-11
    
    for Hv in Hvecs:
        junction.clearLog()
        junction.setLayerExternalFieldDriver("all", AxialDriver(*Hv))
        junction.runSimulation(sim_time, dt, write_freq)


def benchmark_tensor_operations():
    """Benchmark simulation with heavy dipole tensor interactions"""
    from cmtj import Layer, Junction, CVector, constantDriver
    
    Ms1 = 1.0
    Ms2 = 1.0
    damping = 0.01
    Ku = 400e3
    
    # Define dipole tensors for strong interlayer interactions
    dipole_tensor = [
        CVector(0.1, 0.05, 0.02),
        CVector(0.05, 0.1, 0.02),
        CVector(0.02, 0.02, 0.05)
    ]
    
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
    
    l1 = Layer(
        "layer1",
        mag=CVector(0.1, 0.1, 0.99),
        anis=CVector(0, 0, 1),
        Ms=Ms1,
        thickness=1.0e-9,
        damping=damping,
        cellSurface=1e-9,
        demagTensor=demag,
    )
    
    l2 = Layer(
        "layer2",
        mag=CVector(0, 0, 1),
        anis=CVector(0, 0, 1),
        Ms=Ms2,
        thickness=1.0e-9,
        damping=damping,
        cellSurface=1e-9,
        demagTensor=demag,
    )
    
    l1.setAnisotropyDriver(constantDriver(Ku))
    l2.setAnisotropyDriver(constantDriver(Ku))
    l1.setTopDipoleTensor(dipole_tensor)
    l2.setBottomDipoleTensor(dipole_tensor)
    
    junction = Junction([l1, l2])
    
    # Run simulation with many integration steps
    sim_time = 10e-9
    dt = 5e-14  # Smaller timestep for more calculations
    write_freq = 1e-11
    
    junction.runSimulation(sim_time, dt, write_freq)


def build_and_install(git_ref: str = None) -> bool:
    """Build and install CMTJ, optionally from a specific git reference"""
    print(f"Building and installing CMTJ{f' from {git_ref}' if git_ref else ''}...")
    
    workspace = Path("/workspace")
    original_ref = None
    
    try:
        if git_ref:
            # Save current git state
            result = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=workspace,
                capture_output=True,
                text=True,
                check=True
            )
            original_ref = result.stdout.strip()
            
            # Checkout the specified reference
            subprocess.run(
                ["git", "checkout", git_ref],
                cwd=workspace,
                capture_output=True,
                check=True
            )
            subprocess.run(
                ["git", "submodule", "update", "--init", "--recursive"],
                cwd=workspace,
                capture_output=True,
                check=True
            )
        
        # Uninstall existing cmtj
        subprocess.run(
            [sys.executable, "-m", "pip", "uninstall", "-y", "cmtj"],
            capture_output=True
        )
        
        # Build and install
        result = subprocess.run(
            [sys.executable, "-m", "pip", "install", "-e", ".", "-v"],
            cwd=workspace,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"Build failed: {result.stderr}")
            return False
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error during build: {e}")
        return False
    finally:
        if original_ref and git_ref:
            # Restore original git state
            subprocess.run(
                ["git", "checkout", original_ref],
                cwd=workspace,
                capture_output=True
            )
            subprocess.run(
                ["git", "submodule", "update", "--init", "--recursive"],
                cwd=workspace,
                capture_output=True
            )


def run_all_benchmarks() -> dict:
    """Run all benchmarks and return results"""
    benchmark_funcs = {
        "single_layer_rk4": benchmark_single_layer_rk4,
        "multi_layer_rk4": benchmark_multi_layer_rk4,
        "dormand_prince": benchmark_dormand_prince,
        "field_sweep": benchmark_field_sweep,
        "tensor_operations": benchmark_tensor_operations,
    }
    
    results = {}
    
    for name, config in BENCHMARKS.items():
        print(f"Running benchmark: {config['desc']}...")
        func = benchmark_funcs[name]
        result = run_benchmark(
            name,
            func,
            iterations=config["iterations"],
            warmup=config["warmup"]
        )
        results[name] = result
        print(f"  Mean: {result['mean']:.4f}s ± {result['std']:.4f}s")
    
    return results


def print_comparison(baseline: dict, optimized: dict):
    """Print comparison table"""
    print("\n" + "="*80)
    print("PERFORMANCE COMPARISON")
    print("="*80)
    print(f"{'Benchmark':<30} {'Baseline (s)':<15} {'Optimized (s)':<15} {'Speedup':<10}")
    print("-"*80)
    
    for name in baseline.keys():
        base_mean = baseline[name]["mean"]
        opt_mean = optimized[name]["mean"]
        speedup = base_mean / opt_mean
        
        print(f"{BENCHMARKS[name]['desc']:<30} "
              f"{base_mean:>7.4f} ± {baseline[name]['std']:>5.4f}  "
              f"{opt_mean:>7.4f} ± {optimized[name]['std']:>5.4f}  "
              f"{speedup:>6.2f}x")
    
    print("-"*80)
    
    # Calculate overall speedup
    total_baseline = sum(r["mean"] for r in baseline.values())
    total_optimized = sum(r["mean"] for r in optimized.values())
    overall_speedup = total_baseline / total_optimized
    
    print(f"{'OVERALL':<30} {total_baseline:>14.4f}  {total_optimized:>14.4f}  {overall_speedup:>6.2f}x")
    print("="*80)


def plot_comparison(baseline: dict, optimized: dict, output_file: str = "benchmark_results.png"):
    """Create visualization of benchmark results"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    names = [BENCHMARKS[k]["desc"] for k in baseline.keys()]
    x = np.arange(len(names))
    width = 0.35
    
    baseline_means = [baseline[k]["mean"] for k in baseline.keys()]
    baseline_stds = [baseline[k]["std"] for k in baseline.keys()]
    optimized_means = [optimized[k]["mean"] for k in optimized.keys()]
    optimized_stds = [optimized[k]["std"] for k in optimized.keys()]
    
    # Plot 1: Execution times
    ax1.bar(x - width/2, baseline_means, width, label='Baseline', 
            yerr=baseline_stds, capsize=5, alpha=0.8)
    ax1.bar(x + width/2, optimized_means, width, label='Optimized',
            yerr=optimized_stds, capsize=5, alpha=0.8)
    
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Benchmark Execution Times')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Speedup factors
    speedups = [baseline_means[i] / optimized_means[i] for i in range(len(names))]
    colors = ['green' if s > 1 else 'red' for s in speedups]
    
    ax2.barh(x, speedups, color=colors, alpha=0.7)
    ax2.axvline(x=1, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax2.set_xlabel('Speedup Factor')
    ax2.set_title('Performance Improvement')
    ax2.set_yticks(x)
    ax2.set_yticklabels(names)
    ax2.grid(axis='x', alpha=0.3)
    
    # Add speedup values as text
    for i, v in enumerate(speedups):
        ax2.text(v + 0.05, i, f'{v:.2f}x', va='center')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")


def save_results(results: dict, filename: str):
    """Save benchmark results to JSON file"""
    # Convert numpy types to Python types for JSON serialization
    serializable_results = {}
    for name, data in results.items():
        serializable_results[name] = {
            "name": data["name"],
            "times": [float(t) for t in data["times"]],
            "mean": float(data["mean"]),
            "std": float(data["std"]),
            "min": float(data["min"]),
            "max": float(data["max"]),
            "median": float(data["median"]),
        }
    
    with open(filename, 'w') as f:
        json.dump(serializable_results, f, indent=2)
    print(f"Results saved to: {filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark CMTJ performance and optionally compare with previous version"
    )
    parser.add_argument(
        "--compare",
        type=str,
        metavar="GIT_REF",
        help="Git reference to compare against (e.g., HEAD~1, commit hash)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="benchmark_results",
        help="Output filename prefix (default: benchmark_results)"
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="Skip rebuilding (use already installed version)"
    )
    
    args = parser.parse_args()
    
    if args.compare:
        print("="*80)
        print("COMPARISON MODE: Benchmarking baseline vs optimized")
        print("="*80)
        
        # Build and benchmark baseline
        if not args.no_build:
            if not build_and_install(args.compare):
                print("Failed to build baseline version")
                return 1
        
        print("\n--- Running BASELINE benchmarks ---")
        baseline_results = run_all_benchmarks()
        save_results(baseline_results, f"{args.output}_baseline.json")
        
        # Build and benchmark optimized
        if not args.no_build:
            if not build_and_install():
                print("Failed to build optimized version")
                return 1
        
        print("\n--- Running OPTIMIZED benchmarks ---")
        optimized_results = run_all_benchmarks()
        save_results(optimized_results, f"{args.output}_optimized.json")
        
        # Compare and visualize
        print_comparison(baseline_results, optimized_results)
        plot_comparison(baseline_results, optimized_results, f"{args.output}.png")
        
    else:
        print("="*80)
        print("BENCHMARK MODE: Testing current version")
        print("="*80)
        
        results = run_all_benchmarks()
        save_results(results, f"{args.output}.json")
        
        print("\nBenchmark complete! To compare with a previous version, run:")
        print(f"  python {sys.argv[0]} --compare HEAD~1")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
