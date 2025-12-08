# CMTJ Benchmark Results Table

!Recommendation!: 
Install from source, and compile with the following flags:

```bash
export LDFLAGS="-O3"
export CXXFLAGS="-O3 -march=native -ffast-math"
```


| Version | Single layer RK4 (perf) | Single layer RK4 (quick) | Multi-layer RK4 (perf) | Multi-layer RK4 (quick) | Dormand-Prince (perf) | Dormand-Prince (quick) | Field sweep (50 steps) (perf) | Field sweep (50 steps) (quick) | Tensor operations (perf) | Tensor operations (quick) | Total benchmark time (perf) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.11.0 (-O3 -march=native -ffast-math) | 0.4793 s ± 0.0271 s | 0.0573 s ± 0.0146 s | 0.3105 s ± 0.0310 s | 0.0375 s ± 0.0141 s | 0.0001 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.2400 s ± 0.0545 s | 1.2081 s ± 0.0326 s | 2.0131 s ± 0.1256 s | 0.1989 s ± 0.0159 s | 4.0430 s |
| 1.11.0 (-O3 -march=native) | 0.5307 s ± 0.0169 s | 0.0846 s ± 0.0321 s | 0.3361 s ± 0.0318 s | 0.0326 s ± 0.0032 s | 0.0001 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.4698 s ± 0.1361 s | 1.3227 s ± 0.0342 s | 2.2090 s ± 0.0996 s | 0.2101 s ± 0.0062 s | 4.5457 s |
| 1.11.0 (02) | 0.7935 s ± 0.0907 s | 0.0790 s ± 0.0230 s | 0.4096 s ± 0.0339 s | 0.0379 s ± 0.0009 s | 0.0002 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.9164 s ± 0.0204 s | 1.7432 s ± 0.1066 s | 2.5543 s ± 0.0314 s | 0.2610 s ± 0.0139 s | 5.6739 s |
| 1.11.0 (Normal) | 0.5528 s ± 0.0273 s | 0.0594 s ± 0.0079 s | 0.3410 s ± 0.0256 s | 0.0357 s ± 0.0062 s | 0.0001 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.3679 s ± 0.0267 s | 1.3619 s ± 0.0665 s | 2.2183 s ± 0.0580 s | 0.2174 s ± 0.0102 s | 4.4801 s |
| 1.10.0 | 0.7144 s ± 0.0544 s | 0.2049 s ± 0.0741 s | 0.4148 s ± 0.0266 s | 0.1206 s ± 0.0303 s | 0.0001 s ± 0.0000 s | 0.0002 s ± 0.0007 s | 1.6580 s ± 0.0510 s | 1.6958 s ± 0.0452 s | 3.0100 s ± 0.3041 s | 0.2948 s ± 0.0182 s | 5.7974 s |

