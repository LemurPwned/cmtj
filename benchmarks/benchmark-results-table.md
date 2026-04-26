# CMTJ Benchmark Results Table

!Recommendation!: 
Install from source, and compile with the following flags:

```bash
export LDFLAGS="-O3"
export CXXFLAGS="-O3 -march=native -ffast-math"
```

Machine notes:

- `Local ARM64 Linux`: `aarch64`, Ubuntu 26.04, Cortex-X925/A725
- `GitHub Codespaces`: historical measurements taken in Codespaces
- `Not recorded`: machine information was not captured for that run

Relative boost notes:

- `Relative boost vs local 1.7.0` is computed from `Total benchmark time (quick)`
- `Relative boost vs previous local` compares each local row with the next older
  local version in the table
- Both are only shown for rows measured on `Local ARM64 Linux`
- Baseline formula: `(v1.7.0_quick_total / row_quick_total - 1) * 100%`
- Previous-version formula: `(previous_local_quick_total / row_quick_total - 1) * 100%`

| Version | Machine | Single layer RK4 (perf) | Single layer RK4 (quick) | Multi-layer RK4 (perf) | Multi-layer RK4 (quick) | Dormand-Prince (perf) | Dormand-Prince (quick) | Field sweep (50 steps) (perf) | Field sweep (50 steps) (quick) | Tensor operations (perf) | Tensor operations (quick) | Total benchmark time (quick) | Relative boost vs local 1.7.0 | Relative boost vs previous local |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.13.0 | Local ARM64 Linux | 0.0134 s ± 0.0001 s | 0.1320 s ± 0.0001 s | 0.0082 s ± 0.0000 s | 0.0834 s ± 0.0001 s | 0.0000 s ± 0.0000 s | 0.0001 s ± 0.0000 s | 0.3342 s ± 0.0003 s | 0.3305 s ± 0.0002 s | 0.0529 s ± 0.0001 s | 0.5315 s ± 0.0013 s | 1.0776 s | +180.3% | +129.3% |
| master (dd67d27) | Local ARM64 Linux | 0.0303 s ± 0.0002 s | 0.3039 s ± 0.0004 s | 0.0183 s ± 0.0001 s | 0.1868 s ± 0.0004 s | 0.0000 s ± 0.0000 s | 0.0001 s ± 0.0000 s | 0.7473 s ± 0.0025 s | 0.7517 s ± 0.0012 s | 0.1221 s ± 0.0003 s | 1.2290 s ± 0.0011 s | 2.4715 s | +22.2% | +19.4% |
| 1.10.0 | Local ARM64 Linux | 0.0370 s ± 0.0002 s | 0.3731 s ± 0.0008 s | 0.0217 s ± 0.0000 s | 0.2121 s ± 0.0002 s | 0.0000 s ± 0.0000 s | 0.0001 s ± 0.0000 s | 0.8740 s ± 0.0007 s | 0.8788 s ± 0.0010 s | 0.1488 s ± 0.0003 s | 1.4873 s ± 0.0083 s | 2.9514 s | +2.3% | +0.8% |
| 1.9.1 | Local ARM64 Linux | 0.0367 s ± 0.0003 s | 0.3696 s ± 0.0082 s | 0.0217 s ± 0.0003 s | 0.2129 s ± 0.0007 s | 0.0000 s ± 0.0000 s | 0.0001 s ± 0.0000 s | 0.8613 s ± 0.0017 s | 0.8618 s ± 0.0017 s | 0.1479 s ± 0.0011 s | 1.5300 s ± 0.0225 s | 2.9742 s | +1.6% | +0.3% |
| 1.8.0 | Local ARM64 Linux | 0.0369 s ± 0.0003 s | 0.3730 s ± 0.0004 s | 0.0219 s ± 0.0000 s | 0.2172 s ± 0.0009 s | 0.0000 s ± 0.0000 s | 0.0001 s ± 0.0000 s | 0.8744 s ± 0.0197 s | 0.8792 s ± 0.0003 s | 0.1488 s ± 0.0001 s | 1.5135 s ± 0.0029 s | 2.9829 s | +1.3% | +1.3% |
| 1.7.0 | Local ARM64 Linux | 0.0378 s ± 0.0001 s | 0.3776 s ± 0.0007 s | 0.0225 s ± 0.0001 s | 0.2189 s ± 0.0002 s | 0.0000 s ± 0.0000 s | 0.0001 s ± 0.0000 s | 0.8913 s ± 0.0014 s | 0.8921 s ± 0.0009 s | 0.1527 s ± 0.0002 s | 1.5319 s ± 0.0012 s | 3.0205 s | baseline | baseline |
| 1.11.0 (-O3 -march=native -ffast-math) | GitHub Codespaces | 0.4793 s ± 0.0271 s | 0.0573 s ± 0.0146 s | 0.3105 s ± 0.0310 s | 0.0375 s ± 0.0141 s | 0.0001 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.2400 s ± 0.0545 s | 1.2081 s ± 0.0326 s | 2.0131 s ± 0.1256 s | 0.1989 s ± 0.0159 s | 4.0430 s | n/a | n/a |
| 1.11.0 (-O3 -march=native) | GitHub Codespaces | 0.5307 s ± 0.0169 s | 0.0846 s ± 0.0321 s | 0.3361 s ± 0.0318 s | 0.0326 s ± 0.0032 s | 0.0001 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.4698 s ± 0.1361 s | 1.3227 s ± 0.0342 s | 2.2090 s ± 0.0996 s | 0.2101 s ± 0.0062 s | 4.5457 s | n/a | n/a |
| 1.11.0 (02) | GitHub Codespaces | 0.7935 s ± 0.0907 s | 0.0790 s ± 0.0230 s | 0.4096 s ± 0.0339 s | 0.0379 s ± 0.0009 s | 0.0002 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.9164 s ± 0.0204 s | 1.7432 s ± 0.1066 s | 2.5543 s ± 0.0314 s | 0.2610 s ± 0.0139 s | 5.6739 s | n/a | n/a |
| 1.11.0 (Normal) | GitHub Codespaces | 0.5528 s ± 0.0273 s | 0.0594 s ± 0.0079 s | 0.3410 s ± 0.0256 s | 0.0357 s ± 0.0062 s | 0.0001 s ± 0.0000 s | 0.0000 s ± 0.0000 s | 1.3679 s ± 0.0267 s | 1.3619 s ± 0.0665 s | 2.2183 s ± 0.0580 s | 0.2174 s ± 0.0102 s | 4.4801 s | n/a | n/a |
| 1.10.0 (historical) | Not recorded | 0.7144 s ± 0.0544 s | 0.2049 s ± 0.0741 s | 0.4148 s ± 0.0266 s | 0.1206 s ± 0.0303 s | 0.0001 s ± 0.0000 s | 0.0002 s ± 0.0007 s | 1.6580 s ± 0.0510 s | 1.6958 s ± 0.0452 s | 3.0100 s ± 0.3041 s | 0.2948 s ± 0.0182 s | 5.7974 s | n/a | n/a |
