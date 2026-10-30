#!/usr/bin/env python
"""Numerical sanity check of every visible GPU against the CPU (fp64 matmul + int64 bit ops).

Motivation: on 2026-10-29 GPU 0 of this node returned wrong results (matmul rel. error ~1e-2, corrupted
int64 copies) and burrow's GPU kernels hit device-side asserts on it, while GPUs 1-3 were correct.
Run before any GPU job; exit code 1 if any GPU fails.
"""
import sys

import torch

N = 4096            # matrix size
N_INT = 10_000_000  # int64 elements
TOL = 1e-10         # relative tolerance for the fp64 matmul sum

torch.manual_seed(0)
a = torch.randn(N, N, dtype=torch.float64)
ref = (a @ a).sum().item()
x = torch.randint(0, 1 << 40, (N_INT,), dtype=torch.int64)
bad = []
for d in range(torch.cuda.device_count()):
    g = a.to(f"cuda:{d}")
    rel = abs((g @ g).sum().item() - ref) / abs(ref)
    ok_int = bool(((x.to(f"cuda:{d}") & 0xFFFFFFFF).cpu() == (x & 0xFFFFFFFF)).all())
    ok = rel < TOL and ok_int
    print(f"cuda:{d} {torch.cuda.get_device_name(d)} matmul_relerr={rel:.2e} int64_ok={ok_int} -> {'OK' if ok else 'FAIL'}")
    if not ok:
        bad.append(d)
if bad:
    print("FAILING GPUs:", bad, "-> exclude them via CUDA_VISIBLE_DEVICES")
sys.exit(1 if bad else 0)
