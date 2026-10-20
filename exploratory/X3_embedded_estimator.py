#!/usr/bin/env python3
"""EXPLORATORY (non-rigorous, float64): would adaptive step control see the cycle?

Uses the embedded first-order solution of MPRK22 used in practice
(Kopecz & Meister 2018; PositiveIntegrators.jl): sigma = y^(1-1/alpha) u^(1/alpha),
error estimate y+ - sigma, scaled RMS norm with abstol/reltol (defaults of
OrdinaryDiffEq.jl: 1e-6 / 1e-3).  Reports the scaled error on the steps of the
certified two-cycles and runs a standard adaptive controller from the
counterexample's initial state.  Writes X3_embedded_estimator.json.
"""
import json
from pathlib import Path
import numpy as np
from pds_schemes import mprk22, linear_P

HERE = Path(__file__).resolve().parent


def scaled_err(y, yp, sig, atol, rtol):
    sc = atol + rtol * np.maximum(np.abs(y), np.abs(yp))
    return float(np.sqrt(np.mean(((yp - sig) / sc) ** 2)))


def cycle_points(Pf, y0, h, alpha, n=4000):
    y = y0.copy()
    for _ in range(n):
        y = mprk22(Pf, y, h, alpha)
    return y, mprk22(Pf, y, h, alpha)


def adaptive(Pf, y0, T, alpha, dt0, atol=1e-6, rtol=1e-3):
    t, y, dt = 0.0, y0.copy(), dt0
    acc = rej = 0
    max_dt = 0.0
    while t < T and acc + rej < 100000:
        dt = min(dt, T - t)
        yp, sig = mprk22(Pf, y, dt, alpha, embedded=True)
        e = scaled_err(y, yp, sig, atol, rtol)
        fac = 0.9 * max(e, 1e-10) ** -0.5
        if e <= 1:
            t, y = t + dt, yp
            acc += 1
            max_dt = max(max_dt, dt)
            dt *= min(5.0, max(0.2, fac))
        else:
            rej += 1
            dt *= min(1.0, max(0.2, fac))
    return {"accepted": acc, "rejected": rej, "max_dt": max_dt, "final_state": y.tolist(), "t_end": t}


def main():
    out = {}
    alpha = 4 / 7
    systems = {
        "two_state_witness": (np.array([[-9990.0, 10], [9990, -10]]), np.array([0.5, 0.5]),
                              np.array([0.001, 0.999])),
        "three_species_chain": (np.array([[-10.0, 10000, 0], [10, -20000, 1], [0, 10000, -1]]),
                                np.full(3, 1 / 3), np.array([1000, 1, 10000]) / 11001),
    }
    for name, (A, y0, pi) in systems.items():
        Pf = linear_P(A)
        ya, yb = cycle_points(Pf, y0, 1.0, alpha)
        rec = {"cycle_points": [ya.tolist(), yb.tolist()]}
        for tol in [(1e-6, 1e-3), (1e-3, 1e-1), (1e-2, 1.0)]:
            errs = []
            for y in (ya, yb):
                yp, sig = mprk22(Pf, y, 1.0, alpha, embedded=True)
                errs.append(scaled_err(y, yp, sig, *tol))
            rec["scaled_error_on_cycle_steps_atol_%g_rtol_%g" % tol] = errs
        yp, sig = mprk22(Pf, y0, 1.0, alpha, embedded=True)
        rec["scaled_error_first_step_default_tol"] = scaled_err(y0, yp, sig, 1e-6, 1e-3)
        for dt0 in [1.0, 1e-4]:
            r = adaptive(Pf, y0, 1e6, alpha, dt0)
            r["final_error_vs_equilibrium"] = float(np.max(np.abs(np.array(r["final_state"]) - pi)))
            rec["adaptive_run_T_1e6_dt0_%g" % dt0] = r
        out[name] = rec
        print(name, json.dumps({k: v for k, v in rec.items() if "scaled" in k}, indent=1))
        print("  adaptive:", {k: (v["accepted"], v["rejected"], v["max_dt"], v["final_error_vs_equilibrium"])
                              for k, v in rec.items() if k.startswith("adaptive")})
    (HERE / "X3_embedded_estimator.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
