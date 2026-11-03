"""
Analytical vs numerical response of single-degree-of-freedom vibrating systems.

Mechanical Vibrations (EML 6223), Florida Atlantic University, Spring 2025, Project 1.
Problem parameters follow the course's serial-number scheme with S = 3.

Each problem is solved two ways:
  * exactly, in closed form (piecewise for step, pulse and ramp forcing)
  * numerically, with a hand-written fixed-step 4th-order Runge-Kutta integrator
and the two are compared. Run:  python vibrations.py   (plots go to plots/, numbers to results.md)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = 3
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plots")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
lines = []


# ---------------------------------------------------------------- numerical
def rk4(m, c, k, F, x0, v0, t_end, dt=1e-3):
    """Classic 4th-order Runge-Kutta for m x'' + c x' + k x = F(t).

    Step, pulse and ramp forces jump or kink at times that sit on the 1 ms grid. So that RK4 never
    samples the force on the far side of a jump, each step uses the force just after its start and
    just before its end (one-sided limits); between jumps this is the ordinary RK4 scheme."""
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    y = np.zeros((n + 1, 2))
    y[0] = x0, v0
    eps = 1e-9 * dt
    f = lambda Fv, yi: np.array([yi[1], (Fv - c * yi[1] - k * yi[0]) / m])
    for i in range(n):
        ti, yi = t[i], y[i]
        k1 = f(F(ti + eps), yi)
        k2 = f(F(ti + dt / 2), yi + dt / 2 * k1)
        k3 = f(F(ti + dt / 2), yi + dt / 2 * k2)
        k4 = f(F(ti + dt - eps), yi + dt * k3)
        y[i + 1] = yi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return t, y[:, 0]


# ---------------------------------------------------------------- exact
def free(m, c, k, xh0, vh0, tau):
    """Homogeneous response with xh(0) = xh0, xh'(0) = vh0 (any damping level)."""
    wn = np.sqrt(k / m)
    z = c / (2 * np.sqrt(m * k))
    if z < 1:
        wd = wn * np.sqrt(1 - z ** 2)
        s = z * wn
        return np.exp(-s * tau) * (xh0 * np.cos(wd * tau) + (vh0 + s * xh0) / wd * np.sin(wd * tau))
    if z > 1:
        r1 = -z * wn + wn * np.sqrt(z ** 2 - 1)
        r2 = -z * wn - wn * np.sqrt(z ** 2 - 1)
        C2 = (vh0 - r1 * xh0) / (r2 - r1)
        return (xh0 - C2) * np.exp(r1 * tau) + C2 * np.exp(r2 * tau)
    return (xh0 + (vh0 + wn * xh0) * tau) * np.exp(-wn * tau)


def free_v(m, c, k, xh0, vh0, tau, h=1e-7):
    return (free(m, c, k, xh0, vh0, tau + h) - free(m, c, k, xh0, vh0, tau - h)) / (2 * h)


def exact_piecewise_linear(m, c, k, segments, x0, v0, t):
    """Exact response to piecewise-linear forcing.
    segments: list of (t_start, a, b) meaning F = a + b*(t - t_start) from t_start on."""
    x = np.zeros_like(t)
    xs, vs = x0, v0
    bounds = [s[0] for s in segments] + [np.inf]
    for i, (ts, a, b) in enumerate(segments):
        te = bounds[i + 1]
        B = b / k                      # particular solution xp = A + B*tau
        A = (a - c * B) / k
        mask = (t >= ts) & (t < te)
        tau = t[mask] - ts
        x[mask] = A + B * tau + free(m, c, k, xs - A, vs - B, tau)
        if np.isfinite(te):            # state at the end of this segment
            T = te - ts
            xs, vs = A + B * T + free(m, c, k, xs - A, vs - B, T), B + free_v(m, c, k, xs - A, vs - B, T)
    return x


def exact_harmonic(m, c, k, F0, w, x0, v0, t):
    """Exact response to F0 cos(w t): steady state plus transient fitted to the initial conditions."""
    X = F0 / np.sqrt((k - m * w ** 2) ** 2 + (c * w) ** 2)
    phi = np.arctan2(c * w, k - m * w ** 2)
    xp = X * np.cos(w * t - phi)
    xp0, vp0 = X * np.cos(-phi), -X * w * np.sin(-phi)
    return xp + free(m, c, k, x0 - xp0, v0 - vp0, t), X, phi


def compare(name, t, xe, xn, title, ylab="x [m]", extra=None):
    err = np.max(np.abs(xe - xn))
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax[0].plot(t, xe, "k-", lw=2.2, label="Exact (analytical)")
    ax[0].plot(t[::40], xn[::40], "o", ms=3.5, color="#d9541a", label="Runge-Kutta 4 (dt = 1 ms)")
    if extra:
        extra(ax[0])
    ax[0].set_ylabel(ylab); ax[0].set_title(title); ax[0].legend(loc="upper right", fontsize=9)
    ax[1].plot(t, xe - xn, color="#3366cc"); ax[1].set_ylabel("exact − RK4 [m]"); ax[1].set_xlabel("t [s]")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), dpi=130); plt.close(fig)
    return err


# ================================================================ Problem 1
m, k, c = 100.0, S * 100 + 650.0, S * 10 + 85.0
t = np.arange(0, 15 + 1e-9, 1e-3)
xe, X, phi = exact_harmonic(m, c, k, 150.0, 10.0, 0.01, 0.1, t)
_, xn = rk4(m, c, k, lambda s: 150 * np.cos(10 * s), 0.01, 0.1, 15)
e1 = compare("p1_harmonic", t, xe, xn, f"Problem 1: harmonic forcing, m = {m:.0f} kg, k = {k:.0f} N/m, c = {c:.0f} kg/s")
lines.append(f"| 1 | Harmonic forcing 150 cos(10t) | steady-state amplitude {X*1000:.2f} mm, phase {phi:.3f} rad; peak displacement {np.max(np.abs(xe))*1000:.1f} mm | {e1:.1e} m |")

# ================================================================ Problem 2
m, k, c0, F0, w = 100.0, S * 10 + 255.0, 10.0, 150.0, 5.0
limit, design_limit = S * 0.005 + 0.15, 0.2
t = np.arange(0, 40 + 1e-9, 1e-3)
xe2, X2, _ = exact_harmonic(m, c0, k, F0, w, 0.01, 0.5, t)
_, xn2 = rk4(m, c0, k, lambda s: F0 * np.cos(w * s), 0.01, 0.5, 40)
peak = np.max(np.abs(xe2))
def peak_for(cc):
    return np.max(np.abs(exact_harmonic(m, cc, k, F0, w, 0.01, 0.5, t)[0]))
def min_damping(target):
    lo, hi = c0, 2000.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if peak_for(mid) < target else (mid, hi)
    return hi
c_lim, c_des = min_damping(limit), min_damping(design_limit)
xe2b = exact_harmonic(m, c_des, k, F0, w, 0.01, 0.5, t)[0]
def p2extra(a):
    a.plot(t, xe2b, "-", color="#1f9a3a", lw=1.6, label=f"with c = {c_des:.1f} kg/s")
    for y in (limit, -limit):
        a.axhline(y, ls="--", color="crimson", lw=1)
    for y in (design_limit, -design_limit):
        a.axhline(y, ls=":", color="grey", lw=1)
    a.text(t[0], limit, f" spec limit ±{limit:.3f} m", va="bottom", ha="left", color="crimson", fontsize=9)
    a.text(t[0], -design_limit, f" redesign target ±{design_limit:.1f} m", va="top", ha="left", color="grey", fontsize=9)
e2 = compare("p2_amplitude_limit", t, xe2, xn2, f"Problem 2: electronics module, c = {c0:.0f} kg/s, limit {limit:.3f} m", extra=p2extra)
lines.append(f"| 2 | Amplitude limit check | steady state {X2*1000:.1f} mm, but transient peak {peak*1000:.0f} mm > {limit*1000:.0f} mm limit; "
             f"min damping {c_lim:.1f} kg/s for {limit:.3f} m, {c_des:.1f} kg/s for {design_limit:.1f} m | {e2:.1e} m |")

# ================================================================ Problems 3 and 4 (step and square pulse)
k, z, wn, F0 = S * 4 + 925.0, S * 0.01 + 0.05, 3.16, 30.0
m = k / wn ** 2
c = 2 * z * np.sqrt(m * k)
t0 = float(S)
t = np.arange(0, 20 + 1e-9, 1e-3)
xe3 = exact_piecewise_linear(m, c, k, [(0.0, 0.0, 0.0), (t0, F0, 0.0)], 0, 0, t)
_, xn3 = rk4(m, c, k, lambda s: F0 if s >= t0 else 0.0, 0, 0, 20)
e3 = compare("p3_step", t, xe3, xn3, f"Problem 3: step of {F0:.0f} N at t0 = {t0:.0f} s (zeta = {z:.2f}, wn = {wn} rad/s)")
xs = F0 / k
lines.append(f"| 3 | Step input at t0 = {t0:.0f} s | static deflection {xs*1000:.2f} mm, peak {np.max(xe3)*1000:.2f} mm "
             f"(overshoot {100*(np.max(xe3)/xs-1):.0f}%, theory {100*np.exp(-np.pi*z/np.sqrt(1-z*z)):.0f}%) | {e3:.1e} m |")

res4 = []
for t1, tag in ((1.5, "a"), (100.0, "b")):
    T = 20.0 if t1 < 10 else 120.0
    tt = np.arange(0, T + 1e-9, 1e-3)
    xe4 = exact_piecewise_linear(m, c, k, [(0.0, F0, 0.0), (t1, 0.0, 0.0)], 0, 0, tt)
    _, xn4 = rk4(m, c, k, lambda s, t1=t1: F0 if s < t1 else 0.0, 0, 0, T)
    e4 = compare(f"p4{tag}_pulse_{t1:g}s", tt, xe4, xn4, f"Problem 4: square pulse {F0:.0f} N lasting t1 = {t1:g} s")
    res4.append((t1, np.max(np.abs(xe4)), e4))
lines.append("| 4 | Square pulse | " + "; ".join(f"t1 = {a:g} s: peak {b*1000:.2f} mm" for a, b, _ in res4)
             + f" | {max(r[2] for r in res4):.1e} m |")

# ================================================================ Problem 6 (triangular pulse, overdamped)
m, c, k = float(S), 150.0 + S, 750.0
F0, t1, t2 = S + 2000.0, 3.0, 6.0
t = np.arange(0, 10 + 1e-9, 1e-3)
segs = [(0.0, 0.0, F0 / t1), (t1, F0, -F0 / (t2 - t1)), (t2, 0.0, 0.0)]
xe6 = exact_piecewise_linear(m, c, k, segs, 0, 0, t)
Ff = lambda s: F0 * s / t1 if s < t1 else (F0 * (t2 - s) / (t2 - t1) if s < t2 else 0.0)
_, xn6 = rk4(m, c, k, Ff, 0, 0, 10)
e6 = compare("p6_triangular_pulse", t, xe6, xn6, f"Problem 6: triangular pulse, peak {F0:.0f} N (m = {m:.0f} kg, c = {c:.0f} N s/m, k = {k:.0f} N/m)")
z6 = c / (2 * np.sqrt(m * k))
lines.append(f"| 6 | Triangular pulse (overdamped, zeta = {z6:.2f}) | peak {np.max(xe6):.3f} m at t = {t[np.argmax(xe6)]:.2f} s "
             f"(quasi-static F0/k = {F0/k:.3f} m) | {e6:.1e} m |")

with open(os.path.join(HERE, "results.md"), "w") as fh:
    fh.write("| Problem | Case | Result (exact solution) | Max abs. difference, exact vs RK4 |\n| --- | --- | --- | --- |\n")
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
