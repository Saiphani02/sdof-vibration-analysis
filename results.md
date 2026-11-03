| Problem | Case | Result (exact solution) | Max abs. difference, exact vs RK4 |
| --- | --- | --- | --- |
| 1 | Harmonic forcing 150 cos(10t) | steady-state amplitude 16.44 mm, phase 3.015 rad; peak displacement 50.9 mm | 4.4e-13 m |
| 2 | Amplitude limit check | steady state 67.7 mm, but transient peak 355 mm > 165 mm limit; min damping 349.5 kg/s for 0.165 m, 239.7 kg/s for 0.2 m | 3.6e-13 m |
| 3 | Step input at t0 = 3 s | static deflection 32.02 mm, peak 56.90 mm (overshoot 78%, theory 78%) | 1.2e-13 m |
| 4 | Square pulse | t1 = 1.5 s: peak 56.90 mm; t1 = 100 s: peak 56.90 mm | 2.2e-11 m |
| 6 | Triangular pulse (overdamped, zeta = 1.61) | peak 2.557 m at t = 3.15 s (quasi-static F0/k = 2.671 m) | 7.3e-11 m |
