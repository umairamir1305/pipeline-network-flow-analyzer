# Pipeline Network Flow Analyzer — MATLAB

A multi-physics fluid mechanics simulation suite for a 7-pipe, 2-loop water distribution network. Implements industry-standard methods including the **Hardy-Cross algorithm**, **Colebrook-White friction factors**, **pump curve intersection**, **NPSH cavitation assessment**, and **boundary layer development analysis**.

---

## What This Project Does

Given a pipe network with known geometry and demands, this suite:

1. **Solves flow distribution** across all pipes simultaneously using the Hardy-Cross iterative method
2. **Computes friction factors** using the implicit Colebrook-White equation (industry standard — used in every real pipe network software)
3. **Generates a full Moody diagram** and overlays the actual operating points of every pipe
4. **Selects between three pump options**, finds each operating point analytically, and flags cavitation risk using NPSH analysis
5. **Analyses boundary layer development** in the pipe entry region and plots velocity profile evolution
6. **Runs a parametric study** showing how changing the main pipe diameter affects Re, friction factor, head loss, and velocity

---

## Network Layout

```
    [1]────P1────[2]────P2────[3]
     │            │            │
     P6           P4           P3
     │            │            │
    [5]────P5────[4]───────────┘

7 Pipes | 5 Nodes | 2 Loops
Fluid: Water at 20°C
```

| Pipe | Route | L (m) | D (mm) | e (mm) |
|------|-------|-------|--------|--------|
| P1 | Node 1→2 | 500 | 250 | 0.15 |
| P2 | Node 2→3 | 400 | 200 | 0.15 |
| P3 | Node 3→4 | 350 | 150 | 0.26 |
| P4 | Node 2→4 | 300 | 200 | 0.15 |
| P5 | Node 4→5 | 450 | 200 | 0.15 |
| P6 | Node 1→5 | 380 | 180 | 0.26 |
| P7 | Node 2→4 | 320 | 150 | 0.26 |

---

## Scripts

| File | Description | Key Methods |
|------|-------------|-------------|
| `script1_hardy_cross_network.m` | Network flow solver | Hardy-Cross, Colebrook-White, Darcy-Weisbach |
| `script2_moody_diagram.m` | Full Moody diagram | Prandtl smooth-pipe law, Colebrook, operating points |
| `script3_pump_cavitation.m` | Pump selection & cavitation | System/pump curve intersection, NPSH, fzero |
| `script4_boundary_layer_parametric.m` | BL development + parametric | Entry length, velocity profiles, diameter sweep |
| `RUN_ALL.m` | Master runner | Runs all scripts in sequence |

---

## How to Run

```matlab
% In MATLAB Command Window:
cd('path/to/pipeline_network')
RUN_ALL
```

Requires: **MATLAB Optimization Toolbox** (for `fsolve` and `fzero`) — available via NUST campus license.

Scripts can also be run individually. Script 1 must run first as it generates `network_results.mat` used by Scripts 2–4.

---

## Theory

### Hardy-Cross Method
Iteratively corrects assumed flow rates in each loop until head loss balance is satisfied:

```
ΔQ = -ΣhF / Σ(dhF/dQ)

where hF = f·(L/D)·(V²/2g)   [Darcy-Weisbach]
```

### Colebrook-White Equation (implicit)
```
1/√f = -2·log₁₀(ε/3.7D + 2.51/(Re·√f))
```
Solved iteratively using `fsolve` at each pipe for each iteration.

### Pump Operating Point
Intersection of pump curve H(Q) and system curve:
```
H_system(Q) = H_static + R·Q²
H_pump(Q)   = H₀ - a·Q - b·Q²
Operating point: H_system = H_pump  →  solved via fzero
```

### Cavitation — NPSH
```
NPSH_available = (p_atm - p_vapour)/(ρg) + Vs²/2g - z_s
Cavitation if: NPSH_available < NPSH_required
```

### Boundary Layer Entry Length
```
Laminar:    L_h = 0.06·Re·D
Turbulent:  L_h = 4.4·Re^(1/6)·D
```

---

## Output Plots (~12 total)

- Converged flow rates per pipe (colour-coded by regime)
- Reynolds numbers per pipe with laminar/turbulent boundaries
- Head losses per pipe (Darcy-Weisbach)
- Hardy-Cross convergence history (log scale)
- Full Moody diagram with network operating points overlaid
- Pump curves vs system curve (3 pumps, operating points marked)
- NPSH available vs required — cavitation onset marked per pump
- Velocity profile development (5 stations along entry length)
- Wall shear stress distribution along pipe entry
- 4-panel parametric study: Re, f, hf, V vs pipe diameter

---

## Skills Demonstrated

- Numerical methods — iterative solution of nonlinear equation systems
- Fluid mechanics — pipe flow, friction factors, head loss, flow regimes
- Turbomachinery — pump curves, system curves, operating point analysis
- Cavitation theory — NPSH calculation and risk assessment
- Boundary layer theory — entry length, developing velocity profiles
- MATLAB — `fsolve`, `fzero`, data structures, multi-figure plotting
- Engineering judgment — parametric design studies and interpretation

---

## Author

**Umair**
B.E. Mechanical Engineering — 4th Semester, NUST Karachi

---

## References

1. Munson, Young, Okiishi — *Fundamentals of Fluid Mechanics*
2. White, F.M. — *Fluid Mechanics*, 8th Edition
3. Moody, L.F. — *Friction Factors for Pipe Flow*, ASME Trans. 1944
4. Hardy Cross — *Analysis of Flow in Networks of Conduits*, 1936
