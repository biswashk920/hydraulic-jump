# Theory guide

This guide explains every calculation in the program, in the order the program does them. If you can follow this page you can follow the code: each section names the module that implements it.

Symbols: $Q$ discharge (m³/s), $y$ depth, $A$ flow area, $T$ top width, $P$ wetted perimeter, $R = A/P$ hydraulic radius, $V = Q/A$ mean velocity, $S_0$ bed slope (positive when the bed falls downstream), $S_f$ friction slope, $n$ Manning's roughness, $g = 9.81$ m/s².

## 1. Froude number and the two kinds of flow  (`hydraulics.froude`)

$$Fr = \frac{V}{\sqrt{gA/T}}$$

$\sqrt{gA/T}$ is the speed of a small surface wave. If $Fr < 1$ the flow is **subcritical** (calm, deep; waves can travel upstream, so the flow "feels" what is downstream). If $Fr > 1$ it is **supercritical** (fast, shallow; waves are swept downstream, so it is controlled from upstream). At $Fr = 1$ the flow is **critical**.

*Consequence used everywhere in the program:* subcritical profiles are calculated **upstream** from a downstream control, supercritical profiles **downstream** from an upstream control.

## 2. Specific energy and critical depth  (`hydraulics.energy`, `critical_depth`, `alternate_depth`)

Energy per unit weight measured from the bed:

$$E = y + \frac{V^2}{2g} = y + \frac{Q^2}{2gA^2}$$

For a fixed $Q$, plot $E$ against $y$ (Figure 4 in the web app). The curve has a nose: a **minimum energy** $E_{min}$. Setting $dE/dy = 0$ gives $Fr = 1$, i.e. **critical depth** $y_c$:

$$\frac{Q^2 T}{g A^3} = 1 \qquad\text{(rectangle: } y_c = (q^2/g)^{1/3},\; E_{min} = \tfrac32 y_c,\; q = Q/b\text{)}$$

Any $E > E_{min}$ is carried at **two** depths, the *alternate depths* — one subcritical, one supercritical (the vertical line in Figure 4). The upper branch leans along the line $E = y$ because for large depth the velocity head is negligible.

## 3. Uniform flow and slope classes  (`hydraulics.normal_depth`, `profile._build_reaches`)

Manning's equation gives the depth of steady uniform flow, the **normal depth** $y_n$:

$$Q = \frac{1}{n} A R^{2/3} S_0^{1/2}$$

It is found with a root-finder (`brentq`). It exists only if $S_0 > 0$. Comparing $y_n$ with $y_c$ classifies the reach:

| Class | Condition | Uniform flow is |
|---|---|---|
| Mild (M) | $y_n > y_c$ | subcritical |
| Steep (S) | $y_n < y_c$ | supercritical |
| Critical (C) | $y_n = y_c$ (within 1 %) | critical (unstable) |
| Horizontal (H) | $S_0 = 0$ | no normal depth |
| Adverse (A) | $S_0 < 0$ | no normal depth |

A *circular* pipe carries its maximum flow at about 94 % of the diameter, so for a given slope there is a largest $Q$ it can carry partly full; the program reports this instead of failing.

## 4. Gradually varied flow  (`reference.py`, used for validation)

When the depth is not the normal depth, it changes slowly along the channel. Differentiating $H = z + y + V^2/2g$ and using $dH/dx = -S_f$ gives

$$\frac{dy}{dx} = \frac{S_0 - S_f}{1 - Fr^2}, \qquad S_f = \frac{n^2 Q^2}{A^2 R^{4/3}}$$

The sign of numerator and denominator decides the shape. The surface curves are named by **slope letter + zone number**, where zone 1 is above both $y_n$ and $y_c$, zone 2 is between them and zone 3 is below both:

| Slope | Zone 1 | Zone 2 | Zone 3 |
|---|---|---|---|
| Mild | M1 backwater | M2 drawdown | M3 |
| Steep | S1 | S2 | S3 |
| Critical | C1 | none | C3 |
| Horizontal | none | H2 | H3 |
| Adverse | none | A2 | A3 |

Far from a control the profile approaches $y_n$ (the numerator goes to zero); near critical depth it becomes steep (the denominator goes to zero), which is where the method needs care.

## 5. The standard step method  (`step.py`)

Between a known station 1 and an unknown station 2 a distance $\Delta x$ apart, the energy equation is

$$z_1 + y_1 + \frac{V_1^2}{2g} = z_2 + y_2 + \frac{V_2^2}{2g} + \bar S_f\,\Delta x, \qquad \bar S_f = \tfrac12 (S_{f1} + S_{f2})$$

with the unknown $y_2$ inside $V_2$ and $S_{f2}$, so it is solved by iteration. The program uses `brentq` on the correct **branch** of the energy curve: depths above $y_c$ for subcritical flow, below $y_c$ for supercritical flow. This always converges when a solution exists.

* If no solution exists on that branch, the energy is below the minimum for that section: the flow cannot stay on this branch. The step is first split into smaller sub-steps; if that fails, the profile has reached critical depth and a warning is given.
* The method is second-order accurate: halving $\Delta x$ cuts the error about four times (shown in `VALIDATION.md`).

**Reach junctions.** The bed is continuous and energy is conserved across a junction (a zero-length step), so a change of cross-section changes the depth there; a change of slope or roughness only changes how the profile continues. Local contraction or expansion losses are not modelled.

## 6. Controls and how the profile is assembled  (`profile._controls`, `_sweeps`, `_assemble`)

A **control** is a place where the depth is known. The program uses:

* your **upstream** condition (a sluice gate gives $y = C_c a$ with $C_c \approx 0.61$; a lake inlet gives critical depth; or a known depth/level);
* your **downstream** condition (a tailwater depth or level; a weir; a free overfall giving $y_c$; or normal depth);
* **critical depth at a break from a flatter to a steeper reach** (for example mild → steep): the flow passes through critical depth there.

From each control the program sweeps with the standard step method in the direction that is physically allowed. This gives, in general, a supercritical profile and a subcritical profile that overlap. The real water surface follows the supercritical one until a **hydraulic jump** moves it onto the subcritical one.

## 7. The hydraulic jump  (`jump.py`, `profile._locate`)

### 7.1 Momentum and the sequent depth

Across a jump energy is lost in turbulence, so energy is *not* conserved — but momentum is (horizontal bed, no outside forces). The **momentum function** is

$$M(y) = \frac{Q^2}{gA} + A\bar y$$

where $\bar y$ is the depth of the centroid of the flow area ($A\bar y = \int_0^y (y-h)\,T(h)\,dh$, computed in closed form for each section type). The **sequent (conjugate) depth** $y_2$ satisfies $M(y_2) = M(y_1)$, solved with `brentq` for any section. For a rectangle this reduces to the Bélanger equation

$$\frac{y_2}{y_1} = \tfrac12\left(\sqrt{1 + 8 Fr_1^2} - 1\right)$$

### 7.2 Energy loss and type

$E_1 - E_2$ from the two specific energies; for a rectangle this equals $\dfrac{(y_2 - y_1)^3}{4 y_1 y_2}$. The jump type depends on the upstream Froude number: undular ($Fr_1 < 1.7$), weak ($< 2.5$), oscillating ($< 4.5$), steady ($< 9$), strong ($> 9$) (commonly quoted limits).

### 7.3 Where the jump forms  (Figure 3)

Compute both profiles, then compare their momentum functions along the channel. Going downstream, the **toe** of the jump is where the momentum of the supercritical flow falls to that of the subcritical (tailwater-controlled) flow, i.e. where the depth of the tailwater profile equals the sequent depth of the supercritical depth. The position is refined between stations with `brentq`.

* Tailwater too high at the start: the jump is **drowned** (pushed up to the gate).
* Tailwater too low everywhere: the jump is **swept out** and the flow stays supercritical to the end.

### 7.4 The jump has a length

A jump is not a step. The water needs a distance to rise from $y_1$ to $y_2$; this is the length of the **roller**. The program uses the empirical rule

$$L_j \approx 6.1\, y_2$$

(roughly "six times the sequent depth"; an approximation that depends on $Fr_1$ and is written from memory, so check it against your textbook). The toe is at the position found above; the roller ends at $x_{toe} + L_j$, where the flow continues on the subcritical profile. Inside the roller the program draws an idealised S-shaped surface and a smoothly falling energy line between the computed end values. **Only $y_1$, $y_2$, $L_j$ and the energy loss are computed; the shape and the foam in the 3D view are illustrations.**

### 7.5 After the jump

The flow is subcritical and follows the profile controlled by the tailwater: if the tailwater equals the normal depth you get uniform flow (as in the first example); if it is higher you get an M1 backwater rising towards it; if it is lower, an M2 drawdown.

## 8. Assumptions and limitations

Steady, one-dimensional, gradually varied flow; prismatic reaches; hydrostatic pressure; Manning friction; no sediment, no side inflow, no local losses at junctions (no choking or constriction control at section changes); the jump uses the horizontal-floor momentum equation even on a sloping bed; jump length is empirical; circular pipes are limited to 99 % of the diameter (no pressurised flow); a critical-slope reach is treated like a steep reach for control purposes; the free overfall is simplified to critical depth at the last station (the true brink depth is lower).

## 9. Further reading

General references for the whole subject: V. T. Chow, *Open-Channel Hydraulics*; M. H. Chaudhry, *Open-Channel Flow*; F. M. Henderson, *Open Channel Flow*. Use them to compare with your own worked examples (see `examples/TEXTBOOK_EXAMPLE.md`).
