# Adding a worked example from your course textbook

The tests in this repository check the maths against formulas that can be verified analytically.
They do **not** compare against published textbook numbers, because none were entered from memory.
This folder is the place to add one yourself (for example from Chow, *Open-Channel Hydraulics*, or
Chaudhry, *Open-Channel Flow*).

## Steps

1. Copy `textbook_example_template.json` to a new file, e.g. `chow_example_X.json`.
2. Replace Q, the reaches (length, slope, n, section) and the boundary conditions with the textbook problem.
   Textbooks often use US units: convert to SI first (m, m3/s). Manning's n is the same number;
   if the book uses the 1.486 factor, the SI form used here is equivalent.
3. Run it, with a small step so step-size error is not the issue:

       python -m gvf examples/chow_example_X.json --step 1 --egl --out examples/results

4. Write the book's answer next to ours in the table below and compare.

## What to look for

| Quantity | Where to find it here | Expected agreement |
|---|---|---|
| Normal depth, critical depth | Summary / top of the printout | 3 significant digits |
| Profile type (M1, M2, S2 ...) | Summary line for each reach | Must match exactly |
| Depth at the end of the reach / at a given distance | `*_stations.csv` | Within about 1-2 % (the book may use a coarser hand calculation with few steps) |
| Sequent depth y2 and energy loss | Summary (jump line) | Within 1 %, apart from rounding in the book |
| Jump location | Summary (jump line) | Within about one jump length; the book may use a different approach |
| Jump length | Summary (jump line) | Expect larger differences: this tool uses L = 6.1 y2, an approximate rule; the book may use another formula |

If a number disagrees by more than that:

* check units (feet vs metres, and whether the book's n and slope are the same),
* check whether the book assumes a different boundary condition (for example critical depth at a free overfall),
* try halving the step size,
* check whether the book neglects friction in a short reach.

## Your comparison (fill in)

| Quantity | Textbook | This tool | Difference |
|---|---|---|---|
| | | | |
