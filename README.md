# 🧫 Bacterial Growth Simulator

An interactive, browser-based pharmacodynamics simulator modeling bacterial 
growth and antibiotic suppression using live ODE (RK4) solving.

🔗 **Live Demo:** [Launch Simulator](https://jobinjohn24.github.io/Bacterial_Growth_Simulation/)

---

## Features

- **Live ODE Simulation** — RK4 solver runs in the browser; every slider change instantly recomputes all charts
- **Time Course Panel** — Overlaid log₁₀(N) curves for 5 concentrations + custom value with hover tooltips
- **Dose–Response Curve** — Sigmoid AUC curve with live EC₅₀ reference line
- **Bacterial Animation** — Canvas particles that grow/shrink/change color with antibiotic pressure
- **Run Comparison** — Save and compare up to 4 experimental runs (table or bar chart view)

## Parameters You Can Tweak

| Parameter | Description |
|---|---|
| Growth rate `r` | Bacterial replication rate (hr⁻¹) |
| Carrying capacity `K` | Max population size |
| Max kill rate `kmax` | Antibiotic potency (hr⁻¹) |
| EC₅₀ | Concentration for 50% kill effect |
| Hill coefficient `h` | Steepness of dose–response curve |
| Dosing schedule | Constant vs. pulse dosing |

## How to Use

1. Adjust sliders in the **Model Parameters** panel
2. Watch all charts update in real time
3. Click **Run 1–4** to save experimental snapshots
4. Click **LOAD** on any saved run to restore its parameters
5. Toggle **Light/Dark** theme and **Lines/Areas** chart style in the Tweaks panel

## Built With

- Vanilla JavaScript (no frameworks)
- HTML5 Canvas
- RK4 numerical ODE solver
