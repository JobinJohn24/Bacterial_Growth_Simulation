# 🦠 Bacterial Growth Simulation with Dosing Strategies 🦠

## 📋 Overview

This project sets up a compartmental ODE model. It aims to simulate the way bacterial populations shift and grow when facing antibiotic pressure. The model builds in suppression effects that depend on the dose given. It also looks into various dosing approaches by running four organized simulations.

**Key Features:**
- ODE-based bacterial growth model with antibiotic pharmacodynamics (PK/PD)
- Configurable dosing strategies (constant and pulsed)
- Parameter fitting via least-squares optimization
- Publication-quality visualizations
- Modular, documented R code architecture

## 🛠️ Project Structure

```
.
├── src/                      # R source code
│   ├── main.R               # Pipeline orchestration
│   ├── model.R              # ODE definition
│   ├── simulate.R           # Data generation
│   ├── fit.R                # Parameter fitting
│   ├── visualize.R          # Plotting
│   └── io.R                 # Config & I/O
├── config/
│   └── sim_params.yaml      # Parameters
├── figures/                 # Output plots (8 files)
├── results/                 # Output metrics (4 CSV files)
├── data/                    # User data (optional)
└── README.md
```

---

## 📋 The Model

**ODE Equation:**
```
dN/dt = r·N·(1 - N/K) - k_max · [A(t)^h / (EC50^h + A(t)^h)] · N
```

**Parameters:**
- **N(t)**: Bacterial density (cells)
- **r**: Growth rate (hr⁻¹)
- **K**: Carrying capacity (cells)
- **A(t)**: Antibiotic concentration (μg/mL, time-dependent)
- **k_max**: Maximum kill rate (hr⁻¹)
- **EC50**: Half-maximal effect concentration (μg/mL)
- **h**: Hill coefficient (sigmoidal slope)

**Dosing Modes:**
- **Constant**: A(t) = A₀ throughout
- **Pulse**: Square pulses every N hours for M-hour duration

---

## Results: 4-Run Comparison

| Run | Configuration | EC50 (μg/mL) | Hill (h) | Best Conc | AUC | Key Finding |
|-----|---|---|---|---|---|---|
| **1** | Constant, fit EC50+h | 2.43 | 1.17 | 20 | 1.06e+06 | **Baseline reference** |
| **2** | Pulse (12h), fit EC50+h | 2.43 | 1.17 | 20 | 1.06e+06 | **Clinical relevance** |
| **3** | Constant, fit EC50 only | 2.00 | 1.50* | 20 | 2.31e+10 | **Parameter sensitivity** |
| **4** | Extended [0-50], fit EC50+h | 2.00 | 1.50 | 50 | 9.23e+05 | **Saturation effects** |

*Hill coefficient fixed at 1.50

---

## 🏃‍♂️ Run Descriptions

### **Run 1: Baseline Analysis**
**Configuration:** Constant dosing [0–20 μg/mL], both parameters fitted

**Finding:** EC50 = 2.43 μg/mL; Hill = 1.17  
Establishes standard dose-response curve. Lower Hill suggests weak cooperativity.

**Interpretation:** Antibiotic effect is relatively linear rather than sigmoidal.

![Run 1 Timecourses](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run1/run1_timecourses_constant%20copy.png)
![Run 1 Dose Response](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run1/run1_dose_response%20copy.png)

---

### **Run 2: Pulse Dosing Strategy**
**Configuration:** Pulsed dosing (12h interval, 2h duration), both parameters fitted

**Finding:** EC50 = 2.43 μg/mL; Hill = 1.17 (identical to Run 1)  
AUC remains consistent despite intermittent dosing.

**Interpretation:** **Clinical relevance**—simulates realistic antibiotic therapy. Bacterial rebound between doses visible in time-courses, explaining resistance development risk.

![Run 2 Timecourses](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run2/run2_timecourses_pulse.png)
![Run 2 Dose Response](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run2/run2_dose_response.png)

---

### **Run 3: Parameter Sensitivity**
**Configuration:** Constant dosing, EC50 fitted only (Hill fixed at 1.50)

**Finding:** EC50 = 2.00 μg/mL; AUC = 2.31e+10 (much higher)  
EC50 shifts when Hill is constrained.

**Interpretation:** Demonstrates parameter interdependence. When prior knowledge constrains Hill, EC50 adjusts to compensate. Shows importance of data quality for reliable parameter estimation.

![Run 3 Timecourses](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run3/run3_timecourses_constant.png)
![Run 3 Dose Response](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run3/run3_dose_response.png)

---

### **Run 4: Extended Dose Range & Saturation**
**Configuration:** Extended range [0–50 μg/mL], both parameters fitted

**Finding:** EC50 = 2.00 μg/mL; Best concentration shifts to 50 μg/mL; AUC = 9.23e+05  
Model saturates at high doses—diminishing returns beyond 20 μg/mL.

**Interpretation:** Identifies **safe and effective dosing windows**. High doses don't provide additional benefit, suggesting toxicity thresholds should be considered in clinical practice.

![Run 4 Timecourses](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run4/run4_timecourses_constant.png)
![Run 4 Dose Response](https://github.com/JobinJohn24/Bacterial_Growth/blob/main/figures/run4/run4_dose_response.png)

---

## 📊 Understanding Outputs

### Time-Course Plots
- **X-axis:** Time (hours)
- **Y-axis:** Bacterial density (CFU/mL, log scale)
- **Panels:** One per antibiotic concentration
- **Pattern:** Higher concentrations show stronger suppression

### Dose-Response Curves
- **X-axis:** Antibiotic concentration (μg/mL)
- **Y-axis:** AUC (area under curve — cumulative bacterial exposure )
- **Red line:** EC50 mark (half-maximal effect)
- **Sigmoid curve:** S-shaped dose-response

### Metrics CSV
| Column | Meaning |
|--------|---------|
| AUC | Total bacterial burden (lower = better) |
| max_density | Peak population size |
| time_to_50pct_reduction | Hours to 50% suppression |
| estimated_EC50 | Fitted potency |
| estimated_h | Fitted cooperativity |

## 👨‍💻 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/jobinjohn24/bacterial-growth-simulation.git
cd bacterial-growth-simulation

# Install dependencies in R
install.packages(c("yaml", "deSolve", "ggplot2"))
```

### Run Simulation

```bash
Rscript src/main.R
```

**Output:** Plots in `figures/` and metrics in `results/`

## Usage

### Change Antibiotic Concentrations

Edit `config/sim_params.yaml`:
```yaml
A_levels: [0, 0.5, 1.5, 3, 5]  # Custom range
```

### Switch Dosing Schedule

```yaml
dose_schedule: "pulse"
pulse_interval: 8      # Every 8 hours
pulse_duration: 1      # For 1 hour
```

### Use Your Own Data

Place CSV in `data/growth.csv` with columns:
```
time, concentration, density
```

## 👨‍🔬 Customization

### Fit Different Parameters

```yaml
fit_free: ["EC50"]              # Fit only EC50
fit_free: ["r", "K"]            # Fit growth parameters
fit_free: ["EC50", "h"]         # Fit potency & cooperativity
```

### Extend Simulation Time

```yaml
tmax: 72                          # 72 hours instead of 48
```

### Adjust Noise Level

```yaml
noise_sd: 1.0e7                   # Increase measurement noise
```
