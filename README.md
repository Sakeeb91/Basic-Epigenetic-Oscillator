# Basic Epigenetic Oscillator Simulator

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/USERNAME/Basic-Epigenetic-Oscillator/blob/main/notebooks/epigenetic_oscillator.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/USERNAME/Basic-Epigenetic-Oscillator/main?filepath=notebooks%2Fepigenetic_oscillator.ipynb)

This project implements a simulator for the fundamental single negative feedback loop (mRNA-Protein-Metabolite) described by Goodwin (Equations 14 & 15). The simulator models the core dynamics of epigenetic oscillations and allows for visualization and analysis of the system's behavior.

## Run Online

You can run this simulator directly in your browser without installing anything:

1. **Google Colab**: Click the "Open In Colab" badge above to run the notebook in Google Colab
2. **Binder**: Click the "Binder" badge to run in a Binder environment
3. **GitHub Codespaces**: Click the "Code" button on the GitHub repository and select "Open with Codespaces"

## Background

Goodwin's equations describe a simple biochemical oscillator based on negative feedback. The model consists of two state variables:
- X (mRNA) - repressed by Y
- Y (protein) - produced from X

The system exhibits sustained oscillations under certain parameter conditions, which is a fundamental property of many biological clock mechanisms.

## Mathematical Model

### Goodwin Model (Equation 14)
```
dX/dt = a_i / (A_i + k_i*Y) - b_i
dY/dt = α_i*X - β_i*Y
```

### Talandic Energy G (Equation 15)
```
G(X, Y) = (α_i/2)*X^2 - β_i*X + b_i*Y + (a_i/k_i)*log(A_i + k_i*Y)
```

## Features

- Implementation of the Goodwin oscillator model (Equation 14)
- Numerical integration of the system of ordinary differential equations (ODEs)
- Calculation and verification of the Talandic Energy G, a conserved quantity in the idealized model
- Visualization tools for time series, phase portraits, and energy conservation
- Parameter studies to explore the effects of different parameters on oscillation dynamics

## Sample Outputs

The simulator generates several plots to help understand the system dynamics:

### Time Series Plot
![Time Series Plot](https://github.com/USERNAME/Basic-Epigenetic-Oscillator/raw/main/sample_outputs/time_series.png)

### Phase Portrait
![Phase Portrait](https://github.com/USERNAME/Basic-Epigenetic-Oscillator/raw/main/sample_outputs/phase_portrait.png)

### Talandic Energy Conservation
![Talandic Energy](https://github.com/USERNAME/Basic-Epigenetic-Oscillator/raw/main/sample_outputs/talandic_energy.png)

## Local Installation

If you prefer to run the code locally:

1. Clone the repository:
```bash
git clone https://github.com/USERNAME/Basic-Epigenetic-Oscillator.git
cd Basic-Epigenetic-Oscillator
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Run the simulator:
```bash
python src/epigenetic_oscillator.py
```

## Project Structure

```
Basic-Epigenetic-Oscillator/
├── src/
│   └── epigenetic_oscillator.py  # Main simulation code
├── notebooks/
│   └── epigenetic_oscillator.ipynb  # Interactive notebook
├── sample_outputs/
│   ├── time_series.png
│   ├── phase_portrait.png
│   └── talandic_energy.png
├── requirements.txt
└── README.md
```

## Customization

You can modify the parameters in the notebook or Python script to explore different oscillator behaviors:

- Hill coefficient (k_i) - controls the strength of repression
- Degradation rates (b_i, β_i) - affect oscillation period and amplitude
- Production rates (a_i, α_i) - affect the steady-state levels and oscillation amplitude

## References

Goodwin, B.C. (1965). Oscillatory behavior in enzymatic control processes. Advances in Enzyme Regulation, 3, 425-438. 