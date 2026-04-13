# Masterarbeit
<div align="center">

# ⚡ Nonlinear Dynamics Discovery in Magnetoelectric Sensors

**SINDy + Neural Networks for Data-Driven Modeling of Source Nonlinearity**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-99.6%25-3776AB?logo=python&logoColor=white)](.)
[![MATLAB](https://img.shields.io/badge/MATLAB-0.4%25-orange)](.)
[![Kiel University](https://img.shields.io/badge/Kiel%20University-CAU-003366)](https://www.uni-kiel.de)
[![AIMSE 2023](https://img.shields.io/badge/Presented-AIMSE%202023-brightgreen)](.)

<br>

*Master's thesis at Christian-Albrechts-Universität zu Kiel, Faculty of Engineering*
*Chair for Multicomponent Materials · Department of Material Science*

**Md Saidul Islam** · Supervised by Mohammad Sadeghi, Franz Faupel, Stephan Wulfinghoff

</div>

---

## The Problem

Magnetoelectric (ME) composite sensors are promising candidates for biomagnetic signal detection, but their response becomes **nonlinear** under excitation — even when the actuator operates in its linear regime. The source of this nonlinearity? The experimental setup itself: amplifiers, interfaces, and readout electronics introduce **unwanted harmonics** that corrupt measurements.

Traditional physics-based modeling struggles here because the nonlinearity isn't in the physics you expect — it's in the instrumentation chain.

## The Approach

This thesis asks: **can we discover the governing equations directly from measurement data?**

Two complementary strategies are employed:

**1. Sparse Identification of Nonlinear Dynamics (SINDy)** — A data-driven method that constructs a library of candidate functions (polynomials, trigonometric terms) from measured time-series data, then uses sparse regression to identify which terms actually govern the system. The result is an interpretable, closed-form differential equation — not a black box.

**2. Artificial Neural Networks (MLP & LSTM)** — MultiLayer Perceptron and Long Short-Term Memory networks serve as benchmarks, capturing the same nonlinear input-output relationships for comparison against SINDy's discovered equations.

### Readout Modes

Both methods are evaluated on two distinct sensor readout configurations:

| Readout | Excitation | What's Measured |
|---|---|---|
| **Voltage** | 100 mA AC through spiral coil | Piezoelectric voltage (mV) |
| **Displacement** | 449 mA saturation current | Magnetostrictive layer dimension change via laser Doppler |

## Key Findings

**On model interpretability:** SINDy recovers non-dimensional governing equations with explicit terms for stiffness, damping, excitation frequency, and higher harmonics — yielding physical insight that neural networks cannot provide.

**On data quality:** The voltage dataset is submerged in background noise from the power amplifier, restricting all models (SINDy and ANN alike) from achieving perfect prediction. The displacement dataset, recorded with a power amplifier that provides a much cleaner signal, enables all models to reach ~100% predictability.

**On nonlinear behavior:**
- Voltage data reveals a **hardening damped oscillator** with dominant higher harmonics, nonlinear stiffness, and a **supercritical Hopf bifurcation** via self-excited oscillation
- Displacement data shows a **linear damped oscillator** with negligible higher harmonics and linear stiffness/damping — confirming the Duffing equation framework

**On SINDy vs. ANN:** SINDy achieves comparable accuracy to neural networks while producing human-readable equations. The damping ratio in the noisy voltage equations is ~1000× higher than in the clean displacement equations — a quantitative fingerprint of the instrumentation noise.

## Repository Structure

```
Masterarbeit/
├── data/
│   └── Pre-Processing/       # Raw & pre-processed ME sensor datasets
├── notebooks/                 # Jupyter notebooks for exploration & visualization
├── src/                       # Core Python source code
│   ├── SINDy implementation   # Sparse regression, library construction, optimization
│   ├── ANN models             # MLP & LSTM regressors
│   └── Analysis               # Phase portraits, FFT, force-deflection curves
├── summary_thesis/            # Condensed thesis results & figures
├── LICENSE                    # MIT
└── README.md
```

## Methods in Detail

### SINDy Pipeline

```
Load Libraries → Prepare Data → Compute Numerical Derivatives
    → Build Candidate Function Library Θ(X)
        → Sparse Regression (find coefficient vector Ξ)
            → Evaluate Model (training + validation)
                → Tune Hyperparameters (α threshold, window length)
                    → Extract Governing Equation
```

The optimization proceeds by selecting the sparsification threshold α that minimizes RMSE on held-out data. Two key hyperparameters — the sparsity threshold and the window length for numerical differentiation — are swept and the combination yielding lowest error on both recorded data and derivatives is chosen.

### Non-Dimensional Model Equations

**Voltage readout:**
```
Ü = −38.33 U + (2 × 0.2992)U̇ − 1.428 U³ − 0.019 U̇³
    − 0.0222 Cos(2τ) + 0.0223 Cos(2τ)
    − 0.0035 Cos(3τ) − 0.0018 Cos(4τ) + 0.0063 Cos(5τ)
```

**Displacement readout:**
```
Ü = −39.4905 U − (2 × 0.0002)U̇ + 0.0040 Cos(τ)
    − 0.0007 Cos(2τ) − 0.0046 Cos(3τ)
    − 0.0060 Cos(4τ) − 0.0043 Cos(5τ)
```

### Physical Interpretation

The discovered equations map directly onto the **Duffing oscillator** framework, enabling extraction of:

- **Stiffness** (linear & nonlinear) → Force-deflection curves
- **Damping** (linear & nonlinear) → Phase portraits
- **Excitation frequency contributions** → Higher harmonic content
- **Bifurcation behavior** → Supercritical Hopf bifurcation in voltage data

## Tech Stack

| Tool | Purpose |
|---|---|
| **Python** (NumPy, SciPy) | Core computation, SINDy implementation |
| **scikit-learn** | MLP Regressor |
| **TensorFlow / Keras** | LSTM networks |
| **Matplotlib** | Visualization (phase portraits, FFT, force-deflection) |
| **MATLAB** | Signal pre-processing & DSP |

## References

1. S. L. Brunton, J. L. Proctor, J. N. Kutz, *Proceedings of the National Academy of Sciences*, (2016). [DOI:10.1073/pnas.1517384113](https://doi.org/10.1073/pnas.1517384113) — *Discovering governing equations from data by sparse identification of nonlinear dynamical systems.*

2. H. Yabukami, Wiley, (2011). [DOI:10.1002/9780470977859.ch3](https://doi.org/10.1002/9780470977859.ch3) — *Free Vibrations of a Duffing Oscillation with Viscous Damping.*

3. M. J. Brennan, I. Kovacic, Wiley, (2011), Chapter 2. — *Examples of Physical Systems Described by the Duffing Equation.*

## Acknowledgements

Henrik Wolffrom (Chair for High Frequency Engineering, Department of Electrical and Information Engineering, CAU) for designing the spiral coil and charge amplifier. Hanan Leweiz (Inorganic Functional Materials, Department of Materials Science, CAU) for fabrication of the ME sensors.

---

<div align="center">

**Presented at AIMSE 2023** · Christian-Albrechts-Universität zu Kiel

📬 Contact: sayeed.shahriar@gmail.com

</div>

