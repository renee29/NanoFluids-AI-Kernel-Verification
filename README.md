# ICON₁₀: Non-local Kernel Verification Suite

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-green.svg)
![Method: Asymptotic Analysis](https://img.shields.io/badge/Method-Asymptotic_Analysis-purple.svg)
![Status: Validated](https://img.shields.io/badge/Status-Validated-brightgreen.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18116082.svg)](https://doi.org/10.5281/zenodo.18116082)

Analytical verification of asymptotic limits and coercivity bounds for non-local Gaussian convolution operators on bounded domains.

## Overview

This repository provides rigorous symbolic and numerical verification of how domain truncation induces eigenvalue splitting in isotropic convolution kernels. For a Gaussian kernel on a slit domain, boundary effects break symmetry and generate an effective anisotropic operator with distinct eigenvalues.

![Anisotropy Analysis Result](kernel_limit_analysis.png)

### Key Result

For the z-marginal of a 3D isotropic Gaussian, `K(z) = (2πξ²)^(-1/2) exp(-z²/(2ξ²))`, integrated over `z ∈ [-L/2, L/2]`:

```
λ_∥ = 1                        (unbounded in-plane integration)
λ_⊥ = erf(L / (2√2 ξ))         (truncated out-of-plane integration)
```

where:

- `L` is the confinement length [nm]
- `ξ` is the molecular correlation length [nm]
- `erf` is the error function

**Asymptotic limits verified:**

- Local recovery: `lim_{ξ→0} λ_⊥ = 1` (kernel → δ, identity operator)
- Bulk limit: `lim_{L→∞} λ_⊥ = 1` (full integration)
- Strong confinement: `lim_{L→0} λ_⊥ = 0` (vanishing measure)

### Coercivity and Well-posedness

The effective tensor `K_eff = diag(1, 1, λ_⊥)` admits Cholesky factorisation `K_eff = LL^T` with `L = diag(1, 1, √λ_⊥)`, verifying:

- **Symmetry**: `K = K^T` (Onsager reciprocity)
- **Positive definiteness**: `x^T K x > 0` for all `x ≠ 0`
- **Coercivity bound**: `λ_K := min(λ_∥, λ_⊥) = λ_⊥ > 0` for all `L > 0`

For the bilinear form `a(u,v) = ∫_Ω (K_eff ∇u)·∇v dx`, Gårding's inequality yields:

```
a(u,u) ≥ λ_K ‖∇u‖²_{L²}
```

guaranteeing strict ellipticity and unique weak solution in `H¹(Ω)` via Lax-Milgram.

### Quantitative Benchmark

The anisotropy ratio `λ_∥/λ_⊥ = [erf(L/(2√2ξ))]⁻¹` yields ~1.1 for `L = 1 nm`, `ξ = 0.3 nm`, compared to experimental observations of ~500 in confined systems. This ~400× discrepancy demonstrates that Gaussian kernels capture qualitative symmetry breaking but require refinement for quantitative accuracy.

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Clone the repository

```bash
git clone https://github.com/renee29/NanoFluids-AI-Kernel-Verification.git
cd NanoFluids-AI-Kernel-Verification
```

### Install dependencies

```bash
pip install -r requirements.txt
```

## Usage

Run the complete verification pipeline:

```bash
python nonlocal_kernel_asymptotics.py
```

### Output

The script performs five verification steps:

1. **Symbolic Verification**: Derives `λ_⊥` analytically via SymPy; verifies three asymptotic limits
2. **Numerical Validation**: Cross-validates against quadrature (error < 10⁻¹⁵)
3. **Cholesky Verification**: Confirms `K_eff = LL^T` with reconstruction error < 10⁻¹⁵
4. **Coercivity Analysis**: Verifies `λ_K > 0` and Gårding inequality for Lax-Milgram
5. **Publication Figure**: Generates `kernel_limit_analysis.png` with two panels:
   - (a) Eigenvalue evolution λ_∥ and λ_⊥ vs L/ξ
   - (b) Anisotropy ratio in log scale with experimental comparison

## Project Structure

```
NanoFluids-AI-Kernel-Verification/
├── nonlocal_kernel_asymptotics.py    # Main verification script
├── README.md                         # This file
├── LICENSE                           # MIT License
├── requirements.txt                  # Python dependencies
├── CITATION.cff                      # Citation metadata
└── .gitignore                        # Git ignore rules
```

## Citation

If you use this code in your research, please cite:

```bibtex
@software{nanofluids_ai_kernel_verification_2025,
  author       = {Fabregas, R.},
  title        = {ICON₁₀: Non-local Kernel Verification Suite},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.18116082},
  url          = {https://github.com/renee29/NanoFluids-AI-Kernel-Verification}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaboration inquiries, please open a GitHub issue or contact:

**R. Fábregas** — rfabregas@ugr.es

---

**Project Status**: Initial release (v1.0.0) - December 2025
