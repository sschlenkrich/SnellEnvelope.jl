# SnellEnvelope

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sschlenkrich.github.io/SnellEnvelope.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sschlenkrich.github.io/SnellEnvelope.jl/dev/)
[![Build Status](https://github.com/sschlenkrich/SnellEnvelope.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sschlenkrich/SnellEnvelope.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sschlenkrich/SnellEnvelope.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sschlenkrich/SnellEnvelope.jl)


This repository provides the implementation and case study results for the paper *Bermudan Switch Value Approximation for Non-Standard Swaptions*.

We calculate the Snell envelope $\left(H_k\right)_k$ via terminal condition $H_n=0$ and backward induction
$$
  H_k = \mathbb{E_k}\left[ \max\left( H_{k+1}, U_{k+1} \right) \right].
$$

Underlying payoffs $U_k$ are of the form
$$
  U_k = A_k \cdot \left[ \left(S_k + x_k\right) - K \right]
  \qquad
  \text{(normal model variant)}
$$
or
$$
  U_k = A_k \cdot \left[ S_k \exp\left(x_k\right) - K \right].
  \qquad
  \text{(lognormal model variant)}
$$


## Getting Startet

Clone the repository.

```bash
git clone https://github.com/sschlenkrich/SnellEnvelope.jl
```

Switch to repository root.

```bash
cd SnellEnvelope.jl
```

Instantiate Julia environment.

```bash
julia --project -e "using Pkg; Pkg.instantiate()"
```

Run unit tests.

```bash
julia --project -e "using Pkg; Pkg.test()"
```

Run start an interactive Julia session via `julia --project` and run example via

```julia
include("examples/bullet_bermudan.jl")
```
