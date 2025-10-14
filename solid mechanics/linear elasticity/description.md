# Solving the equations of linear elasticity
$$
\newcommand{\vect}[1]{\boldsymbol{\mathbf{#1}}}
$$

## Governing equations

Consider the open bounded Lipschitz domain $\Omega\subset\mathbb{R}^3$ whose closure $\overline{\Omega}$ represents the reference configuration of an elastic body. We wish to solve for the displacement of the body $\vect{u}(\mathbf{x}): \Omega \to \mathbb{R}^3$ given a body force $\vect{f}(\mathbf{x}): \Omega \to \mathbb{R}^3$ acting on the body. The momentum equation for the system is given by
$$

    -\vect{\nabla} \cdot \vect{\sigma} = \vect{f},
$$
where $\vect{\sigma}(\mathbf{x}): \Omega \to \mathbb{R}^{3\times 3}_{\text{sym}}$ is the Cauchy stress tensor. For a homogeneous and isotropic material, the stress strain relationship is given by
$$
    \vect{\sigma}(\vect{u}) = \lambda \text{tr}(\vect{\varepsilon(\vect{u})})\vect{I} + 2\mu \vect{\varepsilon}(\vect{u}),
$$
where $\lambda$ and $\mu$ are the Lamé parameters and $\vect{\varepsilon}(\vect{u})$ is the infinitesimal strain tensor defined as
$$
    \vect{\varepsilon}(\vect{u}) = \text{sym}(\vect{\nabla u}) = \frac{1}{2}\left(\vect{\nabla}\vect{u} + (\vect{\nabla}\vect{u})^T\right),
$$
assuming small displacements. To close the system we need to impose boundary conditions.
We decompose the boundary into two disjoint parts of non-zero measure such that
$\partial \Omega = \Gamma_D \cup \Gamma_N$. Then we set
$$
    \vect{u} = 0 \quad \text{on } \Gamma_D,
$$
$$
    \vect{\sigma} \cdot \vect{n} = \vect{t} \quad \text{on } \Gamma_N,
$$
where $\vect{t}$ is a prescribed traction vector.

# Variational formulation
To derive the variational form of problem we take the inner product of the momentum equation with a test function $\vect{v} \in [H^1(\Omega)]^3$ and integrate over the domain $\Omega$ to obtain
$$
    -\int_{\Omega} (\vect{\nabla} \cdot \vect{\sigma})\cdot\vect{v} \ \mathrm{d}V = \int_{\Omega} \vect{f} \cdot \vect{v} \, \mathrm{d}V.
$$
Integrating the left hand side by parts we get
$$
    \int_{\Omega} \vect{\sigma} : \vect{\nabla v} \, \mathrm{d}V = \int_{\Omega} \vect{f}\cdot\vect{v} \, \mathrm{d}V
    + \int_{\Gamma_N} \vect{t} \cdot \vect{v} \, \mathrm{d}S.
$$
Finally, we not that the inner product of a symmetric tensor and an antisymmetric tensor is zero
we can replace $\vect{\nabla v}$ with $\vect{\varepsilon(v)}$ to obtain the final variational form
$$
    \int_{\Omega} \vect{\sigma} : \vect{\varepsilon(v)} \, \mathrm{d}V = \int_{\Omega} \vect{f}\cdot\vect{v} \, \mathrm{d}V
    + \int_{\Gamma_N} \vect{t} \cdot \vect{v} \, \mathrm{d}S.
$$

## Bending of a cantilever beam

As an example we consider a beam $\Omega = (0,L)\times(0,W)\times(0,W)$ subject to a body force
$$
    \vect{f}(\mathbf{x}) = (0,0,-\rho g).
$$
The boundary is traction free.