---
author:
  - LemurPwned
date: August 2025
title: LLG linearisation
---

# Introduction

In certain systems that exhibit stiffness or known stable states we can omit computing the equilibrium and compute the oscillation modes directly thanks to the linearisation procedure.

# Short derivation

Here we assume a 2-layer FM system that is stable along $z$ axis with $m_z = \pm 1$.

We start out with the effective field:

$$
\mathbf{H}_j = \frac{1}{\mu_0\mathrm{M}_\mathrm{s} t}\nabla_\mathbf{m} \mathcal{F} \tag{1}
$$

where $t$ is the thickness, $\mathrm{M}_\mathrm{s}$ is the magnetisation saturation and $\mathcal{F}$ is the energy expression. For details on how to construct $F$ you can check out [Smit-Beljers model description](sb-model.md). Index $j$ denotes the layer.

The stable oscillation LLG form is given by:

$$
\mathrm{d}_t \mathbf{m}_j = - \mu_0 \gamma \mathbf{m}_j \times \mathbf{H}_j
\tag{2}
$$

In the linearisation procedure for stable $z$ axis we need to assume that $m_x << 1$ and $m_y << 1$, meaning that non $z$ components experience small variation. In effect, this is a picture of a tiny conical oscillation about the $z$-axis. This small precession assumption allows us to drop second order approximations later from expanding $\mathbf{m}_j$.

Therefore, we linearise both $\mathbf{m}_j$ and $\mathbf{H}_j$ as follows:

$$
\mathbf{m}_j = \mathbf{m}_{j,0} + \delta \mathbf{m}_j
$$

and

$$
\mathbf{H}_j = \mathbf{H}_{j,0} + \delta \mathbf{H}_j
$$

and we have

$$
 \delta \mathbf{m}_j = (m_{x,j}, m_{y,j}, 0) \\ \mathbf{m}_{j,0} = (0, 0, s)
$$

with $s = \pm 1$.

To expand the $\mathrm{d}_t \mathbf{m}_j$ we can compute the Jacobian matrix. So, for a single layer, again, we have, using a differential of a cross product rule at equilibrium:

$$
  \partial\,(\mathbf{m}_j \times \mathbf{H}_j) =  (
    \,\partial\mathbf{m}_j \times \mathbf{H}_j +  \mathbf{m}_{j,0} \times \left(\frac{\partial\mathbf{H}_{j,0}}{\partial\mathbf{m_j}}\right)_0 \partial \mathbf{m}_j
) \tag{3}
$$

We can define the corresponding Jacobian **operator** with the help of 3x3 $\mathcal{E}(x)$ matrix which means "cross with x" and later substitutions which cause the higher olders to be neglected.

$$
-\mu_0\gamma J_\mathbf{m} = -\mu_0\gamma (\mathcal{E}(\mathbf{H_0}) + \mathcal{E}(\mathbf{m_0})\left(\partial \mathbf{H}/\partial\mathbf{m}\right)_{\mathbf{m}_0})
\tag{4}
$$

so when we do

$$
-\mu_0\gamma J_\mathbf{m} \mathbf{m}
\tag{5}
$$

it gives us the expression $(3)$ in matrix form.
In practice, we compute, per each layer, the RHS of the above LLG equation ($(2)$), and create a Jacobian matrix from that RHS side, component-wise:

$$
- \mu_0 \gamma J_\mathbf{m}\begin{pmatrix}
\mathbf{m}_1 \times \mathbf{H}_1 \\
\mathbf{m}_2 \times \mathbf{H}_2
\end{pmatrix}  =
- \mu_0 \gamma J_\mathbf{m}\begin{pmatrix}
    f_1(\mathbf{m}_1, \mathbf{H}_1) \\
    f_2(\mathbf{m}_2, \mathbf{H}_2)
\end{pmatrix}
= \\
- \mu_0 \gamma J_\mathbf{m} \begin{pmatrix}
    f_{x,1}(\mathbf{m}_1, \mathbf{H}_1) \\
    f_{y,1}(\mathbf{m}_1, \mathbf{H}_1) \\
    f_{z,1}(\mathbf{m}_1, \mathbf{H}_1) \\
    f_{x,2}(\mathbf{m}_2, \mathbf{H}_2) \\
    f_{y,2}(\mathbf{m}_2, \mathbf{H}_2) \\
    f_{z,2}(\mathbf{m}_2, \mathbf{H}_2) \\
\end{pmatrix}
=
- \mu_0 \gamma\begin{pmatrix}
    \frac{d}{dm_{x,1}} f_{x,1} && \frac{d}{dm_{y,1}} f_{x,1} && \frac{d}{dm_{z,1}} f_{x,1} \\
    ... && ... && ...\\
    \frac{d}{dm_{x,2}} f_{z,2} && \frac{d}{dm_{y,2}} f_{z,2} && \frac{d}{dm_{z,2}} f_{z,2} \\
\end{pmatrix} = \mu_0\gamma \mathcal{J}
\tag{6}
$$

which in fact is a 6x6 matrix because $f_1$ and $f_2$ are expanded per each of 3 cross product components.
To make that matrix follow the linearisation approximations we the following subsitutions: $m_{x, i} = 0, m_{y, i}= 0$. This causes the entries of $\mathcal{J}$ non-dependent solely on $m_z$ vanish (the entries, not the corresponding rows and columns). Crucially, field terms in $\mathbf{H}_1$ depend on both $\mathbf{m}_1$ and $\mathbf{m}_2$ due to presence of coupling terms.

The substitution allows us to drop the rank of the matrix $\mathcal{J}$ from 6 to 4, because 3rd and 6th row and column will disappear since the $z$ component of $m \times H$ cross product does not contain dependency on $m_z$ (due to the nature of said cross product), hence they surely vanish. This new matrix we call $\mathcal{J}_0$.

Up to the first order, we make the following substiution: $\mathrm{d}_t \mathbf{m} = i \omega \mathbf{m}$, hence after collecting terms with:

$$
\mathcal{C} = i \omega \mathbf{m} + \mu_0\gamma \mathcal{J}_0
=
\rightarrow i\omega
\begin{pmatrix}
m_{x,1} \\
m_{y,1} \\
m_{x,2} \\
m_{y,2} \\
\end{pmatrix} +  \mu_0\gamma \mathcal{J}_0
\begin{pmatrix}
m_{x,1} \\
m_{y,1} \\
m_{x,2} \\
m_{y,2} \\
\end{pmatrix}
\tag{7}
$$

where $\mathbf{m}$ contains the leftover $x, y$ components of $m$. Finally, we substute $m_z \pm 1$ depending on the state in $J_0$, leading to fully determined matrix. In $(7)$ we applied the operator $(6)$ to $\mathbf{m}$.
The only thing left is to compute the eigenvalues of the characteristic matrix $\mathcal{C}$ which would be the frequencies of the system.

# TL;DR

Steps we do in code:

1. Magnetisation vector of all layers is defined as $\mathbf{M} = (m_{x,1}, m_{y,1}, m_{z,1}, m_{x,2}, m_{y,2}, m_{z,2})^T$
2. We compute the effective field per layer as:
   $$
   \mathbf{H}_i = -\nabla_{\mathbf{m}_i} \mathcal{F}_i
   $$
3. LLG given by (note how the fields will depend on entire $\mathbf{M}$, not just its layers $m_i$):
   $$
   \frac{d\mathbf{m}}{dt} = -\mu_0 \gamma \begin{pmatrix} \mathbf{m}_1 \times \mathbf{H}_1(\mathbf{M}) \\ \mathbf{m}_2 \times \mathbf{H}_2(\mathbf{M}) \end{pmatrix} = \mathbf{T}(\mathbf{M})
   $$
4. We take Jacobian $\mathcal{J} = \partial \mathbf{T}/\partial \mathbf{M}$
5. We substitute $m_{i, x} = m_{i, y} = 0$ and $m_{i, z} = \pm 1$ in $\mathcal{J}$ to make $\mathcal{J}_0$
6. We drop columns and rows corresponding to $m_z$ component (it's not moving and it's zero). This is due first order approximation $\delta m_{z,i} = 0$
7. Characteristic equation: $\mathcal{C} = i\omega I - \mathcal{J}_0$
8. Compute the determinant
