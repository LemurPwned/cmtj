---
author:
  - LemurPwned
date: June 2021
title: Macromagnetic models
---

# Introduction

In this section we will walk through the basics of the LLG equation and the transformation to the LL-form of the LLG equation.
The LLG equation is a "proper" equation in physical sense (derived from actual mechanics), whereas LL is simply a transformed approximation of the LLG equation that
allows us to solve it numerically (no implicit term on $\frac{dm}{dt}$).

# Landau Lifshitz form of Landau Lifshitz-Gilbert equation

Standard form of the LLG-SOT equation:

$$\frac{d\textbf{m}}{dt} = -\gamma \textbf{m} \times \textbf{H}_{\mathrm{eff}} + \alpha \textbf{m}\times \frac{d\textbf{m}}{dt} -\gamma|\tau_{fl}|\textbf{m} \times \phi-\gamma|\tau_{dl}|\textbf{m}\times\textbf{m}\times \phi$$

Multiply that
equation $\times m$:

$$
\begin{gathered}
       \textbf{m} \times \frac{d\textbf{m}}{dt} = -\gamma  \textbf{m} \times\textbf{m} \times \textbf{H}_{\mathrm{eff}} + \alpha \textbf{m} \times \textbf{m}\times \frac{d\textbf{m}}{dt}  -   \gamma|\tau_{fl}| \textbf{m} \times \textbf{m} \times \phi-   \gamma|\tau_{dl}|\textbf{m} \times \textbf{m}\times\textbf{m}\times \phi \\
   \textbf{m} \times \frac{d\textbf{m}}{dt} = -\gamma  \textbf{m} \times\textbf{m} \times \textbf{H}_{\mathrm{eff}} - \alpha \frac{d\textbf{m}}{dt}  -  \gamma|\tau_{fl}| (\textbf{m}(\textbf{m}\cdot\phi) - \phi) +   \gamma|\tau_{dl}|\textbf{m}\times \phi\end{gathered}
$$

Substitute RHS of the derived equation above into $\textbf{m}\times \frac{d\textbf{m}}{dt}$:

$$
\begin{gathered}
  \frac{d\textbf{m}}{dt} = -\gamma \textbf{m} \times \textbf{H}_{\mathrm{eff}} + \alpha [-\gamma  \textbf{m} \times\textbf{m} \times \textbf{H}_{\mathrm{eff}} - \alpha \frac{d\textbf{m}}{dt}  -  \gamma|\tau_{fl}| (\textbf{m}(\textbf{m}\cdot\phi) - \phi) +   \gamma|\tau_{dl}|\textbf{m}\times \phi] \\ -
   \gamma|\tau_{fl}|\textbf{m} \times \phi \\ -
    \gamma|\tau_{dl}|\textbf{m}\times\textbf{m}\times \phi \\
     \frac{d\textbf{m}}{dt}(1 + \alpha^2)  = -\gamma \textbf{m} \times \textbf{H}_{\mathrm{eff}} - \alpha\gamma\textbf{m}\times\textbf{m}\times\textbf{H}_{\mathrm{eff}} \\
     - \gamma|\tau_{fl}|[\textbf{m} \times \phi  + \alpha (\textbf{m}(\textbf{m}\cdot\phi) - \phi) ] \\
     -  \gamma|\tau_{dl}|[\textbf{m}\times\textbf{m}\times \phi - \alpha\textbf{m}\times \phi]\end{gathered}
$$

Rearranging the terms gives:

$$
\begin{aligned}
    \frac{d\textbf{m}}{dt} = \frac{-\gamma}{1 + \alpha^2}[\textbf{m} \times \textbf{H}_{\mathrm{eff}} + \alpha\textbf{m}\times\textbf{m}\times\textbf{H}_{\mathrm{eff}} &\\
     + |\tau_{fl}|[\textbf{m} \times \phi  + \alpha (\textbf{m}(\textbf{m}\cdot\phi) - \phi) ] &\\
     +  |\tau_{dl}|[\textbf{m}\times\textbf{m}\times \phi - \alpha\textbf{m}\times \phi]]\end{aligned}
$$

In this form, $\gamma$ is the gyromagnetic ratio and is equal to $\gamma \approx 2.2e5 \frac{m}{As}$.
The last part can be re-arranged to:

$$
\begin{aligned}
(...) + |\tau_{fl}|[\textbf{m} \times \phi  + \alpha (\textbf{m}(\textbf{m}\cdot\phi) - \phi) ] + |\tau_{dl}|[\textbf{m}\times\textbf{m}\times\phi - \alpha\textbf{m}\times \phi]]&\\ = \textbf{m} \times \phi(|\tau_{fl}| - \alpha|\tau_{dl}|) + \textbf{m}\times\textbf{m}\times\phi(|\tau_{dl}| + \alpha|\tau_{fl}|)\end{aligned}
$$

What is evident in this form of LL form of the LLG equation is the
mixing of the torques with damping as the scaling factor (the field-like
term for instance now becomes $|\tau_{fl}| - \alpha|\tau_{dl}|$). Proper
LL-form of the LLGS equation is:

$$
(...)+ \frac{-\gamma}{1 + \alpha^2}[
    \tau'_1 \textbf{m}\times(\textbf{m}\times \phi)
    + \tau'_2 \textbf{m}\times \phi]
$$

## STT interaction

The origin of STT is different, thus, we use a different set of
quantities:

$$
\begin{gathered}
    |H_{dl}| = \beta |H_{fl}|, \quad \beta \in [0, 1] \\
    a_j = \frac{\hbar j_e}{eM_s t_{FM}} \\
    \eta = \frac{\sigma \Lambda^2}{\Lambda^2 + 1 + (\Lambda^2 -1)\textbf{m}\cdot\sigma}\end{gathered}
$$

Given those two terms above, the non-torque part remains the same:

$$... +  a_j\eta\beta\textbf{m} \times \sigma + a_j\eta\textbf{m}\times\textbf{m}\times \sigma$$

Then, equation only changes the coefficients scaling the
damping-like and field-like torques.

The LL form of the STT equation is:

$$
\begin{aligned}
\frac{d\textbf{m}}{dt} = \frac{-\gamma}{1 + \alpha^2}[\textbf{m} \times \textbf{H}_{\mathrm{eff}} + \alpha\textbf{m}\times\textbf{m}\times\textbf{H}_{\mathrm{eff}} &\\
+ (-a_j (\mathbf{m}\times\mathbf{p})  + a_j\beta (\mathbf{m}\times\mathbf{m}\times\mathbf{p})]
\end{aligned}
$$

where
$$a_j =  \gamma_0 \eta \frac{\hbar j}{e M_\mathrm{s} t_\mathrm{FM}}$$

# Stochastic LLGS

## Stratonovich formulation of the s-LLGS SDE

A stochastic formulation of LLGS will take the form of a Stratonovich
SDE:

$$\mathrm{d}X_t = f(X_t, t)dt + g(X_t, t)\circ \mathrm{d}W_t\sqrt{\Delta t}$$

where $f(X_t, t)$ is the deterministic part of the
equation and $g(X_t, t)$ is the stochastic part of the equation.
$\mathrm{d}W$ is \"derivative-like\" of the Brownian motion. The symbol
$\circ$ denotes the Stratonovich product which distinguishes it from
Ito's SDE. By assuming that the effective field contains thermal
fluctuations
$\mathbf{H}_{\mathrm{eff}} \rightarrow \mathbf{H}_{\mathrm{eff}} + \mathbf{H}_{\mathrm{T}}$
we transform the standard LLGS equation into the form that fits
Stratonovich SDE. The
thermal fluctuations have zero mean and a preset standard deviation:

$$\sigma(t) = \sqrt{\frac{2\alpha k_bT(t)}{M_sV\gamma_0}}$$

where $V$ is the volume of the cell (layer), and $k_bT(t)$ is the thermal energy of the system. As a
result, $\sigma(t)$ should be dimensionless. Finally, we set
$\mathbf{f}(\mathbf{m}_t, t)$ to LL form where $\mathbf{H}_{\mathrm{eff}}$ contains no
stochastic (thermal) parts and the $g$, the stochastic part, to the
following:

$$
\mathbf{g}(\mathbf{m}_t, t)\circ\mathrm{d}W  =
    - \frac{\sigma\gamma}{1+\alpha^2}[\mathbf{m}\times\mathrm{d}W + \alpha\mathbf{m}\times(\mathbf{m}\times\mathrm{d}W)]
$$

with $\mathrm{d}W \in \mathbf{R}^3 \sim \sqrt{t}\mathcal{N}(0, 1)$, a
multinomial Gaussian distributed random vector (here we make a
transition from $W$ being a generalised Brownian process to a Wiener
process). For numerical solutions, we have have $\Delta W$ instead of $\mathrm{d}W$.
$\Delta W(t) = W(t + \Delta t) - W(t)$, where the stochastic vector is being drrawn from a normal distriubtion, with zero mean and unit variance: $\xi_t \in \mathbf{R}^3 \sim \mathcal{N}(0, 1)$.
The form above follows from the distributive
properties of cross-product over addition. Furthermore, there is some
evidence that the second term in that equation should be skipped if the noise is
sufficiently small which seems to be the case for up to room temperature
experiments.

## Numerical solutions

We generally solve the stochastic model by either with Euler-Heun or Heun method.

### Euler-Heun method

This is in fact first order-method in the limit of 0 K.
Euler-Heun method is suitable for Stratonovich SDEs as Euler-Maruyama
can only be applied to Ito's SDEs. The update of the step is:

$$Y_{n+1} = Y_n + f_n \Delta t + \frac{1}{2}[g_n + g_n(\hat{Y}_n)]\Delta W_n\sqrt{\Delta t}$$

where $\hat{Y}_n = Y_n + g_n\Delta W_n$. Contrary to the Milstein
method, it is easier to the user the Euler-Heun due to the lack of
quadratic terms of $\Delta W_n$. The cost is in the convergence order
which is 0.5 for strong convergence and 1 for weak convergence. For the
solution, we substitute $Y_n = \mathbf{m_t}$,
$f_n = \mathbf{f}_n(\mathbf{m_t}, t)$,
$g_n= \mathbf{g}_n(\mathbf{m_t}, t)$.

### Heun method

Now preferred method to solve stochastic form of the LLG equation is the Heun method. It introduces second order correction to the non-stochastic part as well and
therefore is deemed a better method.

$$Y_{n+1} = Y_n + \frac{1}{2}\left[f_n(\hat{Y}_{n+1}, t_{n+1}) + f_n(Y_n, t_n)\right] + \frac{1}{2}\left[g_n(\hat{Y}_{n+1}, g_{n+1}) + g_n(Y_n, t_n)\right]\Delta W_n$$

where $\hat{Y}_{n+1} = Y_n + f_n(Y_n, t_n)\Delta t + g_n(Y_n, t_n)\Delta W_n\sqrt{\Delta t}$.

### References

[Numerical Integration of SDEs: A Short Tutorial,
Thomas Schaffter, January 19, 2010](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwilpP-T5p_yAhXjAxAIHZosBBgQFnoECAgQAw&url=https%3A%2F%2Finfoscience.epfl.ch%2Frecord%2F143450%2Ffiles%2Fsde_tutorial.pdf&usg=AOvVaw1VNG29Y2knOPBB3Hic2QvU)
