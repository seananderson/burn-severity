---
title: "ZOIB methods"
author: ''
csl: ecological-applications.csl
output:
  pdf_document: default
  word_document: default
  html_document:
    df_print: paged
bibliography: refs.bib
---

We fit zero-one inflated beta (ZOIB) hierarchical regression models [@ospina2012; @liu2015] ... Following the notation of @liu2015, ZOIB models assume that the response data $y$ for observation $i$ follow a piecewise distribution such that:

$$
f(y_i) =
\begin{cases}
p_i  & \text{if $y_i = 0$,} \\
(1-p_i)q_i & \text{if $y_i = 1$,}\\
(1-p_i)(1-q_i)\mathrm{Beta}(a_i, b_i) & \text{if $y_i \in (0,1)$,}
\end{cases}
$$

where $p_i$ represents the probability $\mathrm{Pr}(y_i = 0)$, $q_i$ represents the conditional probability $\mathrm{Pr}(y_i = 1 \mid y_i \ne 0)$, and $a_i$ and $b_i$ represent the beta distribution shape parameters for $y_i \in (0, 1)$. We can combine these components to derive the unconditional estimate of the response $E(y_i)$ [@liu2015] as:

$$
E(y_i) = \left(1 - p_i \right) \left( q_i + \left( 1 - q_i \right) \mu_i^{(0,1)} \right).
$$

In versions of the model where we did not have any zeros, we fit a one-inflated beta (OIB) model, which is similar to the ZOIB model except it omits $p_i$:

$$
E(y_i) = q_i + \left( 1 - q_i \right) \mu_i^{(0,1)}.
$$

Fitting the piecewise distribution requires fitting the following component models in wich we use the superscripts $p$, $q$, and $r$ to represent related data and parameters across the three models. First, we fit a model to the zeros vs. non-zeros:

$$
\begin{aligned}
y_i^p &\sim \mathrm{Binomial}(p_i)\\
  p_i &= \mathrm{logit}^{-1} \left( \alpha^p + \alpha^{p}_{j[i]} + 
  \pmb{X}^p_i \pmb{\beta}^p \right),
\end{aligned}
$$

where $y_i^p$ is a series of zeros and ones with a one representing that the original data $y_i = 0$. The symbol $\pmb{X}_i$ represents a vector of predictors and $\pmb{B}$ a vector of slope coefficients. The parameter $\alpha$ represents an intercept and $\alpha_{j[i]}$ represents a fire-specific intercept (indexed by $j$) that is allowed to vary and is constrained by a normal distribution (i.e. a "random intercept"):

$$
\alpha^p_j \sim \mathrm{Normal}(0, \sigma^2_{p}).
$$

Second, we fit a model to the ones vs. non-ones for all cases that were not $y_i = 0$, with $y_i^q$ being a series of ones and zeros with one representing that the original data $y_i = 1 \mid y_i \ne 0$:

$$
\begin{aligned}
y_i^q &\sim \mathrm{Binomial}(q_i)\\
  q_i &= \mathrm{logit}^{-1} \left( \alpha^q + \alpha^{q}_{j[i]} + 
  \pmb{X}^q_i \pmb{\beta}^q \right)\\
  \alpha^\mathrm{q}_j &\sim \mathrm{Normal}(0, \sigma^2_{q}).
\end{aligned}
$$

Third, we fit a model to the proportional data for all cases of $y_i$ that were not exactly zero or one:

$$
\begin{aligned}
y_i^{r} &\sim \mathrm{Beta}(a_i, b_i)\\
a_i &= \phi \mu_i\\
b_i &= \phi (1 - \mu_i)\\
\mu_i &= \mathrm{logit}^{-1} \left( \alpha^{r} + 
  \alpha^r_{j[i]} + 
  \pmb{X}_i^r \pmb{\beta}^r \right)\\
\alpha^r_j &\sim \mathrm{Normal}(0, \sigma^2_{r}),
\end{aligned}
$$

where $y_i^{r}$ represents a proportional value between zero and one. We rearranged the beta distribution so that the two shape parameters were represented by a dispersion parameter $\phi$ and a mean parameter $\mu_i$ [@ferrari2004].

We fit the ZOIB model with Stan [@carpenter2017] and rstan [@rstan2018] for the statistical software R [@R2018]. We placed weakly informative priors of Normal(0, 5) on the intercepts $\alpha$, Normal(0, 2) on the slope parameters $\pmb{\beta}$, Half-t(3, 0, 25) on the dispersion parameter $\phi$ (i.e. a Student-t distribution with degrees of freedom 3 and scale 25 for values $>$ 0), and Half-t(3, 0, 2) priors on the $\sigma$ parameters constraining the spread of the random intercepts. For each model, we sampled from the posterior with 1000 iteration across four chains and discarded the first half of each chain as a warm-up. We ensured consistency with chain convergence by ensuring $\widehat{R}$ (the potential scale reduction factor) was less than $1.05$ and the minimum effective sample size $n_\mathrm{eff}$ was greater than $100$ for all parameters [@gelman2014].

Stuff to fit in:

- We standardized the predictors by subtracting their mean and dividing them by two times their standard deviations [@gelman2008c]. 

- We removed any satellite predictor data with values greater than 2500.

- We calculated area under the receiver operating curves (AUC) across thresholds as...

# References
