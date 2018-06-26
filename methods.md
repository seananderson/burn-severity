---
title: "ZOIB methods"
author: ''
csl: ecological-applications.csl
output:
  word_document: default
  html_document:
    df_print: paged
  pdf_document: default
bibliography: refs.bib
---

We required a statistical model to predict the proportion of each burn severity
field measure based on satellite burn severity metrics and other possible
biophysical covariates. These proportions of field burn severity could range
from zero (unburned) to one (completely burned). Beta regression, in which
a logit link is coupled with a beta distribution observation model [e.g.
@ferrari2004], would be an obvious choice to model the proportions between zero
and one, but cannot account for the zeros or ones themselves. Therefore, we fit
zero-one inflated beta (ZOIB) hierarchical regression models [@ospina2012;
@liu2015], which allow for zeros, ones, and continuous proportions between
these bounds.

Following the notation of @liu2015, ZOIB models assume that the response data
(here, proportion burned) $y$ for observation $i$ follow a piecewise
distribution such that:

$$
f(y_i) =
\begin{cases}
p_i  & \text{if $y_i = 0$,} \\
(1-p_i)q_i & \text{if $y_i = 1$,}\\
(1-p_i)(1-q_i)\mathrm{Beta}(a_i, b_i) & \text{if $y_i \in (0,1)$,}
\end{cases}
$$

where $p_i$ represents the probability $\mathrm{Pr}(y_i = 0)$, $q_i$ represents
the conditional probability $\mathrm{Pr}(y_i = 1 \mid y_i \ne 0)$, and $a_i$
and $b_i$ represent the beta distribution shape parameters for $y_i \in (0,
1)$. We can combine these components to derive the unconditional estimate of
the response $E(y_i)$ [@liu2015] as:

$$
E(y_i) = \left(1 - p_i \right) \left( q_i + \left( 1 - q_i \right) \mu_i^{(0,1)} \right).
$$

In versions of the model where we did not have any zeros, we fit a one-inflated
beta (OIB) model, which is similar to the ZOIB model except it omits $p_i$:

$$
E(y_i) = q_i + \left( 1 - q_i \right) \mu_i^{(0,1)}.
$$

Fitting the piecewise distribution requires fitting the following component
models in which we use the superscripts $p$, $q$, and $r$ represent related
data and parameters across the three models. First, we fit a model to the zeros
vs. non-zeros:

$$
\begin{aligned}
y_i^p &\sim \mathrm{Binomial}(p_i)\\
  p_i &= \mathrm{logit}^{-1} \left( \alpha^p + \alpha^{p}_{j[i]} +
  \pmb{X}^p_i \pmb{\beta}^p \right),
\end{aligned}
$$

where $y_i^p$ is a series of zeros and ones with a one representing that the
original data $y_i = 0$. The symbol $\pmb{X}_i$ represents a vector of
predictors and $\pmb{B}$ a vector of slope coefficients. The parameter $\alpha$
represents an intercept and $\alpha_{j[i]}$ represents a fire-specific
intercept (indexed by $j$) that is allowed to vary and is constrained by
a normal distribution in logit space (i.e. a "random intercept"):

$$
\alpha^p_j \sim \mathrm{Normal}(0, \sigma^2_{p}).
$$

Second, we fit a model to the ones vs. non-ones for all cases that were not
$y_i = 0$, with $y_i^q$ representing a series of ones and zeros with one
representing that the original data $y_i = 1 \mid y_i \ne 0$:

$$
\begin{aligned}
y_i^q &\sim \mathrm{Binomial}(q_i)\\
  q_i &= \mathrm{logit}^{-1} \left( \alpha^q + \alpha^{q}_{j[i]} +
  \pmb{X}^q_i \pmb{\beta}^q \right)\\
  \alpha^\mathrm{q}_j &\sim \mathrm{Normal}(0, \sigma^2_{q}).
\end{aligned}
$$

Third, we fit a model to the proportional data for all cases of $y_i$ that were
not exactly zero or one:

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

where $y_i^{r}$ represents a proportional value between zero and one. We
rearranged the beta distribution so that the two shape parameters were
represented by a dispersion parameter $\phi$ and a mean parameter $\mu_i$
[@ferrari2004].

We fit the ZOIB models with Stan [@carpenter2017] and rstan [@rstan2018] for
the statistical software R [@R2018]. After standardizing the predictors by
subtracting their means and dividing them by the two times their standard
deviations [@gelman2008c], we placed weakly informative priors of Normal(0, 5)
on the intercepts $\alpha$, Normal(0, 2) on the slope parameters $\pmb{\beta}$,
Half-t(3, 0, 25) on the dispersion parameter $\phi$ (i.e. a Student-t
distribution with degrees of freedom 3 and scale 25 for values $>$ 0), and
Half-t(3, 0, 2) priors on the $\sigma$ parameters constraining the spread of
the random intercepts. For each model, we sampled from the posterior with 1000
iteration across four chains and discarded the first half of each chain as
a warm-up. We ensured consistency with chain convergence by ensuring
$\widehat{R}$ (the potential scale reduction factor) was less than $1.05$ and
the minimum effective sample size $n_\mathrm{eff}$ was greater than $100$ for
all parameters [@gelman2014].

For Question 1, we fit ZOIB models testing the bivariate relationship between
the satellite and field measures of burn severity. For Question 2, we added the
main effects of commonly available biophysical variables and their interactions
with each satellite metric of burn severity. For Question 3, we instead
included less commonly available biophysical variables and their interactions
with satellite metrics. For Questions 2 and 3, our main focus of inference was
on the interaction coefficients rather than the main effects. A positive
interaction, for example, would indicate that as the value of the biophysical
term in the model (e.g. latitude) increased, the slope between the satellite
and field measures of burn severity became more steep. Since there are three
components to the ZOIB model, there are three versions of each interaction,
which can be visualized in aggregate by combining the component models into
their unconditional expectation $E(y_i)$ at various values of the biophysical
variables.

We evaluated model fit via the area under receiver operating curves (AUC)
across a sequence of thresholds. Calculating the AUC involves dichotomizing the
field response proportions into a series of ones and zeros depending on whether
they are above or below a specific proportion threshold, and then calculating
the probability that the model would correctly rank a randomly chosen
observation now coded as a one as as more burned than a less burned observation
now coded as a zero. We evaluated AUC at thresholds of 0.050, 0.275, 0.500,
0.725, and 0.950. For simplicity, we performed these calculations on the median
posterior predictions instead of across all samples from the posterior. We
initially tried a cross-validated evaluation of model performance for Question
1 where we successively left out each fire from model building and predicted
the omitted data. However, we found this made little qualitative difference to
the results since no one fire substantially drove the model. Therefore, we
chose to present the results without cross validation for simplicity.

**Perhaps in earlier paragraphs about the data: We removed any satellite predictor data with values greater than 2500.**

# References
