# The Origin of the Population algorithm

In this section, I describe structured population modelling and its links with the canonical Matrix Population Models (MPM) and the Degree Day (DD) approach. 
 the dynamically-structured matrix population model step-by-step starting from the canonical Matrix Population Model (MPM) and step-by-step arriving at the sPop model and its more contemporary derivatives. 

## The Matrix Population Model (MPM)

The matrix population modelling approach involves matrix algebra to project the state of a population to the next time point (Caswell2001). According to this, a population is structured into a set of stages and/or within-stage age classes, and a projection matrix is constructed to describe the expected change in each stage/class during a time step.

### Age-structured population dynamics (Leslie matrix)

$$
  M(n,\tau) = \left[
  \begin{array}{cccccc}
    f_0     &            f_1 & f_2         & \cdots & f_{n-2}   & f_{n-1} \\
    s_0      &         0 & 0 & \cdots & 0 & 0\\
    0        &           s_1  &        0 & \cdots & 0 & 0\\
    \vdots & \vdots     & \vdots    & \ddots & \vdots           & \vdots\\
    0        &            0  &        0 & \cdots & 0 & 0\\
    0        &            0  &        0 & \cdots & s_{n-2} & 0\\
  \end{array}
\right],
$$

### Stage-structured population dynamics (Lefkovitch matrix)

$$
  M(n,\tau) = \left[
  \begin{array}{cccccc}
    P_0     &            F_1 & F_2         & \cdots & F_{n-2}   & F_{n-1}\\
    G_0        &         P_1 & 0 & \cdots & 0 & 0\\
    0        &           G_1  &        P_2 & \cdots & 0 & 0\\
    \vdots & \vdots     & \vdots    & \ddots & \vdots           & \vdots\\
    0        &            0  &        0 & \cdots & P_{n-2} & 0\\
    0        &            0  &        0 & \cdots & G_{n-2} & P_{n-1}\\
  \end{array}
\right],
$$

$$ P = \frac{q-q^n}{1-q^n} $$

$$ G = \frac{(1-q)q^n}{1-q^n} $$

### Degree-day (DD) development

This one shows how DD can be used to represent lifetime instead of development:

Gómez NN, Venette RC, Gould JR, Winograd DF. A unified degree day model describes survivorship of Copitarsia corruda Pogue & Simmons (Lepidoptera: Noctuidae) at different constant temperatures. Bull Entomol Res. 2009 Feb;99(1):65-72. doi: 10.1017/S0007485308006111. Epub 2008 Nov 12. PMID: 19006579.


DD process assumes that there is a lower temperature threshold for development, and each day spent under each degree above this threshold accummulates heat contributing to progress development. For instance, 2 days spent under the average condition of 12<sup>o</sup>C, when a threshold of 10<sup>o</sup>C is assumed, 4 DD worth of development accumulates.

At the end of the process, when a certain amount of degree-days accumulate, all which attains this state completes development and moves to the next stage simultaneously and instantly.

Since DD accumulation is on a continous domain, it cannot be readily represented using MPM. Instead, the state of development by the time $t$ can be represented as

$$
    D_0^t = \sum_{i=0}^{t}DD_t, \text{where}\, DD_0 = 0,
$$

where $D_0^t$ is the total amount of degree day accumulated from day 0 until $t$ and $DD_t$ is the degree day accumulated specifically at time $t$.

The equation implicitly ignores mortality, as it deals with the state of development, however, one could describe life processes, such as mortality or fecundity, based on the DD accumulated,

$$
    N_{t+1} = P(D_0^t)N_t + F(D_0^t),
$$

where $N_t$ represents population size at time $t$, and $P$ and $F$ are functions of the DD accumulated by $t$, representing survival and fecundity, respectively. 

This equation assumes that all individuals have been subjected to the same conditions from $t=0$. But the new individuals, generated due to fecundity, would not accumulate the same DD as the old ones.

This requires an extra dimension in the matrix to indicate the time of entry of each individual.

If we assume that this is a closed system and that all individuals were the same age initially, 

Up to this point, no population structure is assumed, however, one would not be restricted to do so. The accumulation of DD is independent of the population structure. If $N$ is composed of individuals of different ages, the progression of DD will be apparent as one moves from the youngest to the oldest individual.

$$
  M(n,\tau) = \left[
  \begin{array}{ccccccc}
    0       &         F(DD_t) & F(DD_t) & \cdots & F(DD_t)   & F(DD_t) & 0\\
    P(DD_t) &         0 & 0 & \cdots & 0 & 0 & 0\\
    0        &        P(DD_t) &        0 & \cdots & 0 & 0 & 0\\
    \vdots & \vdots     & \vdots    & \ddots & \vdots & \vdots & \vdots\\
    0        &            0  &        0 & \cdots & 0 & 0 & 0\\
    0        &            0  &        0 & \cdots & P(DD_t) & 0 & 0\\
    0        &            0  &        0 & \cdots & 0 & P(DD_t) & 1\\
  \end{array}
\right],
$$

$$
  M(n,\tau) = \left[
  \begin{array}{ccccccc}
    0       &         F(DD_1) & F(DD_2) & \cdots & F(DD_{n-1})   & F(DD_n) & 0\\
    P(DD_0) &         0 & 0 & \cdots & 0 & 0 & 0\\
    0        &        P(DD_1) &        0 & \cdots & 0 & 0 & 0\\
    \vdots & \vdots     & \vdots    & \ddots & \vdots & \vdots & \vdots\\
    0        &            0  &        0 & \cdots & 0 & 0 & 0\\
    0        &            0  &        0 & \cdots & P(DD_{n-1}) & 0 & 0\\
    0        &            0  &        0 & \cdots & 0 & P(DD_n) & 1\\
  \end{array}
\right],
$$

Why don't we get rid of the age-dimension and keep track of the DD accumulated? So, the dimensions of the matrix will be DD accumulated, but this is a continuous value.

Instead of DD, what if we keep track of a discrete quantity that accumulates at a rate depending on heat accummulation?

Discretisation is inevitable for numerical simulations. So, why not do it now?

When discretised, the equation looks quite similar to the accumulative process.

$$
  M(n,\tau) = \left[
  \begin{array}{llllll}
    p\,f(0) & F_1 & F_2 & \cdots & F_{k-1} & F_k\\
    p\,f(1) &  p\,f(0) & 0 & \cdots & 0 & 0\\
    p\,f(2) &  p\,f(1) & p\,f(0) & \cdots & 0 & 0\\
    \vdots & \vdots     & \vdots    & \ddots & \vdots           & \vdots\\
    p\,f(k-1) & p\,f(k-2) & p\,f(k-3) & \cdots & p\,f(0) & 0\\
    p\,f(k)   & p\,f(k-1) & p\,f(k-2) & \cdots & p\,f(1) & p\,f(0)\\
  \end{array}
\right],
$$

where M represents the DD at time $t$, $f(d)$ is the probability/rate of accumulating $d$ degree days in one iteration, from $t-1$ to $t$, and $k$ is the integer analogue of the total amount of degree days required for development. 


### Development with renewal processes

Here, we derive the projection matrix for a special case of accumulative development with the deterministic assumption. By fixing the number of pseudo-stages to $k$, we eliminate the need to employ the development indicator $q$, and structure the population into $k$ pseudo-stages, $n_0 \dots n_{k-1}$, and $1$ additional stage, $u$, to receive the individuals completing development. By doing so, we map the development stages onto the projection matrix, and proceed to account for transitions from one pseudo-stage to multiple others.

To calculate the expected number of pseudo-stages one individual accumulates in one step, we use the cumulative density function describing the probability of accumulating $i$ random renewal events in a single time step, $F(i) = F(i,\theta)$ (see Methods). For instance, an individual stays in a pseudo-stage with rate $F(0)$, progresses to the next with rate $F(1)-F(0)$, and skips one to reach the second pseudo-stage with rate $F(2)-F(1)$. Consequently, the projection matrix, $M(n,\tau)$, can be written as

$$
  M(n,\tau) = \left[
  \begin{array}{ccccccc}
    F(0)      & 0 & 0 & \cdots & 0 & 0 & 0\\
    F(1)-F(0) &         F(0) & 0 & \cdots & 0 & 0 & 0\\
    F(2)-F(1) &    F(1)-F(0) & F(0) & \cdots & 0 & 0 & 0\\
    \vdots & \vdots     & \vdots    & \ddots & \vdots & \vdots & \vdots\\
    F(k-1)-F(k-2) & F(k-2)-F(k-1) & F(k-3)-F(k-2) & \cdots & F(1)-F(0) & F(0) & 0\\
    1-F(k-1)      & 1-F(k-2)      & 1-F(k-3) & \cdots & 1-F(1) & 1-F(0) & 1\\
  \end{array}
\right],
$$

where $n = [n_0,n_1,\cdots,n_{k-1},u]^T$, $n_i$ is the $i^{th}$ pseudo-stage, and $u$ is the subsequent development stage.
