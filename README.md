# Indometh Models

## Introduction
One-compartment model is the conceptual foundation of PK. It introduces the core parameters used to describe how a drug moves through the body. This model was developed for drugs that distribute into tissues almost instantly (relative to how slowly they are eliminated) (NLM, 2023). In this project, the Indometh dataset represents intravenous (IV) administration. After an IV bolus, where there is no absorption phase because the drug is given directly into the blood. This model can be used for companies to simulate scenarios to see how a drug will behave across diverse populations (Weiss, 2023).

## Objective
The objective of this project is to perform a comparative analysis of single-compartment and two-compartment pharmacokinetic models (Pinheiro, 2026). By evaluating key performance metrics, this study aims to determine the most mathematically robust and physiologically adequate model for Indometh.

## Development
The function nlme() in the library(nlme)performs nonlinear regression on plasma concentration versus time to estimate the elimination rate constant (k_el), elimination half-life (t1/2), initial concentration (C0), apparent volume of distribution (Vd), and clearance (CL) to calculate the population-level averages (fixed effects) and individual-level deviations (random effects). 

### 1-Compartment model
Due to the sample size and data density limitations inherent to the Indometh dataset, a parsimonious modeling strategy was adopted. The structural model was restricted to minimal essential parameters to ensure numerical stability, prevent overfitting, and achieve robust model convergence within the nlme framework:

First-order kinetic equation:

$C(t) = C_0 \cdot e^{-k_{el} \cdot t}$.

Parameters:

- C(t): plasma concentration at time “t”.
- C0: is the initial concentration.
- -kel: is the elimination rate constant. Negative because the drug is removed from the body per unit of time ($hr^{-1}$).

Also, the Pharmacokinetic parameters will be investigated:

- Volume of Distribution (Vd): Relates the amount of drug in the body to the concentration measured in the blood. Vd = CL / k_el
- Clearance (CL): The volume of blood cleared of the drug per unit of time. CL = Dose / AUC
- Elimination Half-Life (t1/2): The time required for the plasma concentration to decrease by 50%. t1/2 = ln(2) / k_el

Figure 1. Image done by Gemini AI to explain the First-order kinetic equation and One-compartment model

<img width="929" height="498" alt="Screenshot 2026-06-10 at 7 32 09 PM" src="https://github.com/user-attachments/assets/90bcde49-20b3-4c11-885b-e95f8dd62013" />

### 2-Compartment model

To build a custom two-compartment intravenous (IV) bolus model on the original linear scale, there are two primary mathematical approaches: the Macro-constant equations (geometric/exponential components) and the Micro-constant equations (direct transfer rates between compartments). The Micro-constant approach is highly preferred because the parameters translate directly to physical biological properties (Pinheiro & Bates, 2000).

In a two-compartment model, the body is divided into:
- The Central Compartment (V1): Highly perfused tissues (blood, organs) where the drug is injected and eliminated.
- The Peripheral Compartment (V2): Deeper tissues where the drug distributes back and forth.

The micro-constant approach is a system of linear differential equations that maps the physical movement of the drug between the central and peripheral spaces, and the calculations down into a step-by-step matrix sequence (Pinheiro & Bates, 2000).
The micro-constant equation defines drug disposition using clearance and volumes directly. The concentration in the central compartment over time is expressed as a bi-exponential decline:

**Step 1 - The Fundamental Rate Matrix (K) **
It describes the fractional distribution of the drug mass over time. Find the constants:

$k_{10} = \frac{CL}{V_1}$, $k_{12} = \frac{Q}{V_1}$, and $k_{21} = \frac{Q}{V_2}$.

- k10 is the elimination rate constant from the central compartment.
- k12 is the transfer rate constant from the central to the peripheral compartment.
- k21 is the transfer rate constant from the peripheral back to the central compartment.

**Step 2 - Solving for the Hybrid Exponents (alpha and beta)**
To find the the steepness of the curves (inclination), the algebraic solution resolves the quadratic characteristic equation of the rate matrix:

$$
\begin{aligned}
\alpha + \beta &= k_{10} + k_{12} + k_{21} \\
\alpha \cdot \beta &= k_{10} \cdot k_{21}
\end{aligned}
$$

Using the quadratic formula, the literal algebraic values for alpha and beta are extracted:

$$
\begin{aligned}
\alpha &= \frac{(k_{10} + k_{12} + k_{21}) + \sqrt{(k_{10} + k_{12} + k_{21})^2 - 4 \cdot k_{10} \cdot k_{21}}}{2} \\
\beta &= \frac{(k_{10} + k_{12} + k_{21}) - \sqrt{(k_{10} + k_{12} + k_{21})^2 - 4 \cdot k_{10} \cdot k_{21}}}{2}
\end{aligned}
$$

**Step 3 - The Structural Concentration Equation**
Once you find alpha and beta, the final concentration in the central compartment, C(t), back to the initial dose (D) and the central volume (V1):

$$
C(t) = \frac{D}{V_1 \cdot (\alpha - \beta)} \cdot \left[ (\alpha - k_{21}) \cdot e^{-\alpha \cdot t} + (k_{21} - \beta) \cdot e^{-\beta \cdot t} \right]
$$

The final equation translates the population's biological parameters (Dose and Central Volume) into a time-dependent curve, accurately predicting the real-world plasma concentration at any given post-dose interval.

Figure 2. Image done by Gemini AI to explain the system of linear differential equations and Two-compartment model

<img width="932" height="504" alt="Screenshot 2026-06-11 at 10 57 45 AM" src="https://github.com/user-attachments/assets/28dc76a9-f7fc-4efc-ac2d-dc9562c37a64" />


# Results

## Model 1 (1-Compartment) plot

<img width="966" height="464" alt="Screenshot 2026-06-10 at 7 33 27 PM" src="https://github.com/user-attachments/assets/a55e0fbf-dd18-41e8-95c6-b7cc87a6f614" />

## Model 2 (2-Compartment) plot

<img width="1200" height="564" alt="Screenshot 2026-06-11 at 11 35 43 AM" src="https://github.com/user-attachments/assets/c119fd59-ab22-442b-a2e5-7a0ca5cad96f" />


### How much does my model predict correctly?
| Accuracy Tier | Model 1 (1-Compartment) | Model 2 (2-Compartment) |
| :--- | :---: | :---: |
| **Highly Accurate (<10% error)** | 18 (27.27%) | 31 (46.97%)|
| **Acceptable (10-20% error)** | 9 (13.64%) | 11(16.67%) |
| **Marginal (20-30% error)** | 1 (1.52%) | 5(7.58%) |
| **Inaccurate (>30% error)** | 38 (57.58%) | 19(28.79%) |

**Plot 1-Compartment model**
<img width="992" height="521" alt="Screenshot 2026-06-10 at 7 36 07 PM" src="https://github.com/user-attachments/assets/41ccbc13-813b-4641-8566-23d3f5a55f2b" />

**Plot 2-Compartment model**
<img width="2232" height="1418" alt="image" src="https://github.com/user-attachments/assets/b3935b67-a8cc-4e4b-9eef-a3f41b9f8455" />



### Metrics
| Metric | Model 1 (1-Compartment) | Model 2 (2-Compartment) |
| :--- | :---: | :---: |
| **AIC** | -41.65 | -97.733 |
| **BIC** | -28.51 | -80.21 |
| **LogLik** | 26.82 | 56.86|
| **Residuals** | 0.138 | 0.079 |

Form model 2, fixed effects establish the baseline human baseline parameters. Random effects represent the average customer. 
 
- Central Volume (V1 = 7.48L)
- Peripheral Volume (V2 = 7.68L)
- Clearance (CL = 9.16L/h), individual variation (StdDev) between patients is 2.28L/h.
- Intercompartmental Clearance (Q = 6.45L/h)
- Correlation (0.628)
- Residual Error = $0.079$

The volumes can be used to map the initial dilution of the drug. As the volume results are similar, we can conclude that Indomethacin divides evenly between the bloodstream and body tissues at equilibrium - safely calculate starting doses for human clinical protocols. [Preclinical to Clinical Scaling]

Because V1 is small, a rapid 25 mg IV bolus will cause an immediate concentration spike in the blood. However, because Q is higher (6.45 L/h), the drug moves very quickly out of the blood and into the tissues.
This movement says that once the drug moves into the deeper tissues (where the therapeutic targets for inflammation or pain usually reside), it binds (connect) and lingers(stay) there, supporting extended dosing intervals.

The clearance standard deviation means some patients clear the drug nearly 25%* faster or slower than average. The correlation shows that patients with larger central volumes(V1) also tend to clear the drug faster (this is the exact justification needed to investigate body weight or kidney function as a covariate, which eventually dictates personalized dosing instructions on the final product label).
Residual error low, permits use Monte Carlo simulation.
*((Random effect SD / Fixed effect CL) *100)


Moving from a 1-compartment model (Model 1) to a 2-compartment model (Model 2) optimizes the statistical fit and biological accuracy.

Decreased Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC). Also, the log-likelihood (LogLik) increased.

The residual variance dropped by nearly half, which means the distance between the model's predictions and the actual patient plasma concentrations has shrunk.

Transitioning to a micro-constant 2-compartment approach successfully handles this structural complexity while yielding much more precise, physiologically accurate parameter estimates for clearance and volume.


## Practical aplication

- The optimization for V1 (Central Volume), V2 (Peripheral Volume), and Q (Intercompartmental Clearance) by this model permit accurately simulate human plasma concentration curves before the first dose is ever given to a human volunteer.


# Source

Manlapaz, P. A. C. (2026). pkpd.Release: Model fitting and simulation for drug release kinetics and PK/PD (Version 0.1.0) [Computer software]. In: One-Compartment IV Bolus Pharmacokinetic Model (Linear). R CRAN. https://cran.r-project.org/web/packages/pkpd.Release/pkpd.Release.pdf

National Library of Medicine. (2023). Indomethacin. In StatPearls. NCBI Bookshelf. https://www.ncbi.nlm.nih.gov/books/NBK557744/

Owen, J. S., & Fiedler-Kelly, J. (2014). Introduction to Population Pharmacokinetic/Pharmacodynamic Analysis with Nonlinear Mixed Effects Models. Chapter 11: Simulation Basis. John Wiley & Sons. https://onlinelibrary.wiley.com/doi/book/10.1002/9781118784860 

Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team.(2026). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1. Available from: https://CRAN.R-project.org/package=nlme

Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications. 4th ed. Philadelphia, PA: Lippincott Williams & Wilkins; 2011. (Chapter 8: Distribution Kinetics). https://downloads.lww.com/wolterskluwer_vitalstream_com/sample-content/9780781750097_rowland/samples/frontmatter.pdf 

Weiss M. (2023). Is the One-Compartment Model with First Order Absorption a Useful Approximation?. Pharmaceutical research, 40(9), 2147–2153. https://doi.org/10.1007/s11095-023-03582-1 





