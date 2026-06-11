# Indometh Models

# Introduction
One-compartment model is the conceptual foundation of PK. It introduces the core parameters used to describe how a drug moves through the body. This model was developed for drugs that distribute into tissues almost instantly (relative to how slowly they are eliminated) (NLM, 2023). In this project, the Indometh dataset represents intravenous (IV) administration. After an IV bolus, where there is no absorption phase because the drug is given directly into the blood. This model can be used for companies to simulate scenarios to see how a drug will behave across diverse populations (Weiss, 2023).

# Objective
The objective of this project is to perform a comparative analysis of single-compartment and two-compartment pharmacokinetic models (Pinheiro, 2026). By evaluating key performance metrics, this study aims to determine the most mathematically robust and physiologically adequate model for Indometh.

# Development
The function nlme() in the library(nlme)performs nonlinear regression on plasma concentration versus time to estimate the elimination rate constant (k_el), elimination half-life (t1/2), initial concentration (C0), apparent volume of distribution (Vd), and clearance (CL) to calculate the population-level averages (fixed effects) and individual-level deviations (random effects). 

## 1-Compartment model
Due to the sample size and data density limitations inherent to the Indometh dataset, a parsimonious modeling strategy was adopted. The structural model was restricted to minimal essential parameters to ensure numerical stability, prevent overfitting, and achieve robust model convergence within the nlme framework:

First-order kinetic equation:

$C(t) = C_0 \cdot e^{-k_{el} \cdot t}$.

Parameters:

- C(t): plasma concentration at time “t”.
- C0: is the initial concentration.
- k_el: is the elimination rate constant. Negative.

Also, the Pharmacokinetic parameters will be investigated:

- Volume of Distribution (Vd): Relates the amount of drug in the body to the concentration measured in the blood. Vd = CL / k_el
- Clearance (CL): The volume of blood cleared of the drug per unit of time. CL = Dose / AUC
- Elimination Half-Life (t1/2): The time required for the plasma concentration to decrease by 50%. t1/2 = ln(2) / k_el

Figure 1. Image done by Gemini AI to explain the First-order kinetic equation and One-compartment model

<img width="929" height="498" alt="Screenshot 2026-06-10 at 7 32 09 PM" src="https://github.com/user-attachments/assets/90bcde49-20b3-4c11-885b-e95f8dd62013" />

## 2-Compartment model

# Results

## Model 1 (1-Compartment) 

<img width="966" height="464" alt="Screenshot 2026-06-10 at 7 33 27 PM" src="https://github.com/user-attachments/assets/a55e0fbf-dd18-41e8-95c6-b7cc87a6f614" />

### How much does my model predict correctly?
| Accuracy Tier | Model 1 (1-Compartment) | Model 2 (2-Compartment) |
| :--- | :---: | :---: |
| **Highly Accurate (<10% error)** | 18 (27.27%) | |
| **Acceptable (10-20% error)** | 9 (13.64%) | |
| **Marginal (20-30% error)** | 1 (1.52%) | |
| **Inaccurate (>30% error)** | 38 (57.58%) | |

**Plot 1-Compartment model**
<img width="992" height="521" alt="Screenshot 2026-06-10 at 7 36 07 PM" src="https://github.com/user-attachments/assets/41ccbc13-813b-4641-8566-23d3f5a55f2b" />




