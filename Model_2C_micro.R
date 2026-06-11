############################# TWO-COMPARTMENT MODEL (NLME 2C)###################
# Load dataset
data("Indometh")

head(Indometh)

# Load library
library(nlme)
library(ggplot2)

# Convert the built-in grouped data into a clean, standalone data frame
my_pk_data <- as.data.frame(Indometh)

# Define the IV Bolus formula
# Note: Indometh provides a 25 mg IV bolus dose to all subjects
calc_2comp_conc <- function(time, CL, V1, Q, V2) {
  dose_mg <- 25
  # STEP 1 - THE FUNDAMENTAL RATE MATRIX (K)
  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2
  
  # STEP 2 - SOLVING FOR THE HYBRID EXPONENTS (alpha and beta)
  sum_k  <- k10 + k12 + k21
  prod_k <- k10 * k21
  disc   <- sqrt(sum_k^2 - 4 * prod_k)
  
  alpha  <- (sum_k + disc) / 2
  beta   <- (sum_k - disc) / 2
  
  # STEP 3 - THE STRUCTURAL CONCENTRATION EQUATION
  A <- (dose_mg / V1) * ((alpha - k21) / (alpha - beta))
  B <- (dose_mg / V1) * ((k21 - beta) / (alpha - beta))
  
  # Render the bi-exponential curve on the original linear scale
  A * exp(-alpha * time) + B * exp(-beta * time)
}

# Reference the helper function inside your nlme formula
pk_formula_2comp <- conc ~ calc_2comp_conc(time, CL, V1, Q, V2)

# NLME model
twocomp_micro_model <- nlme(
  model = pk_formula_2comp,
  data = my_pk_data,
  # We optimize for the 4 physiological parameters
  fixed = CL + V1 + Q + V2 ~ 1,
  # Subject-level random effects (due to sample size constraints, 
  # you can start with random effects on core parameters like CL and V1)
  random = CL + V1 ~ 1 | Subject,
  # Custom formulas require explicit starting values (approximated from Model 1)
  start = c(CL = 12, V1 = 9, Q = 5, V2 = 15),
  control = nlmeControl(maxIter = 100, msMaxIter = 100)
)

# Quality metrics
summary(twocomp_micro_model)


### PLOT----

## Plot the 2-Compartment model----

# Individual predictions (Population trend + Subject variations)
my_pk_data$pred_individual <- predict(twocomp_micro_model)

# Population-only predictions (The average fixed-effect curve)
my_pk_data$pred_population <- predict(twocomp_micro_model, level = 0)

# The True Nonlinear Mixed Model Fit
p2 <- ggplot(data = my_pk_data, aes(x = time, y = conc)) +
  # Actual blood sample measurements
  geom_point(color = "black", size = 1.5) +
  # Model-predicted individual curves (Random Effects)
  geom_line(aes(y = pred_individual, group = Subject), color = "grey80") +
  # Model-predicted population curve (Fixed Effect)
  geom_line(aes(y = pred_population), color = "blue", linewidth = 1) +
  # Aesthetics
  theme_classic() +
  labs(
    title = "Model 2 (2-Compartment)",
    subtitle = "Custom Two-Compartment Micro-Constant NLME Fit",
    x = "Time (h)", 
    y = "Concentration (mg/L)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )
# Visualization
print(p2)


### Accuracy----
# Observed x predict

# Absolute percentage error for each 2-compartment prediction
my_pk_data$percent_error <- abs((my_pk_data$conc - my_pk_data$pred_individual) / my_pk_data$conc) * 100

# Ranges
my_pk_data$Accuracy_Tier <- cut(
  my_pk_data$percent_error,
  breaks = c(0, 10, 20, 30, Inf),
  labels = c("Highly Accurate (<10% error)", 
             "Acceptable (10-20% error)", 
             "Marginal (20-30% error)", 
             "Inaccurate (>30% error)")
)

# Confusion Matrix
accuracy_matrix_2comp <- table(my_pk_data$Accuracy_Tier)
  # Convert to percentages
accuracy_percentage_2comp <- prop.table(accuracy_matrix_2comp) * 100

# Print the matrix results 
print("Number of 2-Compartment Predictions per Tier:")
print(accuracy_matrix_2comp)

print("Percentage of 2-Compartment Predictions per Tier:")
print(round(accuracy_percentage_2comp, 2))

### Plot the accuracy----
ggplot(my_pk_data, aes(x = conc, y = pred_individual)) +
  # Identity line (Perfect agreement where y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "orange", linewidth = 0.8) + 

  # Scatter observed vs predicted, color-coded by accuracy tier
  geom_point(aes(color = Accuracy_Tier), size = 2.5, alpha = 0.8) +
  # 
  theme_classic(base_size = 12) +
  scale_color_brewer(palette = "RdYlGn", direction = -1) + # Nice Red-Yellow-Green scale for tiers
  labs(
    title = "2-Compartment Model Accuracy Matrix",
    subtitle = "Observed vs. Predicted Plasma Concentrations",
    x = "Actual Observed Concentration (mg/L)",
    y = "Model Predicted Concentration (mg/L)",
    color = "Prediction Accuracy"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right",
    legend.direction = "vertical"
  )