# Load necessary libraries
library(bblme)  # For MLE-based Exponential Weibull fitting
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the CSV file
df <- read.csv("../data/MetDistRA21/MD_RSG.csv", header = TRUE)[,-1]

# Function to fit Exponential Weibull using MLE and extract parameters
fit_exponential_weibull_mle <- function(column_data, colname) {
  column_data <- na.omit(column_data)  # Remove NA values
  
  if (length(column_data) < 10) {  # Ensure enough data points
    return(data.frame(shape = NA, scale = NA, alpha = NA, column = colname))
  }

  nll_expweibull <- function(alpha, shape, scale) {
    -sum(dexpweibull(column_data, alpha, shape, scale, log = TRUE))
  }

  # Fit Exponential Weibull using bblme
  mle_fit <- mle2(
    minuslogl = nll_expweibull,
    start = list(alpha = 1.5, shape = 1, scale = median(MD_RTLC_03_RSG[,"RSG_ab"])),
    method = "Nelder-Mead", control = list(trace = 3, maxit = 1e4, reltol = 1e-6)
  )
  
  # Extract estimated parameters
  params <- coef(mle_fit)
  
  return(data.frame(alpha = params["alpha"],
                    shape = params["shape"], 
                    scale = params["scale"], 
                    column = colname))
}

# Store parameter estimates for each column
param_table <- data.frame()

# Percentile matching function for observed vs. predicted
plot_percentile_matching <- function(column_data, colname) {
  column_data <- na.omit(column_data)  # Remove NA values
  
  if (length(column_data) < 10) return(NULL)  # Skip if too few data points
  
  nll_expweibull <- function(alpha, shape, scale) {
    -sum(dexpweibull(column_data, alpha, shape, scale, log = TRUE))
  }
  
  # Fit Exponential Weibull using bblme
  mle_fit <- mle2(
    minuslogl = nll_expweibull,
    start = list(alpha = 1.5, shape = 1, scale = median(MD_RTLC_03_RSG[,"RSG_ab"])),
    method = "Nelder-Mead", control = list(trace = 3, maxit = 1e4, reltol = 1e-6)
  )
  
  summary(mle_fit)
  
  # Extract estimated parameters
  params <- coef(mle_fit)
  
  # Generate observed percentiles
  obs_percentiles <- quantile(column_data, probs = seq(0.01, 0.99, 0.01), na.rm = TRUE)
  
  # Generate predicted percentiles from fitted model
  pred_percentiles <- qexpweibull(seq(0.01, 0.99, 0.01),
                               alpha = params["alpha"],
                               shape = params["shape"], 
                               scale = params["scale"])
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    observed = log10(obs_percentiles),
    predicted = log10(pred_percentiles)
  )
  
  # Generate log-log observed vs. predicted plot
  p <- ggplot(plot_data, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste(colname),
         x = "Log10 Observed",
         y = "Log10 Predicted") +
    theme_minimal()
  
  return(p)
}

# Loop over all columns in the dataset
plots <- list()
for (colname in colnames(df)) {
  column_data <- df[[colname]]
  
  # Fit the model and store parameters
  param_result <- fit_exponential_weibull_mle(column_data, colname)
  param_table <- rbind(param_table, param_result)
  
  # Generate the percentile matching plot
  p <- plot_percentile_matching(column_data, colname)
  if (!is.null(p)) plots[[colname]] <- p
}

# Save the parameter table
write.csv(param_table, "Exponential_Weibull_MLE_Parameters.csv", row.names = FALSE)
print("Parameter estimates saved as Exponential_Weibull_MLE_Parameters.csv")

# Display all percentile matching plots
print(plots)
all_plots <- plot_grid(plotlist = plots)+
  theme(plot.background = element_rect(fill = "white", color = NA))

save_plot("../output/ExpWei_RTLC.png", all_plots, base_width = 20, base_height = 20)
