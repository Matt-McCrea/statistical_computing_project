
# Here is all the code to include in the porject

#----------------------------------------------
#             Load Libraries
#----------------------------------------------
library(tidyverse)
library(lubridate)
library(estimatr)
library(sandwich)
library(lmtest)
library(rsample)
library(purrr)
library(broom)

#----------------------------------------------
#              Data Cleaning
#----------------------------------------------
# Load data
demand_data <- read_csv("SCS_demand_modelling.csv")
temp_hourly <- read_csv("SCS_hourly_temp.csv", col_names = FALSE)

# Clean temp_hourly
temp_hourly <- temp_hourly[-1, ]  # remove header row if present again
colnames(temp_hourly) <- c("timestamp", "temp")

temp_hourly <- temp_hourly %>%
  mutate(
    timestamp = parse_date_time(timestamp, orders = "dmy HM"),
    date = as_date(timestamp),
    hour = hour(timestamp),
    temp = as.numeric(temp)
  )

# Clean demand_data
demand_data <- demand_data %>%
  mutate(
    Date = as_date(Date),
    date = Date)

# Define fixed holiday list
holiday_dates <- c("23-12", "24-12", "25-12", "26-12", 
                   "31-12", "01-01", "01-02", "01-03")

# Mark holidays
demand_data <- demand_data %>%
  mutate(
    holiday_code = format(Date, "%d-%m"),
    is_holiday = holiday_code %in% holiday_dates
  )

# Exclude holidays
demand_data <- demand_data %>% filter(!is_holiday)

# Extract temp at 16:00
temp_16 <- temp_hourly %>%
  filter(hour == 16) %>%
  select(date, temp) %>%
  rename(temp_16 = temp)

# Merge in today & previous day temp
demand_data <- demand_data %>%
  left_join(temp_16, by = "date") %>%
  left_join(
    temp_16 %>%
      rename(temp_16_prior = temp_16) %>%
      mutate(date = date + 1),
    by = "date"
  )

# Create weighted temp variable
demand_data <- demand_data %>%
  mutate(TA = (0.4 * temp_16 + 0.6 * temp_16_prior) / 2)
# Convert Date and extract components
demand_data <- demand_data %>%
  mutate(
    winter_month = case_when(
      month(Date) == 11 ~ 1,  # November
      month(Date) == 12 ~ 2,  # December
      month(Date) == 1  ~ 3,  # January
      month(Date) == 2  ~ 4,  # February
      month(Date) == 3  ~ 5,  # March
      TRUE ~ NA_integer_
    )
  )

# Create a grouping variable for month and weekday/weekend
demand_data <- demand_data %>%
  mutate(
    is_weekend = if_else(wday(Date) %in% c(1, 7), "Weekend", "Weekday"),
    month = month(Date)
  )
demand_data <- demand_data %>%
  mutate(
    month = month(Date),
    is_weekend = if_else(wday(Date) %in% c(1, 7), "Weekend", "Weekday")
  )
demand_data <- demand_data %>%
  arrange(Date) %>%
  mutate(lag_demand = lag(demand_gross, default = first(demand_gross)))





#----------------------------------------------
#              Our Model
#----------------------------------------------

# Base Model
m0 <- lm(demand_gross ~ wind + solar_S + temp + as.factor(wdayindex) + as.factor(monthindex), data = demand_data)

# Initial model
initial_model <-lm(demand_gross ~ TA+ 
                     as.factor(wdayindex) + I(DSN^2)+
                     as.factor(start_year) + month(date),
                   data = demand_data)

# Final Model
model_formula <-lm(demand_gross ~ TA+ lag_demand +  
                     as.factor(wdayindex) + I(DSN^2)+
                     as.factor(start_year) + month(date),
                   data = demand_data)




#----------------------------------------------------
#               Heteroskedasticity Test
#----------------------------------------------------

bptest(good_model)

#----------------------------------------------------
#                Autocorrelation Analysis
#----------------------------------------------------
residuals_model <- residuals(good_model)
par(mfrow=c(1,2))

acf(residuals_model, main = "ACF of residuals")
pacf(residuals_model, main = "PACF of residuals")


#----------------------------------------------------
#             Finding the best Temperature variable
#----------------------------------------------------
# generate_temp_variable() creates a weighted average temperature feature using specified hour ranges for today and the previous day.
# evaluate_model() fits a linear model using the generated temperature feature and returns performance metrics (RMSE, R²).
# search_best_combination() iterates over different hour windows and weights to find the temperature feature that best predicts demand.

temp_data <- read_csv("SCS_hourly_temp.csv", col_names = FALSE)
temp_data <- temp_data[-1, ]

colnames(temp_data) <- c("timestamp", "temp")

temp_data <- temp_data %>%
  rename(datetime = timestamp, temp = temp) %>%
  mutate(
    timestamp = parse_date_time(datetime, orders = "dmy HM"),
    date = as_date(timestamp),
    hour = hour(timestamp),
    temp = as.numeric(temp)
  )

demand_data <- read_csv("SCS_demand_modelling.csv")
demand_data <- demand_data %>%
  rename(date = Date) %>%
  mutate(date = as_date(date))


generate_temp_variable <- function(today_range, prev_range, weight_today = 0.5, weight_prev = 0.5) {
  # Average temp for today window
  today_avg <- temp_data %>%
    filter(hour >= today_range[1], hour <= today_range[2]) %>%
    group_by(date) %>%
    summarise(today_temp = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # Average temp for previous day window, shifted forward by 1 day
  prev_avg <- temp_data %>%
    filter(hour >= prev_range[1], hour <= prev_range[2]) %>%
    mutate(date = date + 1) %>%
    group_by(date) %>%
    summarise(prev_temp = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # Combine and weight
  combined <- today_avg %>%
    inner_join(prev_avg, by = "date") %>%
    mutate(
      custom_temp = weight_today * today_temp + weight_prev * prev_temp
    ) %>%
    select(date, custom_temp)
  
  return(combined)
}
evaluate_model <- function(today_range, prev_range, weight_today = 0.5) {
  weight_prev <- 1 - weight_today
  temp_feat <- generate_temp_variable(today_range, prev_range, weight_today, weight_prev)
  
  model_data <- demand_data %>%
    inner_join(temp_feat, by = "date") %>%
    filter(!is.na(demand_gross), !is.na(custom_temp))
  
  if (nrow(model_data) < 10) return(NULL)
  
  model <- lm(demand_gross ~ custom_temp, data = model_data)
  preds <- predict(model, newdata = model_data)
  rmse_val <- sqrt(mean((model_data$demand_gross - preds)^2))
  
  r_squared <- summary(model)$r.squared
  
  return(list(
    today_start = today_range[1],
    today_end = today_range[2],
    prev_start = prev_range[1],
    prev_end = prev_range[2],
    weight_today = weight_today,
    rmse = rmse_val,
    r_squared = r_squared,
    model = list(model) 
  ))
}

search_best_combination <- function(window_size = 3) {
  # All valid start times between 13 and (20 - window_size + 1)
  start_hours <- 13:(20 - window_size + 1)
  today_windows <- map(start_hours, ~c(.x, .x + window_size - 1))
  prev_windows <- today_windows
  weights <- seq(0.1, 0.9, by = 0.1)
  
  combos <- crossing(
    today_range = today_windows,
    prev_range = prev_windows,
    weight_today = weights
  )
  
  # Evaluate
  results <- map_dfr(1:nrow(combos), function(i) {
    row <- combos[i, ]
    res <- evaluate_model(row$today_range[[1]], row$prev_range[[1]], row$weight_today)
    if (!is.null(res)) return(as_tibble(res))
    return(NULL)
  })
  
  # Return best
  best <- results %>% arrange(rmse) %>% slice(1)
  
  cat("Best combo:\n")
  cat("Today hours:", best$today_start, "-", best$today_end, "\n")
  cat("Prev day hours:", best$prev_start, "-", best$prev_end, "\n")
  cat("Weight (today):", best$weight_today, "\n")
  cat("Lowest RMSE:", round(best$rmse, 5), "\n")
  
  return(best)
}

best_result <- search_best_combination(window_size = 1)
summary(best_result$model[[1]])


#------------------------------------------------
#             Newey -West Standard errors to determine significance
#------------------------------------------------

# Reload all data so that it works nicely
temp_hourly <- read_csv("SCS_hourly_temp.csv") 
demand_data <- read_csv("SCS_demand_modelling.csv")

temp_hourly <- temp_hourly %>%
  mutate(datetime = dmy_hm(Date),     
         date = as_date(datetime),    
         hour = hour(datetime))       
# construct the specification selected by our window tool
temp_16 <- temp_hourly %>%
  filter(hour == 16) %>%
  select(date, temp) %>%
  rename(temp_16 = temp)
demand_data <- demand_data %>%
  mutate(
    Date = as_date(Date),
    date = Date)
demand_data$Date <- as.Date(demand_data$Date)

demand_data <- demand_data %>%
  mutate(date = as_date(Date)) %>%  
  left_join(temp_16, by = "date") %>%  
  left_join(
    temp_16 %>% 
      rename(temp_16_prior = temp_16) %>% 
      mutate(date = date + 1),      
    by = "date"
  )
demand_data <- demand_data %>%
  mutate(TA = (0.4 * temp_16 + 0.6 * temp_16_prior) / 2)

demand_data <- demand_data %>%
  arrange(Date) %>%
  mutate(lag_demand = lag(demand_gross, default = first(demand_gross)))
# define the initial model before cutting variables
model_formula <-lm(demand_gross ~ TA+ lag_demand + solar_S+
                     as.factor(wdayindex) + I(DSN^2)+ DSN+
                     as.factor(start_year) + month(date),
                   data = demand_data)

# Compute Newey-West robust standard errors and run the coefficient test
# we take the ourput of this by printing "print(nw_test)" and get rid of variables not significant at the 5% level
nw_se <- NeweyWest(model_formula, prewhite = TRUE)
nw_test <- coeftest(model_formula, vcov = nw_se)


#------------------------------------------------
#        Cross validation schemes (pulled directly from the Rmd)
#------------------------------------------------

# Make both the LOMO and 80:20 for month below

# For an 80% prediction interval
alpha <- 0.2
# making LOMO redefine months to make sure
data_cv1 <- demand_data %>%
  filter(!is_holiday, month(Date) %in% c(11, 12, 1, 2, 3)) %>%
  mutate(monthindex = factor(case_when(
    month(Date) == 11 ~ "Nov",
    month(Date) == 12 ~ "Dec",
    month(Date) == 1  ~ "Jan",
    month(Date) == 2  ~ "Feb",
    month(Date) == 3  ~ "Mar"
  ), levels = c("Nov", "Dec", "Jan", "Feb", "Mar"))) %>%
  arrange(Date)

# Define model formula, follows from tutorial 7
model_formula1 <- demand_gross ~ TA + lag_demand + wday(Date) + month(Date) +
  I(DSN^2) + as.factor(start_year) + DSN

rmse_vec1 <- numeric(length(levels(data_cv1$monthindex)))
int_score_vec1 <- numeric(length(levels(data_cv1$monthindex)))

for (mi in levels(data_cv1$monthindex)) {
  train_set1 <- data_cv1 %>% 
    filter(monthindex != mi) %>%
    mutate(monthindex = factor(monthindex, levels = levels(data_cv1$monthindex)))
  
  test_set1 <- data_cv1 %>% 
    filter(monthindex == mi) %>%
    mutate(monthindex = factor(monthindex, levels = levels(data_cv1$monthindex)))
  
  if(nrow(test_set1) == 0) next
  
  fit1 <- lm(model_formula1, data = train_set1)
  pred1 <- predict(fit1, newdata = test_set1, interval = "prediction", level = 0.8)
  obs1 <- test_set1$demand_gross
  
  rmse_vec1[which(levels(data_cv1$monthindex) == mi)] <- sqrt(mean((obs1 - pred1[,"fit"])^2, na.rm = TRUE))
  
  int_score_vec1[which(levels(data_cv1$monthindex) == mi)] <- mean(
    (pred1[,"upr"] - pred1[,"lwr"]) + (2/alpha) * (
      ((pred1[,"lwr"] - obs1) * (obs1 < pred1[,"lwr"])) +
        ((obs1 - pred1[,"upr"]) * (obs1 > pred1[,"upr"]))
    ),
    na.rm = TRUE
  )
}

results_method1 <- tibble(
  month = levels(data_cv1$monthindex),
  'RMSE Score LOMO' = rmse_vec1,
  'Interval Score LOMO' = int_score_vec1
)

# 80:20 split
data_cv2 <- demand_data %>%
  filter(!is_holiday, month(Date) %in% c(11, 12, 1, 2, 3)) %>%
  mutate(winter_month = case_when(
    month(Date) == 11 ~ "Nov",
    month(Date) == 12 ~ "Dec",
    month(Date) == 1  ~ "Jan",
    month(Date) == 2  ~ "Feb",
    month(Date) == 3  ~ "Mar",
    TRUE ~ NA_character_
  )) %>%
  mutate(winter_month = factor(winter_month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar"))) %>%
  arrange(Date) %>%
  mutate(lag_demand = lag(demand_gross, default = first(demand_gross)))

# same model again
model_formula2 <- demand_gross ~ TA + lag_demand + as.factor(wdayindex) + month(Date) +
  I(DSN^2) + as.factor(start_year)

rmse_vec2 <- numeric(length(levels(data_cv2$winter_month)))
int_score_vec2 <- numeric(length(levels(data_cv2$winter_month)))

for (mi in levels(data_cv2$winter_month)) {
  train_set2 <- data_cv2 %>% 
    filter(winter_month != mi) %>%
    mutate(winter_month = factor(winter_month, levels = levels(data_cv2$winter_month)))
  
  test_set2 <- data_cv2 %>% 
    filter(winter_month == mi) %>%
    mutate(winter_month = factor(winter_month, levels = levels(data_cv2$winter_month)))
  
  if(nrow(test_set2) == 0) next
  
  fit2 <- lm(model_formula2, data = train_set2)
  pred2 <- predict(fit2, newdata = test_set2, interval = "prediction", level = 0.8)
  obs2 <- test_set2$demand_gross
  
  rmse_vec2[which(levels(data_cv2$winter_month) == mi)] <- sqrt(mean((obs2 - pred2[,"fit"])^2, na.rm = TRUE))
  
  int_score_vec2[which(levels(data_cv2$winter_month) == mi)] <- mean(
    (pred2[,"upr"] - pred2[,"lwr"]) + (2/alpha) * (
      ((pred2[,"lwr"] - obs2) * (obs2 < pred2[,"lwr"])) +
        ((obs2 - pred2[,"upr"]) * (obs2 > pred2[,"upr"]))
    ),
    na.rm = TRUE
  )
}

results_method2 <- tibble(
  month = levels(data_cv2$winter_month),
  'RMSE Score 80:20 Split' = rmse_vec2,
  'Interval Score 80:20 Split' = int_score_vec2
)

# Combine the results from both methods by month
combined_results <- left_join(results_method1, results_method2, by = "month")

kable(combined_results,
      digits = 2,
      caption = "Comparison of Leave-One-Month-Out and 80:20 Cross-Validation Results by Month")

# Repeat for weekend/weekday split

alpha <- 0.2  # 80% prediction interval

# Ensure day type is defined
demand_data <- demand_data %>%
  mutate(
    is_weekend = if_else(wday(Date) %in% c(1, 7), "Weekend", "Weekday"),
    is_weekend = factor(is_weekend, levels = c("Weekday", "Weekend"))
  )

# CV by day type 
cv_results_weekday <- demand_data %>%
  group_by(is_weekend) %>%
  group_modify(~ {
    set.seed(123)
    cv_folds <- vfold_cv(.x, v = 5)
    
    fold_results <- cv_folds %>%
      mutate(
        model = map(splits, ~ lm(model_formula, data = analysis(.x))),
        preds = map2(model, splits, ~ suppressWarnings(
          predict(.x, newdata = assessment(.y), interval = "prediction", level = 0.8)
        )),
        actuals = map(splits, ~ assessment(.x)$demand_gross),
        
        rmse = map2_dbl(preds, actuals, ~ sqrt(mean((.y - .x[,"fit"])^2, na.rm = TRUE))),
        
        interval_score = map2_dbl(preds, actuals, ~ mean(
          (.x[,"upr"] - .x[,"lwr"]) + (2/alpha) * (
            ((.x[,"lwr"] - .y) * (.y < .x[,"lwr"])) +
              ((.y - .x[,"upr"]) * (.y > .x[,"upr"]))
          ),
          na.rm = TRUE
        ))
      )
    
    tibble(
      avg_rmse = mean(fold_results$rmse),
      avg_interval_score = mean(fold_results$interval_score)
    )
  }) %>%
  ungroup()

# Display results as a kable table
kable(cv_results_weekday,
      col.names = c("Day Type", "Average RMSE", "Average Interval Score"),
      digits = 2,
      caption = "Cross-Validation Results by Day Type (Weekday vs Weekend)") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))


#------------------------------------------------
#        Using other winters to predict 2013
#------------------------------------------------
# Exclude 2013 winter from training (i.e. use data from 1991 to 2012)
train_data <- demand_data %>%
  filter(year(Date) != 2013, month(Date) %in% c(11, 12, 1, 2, 3))

# Fit the model on the training data
model_excl_2013 <- lm(demand_gross ~ TA + lag_demand +  
                        as.factor(wdayindex) + I(DSN^2) +
                        as.factor(start_year) + month(date),
                      data = train_data)

## Get unique weather years from training data
sim_years <- sort(unique(train_data$start_year))

# store simulated peak predictions
predictions <- tibble(WeatherYear = sim_years, PredictedPeak = NA_real_)

# Loop over each weather year
# For each weather year, compute the average TA and lag_demand for the winter period.
for(i in seq_along(sim_years)) {
  yr <- sim_years[i]
  weather_yr <- demand_data %>% 
    filter(year(Date) == yr, month(Date) %in% c(11,12,1,2,3))
  
  # Compute average values for TA and lag_demand from that winter
  avg_TA <- mean(weather_yr$TA, na.rm = TRUE)
  avg_lag <- mean(weather_yr$lag_demand, na.rm = TRUE)
  # For DSN, we can use the median value as it's in the middle
  med_DSN <- median(weather_yr$DSN, na.rm = TRUE)
  
  # Create a new data frame for prediction
  new_data <- tibble(
    TA = avg_TA,
    lag_demand = avg_lag,
    wdayindex = 3, # we chose a workday to be baseline as they are higher
    DSN = med_DSN,
    start_year = 2013,  # target winter year effect
    date = as.Date("2013-01-15")  
  )
  
  # Predict peak demand using the fitted model
  pred <- predict(model_excl_2013, newdata = new_data)
  predictions$PredictedPeak[i] <- pred
}

# Obtain the actual peak demand for winter 2013/14 from the full data
actual_2013_14 <- demand_data %>%
  filter((year(date) == 2013 & month(date) %in% c(11,12)) | 
           (year(date) == 2014 & month(date) %in% c(1,2,3))) %>%
  pull(demand_gross) %>%
  max(na.rm = TRUE)

# Identify the highest and lowest predicted peaks
max_pred <- max(predictions$PredictedPeak, na.rm = TRUE)
min_pred <- min(predictions$PredictedPeak, na.rm = TRUE)

#plot the results
ggplot(predictions, aes(x = WeatherYear, y = PredictedPeak)) +
  geom_line(size = 1) +
  geom_point(data = filter(predictions, PredictedPeak %in% c(max_pred, min_pred)),
             aes(x = WeatherYear, y = PredictedPeak), color = "red", size = 3) +
  geom_hline(yintercept = actual_2013_14, linetype = "dashed", color = "blue", size = 1) +
  labs(
    title = "Simulated 2013/14 Peak Electricity Demand Based on Historical Weather",
    x = "Weather Year",
    y = "Predicted Peak Demand (MW)",
    caption = "Red dots denote highest and lowest simulated peaks; dashed line shows actual 2013 peak demand."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  )

peak_data <- data.frame(
  Type = c("Highest simulated peak", "Lowest simulated peak", "Actual peak (2013/14)"),
  Demand_MW = c(max_pred, min_pred, actual_2013_14)
)

kable(peak_data, 
      col.names = c("Demand Type", "Value (MW)"),
      digits = 0,  # Round to whole numbers
      caption = "Simulated vs Actual Peak Electricity Demand (2013/14)") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), 
                full_width = FALSE,
                position = "center") %>%
  row_spec(3, bold = TRUE, background = "#E8F4F8")