library(tidyquant)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(fredr) # Used to set the key globally for tidyquant/quantmod

# Set the API key for the session using the environment variable we created
fredr_set_key(Sys.getenv("FRED_API_KEY"))

# 1. Define variables and start date
tickers <- c("GDPC1", "TB3MS", "GS10", "NFCI")
start_date <- as.Date("1970-01-04")

# 2. Retrieve Data
df_raw <- tq_get(tickers, get = "economic.data", from = start_date)

# 3. Apply Transformations and Frequency Alignment
df_quarterly <- df_raw %>%
  # Convert all dates to the start of the quarter for alignment
  mutate(date = floor_date(date, "quarter")) %>%
  group_by(symbol, date) %>%
  # Aggregate monthly data to quarterly (averaging)
  summarise(price = mean(price, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = symbol, values_from = price) %>%
  mutate(
    recession_flag = as.numeric(GDPC1 / dplyr::lag(GDPC1) - 1 < 0),
    yield_spread = GS10 - TB3MS,
  ) %>%
  filter(complete.cases(.))

# View the final dataset
print(df_quarterly)


# Build recession intervals (contiguous runs of recession_flag == 1)
recessions <- df_quarterly %>%
  arrange(date) %>%
  mutate(
    is_rec = recession_flag == 1,
    grp = cumsum(!is_rec) # increments when not in recession
  ) %>%
  group_by(grp) %>%
  # keep only groups that contain at least one recession observation
  filter(any(is_rec)) %>%
  summarise(
    start = min(date[is_rec]),
    end = max(date[is_rec]),
    .groups = "drop"
  )

print(recessions)

# Plot GDP with shaded recession periods
ggplot(df_quarterly, aes(x = date, y = GDPC1)) +
  geom_rect(
    data = recessions,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(color = "steelblue", linewidth = 0.6) +
  labs(
    x = "Date",
    y = "Real GDP (GDPC1)",
    title = "Real GDP with Recession Periods Highlighted"
  ) +
  theme_minimal()

ggplot(df_quarterly, aes(x = date, y = yield_spread)) +
  geom_rect(
    data = recessions,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(color = "steelblue", linewidth = 0.6) +
  labs(
    x = "Date",
    y = "Real GDP (GDPC1)",
    title = "Real GDP with Recession Periods Highlighted"
  ) +
  theme_minimal()

# Save the final dataset
df_quarterly <- df_quarterly %>%
  select(date, recession_flag, NFCI, yield_spread)
save(df_quarterly, file = 'data/df_quarterly.Rdata')
