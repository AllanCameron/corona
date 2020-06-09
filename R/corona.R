populations <- structure(list(Region = c("Ayrshire and Arran", "Borders",
                                         "Dumfries and Galloway",
"Fife", "Forth Valley", "Grampian", "Greater Glasgow and Clyde",
"Highland", "Lanarkshire", "Lothian", "Orkney", "Shetland", "Tayside",
"Western Isles"), Population = c(369360, 115510, 148860, 373550,
306640, 585700, 1183120, 321700, 661900, 907580, 22270, 22920,
417470, 26720)), row.names = c(NA, -14L), class = "data.frame")

hb_map <- sf::read_sf(paste0("https://raw.githubusercontent.com/mbronsvo/",
                             "covid19_scotland/master/SG_NHS_HealthBoards",
                             "_2019c.geojson"))

get_scotdf <- function()
{
  base_url   <- "https://raw.githubusercontent.com/DataScienceScotland/"
  folder_url <- "COVID-19-Management-Information/master/"
  common_url <- "COVID19 - Daily Management Information - "
  big_url    <- paste0(base_url, folder_url, common_url)

  urls <- list(
    Scotland_Deaths = "Scotland - Deaths.csv",
    Scotland_Tests  = "Scotland - Testing.csv",
    Scotland_IP     = "Scotland - Hospital care.csv",
    Regional_ICU    = "Scottish Health Boards - ICU patients.csv",
    Regional_Cases  = "Scottish Health Boards - Cumulative cases.csv",
    Regional_IP = "Scottish Health Boards - Hospital patients - Confirmed.csv")

  urls <- lapply(urls, function(x) gsub(" ", "%20", paste0(big_url, x)))

  df_list <- lapply(urls, read.csv, stringsAsFactors = FALSE)
  df_list <- lapply(df_list, function(x) {x$Date <- as.POSIXct(x$Date); x})

  df_list$Scotland_Tests %>%
    mutate(Negative = c(0, diff(.[[2]])),
           Positive = c(0, diff(.[[3]])),
           Tests    = c(0, diff(.[[4]])))   %>%
    select(Date, Negative, Positive, Tests) %>%
    filter(seq_along(Date) != 1) %>%
    left_join(df_list$Scotland_Deaths, "Date") %>%
    rename(Deaths = Number.of.COVID.19.confirmed.deaths.registered.to.date) %>%
    mutate(Deaths = c(0, diff(replace(Deaths, is.na(Deaths), 0)))) %>%
    left_join(df_list$Scotland_IP %>%
               select(Date,
                      Inpatients = COVID.19.patients.in.hospital...Confirmed,
                      IP_Tot = COVID.19.patients.in.hospital...Total,
                      ICU = COVID.19.patients.in.ICU...Confirmed,
                      ICU_Tot = COVID.19.patients.in.ICU...Total) %>%
               mutate(Inpatients = ifelse(is.na(Inpatients), IP_Tot, Inpatients),
                      ICU = ifelse(is.na(ICU), ICU_Tot, ICU)) %>%
               select(Date, Inpatients, ICU), "Date") %>%
    mutate(Inpatients = ifelse(is.na(Inpatients), 0, Inpatients),
           ICU = ifelse(is.na(ICU), 0, ICU))
}

healthboard <- function()
{
  base_url   <- "https://raw.githubusercontent.com/DataScienceScotland/"
  folder_url <- "COVID-19-Management-Information/master/"
  common_url <- "COVID19 - Daily Management Information - "
  big_url    <- paste0(base_url, folder_url, common_url)

  urls <- list(
    Scotland_Deaths = "Scotland - Deaths.csv",
    Scotland_Tests = "Scotland - Testing.csv",
    Scotland_IP = "Scotland - Hospital care.csv",
    Regional_ICU = "Scottish Health Boards - ICU patients.csv",
    Regional_Cases = "Scottish Health Boards - Cumulative cases.csv",
    Regional_IP = "Scottish Health Boards - Hospital patients - Confirmed.csv")

  urls <- lapply(urls, function(x) gsub(" ", "%20", paste0(big_url, x)))

  df_list <- lapply(urls, read.csv, stringsAsFactors = FALSE)
  df_list <- lapply(df_list, function(x) {x$Date <- as.POSIXct(x$Date); x})

suppressWarnings(df_list$Regional_Cases %>%
                   mutate_if(is.character, as.numeric)) %>%
  mutate_if(is.numeric, function(x) replace(x, is.na(x), 0)) %>%
  pivot_longer(names(.)[-1]) %>%
  rename(Region = name, Cases = value) %>%
  mutate(Region = gsub("[.]", " ", Region)) %>%
  left_join(suppressWarnings(df_list$Regional_IP %>%
                             mutate_if(is.character, as.numeric)) %>%
              mutate_if(is.numeric, function(x) replace(x, is.na(x), 0)) %>%
              pivot_longer(names(.)[-1]) %>%
              rename(Region = name, Inpatients = value) %>%
              mutate(Region = gsub("[.]", " ", Region)), c("Date", "Region")) %>%
  mutate(Inpatients = ifelse(is.na(Inpatients), 0, Inpatients)) %>%
  left_join(suppressWarnings(df_list$Regional_ICU %>%
                             mutate_if(is.character, as.numeric)) %>%
              mutate_if(is.numeric, function(x) replace(x, is.na(x), 0)) %>%
              pivot_longer(names(.)[-1]) %>%
              rename(Region = name, ICU = value) %>%
              mutate(Region = gsub("[.]", " ", Region)), c("Date", "Region")) %>%
  mutate(ICU = ifelse(is.na(ICU), 0, ICU)) %>%
    left_join(populations, by = "Region") %>%
    rename(per_10000 = Population) %>%
    mutate(per_10000 = 10000 * Cases/per_10000)
}

prediction_plot <- function(scotdf, n = 14, logplot = FALSE)
{
  scotdf$cum_pos <- cumsum(scotdf$Positive)
  scotdf[which(scotdf$Date > as.POSIXct("2020-03-15") &
               scotdf$Date < as.POSIXct("2020-03-25")),] -> earlydf
  earlydf$Day <- 0:8

  earlydf$log_cases <- log(earlydf$cum_pos)
  mod2 <- lm(log_cases ~ Day, earlydf)
  coef <- summary(mod2)$coefficients
  mid <- c(earlydf$cum_pos, round(exp(earlydf$log_cases[9] + (1:n) * coef[2,1])))
  lower <- c(earlydf$cum_pos, round(exp(earlydf$log_cases[9] + (1:n) * (coef[2,1] - 1.96 * coef[2,2]))))
  upper <- c(earlydf$cum_pos, round(exp(earlydf$log_cases[9] + (1:n) * (coef[2,1] + 1.96 * coef[2,2]))))
  dates <- c(earlydf$Date, max(earlydf$Date) + lubridate::days(1:n))
  d <- data.frame(Date = dates, mid, lower, upper)
  all_dates <- scotdf$Date[scotdf$Date > as.POSIXct("2020-03-15")]
  label_dates <- all_dates[(seq_along(all_dates) %% 3) == 1]

  p <- ggplot(d, aes(Date, mid)) +
    geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.5, fill = "red") +
    geom_line() +
    geom_line(aes(y = lower), linetype = 2) +
    geom_line(aes(y = upper), linetype = 2) +
    geom_point(data = scotdf, aes(y = cum_pos)) +
    labs(title = "Predicted (line) versus actual (point) Scottish Covid-19 Cases",
         y = "Confirmed Cases") +
        theme(text = element_text(size = 16),
          legend.position = "none",
          plot.margin = unit(c(20, 20, 20, 20), "pt"),
          axis.title.x = element_text(vjust = -4),
          axis.title.y = element_text(vjust = 4)) +
    scale_x_datetime(limits = c(as.POSIXct("2020-03-16"), max(d$Date)),
                     breaks = label_dates,
                     labels = strftime(label_dates, "%d-%b")) +
    geom_vline(aes(xintercept=as.POSIXct("2020-03-24")), linetype=2)
  if(logplot) print(p + scale_y_log10(limits=c(100,10000))) else print(p)
  return(d)
}

three_day_columns <- function(df, area = "All", start = as.POSIXct("2020-03-14"),
                              stop = lubridate::now())
{
  df %>%
  dplyr::filter(Date >= start & Date <= stop & Region == area) %>%
  mutate(Count = c(0, diff(Cases))) -> df
  df <- df[-1,]
  df %>%
    mutate(trip = (lubridate::yday(Date) - min(lubridate::yday(Date))) %/% 3) %>%
    group_by(trip) %>%
    summarise(Count = sum(Count), Date = min(Date)) -> df
  ggplot(df, aes(Date, Count, fill = Count)) +
    geom_col(alpha = 0.8, colour = "black") +
    scale_fill_continuous(low = "forestgreen", high = "red") +
    labs(y = "New cases per three days",
         x = "3-day period beginning",
         title = paste("New COVID-19 cases in", area, "per three days")) +
    theme(text = element_text(size = 16),
          legend.position = "none",
          plot.margin = unit(c(20, 20, 20, 20), "pt"),
          axis.title.x = element_text(vjust = -4),
          axis.title.y = element_text(vjust = 4)) +
    scale_x_datetime(breaks = df$Date, labels = strftime(df$Date, "%d-%b")) -> p
  print(p)
  return(df)
}

trajectory <- function(df, n = 7)
{
  df %>%
    group_by(Region) %>%
    mutate(new_cases = zoo::rollmean(c(0, diff(Cases)), n, na.pad = T, align = "right")) %>%
    ggplot(aes(Cases, new_cases, colour = Region)) +
    geom_line() + scale_x_log10() + scale_y_log10()
}

week_on_week <- function(scotdf)
{
  scotdf$week_on_week <- c(scotdf$Positive[1:8],
                           (scotdf$Positive / lag(scotdf$Positive, 7))[-(1:8)])
  ggplot(scotdf[-(1:8),], aes(Date, week_on_week, fill = week_on_week)) +
    geom_col() +
    scale_fill_continuous(low = "forestgreen", high = "red") +
    geom_hline(aes(yintercept = 1)) + scale_y_log10()
}

weekly <- function(scotdf)
{
  scotdf %>%
    mutate(Week = lubridate::week(Date + lubridate::days(2)) - 10) %>%
    filter(Week > 0) %>%
    group_by(Week) %>%
    summarize(Cases = sum(Positive), Date = min(Date) + lubridate::days(3) + hours(12)) %>%
    ggplot(aes(Date, Cases, fill = Cases)) +
    geom_rect(aes(xmin = as.POSIXct("2020-03-24"), xmax = lubridate::now(),
          ymin = 0, ymax = Inf),
      inherit.aes = FALSE, fill = "#3355FF06") +
    geom_col() +
    labs(title = "Positive Coronavirus Cases per week in Scotland",
         y = "Positive Tests per week (Monday - Sunday)") +
    scale_fill_continuous(low = "forestgreen", high = "red") +
    theme(text = element_text(size = 16),
          legend.position = "none",
          plot.margin = unit(c(20, 20, 20, 20), "pt"),
          axis.title.x = element_text(vjust = -4),
          axis.title.y = element_text(vjust = 4)) +
    scale_y_continuous(limits = c(0, 3000)) +
    geom_text(inherit.aes = FALSE, aes(x = as.POSIXct("2020-04-01"), y = 2800,
              label = "Lockdown starts"), size = 6, hjust = 0) +
    geom_segment(inherit.aes = FALSE, aes(x = as.POSIXct("2020-03-31"),
                                          xend= as.POSIXct("2020-03-24"),
                                          y = 2800, yend = 2800),
                 arrow = arrow(length = unit(0.15, "inches"), type = "closed"))
}

Regional <- function(df, rolling = 7, per_10000 = FALSE)
{
  if(per_10000) df$Cases <- df$per_10000
  df %>%
    filter(Region != "Orkney" & Region != "Western Isles" & Region != "Shetland") %>%
    group_by(Region) %>%
    mutate(daily = c(0, diff(Cases))) %>%
    mutate(roll = zoo::rollmean(daily, rolling, na.pad = TRUE)) %>%
    ggplot(aes(Date, daily)) +
      geom_col(aes(fill = daily, colour = daily)) +
      geom_line(aes(y = roll)) +
      facet_wrap(.~Region, scales = if(per_10000) "fixed" else "free_y") +
      scale_fill_continuous(low = "forestgreen", high = "red") +
      scale_colour_continuous(low = "forestgreen", high = "red") +
      labs(fill = "Daily\nCases",
           title = "Regional positive Covid tests with 7-day rolling average") +
      theme(axis.title.y = element_blank(), legend.position = "none",
            text = element_text(size = 14))
}

test_plot <- function(scotdf, grad = 0.4)
{
  grad_df <- data.frame(d = rep(c(rep(min(scotdf$Date), 2), Inf, Inf), 40),
                        h = rep(c(0, .5, .5, 0), 40) + rep(0:39/2, each = 4),
                        g = rep(1:40, each = 4),
                        a = rep(seq(grad, 0, length.out = 40), each = 4))
  scotdf %>%
    mutate(Tests = Positive + Negative, Percent = 100 * Positive/Tests) %>%
    filter(seq(nrow(.)) > 1) %>%
    ggplot(aes(Date, Percent, fill = Percent)) +
    geom_col(width = 86400, position = "identity") +
    scale_fill_continuous(low = "forestgreen", high = "red") +
    labs(y = "Percentage of tests that were positive",
         title = "Percentage of Scottish tests for Coronavirus that were positive") +
    theme(legend.position = "none",
          text = element_text(size = 20),
          plot.margin = margin(30, 30, 30, 30, "pt"),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(vjust = 5),
          plot.title  = element_text(vjust = 5)) +
    geom_polygon(data = grad_df, inherit.aes = FALSE,
                 aes(x = d, y = h, group = g, alpha = a), fill = "white") +
    scale_alpha_identity() +
    scale_x_datetime(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 40))
}

lockdown_plot <- function(scotdf, max_date = as.POSIXct("2020-06-15"))
{
ggplot() +
  geom_polygon(aes(x = as.POSIXct(c("2020-03-23", "2020-03-23",
                                    "2020-05-28", "2020-05-28")),
                   y = c(0, 600, 600, 0)), fill = "dodgerblue", alpha = 0.2) +
  geom_col(data = scotdf, mapping = aes(Date, Positive, fill = Positive),
           width = 86400) +
  geom_line(data = scotdf %>% mutate(roll = zoo::rollmean(Positive, 7, fill = NA)),
            size =1, mapping = aes( linetype = "", x = Date, y = roll)) +
  scale_fill_continuous(high = "red", low = "forestgreen", guide = F) +
  scale_x_datetime(limits = c(as.POSIXct("2020-03-01"), max_date)) +
  coord_cartesian(ylim = c(0, 500)) +
  labs(x = "Date", y = "Positive tests per day",
       title = "Daily positive Coronavirus tests in Scotland") +
  scale_linetype(name = "7-day rolling\naverage") +
  geom_text(inherit.aes = FALSE, aes(x = as.POSIXct("2020-04-25"), y = 475,
              label = "Lockdown"), size = 6, hjust = 0.5) +
    geom_segment(inherit.aes = FALSE, aes(x = as.POSIXct("2020-04-10"),
                                          xend= as.POSIXct("2020-03-24"),
                                          y = 475, yend = 475),
                 arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
      geom_segment(inherit.aes = FALSE, aes(x = as.POSIXct("2020-05-10"),
                                          xend= as.POSIXct("2020-05-28"),
                                          y = 475, yend = 475),
                 arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
      theme(legend.position = "none",
          text = element_text(size = 20),
          plot.margin = margin(30, 30, 30, 30, "pt"),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(vjust = 5),
          plot.title  = element_text(vjust = 5))
}


