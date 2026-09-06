# Install packages (run once)
install.packages("forecast")
install.packages("tseries")
install.packages("ggplot2")
install.packages("TSA")

# Load libraries
library(forecast)
library(tseries)
library(ggplot2)
library(TSA)

# Create the Dataset

Year <- 2000:2015

Male <- c(2450,2502,2545,2623,2687,2761,2840,2918,
          2995,3067,3145,3230,3308,3389,3475,3555)

Female <- c(165,171,176,183,187,195,201,210,
            220,227,237,246,253,259,267,276)

Total <- Male + Female

data <- data.frame(Year,Male,Female,Total)

print(data)

# Descriptive Statistics

summary(data)

mean(Male)
mean(Female)
mean(Total)

sd(Male)
sd(Female)
sd(Total)

var(Male)
var(Female)
var(Total)

# Convert to Time Series

male.ts <- ts(Male,start=2000,frequency=1)

female.ts <- ts(Female,start=2000,frequency=1)

total.ts <- ts(Total,start=2000,frequency=1)

# Plot Time Series

plot(total.ts,
     main="Annual Total Inmate Population",
     xlab="Year",
     ylab="Population",
     col="blue",
     lwd=2)

plot(male.ts,
     main="Annual Male Inmate Population",
     xlab="Year",
     ylab="Male Population",
     col="darkgreen",
     lwd=2)

plot(female.ts,
     main="Annual Female Inmate Population",
     xlab="Year",
     ylab="Female Population",
     col="red",
     lwd=2)

# Trend Model

trend <- lm(Total~Year,data=data)

summary(trend)

plot(Year,Total,pch=19)

abline(trend,col="red",lwd=2)

# Augmented Dickey-Fuller Test

adf.test(total.ts)

# First Difference

diff.total <- diff(total.ts)

plot(diff.total,
     main="Differenced Series",
     col="blue")

adf.test(diff.total)


# ACF and PACF

acf(diff.total)

pacf(diff.total)


# Fit ARIMA Model Automatically

model <- auto.arima(total.ts)

summary(model)

# Manual ARIMA(1,1,1)


model2 <- arima(total.ts,order=c(1,1,1))

model2

# Diagnostic Checking

checkresiduals(model)

# Forecast Next Five Years

forecast.values <- forecast(model,h=5)

forecast.values

plot(forecast.values)


# Accuracy Measures

accuracy(model)

# Save Forecast

forecast.table <- data.frame(
  
  Year=2016:2020,
  
  Forecast=round(as.numeric(forecast.values$mean),0)
  
)

forecast.table

write.csv(forecast.table,
          "Forecast.csv",
          row.names=FALSE)
