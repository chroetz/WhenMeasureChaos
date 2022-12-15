loaded <- readRDS("data5-1671101404.RDS")
sds <- loaded$sds
deltaTs <- loaded$deltaTs
msesLst <- loaded$msesLst
projDstLst <- loaded$projDstLst

projDstArray <- unlist(projDstLst)
dim(projDstArray) <- c(dim(projDstLst[[1]]), length(projDstLst))

arithMean <- sapply(msesLst, \(x) mean(x, na.rm = TRUE))
arithVar <- sapply(msesLst, \(x) var(x, na.rm = TRUE))

library(tidyverse)
data <- 
  tibble(
    sd = sds,
    deltaT = deltaTs,
    meanMSE = arithMean,
    varMSE = arithVar,
    relative3sigma = 3*sqrt(arithVar / 1e5) / arithMean,
    lower3sigma = arithMean - 3*sqrt(arithVar / 1e5),
    upper3sigma = arithMean + 3*sqrt(arithVar / 1e5),
    attrDst = apply(projDstArray[3,,], 2, mean, na.rm = TRUE),
    timesStateDst = apply(projDstArray[2,,], 2, mean, na.rm = TRUE),
    timeDst = apply(projDstArray[1,,], 2, mean, na.rm = TRUE)
  )

# 85 is corrupt...
data <- data[-85,]

data <- 
  data |> 
  mutate(sd = factor(sd)) |> 
  arrange(sd, deltaT)
pltMse <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(meanMSE = meanMSE/meanMSE[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = meanMSE, color = sd)) +
  labs(x="t", y="MSE", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0.6, 2.2, by = 0.2))
pltMsePerp <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(attrDst = attrDst/attrDst[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = attrDst, color = sd)) +
  labs(x="t", y="perpendicular MSE", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0.0, 2.2, by = 0.2))
pltAbsParaTime <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(timeDst = timeDst/timeDst[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = timeDst, color = sd)) +
  labs(x="t", y="parallel time distance", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) 
pltMseParaSpace <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(timesStateDst = timesStateDst/timesStateDst[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = timesStateDst, color = sd)) +
  labs(x="t", y="parallel MSE", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) 



