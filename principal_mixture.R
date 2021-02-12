## "Wrong method" v2.0

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(splines)
library(lubridate)
library(furrr)
library(MatchIt)


## Process smoke wave/wildfire data
species <- local({
        load("data/smoke_wave_species.rda") ## 'm2'
        m2 <- rename(m2, smoke = SW111.x, date = date.x, fips = fips.x) %>% 
                mutate(fips = formatC(fips, flag = "0", width = 5))
        names(m2) <- sub("\\.PM2\\.5\\.LC$|\\.PM2\\.5\\.LC.TOR$", "", names(m2))
        cases <- select(m2, Aluminum:Zirconium, -Sulfur)
        m2 %>%
                mutate(include.spec = complete.cases(cases),
                       include.smoke = !is.na(smoke)) %>%
                as_tibble()
})

mb <- local({
        fipsList <- unique(species$fips)
        fileList <- sprintf("MCAPS/%s.rds", fipsList)
        use <- file.exists(fileList)
        fileList <- fileList[use]
        
        ## Health data
        m <- lapply(fileList, function(infile) {
                readRDS(infile) %>%
                        rename(cardio = ALL_CVD, resp = ALL_Resp) %>%
                        dplyr::select(date, denom, cardio, resp) %>%
                        filter(date >= "2004-01-01" & date <= "2009-12-31") %>%
                        group_by(date) %>% 
                        summarize_at(vars(denom:resp), sum) %>% 
                        mutate(dow = factor(weekdays(date)))
        })
        names(m) <- fipsList
        bind_rows(m, .id = "fips")
})

mcaps0 <- left_join(mb, select(species, -dow),
                    by = c("fips", "date")) %>%
        filter(include.smoke) %>%
        mutate(season = factor(quarter(date)),
               year.f = factor(year(date))) %>%
        select(-Sulfur)

## Keep Counties that have at least one smoke day
fips.use <- mcaps0 %>%
        group_by(fips) %>%
        summarize(n.smoke = sum(smoke)) %>%
        filter(n.smoke > 0)
mcaps1 <- mcaps0 %>%
        filter(fips %in% fips.use$fips)

dim(mcaps1)

## What's the difference between mcaps1 and counties w/speciation data?
spec.use <- mcaps1 %>%
        group_by(fips) %>%
        summarize(n.spec = sum(include.spec)) %>%
        mutate(has.spec = as.integer(n.spec > 0))

## Only counties with any speciation data
mcaps2 <- mcaps1 %>%
        inner_join(filter(spec.use, has.spec > 0) %>%
                           select(-n.spec, has.spec),
                   by = "fips")
dim(mcaps2)

mcaps3 <- mcaps2 %>%
        mutate(has.spec_day = complete.cases(select(mcaps2, Aluminum:Zirconium)))


### Use MatchIt package
mcaps.m00 <- local({
        modeldata <- mcaps3 %>%
                select(smoke, fips, date, tmpd, dptp, has.spec_day, total) %>%
                mutate(fips = factor(fips),
                       season = factor(quarter(date)),
                       year = factor(year(date)),
                       month = factor(month(date)),
                       season = factor(quarter(date)))
        use <- complete.cases(modeldata)
        modeldata <- modeldata[use, ]
        
        ## Matching model
        pfit <- matchit(smoke ~ month + year + fips + tmpd + dptp, 
                        data = modeldata, ratio = 4, method = "nearest", 
                        replace = TRUE)
        ## summary(pfit)
        mdata <- match.data(pfit) %>% 
                as_tibble() %>%
                mutate(fips = as.character(fips))
        inner_join(mcaps3, 
                   select(mdata, fips, date, weights), 
                   by = c("fips", "date")) %>%
                mutate(logtotal = ifelse(total > 0, log(total), NA),
                       has.spec_day = as.integer(has.spec_day))
})
        


## Check matching
mcaps.m00 %>%
        #mutate(weights = 1) %>%
        group_by(smoke) %>%
        summarize(tmpd = weighted.mean(tmpd, weights, na.rm=T), 
                  dptp = weighted.mean(dptp, weights, na.rm=T), 
                  fips = weighted.mean(fips == "06037", weights),
                  dt = weighted.mean(date, weights),
                  spec = weighted.mean(has.spec_day, weights),
                  quarter = weighted.mean(quarter(date), weights),
                  total = weighted.mean(total, weights)
        )

mcaps.m <- mcaps.m00 %>%
        filter(has.spec_day > 0)


mcaps.m %>%
        mutate(smoke.f = factor(smoke, labels = c("Non-smoke wave", "smoke wave"))) %>%
        ggplot(aes(smoke.f, log(total))) + 
        xlab("") +
        ylab(expression(log * " " * PM[2.5] * " concentration")) +
        geom_boxplot() + 
        theme_minimal()

################################################################################
## PM Model

sd_x <- 0.75
mcaps.m <- local({
        omega <- 0.2
        wfun <- function(x) {
                dnorm(x, sd = sd_x)
        }
        fit0 <- lm(logtotal ~ fips + dow + ns(tmpd, 3), 
                   data = filter(mcaps.m, smoke == 0))
        mu0 <- predict(fit0)
        eta0 <- summary(fit0)$sigma
        mu01 <- predict(fit0, filter(mcaps.m, smoke == 1))
        pm1 <- filter(mcaps.m, smoke == 1)$logtotal
        eta1 <- summary(lm(pm1 ~ offset(mu01)))$sigma
        pm0 <- filter(mcaps.m, smoke == 0)$logtotal
        tau <- coef(lm(logtotal ~ smoke, mcaps.m))[2]
        pm01 <- mu01 + omega * eta0/eta1 * (pm1 - mu01 - tau)
        #pm01 <- rnorm(length(mu01), mu01 + omega * eta0/eta1 * (pm1 - mu01 - tau),                      eta0^2 * (1 - omega))
        pm10 <- (mu0 + tau) + omega * eta1/eta0 * (pm0 - mu0)
        #pm10 <- rnorm(length(mu0), (mu0 + tau) + omega * eta1/eta0 * (pm0 - mu0), eta1^2 * (1 - omega))
        mcaps.m <- mutate(mcaps.m, pm_weight = 1)
        mcaps.m$pm_weight[mcaps.m$smoke == 0] <- wfun(pm0 - pm10)
        mcaps.m$pm_weight[mcaps.m$smoke == 1] <- wfun(pm1 - pm01)
        mcaps.m <- mcaps.m %>%
                mutate(ppm0 = NA, ppm10 = NA, ppm1 = NA, ppm01 = NA) %>%
                mutate(ppm0 = replace(ppm0, smoke == 0, pm0),
                       ppm10 = replace(ppm10, smoke == 0, pm10),
                       ppm1 = replace(ppm1, smoke == 1, pm1),
                       ppm01 = replace(ppm01, smoke == 1, pm01))
        mcaps.m
})

mcaps.m %>%
        select(ppm0, ppm10) %>%
        rename(pm0 = ppm0, pm1 = ppm10) %>%
        bind_rows(select(mcaps.m, ppm01, ppm1, smoke) %>%
                          rename(pm0 = ppm01, 
                                 pm1 = ppm1),
                  .id = "observed") %>%
        mutate(observed = factor(observed, labels = c("0", "1"))) %>%
        ggplot(aes(pm0, pm1)) + 
        #geom_point(aes(color = observed)) + 
        geom_point(alpha = 1/2) +
        geom_rug(alpha = 1/8, length = unit(0.015, "npc")) +
        geom_abline(intercept = 0, slope = 1) + 
        geom_abline(intercept = sd_x, slope = 1, lty = 2) + 
        geom_abline(intercept = -sd_x, slope = 1, lty = 2) +
        xlab(expression(PM[t](0))) + 
        ylab(expression(PM[t](1))) +
        theme_minimal()

mcaps.m %>%
        mutate(obs0 = ppm10 - ppm0,
               obs1 = ppm1 - ppm01) %>%
        select(obs0, obs1) %>%
        gather(observed, diff) %>%
        mutate(observed = factor(observed)) %>%
        #ggplot(aes(diff, fill = observed)) + 
        ggplot(aes(diff)) + 
        geom_histogram(bins = 10) + 
        geom_vline(xintercept = 0) + 
        geom_vline(xintercept = sd_x, lty = 2) + 
        geom_vline(xintercept = -sd_x, lty = 2) +
        xlab(expression(X[t](1)-X[t](0))) +
        theme_minimal()

## SIR for binary outcome
vnorm <- function(x) {
        x / sqrt(sum(x * x))
}
sir_b <- function(x, z, w = 1) {
        x <- data.matrix(x)
        sw <- sqrt(w)
        x <- sweep(x, 1, sw, "*")
        xd <- sweep(x, 2, colMeans(x), "-")
        R <- cov(xd) %>% 
                chol %>% 
                solve
        mg <- xd %*% R %>%
                data.frame(storm = z) %>%
                group_by(storm) %>%
                summarize_all(mean) %>%
                select(-storm)
        cf <- (mg[2, ] - mg[1, ]) %>% 
                unlist
        drop(vnorm(R %*% cf))
}


x <- mcaps.m %>%
        select(Aluminum:Zirconium) %>% 
        data.matrix
xd <- sweep(x, 2, colMeans(x), "-")
z <- mcaps.m$smoke
wts <- mcaps.m$pm_weight
b <- sir_b(x, z, wts)

## Create principal mixture score
mcaps.m <- mcaps.m %>%
        mutate(zscore = drop(xd %*% b),
               smoke.f = factor(smoke, labels = c("None", "Smoke")))
xstd <- apply(xd, 2, sd)

dim(mcaps.m)

tibble(b = b, elt = names(b)) %>%
        ggplot(aes(elt, b)) + 
        geom_linerange(aes(ymin = 0, ymax = b), lwd = 1) + 
        theme_minimal() + 
        xlab("") + 
        ylab("Coefficient") +
        theme(axis.text.x = element_text(angle = 90))

mcaps.m %>%
        group_by(smoke) %>%
        summarize(zsd = sd(zscore))

lm(zscore ~ smoke, mcaps.m) %>%
        tidy()

###############################################################################
## Modeling

mcaps.m0 <- mcaps.m

get_pred <- function(mcaps.m, ydf = 4) {
        dat <- with(mcaps.m, {
                data.frame(z = smoke,
                           x = zscore,
                           y = resp,
                           fips = fips,
                           dow = dow,
                           tmpd = tmpd,
                           dptp = dptp,
                           date = date,
                           weights = weights,
                           pm_weight = pm_weight,
                           year = factor(year(date)),
                           denom = denom
                )
        }) %>% as_tibble()
        try({
                rho <- 0.8
                fit0 <- lm(x ~ ns(date, 3), 
                           data = filter(dat, z == 0), 
                           weights = weights)
                sigma0 <- summary(fit0)$sigma
                mu1 <- predict(fit0, filter(dat, z == 1))
                x1 <- filter(dat, z == 1)$x
                sigma1 <- summary(lm(x1 ~ offset(mu1)))$sigma
                delta <- lm(x1 ~ offset(mu1), 
                            data = filter(dat, z == 1)) %>%
                        coef()
                x01 <- mu1 + rho * sigma0/sigma1 * (x1 - mu1 - delta)
                diffx <- x1 - x01
                
                fit0 <- glm(y ~ offset(log(denom)) + dow + tmpd + dptp,
                            data = filter(dat, z == 0),
                            family = stats::poisson, 
                            weights = weights,
                            control = glm.control(maxit = 100)) 
                v1 <- predict(fit0, filter(dat, z == 1), type = "link")
                fit <- glm(y ~ offset(v1) 
                           -1 + ns(diffx, ydf, intercept = TRUE),
                           data = filter(dat, z == 1), 
                           weights = weights,
                           family = stats::poisson,
                           control = glm.control(maxit = 100))
                pred <- predict(fit, type = "terms")
                dup <- duplicated(diffx)
                f <- splinefun(diffx[!dup], pred[!dup, 1])
                f(xpts)
        })
}

xpts <- seq(-0.005, 0.008, len = 100)

plan(multisession, workers = 4L)
#plan(sequential)
r <- future_map(seq_len(5000), function(dummy) {
        mcaps.m <- mcaps.m0 %>%
                group_by(smoke) %>%
                sample_n(n(), replace = TRUE) %>%
                ungroup()
        get_pred(mcaps.m, 2 + 1)
}, .progress = TRUE)

pdata <- local({
        u <- !sapply(r, inherits, what = "try-error")
        rr <- do.call("cbind", r[u])
        f <- splinefun(xpts, rowMeans(rr))
        tibble(xpts = xpts,
               fit = rowMeans(rr, na.rm = TRUE),
               lo = apply(rr, 1, quantile, 0.025, na.rm = TRUE),
               hi = apply(rr, 1, quantile, 0.975, na.rm = TRUE)) 
})
print(
        pdata %>%
                ggplot(aes(xpts, exp(fit))) + 
                geom_line() + 
                geom_line(aes(xpts, exp(lo)), lty = 2) + 
                geom_line(aes(xpts, exp(hi)), lty = 2) + 
                geom_hline(yintercept = 1, lty = 3) +
                geom_vline(xintercept = 0, lty = 3) + 
                xlab(expression(X[t](1)-X[t](0))) + 
                ylab(expression("E[" * Y[t](1) * "]/E[" * Y[t](0) * "]")) + 
                theme_minimal()
)





