require(sparseFLMM)

# (skew-)symmetric smooths ---------------------------------------

# generate random surface 
dat1 <- data.frame(arg1 = 1:50)
dat2 <- expand.grid(arg1 = 1:50, arg2 = 1:50)

Bskew <- Predict.matrix(
  smooth.construct( 
    s(arg1, arg2, bs = "symm", xt = list(skew = TRUE)),
    data = dat2, knots = NULL ),
  data = dat2 )
Bsymm <- Predict.matrix(
  smooth.construct( 
    s(arg1, arg2, bs = "symm", xt = list(skew = FALSE)),
    data = dat2, knots = NULL ),
  data = dat2 )

set.seed(934811)
dat2$yskew <- c(Bskew %*% rnorm(ncol(Bskew)))
dat2$ysymm <- c(Bsymm %*% rnorm(ncol(Bsymm)))

# fit sum of skew-symmetric and symmetric parts with corresponding smooths
modpa <- gam( I(yskew + ysymm) ~ s(arg1, arg2, bs = "symm", xt = list(skew = TRUE)) + 
                s(arg1, arg2, bs = "symm", xt = list(skew = FALSE)), data = dat2)
# predict surfaces
preds <- predict(modpa, type = "terms")
dat1 <- as.list(dat1)
dat1$arg2 <- dat1$arg1
dat1$predskew <- matrix(preds[,1], nrow = length(dat1$arg1))
dat1$predsymm <- matrix(preds[,2], nrow = length(dat1$arg1))

cols <- hcl.colors(12, "RdBu")
opar <- par(mfcol = c(2,2))
# symm part (intercept missing)
with(dat1, image(arg1, arg2, predsymm, asp = 1,
                 main = "Symmetric part of y",
                 col = cols))
with(dat1, image(arg1, arg2, asp = 1, 
                 main = "Fit via symm.smooth",
                 matrix(dat2$ysymm, nrow = length(arg1)), 
                                    col = cols))
# skew-symm part
with(dat1, image(arg1, arg2, predskew, asp = 1,
                 main = "Skew-symmetric part of y",
                 col = cols))
with(dat1, image(arg1, arg2, asp = 1, 
                 main = "Fit via symm.smooth",
                 matrix(dat2$yskew, nrow = length(arg1)), 
                 col = cols))
par(opar)

stopifnot(all.equal(dat1$predskew, - t(dat1$predskew)))
stopifnot(all.equal(dat1$predsymm, t(dat1$predsymm)))




# cyclic (skew-)symmetric splines ---------------------------------------

# fit the above example with cyclic smooths
modpac <- gam( I(yskew + ysymm) ~ s(arg1, arg2, bs = "symm", 
                                   xt = list(skew = TRUE, cyclic = TRUE)) + 
                s(arg1, arg2, bs = "symm", xt = list(skew = FALSE, cyclic = TRUE)),
              knots = list(arg1 = c(1, 50), arg2 = c(1,50)), 
              # specify arg range to specify 'wavelength'! 
              data = dat2)
plot(modpac, asp = 1, se = FALSE, pages = 1)

predsc <- predict(modpac, type = "terms")
dat1$predskewc <- matrix(predsc[,1], nrow = length(dat1$arg1))
dat1$predsymmc <- matrix(predsc[,2], nrow = length(dat1$arg1))

# check cyclic margins
opar <- par(mfrow = c(1,2))
with(dat1, matplot(arg1, predsymmc[, c(1,10, 40)], t = "l",
                   main = "symmetric smooth"))
abline(h = dat1$predsymmc[1, c(1,10, 40)], col = "darkgrey")
abline(v = c(1,50), col = "darkgrey")

with(dat1, matplot(arg1, predskewc[, c(1,10, 40)], t = "l",
                   main = "skew-symmetric smooth"))
abline(h = dat1$predskewc[1, c(1,10, 40)], col = "darkgrey")
abline(v = c(1,50), col = "darkgrey")
par(opar)



# 1D point symmetric B-splines --------------------------------------------

# generate toy data
dat <- data.frame( x = 1:100 )
ps_obj <- with(dat, s(x, bs = "ps"))
B <- Predict.matrix(smooth.construct(ps_obj, dat, NULL), dat)               
set.seed(3904)
dat$y <- B %*% rnorm(ncol(B))
plot(dat, t = "l")

# fit skew-symmetric spline
mod0 <- gam( y ~ s(x, bs = "symm", xt = list(skew = TRUE)), 
             knots = list(x = c(0,100)), # specify x range to determine inversion point 
             dat = dat )
lines(dat$x, predict(mod0), col = "cornflowerblue", lty = "dashed")

# or a symmetric spline to first part only
mod1 <- gam( y ~ s(x, bs = "symm"), 
             knots = list(x=c(0,50)), 
                          dat = dat[1:50, ])
lines(dat[1:50, ]$x, predict(mod1), col = "darkred", lty = "dashed")
