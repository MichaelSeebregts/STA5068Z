### Q1 
legendre = function(order, x, sigma = 0)
{
  
  L = 0
  
  for(i in 0:order)
  {
    L = L + (x^i)*choose(order, i)*choose(((order + i - 1)/2), order)
  }
  
  L = (2^order)*L
  
}

xSeq = seq(-1, 1, by = 1/500)

l0 = legendre(0, xSeq)
l1 = legendre(1, xSeq)
l2 = legendre(2, xSeq)
l3 = legendre(3, xSeq)
l4 = legendre(4, xSeq)
l5 = legendre(5, xSeq)

plot(xSeq, l0, type = "l", ylim = c(-1, 1))
lines(xSeq, l1, col = "red")
lines(xSeq, l2, col = "blue")
lines(xSeq, l3, col = "green")
lines(xSeq, l4, col = "pink")
lines(xSeq, l5, col = "lightblue")

set.seed(123)

legenTarget = function(targets, Qf, x, plot = TRUE, sigma = 0, beta = c())
{
  
  if(plot == TRUE)
  {
    plot(x, legendre(Qf, x), type = "l")
  }
  
  for(j in 1:targets)
  {
    f = 0
    if (length(beta) == 0)
    {
      tempBeta = runif(Qf+1, -1, 1)
    }
    else
    {
      tempBeta = beta
    }
    
    for(i in 0:Qf)
    {
      f = f + tempBeta[i+1]*legendre(i, x)
    }
    f = f + rnorm(length(x), 0, sigma)
    
    if (plot == TRUE)
    {
      lines(x, f, col = (j+2))
    }
    
  }
  
  return(invisible(f))
  
  
  
}

polyTarget = function(targets, Qf, x)
{
  plot(x, legendre(Qf, x), type = "l")
  
  for(j in 1:targets)
  {
    
    f = 0
    tempAlpha = runif(Qf + 1, -1, 1)
    for (i in 0:Qf)
    {
      f = f + tempAlpha[i+1]*(x^i)
    }
    
    lines(x, f, col = (j+2))
  }
}

polyTarget(3, 1, xSeq)
legenTarget(3, 1, xSeq)

polyTarget(3, 2, xSeq)
legenTarget(3, 2, xSeq)

polyTarget(3, 3, xSeq)
legenTarget(3, 3, xSeq)

polyTarget(3, 4, xSeq)
legenTarget(3, 4, xSeq)

polyTarget(3, 5, xSeq)
legenTarget(3, 5, xSeq)

###


### Q2

# a

pal  <- colorRampPalette(c("navy","cyan","yellow","red"))(100)
sc   <- function(z, zmin, zmax) pmax(1, pmin(100, round(1 + 99*(z - zmin)/(zmax - zmin))))

target = legendre(10, xSeq, 0)
plot(target ~ xSeq, type = "l")

NSeq = c(20:110)
sigmaSeq = seq(0.2, 1.1, length.out = length(NSeq))

overFitDF = data.frame()

betaTemp = runif(11, -1, 1)

#plot(1,1, type = "n", xlim = c(20, 110), ylim = c(0.2, 1.1))

for (i in 1:length(NSeq))
{
  #overFitRow = c()
  
  for (j in 1:length(sigmaSeq))
  {
    xTemp = runif(NSeq[i], -1, 1)
    yTemp = legenTarget(1, 10, xTemp, FALSE, sigmaSeq[j], betaTemp)
    dat = data.frame(y = yTemp, x = xTemp)
    H2  = lm(y ~ x + I(x^2), data = dat)
    H10 = lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) +
               I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10), data = dat)
    xOut = runif(500, -1, 1)
    y2Pred = predict(H2, newdata = data.frame(x = xOut))
    y10Pred = predict(H10, newdata = data.frame(x = xOut))
    
    yOut = legenTarget(1, 10, xOut, FALSE, 0, betaTemp)
    
    H2Error = mean((yOut - y2Pred)^2)
    H10Error = mean((yOut - y10Pred)^2)
    
    overFitTemp = H10Error - H2Error
    overFitDF[i, j] = overFitTemp
    #points(NSeq[i], sigmaSeq[j], col = pal[ sc(overFitTemp, zmin, zmax) ])
  }
}

zmin <- min(overFitDF, na.rm = TRUE)
zmax <- max(overFitDF, na.rm = TRUE)

plot(1, 1, type = "n", xlim = range(NSeq), ylim = range(sigmaSeq),
     xlab = "N", ylab = "sigma")
for (i in seq_along(NSeq)) {
  for (j in seq_along(sigmaSeq)) {
    points(NSeq[i], sigmaSeq[j],
           col = pal[ sc(overFitDF[i, j], zmin, zmax) ],
           pch = 16, cex = 0.6)
  }
}

# pal <- colorRampPalette(c("navy","cyan","yellow","red"))
# 
# filled.contour(x = NSeq,
#                y = sigmaSeq,
#                z = t(overFitDF),
#                color.palette = pal,
#                xlab = "N", ylab = "sigma",
#                plot.title = title(main = "Overfit heatmap"),
#                nlevels = 50)

# b

Nseq2 = 20:60
Order = 1:40
sig = 0.2

overFitDF2 = data.frame()

#plot(1,1, type = "n", xlim = c(20, 110), ylim = c(0.2, 1.1))

for (i in 1:length(Nseq2))
{
  #overFitRow = c()
  
  for (j in 1:length(Order))
  {
    xTemp = runif(NSeq[i], -1, 1)
    yTemp = legenTarget(1, Order[j], xTemp, FALSE, sig)
    dat = data.frame(y = yTemp, x = xTemp)
    H2 = lm(y ~ poly(x, 2), data = dat)
    H10 = lm(y ~ poly(x, 10), data = dat)
    xOut = runif(500, -1, 1)
    y2Pred = predict(H2, newdata = data.frame(x = xOut))
    y10Pred = predict(H10, newdata = data.frame(x = xOut))
    
    yOut = legenTarget(1, Order[j], xOut, FALSE, 0)
    
    H2Error = mean((yOut - y2Pred)^2)
    H10Error = mean((yOut - y10Pred)^2)
    
    overFitTemp = H10Error - H2Error
    overFitDF2[i, j] = overFitTemp
    #points(NSeq[i], sigmaSeq[j], col = pal[ sc(overFitTemp, zmin, zmax) ])
  }
}

zmin <- min(overFitDF2, na.rm = TRUE)
zmax <- max(overFitDF2, na.rm = TRUE)

plot(1, 1, type = "n", xlim = range(Nseq2), ylim = range(Order),
     xlab = "N", ylab = "sigma")
for (i in seq_along(Nseq2)) {
  for (j in seq_along(Order)) {
    points(Nseq2[i], Order[j],
           col = pal[ sc(overFitDF2[i, j], zmin, zmax) ],
           pch = 16, cex = 0.6)
  }
}

filled.contour(x = Nseq2,
               y = Order,
               z = t(overFitDF2),
               color.palette = pal,
               xlab = "N", ylab = "Order",
               plot.title = title(main = "Overfit heatmap"),
               nlevels = 50)

###


### Q3

set.seed(2023)

underlying = function(x)
{
  
  y = rep(0, length(x))
  for (j in 1:length(x))
  {
    if (x[j] < 0)
    {
      y[j] = abs(x[j] + 1) - 0.5
    }
    
    if (x[j] >= 0)
    {
      y[j] = abs(x[j] - 1) - 0.5
    }
  }
  
  return( y)
}

DGP = function(N, sigma)
{
  x = runif(N, -2, 2)
  e = rnorm(N, 0, sigma)
  y = underlying(x) + e
  return(list(x = x, y = y))
}

xBase = seq(-2, 2, length.out = 500)
runUnderlying = underlying(xBase)
plot(xBase, runUnderlying, type = "l")

hyp1 = function(x, xPred, y, color, plt)
{
  mod = lm(y ~ x)
  
  pred = predict(mod, newdata = data.frame(x = xPred))
  
  if (plt)
  {
    lines(xPred, pred, col = color) 
  }
  return(pred)
}

hyp2 = function(x, xPred, y, color, plt)
{
  mod = lm(y ~ sin(pi*x) + cos(pi*x) - 1)
  pred = predict(mod, newdata = data.frame(x = xPred))
  if (plt)
  {
    lines(xPred, pred, col = color)
  }
  return(pred)
}

dgp1 = DGP(5, 0)
dgp2 = DGP(5, 0)
dgp3 = DGP(5, 0)
dgp4 = DGP(5, 0)
dgp5 = DGP(5, 0)

plot(xBase, runUnderlying, type = "l")
hyp1(dgp1$x, xBase, dgp1$y, "red", TRUE)
points(dgp1$x, dgp1$y, col = "red", type = "p", lwd = 5)
hyp1(dgp2$x, xBase, dgp2$y, "blue", TRUE)
points(dgp2$x, dgp2$y, col = "blue", type = "p", lwd = 5)
hyp1(dgp3$x, xBase, dgp3$y, "green", TRUE)
points(dgp3$x, dgp3$y, col = "green", type = "p", lwd = 5)
hyp1(dgp4$x, xBase, dgp4$y, "pink", TRUE)
points(dgp4$x, dgp4$y, col = "pink", type = "p", lwd = 5)
hyp1(dgp5$x, xBase, dgp5$y, "magenta", TRUE)
points(dgp5$x, dgp5$y, col = "magenta", type = "p", lwd = 5)

plot(xBase, runUnderlying, type = "l")
hyp2(dgp1$x, xBase, dgp1$y, "red", TRUE)
points(dgp1$x, dgp1$y, col = "red", type = "p", lwd = 5)
hyp2(dgp2$x, xBase, dgp2$y, "blue", TRUE)
points(dgp2$x, dgp2$y, col = "blue", type = "p", lwd = 5)
hyp2(dgp3$x, xBase, dgp3$y, "green", TRUE)
points(dgp3$x, dgp3$y, col = "green", type = "p", lwd = 5)
hyp2(dgp4$x, xBase, dgp4$y, "pink", TRUE)
points(dgp4$x, dgp4$y, col = "pink", type = "p", lwd = 5)
hyp2(dgp5$x, xBase, dgp5$y, "magenta", TRUE)
points(dgp5$x, dgp5$y, col = "magenta", type = "p", lwd = 5)


bias_var = function(N, M, hyp, dx, sig)
{
  xx_lat = seq(-2, +2, dx)
  Nx = length(xx_lat)
  
  gBar = xx_lat*0
  
  G_D = matrix(0, M, Nx)
  
  testError = 0
  for (i in 1:M)
  {
    dat = DGP(N, sig)
    x = dat$x
    y = dat$y
    
    if (hyp == 1)
    {
      mod = lm(y ~ x)
      
    }
    if (hyp == 2)
    { 
      mod = lm(y ~ sin(pi*x) + cos(pi*x) - 1)
    }
    g_D = predict(mod, newdata = data.frame(x = xx_lat))
    
    G_D[i, ] = g_D
    gBar = gBar + g_D
    
    datOOS = DGP(N, sig)
    pred_oos = predict(mod, data.frame(x = datOOS$x))
    testErrorD = mean((datOOS$y - pred_oos)^2)
    testError = testError + testErrorD
  }
  
  gBar = gBar/M
  
  phi_x = 1/4
  bias2 = sum((gBar - underlying(xx_lat))[-Nx]^2*phi_x*dx)
  bias2
  
  
  testError = testError/M
  
  ones = matrix(1, M, 1)
  varAtX = colSums((G_D - ones%*%gBar)^2)/M # calculate variance at each point x
  var = sum(varAtX[-Nx]*phi_x*dx) # Usual riemann integral 
  var
  
  both = bias2 + var
  
  return(list("testError" = testError, "bias" = bias2, "variance" = var, "OOSE" = both))
}

biasVarH1 = bias_var(5, 1000, 1, 1/100, 0)
biasVarH2 = bias_var(5, 1000, 2, 1/100, 0)

varyingH1 = data.frame("Bias" = numeric(0), "Var" = numeric(0), "OOSE" = numeric(0))
varyingH2 = data.frame("Bias" = numeric(0), "Var" = numeric(0), "OOSE" = numeric(0))

sigSeq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

for (j in 1:length(sigSeq))
{
  tempH1 = bias_var(5, 1000, 1, 1/100, sigSeq[j])
  tempH2 = bias_var(5, 1000, 2, 1/100, sigSeq[j])
  
  tempH1Res = c(tempH1$bias, tempH1$variance, tempH1$testError)
  tempH2Res = c(tempH2$bias, tempH2$variance, tempH2$testError)
  
  varyingH1 = rbind(varyingH1, tempH1Res)
  varyingH2 = rbind(varyingH2, tempH2Res)
}

colnames(varyingH1) = c("Bias", "Variance", "OOSE")
colnames(varyingH2) = c("Bias", "Variance", "OOSE")

varyingH1 = cbind(sigSeq, varyingH1)
varyingH2 = cbind(sigSeq, varyingH2)

plot(varyingH1$sigSeq, varyingH1$Bias, type = "l", ylim = c(min(varyingH1$Bias, varyingH2$Bias), max(varyingH1$Bias, varyingH2$Bias)), 
     xlab = "Sigma", ylab = "Bias")
lines(varyingH2$sigSeq, varyingH2$Bias, col = "red")

plot(varyingH1$sigSeq, varyingH1$Variance, type = "l", ylim = c(min(varyingH1$Variance, varyingH2$Variance), max(varyingH1$Variance, varyingH2$Variance)), 
     xlab = "Sigma", ylab = "Variance")
lines(varyingH2$sigSeq, varyingH2$Variance, col = "red")

plot(varyingH1$sigSeq, varyingH1$OOSE, type = "l", ylim = c(min(varyingH1$OOSE, varyingH2$OOSE), max(varyingH1$OOSE, varyingH2$OOSE)), 
     xlab = "Sigma", ylab = "Out of Sample Error")
lines(varyingH2$sigSeq, varyingH2$OOSE, col = "red")

###



### Q4 

library(colorspace)

color.gradient = function(x, colors=c('cyan','orange','magenta'), colsteps=50)
{
  colpal = colorRampPalette(colors)
  return( colpal(colsteps)[ findInterval(x, seq(min(x),max(x), length=colsteps)) ] )
}

dat = read.table("Collider_Data_2022.txt", h = TRUE,stringsAsFactors =TRUE)

head(dat)

plot(dat$X1, dat$X2, col = max.col(dat[, 3:5]), type = "p", lwd = 4)

X = cbind(dat[, 1:2])
Y = dat[, 3:5]

softMax = function(x)
{
  denom = sum(exp(x))
  softVec = rep(0, length.out = length(x))
  for(i in 1:length(x))
  {
    softVec[i] = exp(x[i])/denom
  }
  return(softVec)
}

softMaxMatrix = function(X)
{
  softMat = apply(X, 2, softMax)
  return(softMat)
}


neural_net = function(X,Y,theta,m,nu)
{
  # Infer dimensions:
  N = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  
  d = c(p,m,q) 
  
  # Populate weight and bias matrices:
  index = 1:(d[1]*d[2])
  W1    = matrix(theta[index],d[1],d[2])
  index = max(index)+(1:(d[2]*d[3]))
  W2    = matrix(theta[index],d[2],d[3])
  index = max(index)+(1:d[2])
  b1    = matrix(theta[index],d[2],1)
  index =max(index)+(1:d[3])
  b2    = matrix(theta[index],d[3],1)
  
  ones = matrix(1,N,1)
  # Evaluate the updating equation in matrix form
  #for(i in 1:N)
  #{
  A0 = t(X)
  A1 = tanh(t(W1)%*%A0+b1%*%t(ones))
  
  #print(head(t(W2)%*%A1+b2%*%t(ones)))
  A2 = softMaxMatrix(t(W2)%*%A1+b2%*%t(ones))
  #print(head(A2))
  #}
  yHat = A2
  # print(dim(yHat))
  # print(A2)
  
  error = rep(NA,N)    #-(Y*log(pi_hat)+(1-Y)*log(1-pi_hat))
  logA2 = log(A2)
  crossEnt = t(Y)*logA2
  #print(head(crossEnt))
  error = apply(crossEnt, 2, sum)
  # Evaluate an appropriate objective function and return some predictions:
  E1 = (-1)*sum(error)/N
  E2 = E1+nu/N*(sum(W1^2)+sum(W2^2))
  # Return a list of relevant objects:
  return(list(A2 = A2,A1 = A1, E1 = E1, E2 = E2))
}
m   = 10
p     = dim(X)[2]
q     = dim(Y)[2]
npars = p*m+m*q+m+q
theta_rand = runif(npars,-1,+1)
res = neural_net(X,Y,theta_rand,m,0)

set.seed(2022)

M  = 100
x1 = seq(min(dat$X1),max(dat$X1),length = M)
x2 = seq(min(dat$X2), max(dat$X2),length = M)
xx1 = rep(x1,M) 
xx2 = rep(x2, each = M)
points(xx2~xx1, pch = 16,cex = 0.75)
abline(v = x1,h = x2,lty =2)

XX = cbind(xx1,xx2)
YY = matrix(1,M^2,3)


N      = dim(X)[1]
set    = sample(1:N,0.5*N,replace = FALSE)
Xtrain = X[set,,drop = FALSE] 
Ytrain = Y[set,,drop = FALSE] 
Xval   = X[-set,,drop = FALSE]
Yval   = Y[-set,,drop = FALSE]


# nu = 5
obj_pen = function(pars)
{
  res = neural_net(Xtrain,Ytrain,pars,m,nu)
  return(res$E2)
}
# obj_pen(theta_rand)


n_nus     = 10
nu_seq    = exp(seq(-5,0,length = n_nus))
val_error = rep(NA,n_nus)
for(i in 1:n_nus)
{
  nu  = nu_seq[i]
  # Fit the neural network using a standard optimizer in R:
  res_opt = nlm(obj_pen,theta_rand, iterlim = 1000)
  
  res_val      = neural_net(Xval,Yval,res_opt$estimate,m,0)
  val_error[i] = res_val$E1
  
  # Draw a response curve over the 2D input space to see
  # what pattern the neural network predicts
  
  # res_fitted = neural_net(XX,YY,res_opt$estimate,m,0)
  # 
  # plot(XX[,2]~XX[,1], pch = 16, col = color.gradient(max.col(t(res_fitted$A2))))
  # points(X$X2~X$X1, col = max.col(Y), pch = 16)
  # 
  # print(paste0('Validation run ', i,'| nu = ',round(nu,5)))
}

plot(val_error~nu_seq,type = 'b',pch = 16,cex = 0.75)

wh = which.min(val_error)

nu_final = nu_seq[wh]
nu_final_conservative = nu_seq[wh+1]
abline(v = nu_final)
abline(v = nu_final_conservative)
obj_final = function(pars)
{
  res = neural_net(X,Y,pars,m,nu_final_conservative)
  return(res$E2)
}
res_opt_final = nlm(obj_final,theta_rand, iterlim = 1000)


res_fitted = neural_net(XX,YY,res_opt_final$estimate,m,0)

plot(XX[,2]~XX[,1], pch = 16, col = color.gradient(max.col(t(res_fitted$A2))))
points(X$X2~X$X1, col = max.col(Y), pch = 16)


####

