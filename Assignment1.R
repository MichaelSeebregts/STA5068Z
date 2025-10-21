### Q1 
legendre = function(order, x)
{
  
  L = 0
  
  for(i in 0:order)
  {
    L = L + (x^i)*choose(order, i)*choose(((order + i - 1)/2), order)
  }
  
  L = (2^order)*L
  
}

xSeq = seq(-1, 1, by = 1/100)

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

legenTarget = function(targets, Qf, x)
{
  
  plot(x, legendre(Qf, x), type = "l")
  
  for(j in 1:targets)
  {
    f = 0
    tempBeta = runif(Qf+1, -1, 1)
    for(i in 0:Qf)
    {
      f = f + tempBeta[i+1]*legendre(i, x)
    }
    
    lines(x, f, col = (j+2))
  }
  
  
  
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
  
  phi_x = 1/2
  bias2 = sum((gBar - underlying(xx_lat))[-Nx]^2*phi_x*dx)
  bias2
  
  
  testError = testError/M
  
  ones = matrix(1, M, 1)
  varAtX = colSums((G_D - ones%*%gBar)^2)/M # calculate variance at each point x
  var = sum(varAtX[-Nx]*phi_x*dx) # Usual riemann integral 
  var
  
  both = bias2 + var
  
  return(list("testError" = testError, "biasSquared" = bias2, "variance" = var))
}

biasVarH1 = bias_var(5, 1000, 1, 1/100, 0)
biasVarH2 = bias_var(5, 1000, 2, 1/100, 0)

varyingH1 = data.frame()

###



### Q4 
dat = read.table("Collider_Data_2022.txt", h = TRUE,stringsAsFactors =TRUE)

head(dat)

max.col(dat[, 3:5])

plot(dat$X1, dat$X2, col = max.col(dat[, 3:5]), type = "p", lwd = 4)

####

