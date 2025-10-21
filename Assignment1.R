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

underlying = function(x, a, b)
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
  
  return(list(y = y))
}

xBase = seq(-2, 2, length.out = 500)
runUnderlying = underlying(xBase, -2, 2)
plot(xBase, runUnderlying$y, type = "l")

x1 = runif(5, -2, 2)
x2 = runif(5, -2, 2)
x3 = runif(5, -2, 2)
x4 = runif(5, -2, 2)
x5 = runif(5, -2, 2)

run1 = underlying(x1, -2, 2)
run2 = underlying(x2, -2, 2)
run3 = underlying(x3, -2, 2)
run4 = underlying(x4, -2, 2)
run5 = underlying(x5, -2, 2)

hyp1 = function(x, xPred, y, color, plt)
{
  mod = lm(y ~ x)
  
  pred = predict(mod, newdata = data.frame(x = xPred))
  
  if (plt)
  {
    lines(xPred, pred, col = color) 
  }
  return(y)
}

hyp2 = function(x, xPred, y, color, plt)
{
  mod = lm(y ~ sin(pi*x) + cos(pi*x) - 1)
  pred = predict(mod, newdata = data.frame(x = xPred))
  if (plt)
  {
    lines(xPred, pred, col = color)
  }
  return(y)
}
plot(xBase, runUnderlying$y, type = "l")
hyp1(x1, xBase, run1$y, "red", TRUE)
points(x1, run1$y, col = "red", type = "p", lwd = 5)
hyp1(x2, xBase, run2$y, "blue", TRUE)
points(x2, run2$y, col = "blue", type = "p", lwd = 5)
hyp1(x3, xBase, run3$y, "green", TRUE)
points(x3, run3$y, col = "green", type = "p", lwd = 5)
hyp1(x4, xBase, run4$y, "pink", TRUE)
points(x4, run4$y, col = "pink", type = "p", lwd = 5)
hyp1(x5, xBase, run5$y, "magenta", TRUE)
points(x5, run5$y, col = "magenta", type = "p", lwd = 5)

plot(xBase, runUnderlying$y, type = "l")
hyp2(x1, xBase, run1$y, "red", TRUE)
points(x1, run1$y, col = "red", type = "p", lwd = 5)
hyp2(x2, xBase, run2$y, "blue", TRUE)
points(x2, run2$y, col = "blue", type = "p", lwd = 5)
hyp2(x3, xBase, run3$y, "green", TRUE)
points(x3, run3$y, col = "green", type = "p", lwd = 5)
hyp2(x4, xBase, run4$y, "pink", TRUE)
points(x4, run4$y, col = "pink", type = "p", lwd = 5)
hyp2(x5, xBase, run5$y, "magenta", TRUE)
points(x5, run5$y, col = "magenta", type = "p", lwd = 5)


bias_var = function(N, M, hyp, dx)
{
  xx_lat = seq(-2, +2, dx)
  Nx = length(xx_lat)
  
  gBar = xx_lat*0
  
  G_D = matrix(0, M, Nx)
  
  testError = 0
  for (i in 1:M)
  {
    resTemp = underlying(5, -2, 2, TRUE)
    
    if (hyp == 1)
    {
      resModelTemp = hyp1(xx_lat, resTemp$y, "", FALSE)
      g_D = resModelTemp
      
    }
    if (hyp == 2)
    { 
      g_D = hyp2(xx_lat, resTemp$y, "", FALSE)
    }
    print(g_D)
    G_D[i, ] = g_D
    gBar = gBar + g_D
    
    dat_oos = underlying(N, -2, 2, TRUE)
    pred_oos = predict(resModelTemp, data.frame(x = dat_oos$x))
    testErrorD = mean((dat_oos$y - pred_oos)^2)
    testError = testError + testErrorD
  }
  
  gBar = gBar/M
  
  phi_x = 1/2
  bias2 = (gBar - f(xx_lat))[-Nx]^2*phi_x*dx
  bias2
  
  
  testError = testError/M
  
  ones = matrix(1, M, 1)
  varAtX = colSums((G_D - ones%*%g_bar)^2)/M # calculate variance at each point x
  var = sum(varAtX[-Nx]*phi_x*dx) # Usual riemann integral 
  var
  
  both = bias2 + var
  
  return(list("testError" = testError, "biasSquared" = bias2, "variance" = var))
}

bias_var(5, 1000, 1, 1/100)

###



### Q4 
dat = read.table("Collider_Data_2022.txt", h = TRUE,stringsAsFactors =TRUE)

head(dat)

max.col(dat[, 3:5])

plot(dat$X1, dat$X2, col = max.col(dat[, 3:5]), type = "p", lwd = 4)

####

