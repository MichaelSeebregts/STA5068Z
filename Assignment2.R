#### Q1

options(scipen = 999)

rm(list = ls())

draw_objects = function(J)
{  
  r_o       =  runif(J,0.3,1)
  pi_vals   =  runif(J,-1,1)*pi
  o1        =  r_o*cos(pi_vals)     # x-coordinates of objects
  o2        =  r_o*sin(pi_vals)
  return(cbind(o1,o2))
}
draw_starts = function(N = 1)
{
  r_x      =  runif(N,0,0.25)
  pi_x     =  runif(N,-1,1)*pi
  xt       =  matrix(cbind(r_x*cos(pi_x),r_x*sin(pi_x)),nrow = N,byrow = TRUE)
}

set.seed(2024)
J       =  50  
objects = draw_objects(J)
xt      = draw_starts()
delt    = 0.025             # How big are the steps you can take?
rad     = 0.05              # How close to an object before you crash?

# Just a function to plot the game:
plot_game = function(xt,objects,rad,cols = 1)
{
  plot(objects[,2]~objects[,1],type = 'n', ylim = c(-1,1), xlim = c(-1,1))
  symbols(objects[,1],objects[,2],circles = rep(rad,J),ylim = c(-1,1),xlim = c(-1,1),inches = FALSE,add = TRUE,bg = 'seagreen')
  points(xt[,2]~xt[,1], pch = 16, cex = 1,col = cols)
  pi_v = seq(-1,1,1/100)*pi
  y_edge = sin(pi_v)
  x_edge = cos(pi_v)
  lines(y_edge~x_edge)
  lines(c(0.25*y_edge)~c(0.25*x_edge),lty = 2)
}

plot_game(xt,objects,rad,'black')


game_status=function(xt,objects, rad)
{
  sk    = rep(0,dim(xt)[1])
  min_dists = rep(0,dim(xt)[1])
  J         = dim(objects)[1]
  ones      = matrix(1,J,1)
  for(i in 1:dim(xt)[1])
  {
    min_dists[i] = min(sqrt(rowSums((objects-ones%*%xt[i,])^2)))
    
    if (((xt[i, 1])^2 + (xt[i, 2]^2 ) >= 1))
    {
      sk[i] = 1    
    }
    else if (min_dists[i] < rad)
    {
      sk[i] = -1
    }
    else 
    {
      sk[i] = 0
    }
  }
  
  #sk = 1*(xt[,2]>1)-1*((xt[,1]< -1)|(xt[,1]> +1)|(xt[,2]< -1)|(min_dists<rad))
  ret = list(status = sk, minDist = min_dists)
  return(ret)
}


play = function(x0,delt,objects,rad,theta,plt = FALSE,trace = FALSE)
{
  k           =   0 # Count how many steps
  xt          = x0   # Set the initial coordinate(s) for the drone.
  trajectories = NULL
  # Check the game status:
  res_status  = game_status(xt, objects, rad)
  status      = res_status$status
    # Check which games are still active:
  terminal      = (status != 0)
  if(trace)
  {
    trajectories = array(dim = c(dim(xt)[1],dim(xt)[2],101))
    trajectories[,,1] = xt
  }
  if(plt){plot_game(xt,objects,rad,'black')}
  while((any(status==0))&(k<100))
  {
    k = k + 1
      
    # Now, let the model update the position of the pieces:
      
    ct = control(xt, theta)
    xt = xt+ct*delt*cbind(1-terminal,1-terminal)
    
    # Checkk the game status after the positions are updates:
    res_status  = game_status(xt, objects, rad)
    status      = res_status$status
    terminal    = (status != 0)
    if(trace){trajectories[,,k] = xt}
    if(plt){plot_game(xt,objects,rad,c('red','black','green')[status+2])}
  }
  return(list(k = k, status = status,xt= xt,trajectories = trajectories))
}


model = function(X,theta,nodes)
{
  # Infer dimensions:
  N = dim(X)[1]
  p = dim(X)[2]
  q = 2
  dims = c(p,nodes,q)
  
  # Populate weight and bias matrices:
  index = 1:(dims[1]*dims[2])
  W1    = matrix(theta[index],dims[1],dims[2])
  index = max(index)+1:(dims[2]*dims[3])
  W2    = matrix(theta[index],dims[2],dims[3])
  index = max(index)+1:(dims[3]*dims[4])
  W3    = matrix(theta[index],dims[3],dims[4])
  
  index = max(index)+1:(dims[2])
  b1    = matrix(theta[index],dims[2],1)
  index = max(index)+1:(dims[3])
  b2    = matrix(theta[index],dims[3],1)
  index = max(index)+1:(dims[4])
  b3    = matrix(theta[index],dims[4],1)
  
  ones    = matrix(1, 1, N)
  a0      = t(X)
    
    # Evaluate the updating equation in matrix form
  a1 = tanh(t(W1)%*%a0 + b1%*%ones)
  a2 = tanh(t(W2)%*%a1 + b2%*%ones)
  a3 = tanh(t(W3)%*%a2 + b3%*%ones)
  
  # Return a list of relevant objects:
  return(list(a3 = t(a3)))
}

p     = 2
q     = 2
nodes = 3
npars = p*nodes+nodes*nodes+nodes*q+nodes+nodes+q
npars
theta_rand = runif(npars,-1,1)

control = function(xt,pars)
{
  
  res_model = model(xt,pars,rep(nodes,2))
  return(res_model$a3)
}
control(xt,theta_rand)

play_a_game = function(theta)
{
  xt    = draw_starts() 
  res   = play(xt,delt,objects,rad,theta, trace = TRUE)
  score = mean(res$status==1)
  return(score)
}
play_a_game(theta_rand)

library('GA')
obj = play_a_game
GA  = ga(type = 'real-valued',fitness = play_a_game,lower = rep(-10,npars),upper = rep(10,npars), popSize = 100,maxiter = 100,keepBest = TRUE)
plot(GA)


theta_hat = GA@solution[1,]
theta_hat

win = c()
for(i in 1:100)
{

  #objects = draw_objects(J)
  xt_try    = draw_starts()
  plot_game(xt_try,objects,rad,'black')
  res_final = play(xt_try,delt,objects,rad,theta_hat,plt = FALSE,trace = TRUE)
  #print(typeof(res_final$trajectories))
  
  trajs = na.omit(res_final$trajectories)
  lines(trajs[, 2,]~trajs[, 1, ])
  win[i] = (res_final$status==1)
}
mean(win)

winNewObjs = c()
for(i in 1:100)
{
  
  objects = draw_objects(J)
  xt_try    = draw_starts()
  plot_game(xt_try,objects,rad,'black')
  res_final = play(xt_try,delt,objects,rad,theta_hat,plt = FALSE,trace = TRUE)
  #print(typeof(res_final$trajectories))
  
  trajs = na.omit(res_final$trajectories)
  lines(trajs[, 2,]~trajs[, 1, ])
  winNewObjs[i] = (res_final$status==1)
}
mean(winNewObjs)

### c

playAdj = function(x0,delt,objects,rad,theta,plt = FALSE,trace = FALSE)
{
  k           =   0 # Count how many steps
  xt          = x0   # Set the initial coordinate(s) for the drone.
  trajectories = NULL
  # Check the game status:
  res_status  = game_status(xt, objects, rad)
  status      = res_status$status
  minDists = res_status$minDist
  
  inverseDist = 1/(minDists-rad)
  # Check which games are still active:
  terminal      = (status != 0)
  if(trace)
  {
    trajectories = array(dim = c(dim(xt)[1],dim(xt)[2],101))
    trajectories[,,1] = xt
  }
  if(plt){plot_game(xt,objects,rad,'black')}
  while((any(status==0))&(k<100))
  {
    k = k + 1
    
    # Now, let the model update the position of the pieces:
    
    ct = controlAdjusted(xt, theta, inverseDist)
    xt = xt+ct*delt*cbind(1-terminal,1-terminal)
    
    # Checkk the game status after the positions are updates:
    res_status  = game_status(xt, objects, rad)
    status      = res_status$status
    minDists = res_status$minDist
    inverseDist = 1/(minDists-rad)
    terminal    = (status != 0)
    if(trace){trajectories[,,k] = xt}
    if(plt){plot_game(xt,objects,rad,c('red','black','green')[status+2])}
  }
  return(list(k = k, status = status,xt= xt,trajectories = trajectories))
}


controlAdjusted = function(xt,pars, inverseDist)
{
  xt_aug = cbind(xt, inverseDist)
  res_model = model(xt_aug,pars,rep(nodes,2))
  return(res_model$a3)
}

play_a_game_adj = function(theta)
{
  xt    = draw_starts() 
  res   = playAdj(xt,delt,objects,rad,theta, trace = TRUE)
  score = mean(res$status==1)
  return(score)
}

p     = 3
q     = 2
nodes = 3
npars = p*nodes+nodes*nodes+nodes*q+nodes+nodes+q
npars
theta_rand = runif(npars,-1,1)

objAdj = play_a_game_adj
GA_adj  = ga(type = 'real-valued',fitness = play_a_game_adj,lower = rep(-10,npars),upper = rep(10,npars), popSize = 100,maxiter = 100,keepBest = TRUE)
plot(GA_adj)


theta_hat_adj = GA_adj@solution[1,]
theta_hat_adj

win_adj = c()
for(i in 1:100)
{
  
  #objects = draw_objects(J)
  xt_try    = draw_starts()
  plot_game(xt_try,objects,rad,'black')
  res_final = playAdj(xt_try,delt,objects,rad,theta_hat_adj,plt = FALSE,trace = TRUE)
  #print(typeof(res_final$trajectories))
  
  trajs = na.omit(res_final$trajectories)
  lines(trajs[, 2,]~trajs[, 1, ])
  win_adj[i] = (res_final$status==1)
}
mean(win_adj)

winNewObjs_adj = c()
for(i in 1:100)
{
  
  objects = draw_objects(J)
  xt_try    = draw_starts()
  plot_game(xt_try,objects,rad,'black')
  res_final = playAdj(xt_try,delt,objects,rad,theta_hat_adj,plt = FALSE,trace = TRUE)
  #print(typeof(res_final$trajectories))
  
  trajs = na.omit(res_final$trajectories)
  lines(trajs[, 2,]~trajs[, 1, ])
  winNewObjs_adj[i] = (res_final$status==1)
}
mean(winNewObjs_adj)

### d
play_a_game_adj_newObs = function(theta)
{
  xt    = draw_starts() 
  objects = draw_objects(J)
  res   = playAdj(xt,delt,objects,rad,theta, trace = TRUE)
  score = mean(res$status==1)
  return(score)
}


objAdj = play_a_game_adj_newObs
GA_adj_newObs  = ga(type = 'real-valued',fitness = play_a_game_adj_newObs,lower = rep(-10,npars),upper = rep(10,npars), popSize = 100,maxiter = 100,keepBest = TRUE)
plot(GA_adj_newObs)


theta_hat_adj_newObs = GA_adj_newObs@solution[1,]
theta_hat_adj_newObs

win_adj_newObs = c()
for(i in 1:100)
{
  
  #objects = draw_objects(J)
  xt_try    = draw_starts()
  plot_game(xt_try,objects,rad,'black')
  res_final = playAdj(xt_try,delt,objects,rad,theta_hat_adj_newObs,plt = FALSE,trace = TRUE)
  #print(typeof(res_final$trajectories))
  
  trajs = na.omit(res_final$trajectories)
  lines(trajs[, 2,]~trajs[, 1, ])
  win_adj_newObs[i] = (res_final$status==1)
}
mean(win_adj_newObs)

winNewObjs_adj_genNewObs = c()
for(i in 1:100)
{
  
  objects = draw_objects(J)
  xt_try    = draw_starts()
  plot_game(xt_try,objects,rad,'black')
  res_final = playAdj(xt_try,delt,objects,rad,theta_hat_adj_newObs,plt = FALSE,trace = TRUE)
  #print(typeof(res_final$trajectories))
  
  trajs = na.omit(res_final$trajectories)
  lines(trajs[, 2,]~trajs[, 1, ])
  winNewObjs_adj_genNewObs[i] = (res_final$status==1)
}
mean(winNewObjs_adj_genNewObs)

####


#### Q2

rm(list = ls())

S = matrix(0,9,8)
# Row Sums:
S[1:3,1] = 1
S[4:6,2] = 1
S[7:9,3] = 1
# Col Sums:
S[c(1:3)*3-2,4] = 1
S[c(1:3)*3-1,5] = 1
S[c(1:3)*3-0,6] = 1
# Diag Sums
S[c(1,5,9),7] = 1
S[c(3,5,7),8] = 1
S


rho = function(m,S)
{
  m = matrix(m,ncol = 1)
  player1  = any(t(m)%*%S==-3) # X
  player2  = any(t(m)%*%S==+3) # O
  m_1      = (m==-1) #m_
  m_2      = (m==+1) #m+
  tie      = sum((t(m_1)%*%S>0)*(t(m_2)%*%S>0))==8
  winner   = c(-1,0,1)[c(player1,tie,player2)]
  terminal = player1|tie|player2
  ret      = list(terminal = terminal,winner = winner)
  return(ret)
}

# m = as.matrix(c(1,1,-1,0,0,0,0,0,-1))
# matrix(m,3,3,byrow = TRUE)
# rho(m,S)
# 
# m = as.matrix(c(1,1,-1,0,0,-1,0,0,-1)) 
# matrix(m,3,3,byrow = TRUE)
# rho(m,S)


game_tree= function(m,k, currMoves = 0, maxMoves = 9)
{
  g = c()
  game_state = rho(m,S)
  if(game_state$terminal)
  {
    g = c(g,game_state$winner)
    return(g)
  }
  else if(currMoves == maxMoves)
  {
    g = c(g, 2)
    return(g)
  }
  
  else{
    Index = which(m == 0)
    for(i in 1:length(Index))
    {
      x = m 
      x[Index[i]]=k
      g = c(g,game_tree(x,-1*k, currMoves + 1, maxMoves))
    }
    return(g)
  }
}

m = as.matrix(c(0,0,0,0,0,0,0,0,0))



res5 = game_tree(m,1, 0, 5)
n_g5  = 	length(res5)
Xwins5 = sum(res5==-1)
Draws5 = sum(res5== 0)
Owins5 = sum(res5==+1)
Unfinished5 = sum(res5==+2)
c(n_g5,Xwins5,Draws5,Owins5, Unfinished5)

res8 = game_tree(m,1, 0, 8)
n_g8  = 	length(res8)
Xwins8 = sum(res8==-1)
Draws8 = sum(res8== 0)
Owins8 = sum(res8==+1)
Unfinished8 = sum(res8==+2)
c(n_g8,Xwins8,Draws8,Owins8, Unfinished8)

res9 = game_tree(m,1, 0, 9)
n_g9  = 	length(res9)
Xwins9 = sum(res9==-1)
Draws9 = sum(res9== 0)
Owins9 = sum(res9==+1)
Unfinished9 = sum(res9==+2)
CombinedFull = c(n_g9,Xwins9/n_g9,Draws9/n_g9,Owins9/n_g9, Unfinished9)


mcts = function(m, k, alpha)
{
  g = c()
  game_state = rho(m,S)
  
  randAlpha = runif(1, 0, 1)
  
  if(game_state$terminal)
  {
    g = c(g,game_state$winner)
    return(g)
  }
  
  else
  {
    if (randAlpha <= alpha)
    {
      Index = which(m == 0)
      for(i in 1:length(Index))
      {
        x = m 
        x[Index[i]]=k
        g = c(g,mcts(x,-1*k, alpha))
      }
      return(g)
    }
    else
    {
      g = c(g,2)
      return(g)
    }

  }
}

m = as.matrix(c(1,1,0,0,0,0,0,0,-1))

resMCTS = mcts(m,-1, 0.95)
n_gMCTS  = 	length(resMCTS)
XwinsMCTS = sum(resMCTS==-1)
DrawsMCTS = sum(resMCTS== 0)
OwinsMCTS = sum(resMCTS==+1)
UnfinishedMCTS = sum(resMCTS==+2)
c(n_gMCTS,XwinsMCTS,DrawsMCTS,OwinsMCTS, UnfinishedMCTS)

m = as.matrix(c(1,1,-1,0,0,0,0,0,-1))

df90 = data.frame()
df70 = data.frame()

for(i in 1:100)
{
  resMCTS90 = mcts(m,1, 0.9)
  resMCTS70 = mcts(m,1, 0.7)
  
  n_gMCTS90  = 	length(resMCTS90)
  XwinsMCTS90 = sum(resMCTS90==-1)
  DrawsMCTS90 = sum(resMCTS90== 0)
  OwinsMCTS90 = sum(resMCTS90==+1)
  UnfinishedMCTS90 = sum(resMCTS==+2)
  Combined90 = c(n_gMCTS90,XwinsMCTS90/n_gMCTS90,DrawsMCTS90/n_gMCTS90,OwinsMCTS90/n_gMCTS90, UnfinishedMCTS90)
  
  n_gMCTS70  = 	length(resMCTS70)
  XwinsMCTS70 = sum(resMCTS70==-1)
  DrawsMCTS70 = sum(resMCTS70== 0)
  OwinsMCTS70 = sum(resMCTS70==+1)
  UnfinishedMCTS70 = sum(resMCTS70==+2)
  Combined70 = c(n_gMCTS70,XwinsMCTS70/n_gMCTS70,DrawsMCTS70/n_gMCTS70,OwinsMCTS70/n_gMCTS70, UnfinishedMCTS70)
  
  df90 = rbind(df90, Combined90)
  df70 = rbind(df70, Combined70)
}

colnames(df90) = c("Games", "Xwins%", "Ties%", "Owins%", "Incomplete")
colnames(df70) = c("Games", "Xwins%", "Ties%", "Owins%", "Incomplete")

boxplot(df90[, 2:4])
points(c(CombinedFull[2], CombinedFull[3], CombinedFull[4]), pch = "*", col = "red", cex = 3)
boxplot(df70[, 2:4])
points(c(CombinedFull[2], CombinedFull[3], CombinedFull[4]), pch = "*", col = "red", cex = 3)

mcts_minBranch = function(m, k, alpha, currMoves, minBranch)
{
  g = c()
  game_state = rho(m,S)
  
  randAlpha = runif(1, 0, 1)
  
  if(game_state$terminal)
  {
    g = c(g,game_state$winner)
    return(g)
  }
  
  else
  {
    if (currMoves <= minBranch)
    {
      Index = which(m == 0)
      for(i in 1:length(Index))
      {
        x = m 
        x[Index[i]]=k
        g = c(g,mcts_minBranch(x,-1*k, alpha, currMoves + 1, minBranch))
      }
      return(g)
    }
    else
    {
      if (randAlpha <= alpha)
      {
        Index = which(m == 0)
        for(i in 1:length(Index))
        {
          x = m 
          x[Index[i]]=k
          g = c(g,mcts_minBranch(x,-1*k, alpha, currMoves + 1, minBranch))
        }
        return(g)
      }
      else
      {
        g = c(g,2)
        return(g)
      }
    }
    
  }
}

# m = as.matrix(c(1,1,0,0,0,0,0,0,-1))
# 
# resMCTS_Branch = mcts_minBranch(m,-1, 0.95, 3, 9)
# n_gMCTS_Branch  = 	length(resMCTS_Branch)
# XwinsMCTS_Branch = sum(resMCTS_Branch==-1)
# DrawsMCTS_Branch = sum(resMCTS_Branch== 0)
# OwinsMCTS_Branch = sum(resMCTS_Branch==+1)
# UnfinishedMCTS_Branch = sum(resMCTS_Branch==+2)
# c(n_gMCTS_Branch,XwinsMCTS_Branch,DrawsMCTS_Branch,OwinsMCTS_Branch, UnfinishedMCTS_Branch)

m = as.matrix(c(1,1,-1,0,0,0,0,0,-1))

df90_Branch_new = data.frame()
df70_Branch_new = data.frame()

for(i in 1:100)
{
  resMCTS90_new = mcts_minBranch(m,1, 0.9, 4, 7)
  resMCTS70_new = mcts_minBranch(m,1, 0.7, 4, 7)
  
  n_gMCTS90_new  = 	length(resMCTS90_new)
  XwinsMCTS90_new = sum(resMCTS90_new==-1)
  DrawsMCTS90_new = sum(resMCTS90_new== 0)
  OwinsMCTS90_new = sum(resMCTS90_new==+1)
  UnfinishedMCTS90_new = sum(resMCTS90_new==+2)
  Combined90_new = c(n_gMCTS90_new,XwinsMCTS90_new/n_gMCTS90_new,
                     DrawsMCTS90_new/n_gMCTS90_new,OwinsMCTS90_new/n_gMCTS90_new, UnfinishedMCTS90_new)
  
  n_gMCTS70_new  = 	length(resMCTS70_new)
  XwinsMCTS70_new = sum(resMCTS70_new==-1)
  DrawsMCTS70_new = sum(resMCTS70_new== 0)
  OwinsMCTS70_new = sum(resMCTS70_new==+1)
  UnfinishedMCTS70_new = sum(resMCTS70_new==+2)
  Combined70_new = c(n_gMCTS70_new,XwinsMCTS70_new/n_gMCTS70_new,
                     DrawsMCTS70_new/n_gMCTS70_new,OwinsMCTS70_new/n_gMCTS70_new, UnfinishedMCTS70_new)
  
  df90_Branch_new = rbind(df90_Branch_new, Combined90_new)
  df70_Branch_new = rbind(df70_Branch_new, Combined70_new)
}

colnames(df90_Branch_new) = c("Games", "Xwins%", "Ties%", "Owins%", "Incomplete")
colnames(df70_Branch_new) = c("Games", "Xwins%", "Ties%", "Owins%", "Incomplete")

boxplot(df90_Branch_new[, 2:4], ylim = c(0, 1))
points(c(CombinedFull[2], CombinedFull[3], CombinedFull[4]), pch = "*", col = "red", cex = 3)
boxplot(df70_Branch_new[, 2:4], ylim = c(0, 1))
points(c(CombinedFull[2], CombinedFull[3], CombinedFull[4]), pch = "*", col = "red", cex = 3)

####



#### Q3 


library(quadprog)

polKern = function(x1, x2, gm, cf, dg)
{
  return((cf+gm*t(x1)%*%x2)^dg)
}

radialKern = function(x1, x2, gm)
{
  return(exp(-gm*(sum((x1 - x2)^2))))
}

my_svm = function(y, x, kern, cost = 0, softMarg = FALSE, gm = 0, cf = 0, dg = 0, plt = FALSE)
{
  N = dim(x)[1]
  DD = matrix(0,N,N)
  
  if(kern == "none")
  {
    for(i in 1:N)
    {
      for(j in 1:N)
      {
        DD[i,j]  = y[i]*y[j]*(t(x[i,])%*%x[j,])
      }
    }
  }
  else if (kern == "poly")
  {
    
    for(i in 1:N)
    {
      for(j in 1:N)
      {
        KK = polKern(x[i,], x[j,], gm, cf, dg)
        DD[i,j]  = y[i]*y[j]*KK
      }
    }
  }
  else if (kern == "radial")
  {
    for(i in 1:N)
    {
      for(j in 1:N)
      {
        KK = radialKern(x[i,], x[j,], gm)
        DD[i,j]  = y[i]*y[j]*KK
      }
    }
  }
  
  eps = 5e-6
  DD  = DD+eps*diag(N)
  Amat = cbind(y,diag(N)) # y will be on first row of t(Amat)
  bvec = matrix(0,N+1,1)
  d    = matrix(1,N,1)
  
  if (softMarg == TRUE)
  {
    negativeC = (-1)*cost
    Amat = cbind(Amat, -diag(N))
    Cvec = rep(negativeC, N)
    vec0 = rep(0, N+1)
    bvec = matrix(c(vec0, Cvec), ncol = 1, nrow = (2*N + 1))
  }
  print(dim(Amat))
  print(dim(bvec))
  
  res = solve.QP(Dmat = DD,dvec = d,Amat = Amat,bvec = bvec,meq = 1,factorized = FALSE)
  a   = res$solution
  
  if (plt == TRUE)
  {
    plot(a,type= 'h', main = expression(alpha[i]),xlab = expression(i), lwd = 2)
  }
  
  pad.a     = round(a,3)
  wh        = which.max(a)
  
  # T1 = rep(0,N)
  # for(i in 1:N)
  # {
  #   if (kern == "none")
  #   {
  #     T1[i] = sum(a*y*(x[i,]%*%t(x)))
  #     intercept = 1/y[wh]-sum(a*y*(t(x[i,])%*%x))
  #   }
  #   else if (kern == "poly")
  #   {
  #     KK = polKern(x[i,], t(x), gm, cf, dg)
  #     
  #     print(KK)
  #     T1[i] = sum(a*y*KK)
  #     KKNew = polKern(x[wh,], t(x), gm, cf, dg)
  #     intercept = 1/y[wh]-sum(a*y*KKNew)
  #   }
  #   else if (kern == "radial")
  #   {
  #     KK = radialKern(x[i,], t(x), gm)
  #     T1[i] = sum(a*y*KK)
  #     KKNew = radialKern(x[wh,], t(x), gm)
  #     intercept = 1/y[wh]-sum(a*y*KKNew)
  #   }
  # }
  
  
  if (kern == "none")
  {
    ww = t(a*y)%*%X
    intercept = 1/y[wh] - X[wh, ]%*%t(ww)
    
    yhat = sign(X%*%t(ww) + intercept[1])
    
  }
  else if (kern == "poly")
  {
    
    T1 = rep(0,N)
    for(i in 1:N)
    {
      KK = polKern(x[i, ], t(X), gm, cf, dg)
      T1[i] = sum(a*y*(KK))
    }
    print(KK)
    
    KKNew = polKern(x[wh,], t(x), gm, cf, dg)
    intercept = 1/y[wh]-sum(a*y*KKNew)
    
    yhat      = sign(T1+intercept[1])
  }
  else if (kern == "radial")
  {
    T1 = rep(0,N)
    for(i in 1:N)
    {
      T1[i] = sum(a*y*(x[i,]%*%t(x)))
    }
    
    KKNew = radialKern(x[wh,], t(x), gm)
    intercept = 1/y[wh]-sum(a*y*KKNew)
    
    yhat      = sign(T1+intercept[1])
  }
  
  res = list("yhat" = yhat, "padA" = pad.a, "a" = a, "intercept" = intercept)
  return(res)
  
}

PLADat = read.table("PLA Dynamics.txt")

plot(X2 ~ X1, col = (Y + 3), data = PLADat)

# Create data matrix:
X  = cbind(PLADat$X1,PLADat$X2)
Y = PLADat$Y

gm = 2
cf = 1
dg = 2

mySVM = my_svm(Y, X, "radial", 10000, TRUE, gm, cf, dg, plt = TRUE)

yhat = mySVM$yhat
padA = mySVM$padA


par(mfrow = c(2,2))
plot(X2~X1,pch = c(1,16)[(Y+1)/2+1], data = PLADat)
plot(PLADat$X2~PLADat$X1, pch = 16, col = c('grey','blue')[(yhat+1)/2+1])

wh.text = which(padA!=0)
points(PLADat$X2~PLADat$X1, pch = 1, col = c(NA,'red')[(padA>0)+1],cex=2)
text(PLADat$X2[wh.text]~PLADat$X1[wh.text],labels = paste0('n = ',wh.text))



library('e1071')

model  = svm(Y~(X1 + X2), data = PLADat, scale = FALSE,kernel = 'polynomial',degree =dg,gamma = gm,coef0 = cf,cost = 10000)

# Our solution
plot(X2~X1,pch = c(1,16)[(Y+1)/2+1], data = PLADat)
plot(PLADat$X2~PLADat$X1, pch = 16, col = c('grey','blue')[(yhat+1)/2+1])

# Package solution
plot(1, 1, type = "n", xlim = c(-25, 25), ylim = c(-17, 17))
points(model$SV[,2]~model$SV[,1],pch = '+', col = 'blue',cex = 2)
 N = dim(X)[1]
plot(model$coefs~model$index,type = 'h',xlim = c(0,N))

 ####