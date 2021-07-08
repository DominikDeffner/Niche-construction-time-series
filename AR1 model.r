
# Code for double-hierarchical model
library(rethinking)

setwd("C:/Users/dominik_deffner/Documents/GitHub/Niche-construction-time-series")

#Prepare data
d <- read.csv("Temporal Analysis.csv")

#Non-constructed (autonomous) cases
dat <- list(N = length(d$Grad.linear.value),
            N_cat = 3,
            N_Study = length(unique(d$Subset)),
            Grad_obs = d$Grad.linear.value,
            Grad_sd = d$Grad.linear.StErr,
            Study = d$Subset)

dat$Length_Study <- sapply(1:dat$N_Study, function(x) length(which(dat$Study == x)) )

NC <- c()
for (j in 1:dat$N) {
  if (d$NC[j] == "N"){
    NC[j] <- 1
  } else if (d$NC[j] == "M"){
    NC[j] <- 2
  } else {
    NC[j] <- 3
  }
}


dat$NC <- sapply(1:dat$N_Study, function(i) unique(NC[which(d$Subset == i)]))


AR1 <- "

data{

    int N;
    int N_cat;
    int N_Study;
    int Length_Study[N_Study]; 
    int NC[N_Study];
    real Grad_obs[N];
    int Study[N];
    real Grad_sd[N];
    
}

parameters{

    vector[N] Grad_est;         // Estimated gradients
    
    real a_NC[N_cat];
    real beta[N_cat];
    real<lower=0> sigma[N_cat];

}


model{
    int pos; // Counter for position
    int L;   // Local var for length of subsets

    a_NC ~ normal(0,1);
    beta ~ normal(0,1);
    sigma ~ exponential(3);
  
  
  // Observed gradients come from distribution of 'true' gradients and error
  
  Grad_obs ~ normal( Grad_est , Grad_sd );

   // Loop over all studies (we want timeseries nested in studies)
  
   pos = 1;
  
   for (k in 1:N_Study){
    
    L = Length_Study[k];

    segment(Grad_est, pos, Length_Study[k])[2:L] ~ normal( a_NC[NC[k]] + beta[NC[k]] * segment(Grad_est, pos, Length_Study[k])[1:(L - 1)], sigma[NC[k]]);
    pos = pos + Length_Study[k];
  } 
    
}


"


AR1_random <- "

data{

    int N;
    int N_cat;
    int N_Study;
    int Length_Study[N_Study]; 
    int NC[N_Study];
    real Grad_obs[N];
    int Study[N];
    real Grad_sd[N];
    
}

parameters{

    vector[N] Grad_est;         // Estimated gradients
    
    real a_NC[N_cat];
    real beta[N_cat];
    real log_sigma[N_cat];


matrix[3,N_Study] z_Study;               //Matrix of uncorrelated z - values
vector<lower=0>[3] sigma_Study;       //SD of parameters among individuals
cholesky_factor_corr[3] Rho_Study;    // This is the Cholesky factor: if you multiply this matrix and it's transpose you get correlation matrix

}

transformed parameters{

matrix[N_Study,3] v_Study; // Matrix of varying effects for each subset
v_Study = ( diag_pre_multiply( sigma_Study , Rho_Study ) * z_Study )';

}

model{
    int pos; // Counter for position
    int L;   // Local var for length of subsets
    real mu;
    real sigma;
    
    a_NC ~ normal(0,1);
    beta ~ normal(0,1);
    log_sigma ~ normal(-3,1);
  
  
  
   // Varying effects
   to_vector(z_Study) ~ normal(0,1);
   sigma_Study ~ exponential(1);
   Rho_Study ~ lkj_corr_cholesky(4);
  
  
    Grad_obs ~ normal( Grad_est , Grad_sd );


   // Loop over all studies (we want timeseries nested in studies)
  
   pos = 1;
  
   for (k in 1:N_Study){
    
    L = Length_Study[k];

    segment(Grad_est, pos, Length_Study[k])[2:L] ~ normal( (a_NC[NC[k]] + v_Study[k,1]) + (beta[NC[k]] + v_Study[k,2]) * segment(Grad_est, pos, Length_Study[k])[1:(L - 1)], exp(log_sigma[NC[k]] + v_Study[k,3]) );
    pos = pos + Length_Study[k];
  } 
    
}


"


#Run stan models
m <- stan( model_code= AR1 , data=dat , chains=4, iter = 10000, cores=4, control=list(adapt_delta=0.99, max_treedepth=15))

s <- extract.samples(m)




library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values


####
###
# Temporal Dynamics Plot
###
###


graphics.off()
png("Ar1SimRand.png", res = 900, height = 18, width = 12, units = "cm")

seq <- 2:50

par(mfrow = c(dat$N_cat,1),
    oma=c(3,2.5,2,0),
    mar=c(3,2.5,2,0.2))

for (NC in 1:dat$N_cat) {
  
  
  plot(c(1,seq), type = "n", ylim = c(-0.4,0.4), bty = "n")
  
  for (i in sample(1:length(s$a_NC[,NC]), 100)) {
    
    y <- c()
    y[1] <- rnorm(1, mean(s$a_NC[,NC]), mean(exp(s$log_sigma[,NC]))) 
    for (t in seq) {
      mu <- s$a_NC[i,NC] + s$beta[i,NC] * y[t-1]
      sigma <- exp(s$log_sigma[i, NC])
      
      y[t] <- rnorm(1, mu, sigma )
    }
    
    lines(y, col = alpha("darkgrey", alpha = ifelse(NC==1, 0.4,0.4)))
  }
  
  if (NC == 1) mtext("Autonomous", side = 3, line=-1.5, cex = 1.3) 
  if (NC == 2) mtext("Mixed", side = 3, line=-1.5, cex = 1.3) 
  if (NC == 3) mtext("Constructed", side = 3, line=-1.5, cex = 1.3)
  #abline(h=0, lty = 2)
  
}
mtext("Time [years]", side = 1, line=1.5, cex = 1.1, outer = TRUE)
mtext("Selection Gradient", side = 2, line=1.2, cex = 1.1, outer = TRUE)
dev.off()


####
###
# Densities
###
####


graphics.off()
png("Ar1DensRand.png", res = 900, height = 9, width = 18, units = "cm")

par(mfrow = c(1,2),
    oma=c(0,1,0.1,0),
    mar=c(3,3,2,0.5))


NC = 1
dens <- density(s$beta[,NC])
x1 <- min(which(dens$x >= quantile(s$beta[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$beta[,NC], 0.95)))
plot(dens, xlim = c(-1,1), ylim = c(0,7), type="n", ann = FALSE)
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)


NC = 2
dens <- density(s$beta[,NC])
x1 <- min(which(dens$x >= quantile(s$beta[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$beta[,NC], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

NC = 3
dens <- density(s$beta[,NC])
x1 <- min(which(dens$x >= quantile(s$beta[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$beta[,NC], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

legend("topleft", c("Autonomous", "Mixed", "Constructed"), col = c(alpha(col.pal[1],alpha = 1),alpha(col.pal[2],alpha = 1),alpha(col.pal[3],alpha = 1)), cex = 0.8, lty =c(1,2,3), lwd = 2, bty = "n" )
mtext("Density", side = 2, line=3, cex = 1.1)
mtext("Temporal Stability", side = 3, line=0.5, cex = 1.1)

NC = 1
dens <- density(exp(s$log_sigma[,NC]))
x1 <- min(which(dens$x >= quantile(exp(s$log_sigma[,NC]), 0.05)))
x2 <- max(which(dens$x <  quantile(exp(s$log_sigma[,NC]), 0.95)))
plot(dens, xlim = c(0,0.1), ylim = c(0,100), type="n", ann = FALSE)
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)



NC = 2
dens <- density(exp(s$log_sigma[,NC]))
x1 <- min(which(dens$x >= quantile(exp(s$log_sigma[,NC]), 0.05)))
x2 <- max(which(dens$x <  quantile(exp(s$log_sigma[,NC]), 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

NC = 3
dens <- density(exp(s$log_sigma[,NC]))
x1 <- min(which(dens$x >= quantile(exp(s$log_sigma[,NC]), 0.05)))
x2 <- max(which(dens$x <  quantile(exp(s$log_sigma[,NC]), 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
mtext("Residual variation", side = 3, line=0.5, cex = 1.1)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

dev.off()






















####
###
# Temporal Dynamics Plot
###
###


graphics.off()
png("Ar1Sim.png", res = 900, height = 18, width = 12, units = "cm")

seq <- 2:50

par(mfrow = c(dat$N_cat,1),
    oma=c(3,2.5,2,0),
    mar=c(3,2.5,2,0.2))

for (NC in 1:dat$N_cat) {
  
  
  plot(c(1,seq), type = "n", ylim = c(-1,1), bty = "n")
  
  for (i in sample(1:length(s$a_NC[,NC]), 100)) {
    
    y <- c()
    y[1] <- rnorm(1, mean(s$a_NC[,NC]), mean(s$sigma[,NC])) 
    for (t in seq) {
      mu <- s$a_NC[i,NC] + s$beta[i,NC] * y[t-1]
      sigma <- s$sigma[i, NC]
      
      y[t] <- rnorm(1, mu, sigma )
    }
    
    lines(y, col = alpha("darkgrey", alpha = ifelse(NC==1, 0.4,0.4)))
  }
  
  if (NC == 1) mtext("Autonomous", side = 3, line=-1.5, cex = 1.3) 
  if (NC == 2) mtext("Mixed", side = 3, line=-1.5, cex = 1.3) 
  if (NC == 3) mtext("Constructed", side = 3, line=-1.5, cex = 1.3)
  #abline(h=0, lty = 2)
  
}
mtext("Time [years]", side = 1, line=1.5, cex = 1.1, outer = TRUE)
mtext("Selection Gradient", side = 2, line=1.2, cex = 1.1, outer = TRUE)
dev.off()


####
###
# Densities
###
####


graphics.off()
png("Ar1Dens.png", res = 900, height = 9, width = 18, units = "cm")

par(mfrow = c(1,2),
    oma=c(0,1,0.1,0),
    mar=c(3,3,2,0.5))


NC = 1
dens <- density(s$beta[,NC])
x1 <- min(which(dens$x >= quantile(s$beta[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$beta[,NC], 0.95)))
plot(dens, xlim = c(-1,1), ylim = c(0,25), type="n", ann = FALSE)
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)


NC = 2
dens <- density(s$beta[,NC])
x1 <- min(which(dens$x >= quantile(s$beta[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$beta[,NC], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

NC = 3
dens <- density(s$beta[,NC])
x1 <- min(which(dens$x >= quantile(s$beta[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$beta[,NC], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

legend("topleft", c("Autonomous", "Mixed", "Constructed"), col = c(alpha(col.pal[1],alpha = 1),alpha(col.pal[2],alpha = 1),alpha(col.pal[3],alpha = 1)), cex = 0.8, lty =c(1,2,3), lwd = 2, bty = "n" )
mtext("Density", side = 2, line=3, cex = 1.1)
mtext("Temporal Stability", side = 3, line=0.5, cex = 1.1)

NC = 1
dens <- density(s$sigma[,NC])
x1 <- min(which(dens$x >= quantile(s$sigma[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$sigma[,NC], 0.95)))
plot(dens, xlim = c(0,0.20), ylim = c(0,100), type="n", ann = FALSE)
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)



NC = 2
dens <- density(s$sigma[,NC])
x1 <- min(which(dens$x >= quantile(s$sigma[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$sigma[,NC], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

NC = 3
dens <- density(s$sigma[,NC])
x1 <- min(which(dens$x >= quantile(s$sigma[,NC], 0.05)))
x2 <- max(which(dens$x <  quantile(s$sigma[,NC], 0.95)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[NC],alpha = 0.3), border = NA), add=TRUE)
mtext("Residual variation", side = 3, line=0.5, cex = 1.1)
lines(dens, col=alpha(col.pal[NC],alpha = 1), lty = NC, lwd = 2)

dev.off()


