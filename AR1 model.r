
# Code for double-hierarchical model
library(rethinking)

setwd("C:/Users/dominik_deffner/Documents/GitHub/Niche-construction-time-series")

#Prepare data
d <- read.csv("Temporal Analysis.csv")

#Non-constructed (autonomous) cases
dat <- list(N = length(d$Grad.linear.value),
            N_Study = length(unique(d$Subset)),
            Grad_obs = d$Grad.linear.value,
            Grad_sd = d$Grad.linear.StErr,
            Study = d$Subset)

dat$Length_Study <- sapply(1:dat$N_Study, function(x) length(which(dat$Study == x)) )


Combined <- ifelse(d$NC == "N", 1, 2)
dat$NC <- sapply(1:dat$N_Study, function(i) unique(Combined[which(d$Subset == i)]))



AR1 <- "

data{

    int N;
    int N_Study;
    int Length_Study[N_Study]; 
    int NC[N_Study];
    real Grad_obs[N];
    int Study[N];
    real Grad_sd[N];
    
}

parameters{

    vector[N] Grad_est;         // Estimated gradients
    
    
    real b_NC[2];
    real beta[2];
    real<lower=0> sigma[2];
}

model{
    int pos; // Counter for position
    int L;

  
  
  b_NC ~ normal(0,1);
  beta ~ normal(0,1);
  sigma ~ exponential(1);
  
  Grad_obs ~ normal( Grad_est , Grad_sd );


  // Loop over all studies (we want timeseries nested in studies)
  
  pos = 1;
  
  for (k in 1:N_Study){
    
    L = Length_Study[k];

    segment(Grad_est, pos, Length_Study[k])[2:L] ~ normal(b_NC[NC[k]] + beta[NC[k]] * segment(Grad_est, pos, Length_Study[k])[1:(L - 1)], sigma[NC[k]]);
    pos = pos + Length_Study[k];
  } 
    
}


"


#Run stan models
m <- stan( model_code= AR1 , data=dat , chains=1, iter = 3000, cores=1, control=list(adapt_delta=0.9, max_treedepth=15))

s <- extract.samples(m)
