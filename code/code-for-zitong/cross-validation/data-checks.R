
#how eta.hat parameter defined for JAGS model 
eta.hat <- pt.data$bx.pgg[is.na(eta.data)] 

sum(is.na(pt.data$bx.pgg)) #should be 0
sum(is.na(eta.hat)) #should be 0

(length(eta.hat) + n_eta_known) == n #should be true

length(eta.track)==n_mask #should be true

(n_mask + n_eta_known)==length(eta.true) #should be true 


#### if any of the above aren't true, please send me the output including the values for these variables here
#they are all different counts of patients 
n #total sample size
n_mask #number in left out fold
n_eta_known #number with true state known after left-out fold excluded
length(eta.hat) 
length(eta.track)
length(eta.true)