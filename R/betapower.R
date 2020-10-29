#Given mu and SD, calculate beta parameters a and b
calculate_a_b <- function(mu, SD){
  var=SD^2
  a=((mu^2)*(1-mu) - (var*mu))/var
  b=a*(1-mu)/mu
  return(c(a,b))
}

#Given parameters a and b and sample size n, simulate beta numbers
simulate_beta <- function(a,b,n){
  sim <- rbeta(n,a,b)
  return(sim)
}

#Calculate power through N simulations, assuming case and control vectors of beta numbers
##Calculate p-value for simulation
beta_p <- function(case, control){
  #Record p-value for simulation
  data <- data.frame("values"=c(case,control),"type" = c(rep("case",length(case)), rep("control",length(control))))
  return(coef(summary(betareg(values ~ type, data = data)))$mean[2,4])
}
##Calculate power at significance level a or the fraction of N simulation p-values that are below a
beta_power <- function(mu0, sd0, mu1, sd1, n0, n1, N,a,MHT,q){
  a0 <- calculate_a_b(mu0, sd0)[1]
  b0 <- calculate_a_b(mu0, sd0)[2]
  a1 <- calculate_a_b(mu1, sd1)[1]
  b1 <- calculate_a_b(mu1, sd1)[2]
  p_values <- c()
  pb = txtProgressBar(min = 0, max = N, initial = 0, title = "Progress")
  for(i in 1:N){
    #add progress bar
    setTxtProgressBar(pb, i, title="Progress", label= paste0(round((i/N)*100, 0),"% done"))
    p <- beta_p(simulate_beta(a1, b1, n1), simulate_beta(a0, b0, n0))
    p_values <- c(p_values,p)
  }
  close(pb)
  p_adj <- p.adjust(p_values, method = MHT)
  power_p <- sum(p_values < a)/length(p_values)
  power_adj <- sum(p_adj < q)/length(p_adj)
  return(list(power_p, p_values, power_adj, p_adj))
}

#Test multiple mu1s for given sample size
mu1_curve_data <- function(mu0, sd0, mu1_min, mu1_max, mu1_step, sd1, n0, n1, N, a, MHT, q){
  #Get all mu1 values
  mu1s <- seq(mu1_min, mu1_max, mu1_step)
  power_p_est <- c()
  power_adj_est <- c()
  for(i in 1:length(mu1s)){
    print(paste0("Calculating power for mu1=",mu1s[i],"..."))
    pr <- beta_power(mu0=mu0, sd0=sd0, mu1=mu1s[i], sd1=sd1, n0=n0, n1=n1, N=N, a=a, MHT=MHT, q=q)
    power_p_est <- c(power_p_est,pr[[1]])
    power_adj_est <- c(power_adj_est,pr[[3]])
  }
  data <- data.frame("mu1"=mu1s, "power_p"=power_p_est, "power_adj"=power_adj_est)
  data$mu0 <- mu0
  return(data[,c("mu1","mu0","power_p","power_adj")])
}

#Produce dataframe for different case:control ratios
beta_curve <- function(mu0, sd0, mu1_min, mu1_max, mu1_step, sd1, n0s, n1s, N, a, MHT, q){
  if(length(n1s) != length(n0s)){
    stop("The vectors of n1 and n0 sample sizes should be of equal length!")
  }
  print("Calculating power...")
  final_data <- data.frame(matrix(ncol=6))
  colnames(final_data) <- c("mu1","mu0","power_p", "power_adj", "n1","n0")
  for(i in 1:length(n1s)){
    n1 <- n1s[i]
    n0 <- n0s[i]
    print(paste0("Collecting data with case:control ratio of ",n1,":",n0))
    data <- mu1_curve_data(mu0=mu0, sd0=sd0, mu1_min = mu1_min, mu1_max = mu1_max, mu1_step = mu1_step, sd1=sd1, n0=n0, n1=n1, N=N, a=a, MHT=MHT, q=q)
    data$n1 <- n1
    data$n0 <- n0
    final_data <- rbind(final_data,data)
  }
  final_data <- final_data[-1,]
  final_data$effect_size <- final_data$mu1-final_data$mu0
  final_data$sample_size <- factor(paste0(final_data$n1,":",final_data$n0))
  return(final_data)
}

plot_beta_curve <- function(data, a, MHT, q){
  effect_size <- power_p <- sample_size <- power_adj <- NULL
  ggplot(data) + geom_line(aes(x=effect_size,y=power_p,color=sample_size)) +geom_line(aes(x=effect_size, y=power_adj, color=sample_size), linetype="dashed")+ labs(color="Case:Control (n)") + theme_bw() + xlab("Effect size (Case - Control)") + ylab("Power") + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 0.8, 0.9, 1)) + ggtitle(paste0("Solid lines: p-value < ",a, ", Dashed lines:",MHT," < ",q))
}



