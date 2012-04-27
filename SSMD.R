ssmdMLE <- function(pop1, pop2)
{
  factor <- function(pop)
  {
    return ((length(pop) - 1) / length(pop) * sd(pop) * sd(pop))
  }
  beta <- mean(pop1) - mean(pop2)
  beta <- beta / sqrt(factor(pop1) + factor(pop2))
  return (beta)
}

ssmdrobust <- function(pop1, pop2)
{
  beta <- median(pop1) - median(pop2)
  beta <- beta / sqrt((mad(pop1))^2 + (mad(pop2))^2)
  return (beta)
}

ssmdMM <- function(pop1, pop2)
{
  beta <- mean(pop1) - mean(pop2)
  beta <- beta / sqrt((sd(pop1))^2 + (sd(pop2))^2)
  return (beta)  
}