# ABM_v2 rewrites ABM.R to be more robust and efficient
AgentGen1 <- function(nPop1, E0, I0, vaccine = NA){                                           # population size, n_exposed, n_infected
  # Create a population of susceptible agents
  Agent1 <- data.frame(AgentNo = 1:nPop1,
                       State = 'S',
                       Mixing = runif(nPop1,0,1),
                       TimeE = 0,
                       TimeI = 0,
                       Vaccination = 0,
                       stringsAsFactors = FALSE)
  # Initializing exposed and infected agents
  Agent1$State[1:E0] <- "E"
  Agent1$TimeE[1:E0] <- rbinom(E0, 13, 0.5) + 1                                 # Agents can be in exposed state for 14 days (+1 so that rbinom doesn't return 0)
  Agent1$State[(E0+1):(E0 + I0)] <- "I"    
  Agent1$TimeI[(E0+1):(E0 + I0)] <- rbinom(I0, 9, 0.5) + 1                      # Agents can be in infected state for 10 days (+1 so that rbinom doesn't return 0) - infectious period
  # Randomly assign vaccination status based on the specified prevalence
  if (!is.na(vaccine$prevalence)) {
    num_vaccinated <- round(nPop1 * vaccine$prevalence)
    Agent1$Vaccination[sample(1:nPop1, num_vaccinated)] <- 1
  }
  return(Agent1)
}

ABM1 <- function(Agent1, par1, nTime1, vaccine = NA) {                                        # Function of agents, model parameters, and simulation time
  nPop1 <- nrow(Agent1)
  Out1 <- data.frame(S = rep(0, nTime1),                                        # Initializing output df
                     E = rep(0, nTime1),
                     I = rep(0, nTime1),
                     R = rep(0, nTime1),
                     D = rep(0, nTime1))
  for(k in 1:nTime1){
    # Move agents through time
    StateS1 <- (1:nPop1)[Agent1$State == "S"] 
    StateSE1 <- (1:nPop1)[Agent1$State == "S" | Agent1$State == "E"] 
    for(i in StateS1) {
      # Determine if they like to meet others
      Mix1 <- Agent1$Mixing[i]
      # How many agents will they meet (Each will meet up to MaxMix other agents per day)
      Meet1 <- round(Mix1*par1$MaxMix,0) + 1
      # Grab the agents they will meet
      Meet2 <- sample(StateSE1, Meet1, replace = TRUE, prob = Agent1$Mixing[StateSE1])     # Replace = TRUE means they could meet the same agent over and over 
      
      for(j in 1:length(Meet2)) {
        # Grab who they will meet 
        Meet1a <- Agent1[Meet2[j], ]
        # If exposed, change state
        if(Meet1a$State == "E"){
          Urand1 <- runif(1,0,1)
          if(Urand1 < par1$S2E){
          Agent1$State[i] <- "E" 
          }
        }
      }
    }
    # Grab agents who have been exposed and add a day to their exposure time 
    StateE1 <- (1:nPop1)[Agent1$State == "E"] 
    Agent1$TimeE[StateE1] = Agent1$TimeE[StateE1] + 1
    StateE2 <- (1:nPop1)[Agent1$State == "E" & Agent1$TimeE > 14]               # If exposure time is greater than 14, move to recovered
    Agent1$State[StateE2] <- "R"
    # Grab agents who could possibly become sick
    StateE3 <- (1:nPop1)[Agent1$State == "E" & Agent1$TimeE > par1$incubation_period]                # Incubation days: the time it takes for an infection to develop after an agent has been exposed
    for(i in StateE3){
      # Randomly assign whether they get sick or not
      Urand1 <- runif(1,0,1)
      if(any(is.na(vaccine))){
        if(Urand1 < par1$E2I) {
          Agent1$State[i] <- "I"
        }
      }
      else {
        if(Urand1 < par1$E2I - vaccine$E2I_reduction) { 
          Agent1$State[i] <- "I"
        }
      }
    }
    
    # Update how long the agents have been sick
    StateI1 <- (1:nPop1)[Agent1$State == "I"]
    Agent1$TimeI[StateI1] = Agent1$TimeI[StateI1] + 1
    StateI2 <- (1:nPop1)[Agent1$State == "I" & Agent1$TimeI > 14]
    Agent1$State[StateI2] <- "R"
    # Grab agents who have been infected for 10 days or less to see if they die
    StateI3 <- (1:nPop1)[Agent1$State == "I" & Agent1$TimeI < par1$infectious_period + 1]
    Agent1$State[StateI3] <- ifelse(
      runif(length(StateI3),0,1) > par1$I2D, "I", "D")
    
    Out1$S[k] <- length(Agent1$State[Agent1$State == "S"])
    Out1$E[k] <- length(Agent1$State[Agent1$State == "E"])
    Out1$I[k] <- length(Agent1$State[Agent1$State == "I"])
    Out1$R[k] <- length(Agent1$State[Agent1$State == "R"])
    Out1$D[k] <- length(Agent1$State[Agent1$State == "D"])
  }
  return(Out1)
}

# Define a function to calculate moving average
library(zoo)
smooth_lines <- function(data, window_size = 5) {
  for (col in names(data)) {
    data[[col]] <- rollmean(data[[col]], k = window_size, align = "center", fill = NA)
  }
  return(data)
}

# Model Parameters
set.seed(123)
par1 <- data.frame(MaxMix = 7, 
                   S2E = 0.12,  # Probability an agent moves from susceptible to exposed
                   E2I = 0.15,   # Probability an agent moves from exposed to infected - prob_infection
                   I2D = 0.01, # Probability an agent moves from infected to deceased - death_rate
                   incubation_period = 6,
                   infectious_period = 10)  

vaccine_tool <- data.frame(E2I_reduction = 0.05, 
                           prevalence = 0.5)
masking_tool <- data.frame(S2E_reduction = 0.05,
                           prevalence = 0.5)

# Running the Model
Agent1 <- AgentGen1(nPop1 = 1000, E0 = 5, I0 = 2, vaccine = vaccine_tool)
model <- ABM1(Agent1, par1, nTime1 = 60, vaccine = vaccine_tool)

# Plotting
# Smooth the columns using a moving average with a window size of 5
smoothed_output <- smooth_lines(model, window_size = 2)
plot(smoothed_output$S, type = "l", col = "blue",
     xlab = "Time (days)",
     ylab = "Number of Individuals", main = "SEIRD Agent-Based Model",
     ylim = c(0, 1000))
lines(smoothed_output$E, col = "orange")
lines(smoothed_output$I, col = "red")
lines(smoothed_output$R, col = "green")
lines(smoothed_output$D, col = "black")

# Add legend
legend("right", legend = c("Susceptible", "Exposed", "Infected", "Recovered",
                              "Deceased"), col = c("blue", "orange", "red",
                                                   "green", "black"), lty = 1)
print(model)

