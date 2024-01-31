# ABM_v4 implements the model capability for social distancing
# -------------------------------- Functions -----------------------------------
# population size, initial_exposed, initial_infected, vaccine, social distancing
AgentGen1 <- function(nPop, E0, I0, vaccine = NA, social_distancing = NA){
  # Create a population of susceptible agents
  Agent1 <- data.frame(AgentNo = 1:nPop,
                       State = 'S',
                       Mixing = runif(nPop,0,1),
                       TimeE = 0,
                       TimeI = 0,
                       Vaccination = 0,
                       Social_Distancing = 0,
                       stringsAsFactors = FALSE)
  # Initializing exposed and infected agents
  Agent1$State[1:E0] <- "E"
  # Agents can be in exposed state for 14 days (+1 so that rbinom doesn't
  # return 0), randomly assigning exposure time
  Agent1$TimeE[1:E0] <- rbinom(E0, 13, 0.5) + 1
  Agent1$State[(E0+1):(E0 + I0)] <- "I"
  # Agents can be in infected state for 10 days (+1 so that rbinom doesn't
  # return 0) (infectious period) randomly assigning infectious time
  Agent1$TimeI[(E0+1):(E0 + I0)] <- rbinom(I0, 9, 0.5) + 1
  # Randomly assign vaccination/social distancing status based on the specified
  # prevalence
  num_vaccinated <- round(nPop * vaccine$prevalence)
  Agent1$Vaccination[sample(1:nPop, num_vaccinated)] <- 1
  num_social_distancing <- round(nPop * social_distancing$prevalence)
  Agent1$Social_Distancing[sample(1:nPop, num_social_distancing)] <- 1
  return(Agent1)
}

# Function of agents, model parameters, and simulation time
ABM1 <- function(Agent1, par1, nTime1, vaccine = NA, social_distancing = NA) {
  nPop <- nrow(Agent1)
  # Initializing output df, all states start at 0
  output <- data.frame(S = rep(0, nTime1),
                       E = rep(0, nTime1),
                       I = rep(0, nTime1),
                       R = rep(0, nTime1),
                       D = rep(0, nTime1))
  # The model cycles through every day of the simulation
  for(k in 1:nTime1){
    # Move agents through time
    StateS1 <- (1:nPop)[Agent1$State == "S"]
    StateSE1 <- (1:nPop)[Agent1$State == "S" | Agent1$State == "E"]

    # Cycling through each susceptible agent
    for(i in StateS1) {
      # Here, the higher number of exposed/infected individuals there are, the
      # more likely a susceptible agent is to come in contact and become exposed
      # Determine if they like to meet others
      Mix1 <- Agent1$Mixing[i]
      # How many agents will they meet (Each will meet up to MaxMix other agents
      # per day)
      Meet1 <- round(Mix1*par1$MaxMix,0) + 1
      # Grab the agents they will meet: sampling susceptible or exposed
      # individuals of size meet1, using mixing probability
      # Replace = TRUE means they could meet the same agent over and over
      Meet2 <- sample(StateSE1, Meet1, replace = TRUE,
                      prob = Agent1$Mixing[StateSE1])

      # Cycling through the agents each susceptible person will meet
      for(j in 1:length(Meet2)) {
        # Grab who they will meet
        Meet1a <- Agent1[Meet2[j], ]
        if(Meet1a$State == "E"){
          Urand1 <- runif(1,0,1)
          # Apply social distancing if the agent practices it
          if (Agent1$Social_Distancing[i] == 1) {
            if (Urand1 < par1$S2E - social_distancing$S2E_reduction) {
              Agent1$State[i] <- "E"
            }
          } else {
              # If not practicing social distancing, apply regular exposure
              if (Urand1 < par1$S2E) {
                Agent1$State[i] <- "E"
              }
          }
        }
      }
    }
    # Grab agents who have been exposed and add a day to their exposure time
    StateE1 <- (1:nPop)[Agent1$State == "E"]
    Agent1$TimeE[StateE1] = Agent1$TimeE[StateE1] + 1
    # If exposure time is greater than 14, move to recovered
    StateE2 <- (1:nPop)[Agent1$State == "E" & Agent1$TimeE > 14]
    Agent1$State[StateE2] <- "R"
    # Grab agents who could possibly become sick
    # Incubation days: the time it takes for an infection to develop after an
    # agent has been exposed
    StateE3 <- (1:nPop)[Agent1$State == "E" &
                         Agent1$TimeE > par1$incubation_period]
    for(i in StateE3){
      # Randomly assign whether they get sick or not based on vaccination status
      Urand1 <- runif(1,0,1)
      if (Agent1$Vaccination[i] == 1) {
          if (Urand1 < par1$E2I - vaccine$E2I_reduction) {
            Agent1$State[i] <- "I"
          }
      } else {
          if (Urand1 < par1$E2I) {
            Agent1$State[i] <- "I"
          }
      }
    }

    # Update how long the agents have been sick
    StateI1 <- (1:nPop)[Agent1$State == "I"]
    Agent1$TimeI[StateI1] = Agent1$TimeI[StateI1] + 1
    StateI2 <- (1:nPop)[Agent1$State == "I" & Agent1$TimeI > 14]
    Agent1$State[StateI2] <- "R"
    # Grab agents who have been infected for 10 days or less to see if they die
    StateI3 <- (1:nPop)[Agent1$State == "I" &
                         Agent1$TimeI < par1$infectious_period + 1]
    Agent1$State[StateI3] <- ifelse(
      runif(length(StateI3),0,1) > par1$I2D, "I", "D")

    output$S[k] <- length(Agent1$State[Agent1$State == "S"])
    output$E[k] <- length(Agent1$State[Agent1$State == "E"])
    output$I[k] <- length(Agent1$State[Agent1$State == "I"])
    output$R[k] <- length(Agent1$State[Agent1$State == "R"])
    output$D[k] <- length(Agent1$State[Agent1$State == "D"])
  }
  return(output)
}

# Define a function to calculate moving average
library(zoo)
smooth_lines <- function(data, window_size = 5) {
  for (col in names(data)) {
    data[[col]] <- rollmean(data[[col]], k = window_size, align = "center",
                            fill = NA)
  }
  return(data)
}

# --------------------------- Building/Running Model ---------------------------

# Model Parameters
set.seed(123)
par1 <- data.frame(MaxMix = 7,  # Maximum number each agent can meet in a day
                   S2E = 0.10,  # prob_exposure
                   E2I = 0.15,  # prob_infection
                   I2D = 0.01,  # death_rate
                   incubation_period = 6, # Time between when an agent is
                                          # exposed to potentially becoming an
                                          # infectious agent
                   infectious_period = 9) # Time which an agent is infectious

vaccine_tool <- data.frame(E2I_reduction = 0.05,
                           prevalence = 0.0)
social_distancing_tool <- data.frame(S2E_reduction = 0.05,
                           prevalence = 0.0)

# Running the Model
Agent1 <- AgentGen1(nPop = 1000, E0 = 10, I0 = 5, vaccine = vaccine_tool,
                    social_distancing = social_distancing_tool)
model <- ABM1(Agent1, par1, nTime1 = 60, vaccine = vaccine_tool,
              social_distancing = social_distancing_tool)

# --------------------------------- Plotting -----------------------------------

# Plotting
# Smooth the columns using a moving average with a window size of 2
smoothed_output <- smooth_lines(model, window_size = 2)
plot(smoothed_output$S, type = "l", col = "blue",
     xlab = "Time (days)",
     ylab = "Number of Individuals", main = "SEIRD ABM - COVID-19",
     ylim = c(0, 1000))
lines(smoothed_output$E, col = "orange")
lines(smoothed_output$I, col = "red")
lines(smoothed_output$R, col = "green")
lines(smoothed_output$D, col = "black")

# Add legend
legend("right", legend = c("Susceptible", "Exposed", "Infected", "Recovered",
                              "Deceased"), col = c("blue", "orange", "red",
                                                   "green", "black"), lty = 1)


