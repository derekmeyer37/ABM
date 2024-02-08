# ABM_v5 fixes the main issue with exposed agents exposing susceptible agents
# now, only infected agents can expose suceptible agents to COVID-19.

# -------------------------------- Functions -----------------------------------
# population size, initial_exposed, initial_infected, vaccine, social distancing
AgentGen <- function(nPop, E0, I0, vaccine = NA, social_distancing = NA){
  # Create a population of susceptible agents
  Agents <- data.frame(AgentNo = 1:nPop,
                       State = 'S',
                       Mixing = runif(nPop,0,1),
                       TimeE = 0,
                       TimeI = 0,
                       Vaccination = 0,
                       Social_Distancing = 0,
                       stringsAsFactors = FALSE)
  # Initializing exposed and infected agents
  Agents$State[1:E0] <- "E"
  # Agents can be in exposed state for 14 days (+1 so that rbinom doesn't
  # return 0), randomly assigning exposure time
  Agents$TimeE[1:E0] <- rbinom(E0, 13, 0.5) + 1
  Agents$State[(E0+1):(E0 + I0)] <- "I"
  # Agents can be in infected state for 10 days (+1 so that rbinom doesn't
  # return 0) (infectious period) randomly assigning infectious time
  Agents$TimeI[(E0+1):(E0 + I0)] <- rbinom(I0, 9, 0.5) + 1
  # Randomly assign vaccination/social distancing status based on the specified
  # prevalence
  num_vaccinated <- round(nPop * vaccine$prevalence)
  Agents$Vaccination[sample(1:nPop, num_vaccinated)] <- 1
  num_social_distancing <- round(nPop * social_distancing$prevalence)
  Agents$Social_Distancing[sample(1:nPop, num_social_distancing)] <- 1
  return(Agents)
}

# Function of agents, model parameters, and simulation time
ABM <- function(Agents, par1, nTime1, vaccine = NA, social_distancing = NA) {
  nPop <- nrow(Agents)
  # Initializing output df, all states start at 0
  output <- data.frame(S = rep(0, nTime1),
                       E = rep(0, nTime1),
                       I = rep(0, nTime1),
                       R = rep(0, nTime1),
                       D = rep(0, nTime1))
  # The model cycles through every day of the simulation
  for(k in 1:nTime1){
    # Move agents through time
    StateS <- (1:nPop)[Agents$State == "S"]
    StateSI <- (1:nPop)[Agents$State == "S" | Agents$State == "I"]

    # Cycling through each susceptible agent
    for(i in StateS) {
      # Here, the higher number of infected individuals there are, the more 
      # likely a susceptible agent is to come in contact and become exposed
      # Determine if they like to meet others
      Mix1 <- Agents$Mixing[i]
      # How many agents will they meet (Each will meet up to MaxMix other agents
      # per day)

      # [2024-02-07] George:
      # For identifying how many agents they could meet, you could leverage the 
      # binomial distribution. If you want them to meet on average 7 other agents,
      # then, you could do the following:
      #
      # Meet1 <- rbinom(1, size = nPop, prob = par1$MaxMix/nPop)
      #
      # That will yield on average MaxMix agents. This is only a recommendation.
      # You can keep the current approach if you want.

      Meet1 <- round(Mix1*par1$MaxMix,0) + 1
      # Grab the agents they will meet: sampling susceptible or exposed
      # individuals of size meet1, using mixing probability
      # Replace = TRUE means they could meet the same agent over and over
      Meet2 <- sample(StateSI, Meet1, replace = TRUE,
                      prob = Agents$Mixing[StateSI])

      # Cycling through the agents each susceptible person will meet
      for(j in 1:length(Meet2)) {
        # Grab who they will meet
        Meet1a <- Agents[Meet2[j], ]
        
        if(Meet1a$State == "I"){
          Urand1 <- runif(1)
          # Apply social distancing if the agent practices it
          if (Agents$Social_Distancing[i] == 1) {
            if (Urand1 < par1$S2E - social_distancing$S2E_reduction) {
              Agents$State[i] <- "E"
              break
            }
          } else {
              # If not practicing social distancing, apply regular exposure
              if (Urand1 < par1$S2E) {
                Agents$State[i] <- "E"
                break
              }
          }
        }
      }
    }
    # Grab agents who have been exposed and add a day to their exposure time
    StateE1 <- (1:nPop)[Agents$State == "E"]
    Agents$TimeE[StateE1] = Agents$TimeE[StateE1] + 1
    # If exposure time is greater than 14, move to recovered
    StateE2 <- (1:nPop)[Agents$State == "E" & Agents$TimeE > 14]
    Agents$State[StateE2] <- "R"
    # Grab agents who could possibly become sick
    # Incubation days: the time it takes for an infection to develop after an
    # agent has been exposed
    StateE3 <- (1:nPop)[Agents$State == "E" &
                         Agents$TimeE > par1$incubation_period]
    for(i in StateE3){
      # Randomly assign whether they get sick or not based on vaccination status
      Urand1 <- runif(1)
      if (Agents$Vaccination[i] == 1) {
          if (Urand1 < par1$E2I - vaccine$E2I_reduction) {
            Agents$State[i] <- "I"
          }
      } else {
          if (Urand1 < par1$E2I) {
            Agents$State[i] <- "I"
          }
      }
    }

    # Update how long the agents have been sick
    StateI1 <- (1:nPop)[Agents$State == "I"]
    Agents$TimeI[StateI1] = Agents$TimeI[StateI1] + 1

    # [2024-02-07] George
    # Having a hard-stop at 14 days of the infection is OK, but I would also
    # randomize this. What you can do is have a prob. of recovery with mean
    # 14 days. All you need is to do the same randomization you have done for
    # other things. e.g.
    #
    # Agents$State[StateI2] <- ifelse(runif(length(StateI2),0,1) < 1/14, "R", "I")
    #
    # But is all up to you.
    StateI2 <- (1:nPop)[Agents$State == "I" & Agents$TimeI > 14]
    Agents$State[StateI2] <- "R"
    # Grab agents who have been infected for 10 days or less to see if they die
    StateI3 <- (1:nPop)[Agents$State == "I" &
                         Agents$TimeI < par1$infectious_period + 1]
    Agents$State[StateI3] <- ifelse(
      runif(length(StateI3),0,1) > par1$I2D, "I", "D")

    output$S[k] <- length(Agents$State[Agents$State == "S"])
    output$E[k] <- length(Agents$State[Agents$State == "E"])
    output$I[k] <- length(Agents$State[Agents$State == "I"])
    output$R[k] <- length(Agents$State[Agents$State == "R"])
    output$D[k] <- length(Agents$State[Agents$State == "D"])
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
                   S2E = 0.25,  # prob_exposure
                   E2I = 0.15,  # prob_infection
                   I2D = 0.01,  # prob_death
                   incubation_period = 6, # Time between when an agent is
                                          # exposed to potentially becoming an
                                          # infectious agent
                   infectious_period = 9) # Time which an agent is infectious

vaccine_tool <- data.frame(E2I_reduction = 0.05,
                           prevalence = 0.0)
social_distancing_tool <- data.frame(S2E_reduction = 0.05,
                           prevalence = 0.0)

# Running the Model
Agents <- AgentGen(nPop = 1000, E0 = 50, I0 = 5, vaccine = vaccine_tool,
                    social_distancing = social_distancing_tool)
model <- ABM(Agents, par1, nTime1 = 60, vaccine = vaccine_tool,
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
