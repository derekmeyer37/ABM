# ABM_v5 fixes the main issue with exposed agents exposing susceptible agents
# now, only infected agents can expose susceptible agents to COVID-19.
# Other minor fixes

# -------------------------------- Functions -----------------------------------
# population size, initial_exposed, initial_infected, vaccine, social distancing
agent_generator <- function(n_pop, E0, I0, vaccine = NA, social_distancing = NA){
  # Create a population of susceptible agents
  agents <- data.frame(agent = 1:n_pop,
                       state = 'S',
                       mixing = runif(n_pop,0,1),
                       timeE = 0,
                       timeI = 0,
                       vaccination = 0,
                       social_distancing = 0,
                       stringsAsFactors = FALSE)
  # Initializing exposed and infected agents
  agents$state[1:E0] <- "E"
  # agents can be in exposed state for 14 days (+1 so that rbinom doesn't
  # return 0), randomly assigning exposure time
  agents$timeE[1:E0] <- rbinom(E0, 13, 0.5) + 1
  agents$state[(E0+1):(E0 + I0)] <- "I"
  # agents can be in infected state for 10 days (+1 so that rbinom doesn't
  # return 0) (infectious period) randomly assigning infectious time
  agents$timeI[(E0+1):(E0 + I0)] <- rbinom(I0, 9, 0.5) + 1
  # Randomly assign vaccination/social distancing status based on the specified
  # prevalence
  num_vaccinated <- round(n_pop * vaccine$prevalence)
  agents$vaccination[sample(1:n_pop, num_vaccinated)] <- 1
  num_social_distancing <- round(n_pop * social_distancing$prevalence)
  agents$social_distancing[sample(1:n_pop, num_social_distancing)] <- 1
  return(agents)
}

# Function of agents, model parameters, and simulation time
abm <- function(agents, parms, n_days, vaccine = NA, social_distancing = NA) {
  n_pop <- nrow(agents)
  # Initializing output df, all states start at 0
  output <- data.frame(S = rep(0, n_days),
                       E = rep(0, n_days),
                       I = rep(0, n_days),
                       R = rep(0, n_days),
                       D = rep(0, n_days))
  # The model cycles through every day of the simulation
  for(k in 1:n_days){
    # Move agents through time
    stateS <- (1:n_pop)[agents$state == "S"]
    stateSI <- (1:n_pop)[agents$state == "S" | agents$state == "I"]

    # Cycling through each susceptible agent
    for(i in stateS) {
      # Here, the higher number of infected individuals there are, the more 
      # likely a susceptible agent is to come in contact and become exposed
      # Determine if they like to meet others
      mix <- agents$mixing[i]
      # How many agents will they meet (Each will meet up to max_mix other agents
      # per day)

      # [2024-02-07] George:
      # For identifying how many agents they could meet, you could leverage the 
      # binomial distribution. If you want them to meet on average 7 other agents,
      # then, you could do the following:
      #
      # meet1 <- rbinom(1, size = n_pop, prob = parms$max_mix/n_pop)
      #
      # That will yield on average max_mix agents. This is only a recommendation.
      # You can keep the current approach if you want.

      meet1 <- round(mix*parms$max_mix,0) + 1
      # Grab the agents they will meet: sampling susceptible or exposed
      # individuals of size meet1, using mixing probability
      # Replace = TRUE means they could meet the same agent over and over
      meet2 <- sample(stateSI, meet1, replace = TRUE,
                      prob = agents$mixing[stateSI])

      # Cycling through the agents each susceptible person will meet
      for(j in 1:length(meet2)) {
        # Grab who they will meet
        meet1a <- agents[meet2[j], ]
        
        if(meet1a$state == "I"){
          Urand1 <- runif(1)
          # Apply social distancing if the agent practices it
          if (agents$social_distancing[i] == 1) {
            if (Urand1 < parms$S2E - social_distancing$S2E_reduction) {
              agents$state[i] <- "E"
              break
            }
          } else {
              # If not practicing social distancing, apply regular exposure
              if (Urand1 < parms$S2E) {
                agents$state[i] <- "E"
                break
              }
          }
        }
      }
    }
    # Grab agents who have been exposed and add a day to their exposure time
    stateE1 <- (1:n_pop)[agents$state == "E"]
    agents$timeE[stateE1] = agents$timeE[stateE1] + 1
    # If exposure time is greater than 14, move to recovered
    stateE2 <- (1:n_pop)[agents$state == "E" & agents$timeE > 14]
    agents$state[stateE2] <- "R"
    # Grab agents who could possibly become sick
    # Incubation days: the time it takes for an infection to develop after an
    # agent has been exposed
    stateE3 <- (1:n_pop)[agents$state == "E" &
                         agents$timeE > parms$incubation_period]
    for(i in stateE3){
      # Randomly assign whether they get sick or not based on vaccination status
      Urand1 <- runif(1)
      if (agents$vaccination[i] == 1) {
          if (Urand1 < parms$E2I - vaccine$E2I_reduction) {
            agents$state[i] <- "I"
          }
      } else {
          if (Urand1 < parms$E2I) {
            agents$state[i] <- "I"
          }
      }
    }

    # Update how long the agents have been sick
    stateI1 <- (1:n_pop)[agents$state == "I"]
    agents$timeI[stateI1] = agents$timeI[stateI1] + 1

    # [2024-02-07] George
    # Having a hard-stop at 14 days of the infection is OK, but I would also
    # randomize this. What you can do is have a prob. of recovery with mean
    # 14 days. All you need is to do the same randomization you have done for
    # other things. e.g.
    #
    # agents$state[stateI2] <- ifelse(runif(length(stateI2),0,1) < 1/14, "R", "I")
    #
    # But is all up to you.
    stateI2 <- (1:n_pop)[agents$state == "I" & agents$timeI > 14]
    agents$state[stateI2] <- "R"
    # Grab agents who have been infected for 10 days or less to see if they die
    stateI3 <- (1:n_pop)[agents$state == "I" &
                         agents$timeI < parms$infectious_period + 1]
    agents$state[stateI3] <- ifelse(
      runif(length(stateI3),0,1) > parms$I2D, "I", "D")

    output$S[k] <- length(agents$state[agents$state == "S"])
    output$E[k] <- length(agents$state[agents$state == "E"])
    output$I[k] <- length(agents$state[agents$state == "I"])
    output$R[k] <- length(agents$state[agents$state == "R"])
    output$D[k] <- length(agents$state[agents$state == "D"])
  }
  return(output)
}

# Define a function to calculate moving average
library(zoo)
smooth_lines <- function(data, window_size) {
  for (col in names(data)) {
    data[[col]] <- rollmean(data[[col]], k = window_size, align = "center",
                            fill = NA)
  }
  return(data)
}

# --------------------------- Building/Running Model ---------------------------

# Model Parameters
set.seed(123)
parms <- data.frame(max_mix = 7,  # Maximum number each agent can meet in a day
                   S2E = 0.2,  # prob_exposure
                   E2I = 0.15,  # prob_infection
                   I2D = 0.01,  # prob_death
                   incubation_period = 6, # Time between when an agent is
                                          # exposed to potentially becoming an
                                          # infectious agent
                   infectious_period = 9) # Time which an agent is infectious

vaccine_tool <- data.frame(E2I_reduction = 0.07,
                           prevalence = 0.0)
social_distancing_tool <- data.frame(S2E_reduction = 0.05,
                           prevalence = 0.0)

# Running the Model
agents <- agent_generator(n_pop = 1000, E0 = 50, I0 = 5, vaccine = vaccine_tool,
                          social_distancing = social_distancing_tool)
model <- abm(agents, parms, n_days = 60, vaccine = vaccine_tool,
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
