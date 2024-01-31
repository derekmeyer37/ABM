set.seed(123)

# Parameters
population_size <- 5000
initial_infected <- 50        # Prevalence = 0.10
infection_probability <- 0.15 # Probability of becoming infected
exposure_probability <- 0.1  # probability of becoming exposed 
exposure_duration <- 8        # Duration an agent remains exposed before becoming infectious
infection_duration <- 7       # Duration an agent remains infectious before possibility of recovering
recovery_probability <- 0.7   # Probability of recovering after the baseline infection duration
simulation_days <- 60

# Create agents with initial states
agents <- data.frame(
  id = 1:population_size,
  state = rep("susceptible", population_size),
  days_exposed = rep(0, population_size),
  days_infected = rep(0, population_size)
)

# Initialize some agents as infected
initial_infected_ids <- sample(1:population_size, initial_infected)
agents$state[initial_infected_ids] <- "infected"

# Create a data frame to store the simulation results
simulation_results <- data.frame(
  day = 1:simulation_days,
  susceptible = numeric(simulation_days),
  exposed = numeric(simulation_days),
  infected = numeric(simulation_days),
  recovered = numeric(simulation_days)
)

# Simulate the spread of COVID-19 and record daily counts
for (day in 1:simulation_days) {
  for (i in 1:population_size) {
    if (agents$state[i] == "susceptible") {
      # Check for potential exposure
      if (runif(1) <= exposure_probability) {
        agents$state[i] <- "exposed"
      }
    }
    else if (agents$state[i] == "exposed") {
      agents$days_exposed[i] <- agents$days_exposed[i] + 1
      if (agents$days_exposed[i] >= exposure_duration) {
        agents$state[i] <- "infected"
        agents$days_exposed[i] <- 0
      }
    } else if (agents$state[i] == "infected") {
      agents$days_infected[i] <- agents$days_infected[i] + 1
      if (agents$days_infected[i] >= infection_duration) {                       # if the number of days an agent has been infected is greater than 5, 
         if (runif(1) < recovery_probability) {                                  # and if a random draw is less than the recovery probability
          agents$state[i] <- "recovered" 
          agents$days_infected[i] <- 0                                           # the agent becomes recovered or remains infected
         } 
      }
    } 
  }

  # Record daily counts
  simulation_results$susceptible[day] <- sum(agents$state == "susceptible")
  simulation_results$exposed[day] <- sum(agents$state == "exposed")
  simulation_results$infected[day] <- sum(agents$state == "infected")
  simulation_results$recovered[day] <- sum(agents$state == "recovered")
}

library(ggplot2)
ggplot(simulation_results, aes(x = day)) +
  geom_line(aes(y = susceptible, color = "Susceptible")) +
  geom_line(aes(y = exposed, color = "Exposed")) +
  geom_line(aes(y = infected, color = "Infected")) +
  geom_line(aes(y = recovered, color = "Recovered")) +
  labs(x = "Days", y = "Number of Individuals", color = "State") +
  scale_color_manual(values = c("Susceptible" = "blue", "Exposed" = "orange", 
               "Infected" = "red", "Recovered" = "green"),
    limits = c("Susceptible", "Exposed", "Infected", "Recovered")) +
  ggtitle("Susceptible-Exposed-Infected-Recovered ABM") +
  theme_minimal()
  
  
