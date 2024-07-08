#clear workspacerm
rm(list = ls())

# load some libraries

library(tidyverse)
library(reshape2)
library(ggplot2)

#### to run this scripts, you need three data sets:
# a) the sex file, which is basically just a file containing one column with each individuals name, and a second column with their sex
# b) the interaction file, which contains the pant grunts. Basically, just three columns: Winner, Loser, and Date
# c) the presence file. This one is a bit more complicated, it is basically a file with one row for each day, and each individual has it's column, and 0 indicates that the individual was not part of the hierarchy at that point, while 1 indicates it was.
# for me, I have only included individuals into the presence file (and therefore, the hierarchy), if they fulfilled the followin criteria:
# they only get a 1 if they are more than 10 years/ 120 months old (males) or more than 12 years / 144 months old (females). They are also only included if they were in the community for at least 2 years after reaching that age.
# i usually only get a death month and migration month from the genealogy, so individuals get a 1 starting the month they reach the age limit, and a 0 again the month they disappear
# why I have restricted this is because otherwise, individuals who only appear for a short period of time, like natal females who reach maturity and subsequently migrate, mess with everyone elses ranks. If you need clear ranks for younger individuals or newer individuals, feel free to change this

# let me know if any questions arise

# read in data: presence data, sex file and pant grunt data; might still have to prepare it
pres <- read.csv("/Users/gb64/Desktop/PhD-GesturalOrigins/Data - clean/Waibira_pres_day.csv", header = T, stringsAsFactors = T, sep = ",") # read the presence file for your group
subj.to.keep <- colnames(pres)
elo.data <- read.csv("/Users/gb64/Desktop/PhD-GesturalOrigins/Data - clean/Waibira_pgs_20220123.csv", header = T, stringsAsFactors = T, sep = ",") #### read or prepare data that contains the pant grunts, with Winner, Loser and Date
elo.data <- subset(elo.data, Winner %in% subj.to.keep & Loser %in% subj.to.keep)

demo.data <- read.csv("/Users/gb64/Desktop/PhD-GesturalOrigins/Data - clean/Wab_sex.csv", header = T, stringsAsFactors = T, sep = ",") # file containing the names and their sex ("m", "f")
elo.data$Date <- as.Date(elo.data$Date,format ="%d/%m/%Y")

## attach sexes to elo.data and remove females for Issa data
elo.data <- elo.data %>%
  left_join(demo.data, by = c("Winner" = "ID")) %>%
  left_join(demo.data, by = c("Loser" = "ID")) %>%
  filter(Sex.x == "M" & Sex.y == "M") %>%
  select(-Sex.x, -Sex.y) %>% arrange(Date)

# prepare data
pres <- data.frame(pres)
Date <- as.Date(as.character(pres$Date), 
                format = "%d/%m/%Y") # create date for presence file
pres <- pres[, c(intersect(colnames(pres), 
                           c(as.character(elo.data$Loser), 
                             as.character(elo.data$Winner))))] # remove all individuals from presence that never have any pant grunts
pres <- cbind(Date, pres)
pres <- subset(pres, 
               Date >= min(as.Date(elo.data$Date, format = "%d/%m/%Y")) & 
                 Date <= max(as.Date(elo.data$Date, format = "%d/%m/%Y"))) # select only presence for the time period that you actually have data for

### make sure to only use pant grunts of individuals when they were "present" (mainly important if males already have data before they are counted)
elo.data$use <- 0
elo.data$Loser <- as.character(elo.data$Loser)
elo.data$Winner <- as.character(elo.data$Winner)
elo.data$Date <- as.character(elo.data$Date)
pres$Date <- as.character(pres$Date)

##GB added step because previous step introduces duplicates to presence data
pres<-pres[!duplicated(pres),]

#back to original script
for (i in 1:nrow(elo.data)) { # this thing checks if the individuals in the interactions where present on the day. if they were not, the interaction is removed
  if (pres[pres$Date == elo.data$Date[i], elo.data$Winner[i]] == 1 & pres[pres$Date == elo.data$Date[i], elo.data$Loser[i]] == 1) {
    elo.data$use[i] <- 1
  }
}
elo.data <- subset(elo.data, use == 1)
elo.data <- subset(elo.data, Winner != Loser)


#### congratulations, the data should be prepared now. It is time to create the functions used in the Foerster et al Paper. I just copied their scripts. There are three models: one testing Elo with everyone starting on average and estimated k, two testing Elo with everyone starting on bottom and estimated k, and three estimating the start value and k (takes longer)
### As I have argued before, the first two models of the Foerster paper make little sense, so I will leave them in here but won't calculate them
### Okay, the following is the Elo function as defined by Foerster et al. 2016. The function optimises both k (the strength with which each pg changes the rank) and the entry point of each female

elo.model3 <-function(par, IA_data, all_ids, return_likelihood = T)
{
  k <-par[1]
  init_elo <-par[2:length(par)]
  # Initialize output columns
  if(!return_likelihood)IA_data$elo_l_before <- IA_data$elo_w_before <- IA_data$elo_l_after <- IA_data$elo_w_after <-NA
  # Set intitial elo scores
  currentELO <-c(init_elo)
  names(currentELO) <-all_ids
  # Initialize the log likelihood
  L <-0
  # Start loop
  for(i in 1:nrow(IA_data))
  {
    ind1 <- which(names(currentELO)==IA_data$Winner[i])
    ind2 <- which(names(currentELO)==IA_data$Loser[i])
    if (!return_likelihood){
      IA_data$elo_w_before[i]<-currentELO[ind1]
      IA_data$elo_l_before[i]<-currentELO[ind2]
    }
    # calculate predited winning probablity of the winner
    p_win <-1/(1+exp(-.01*(currentELO[ind1]-currentELO[ind2])))
    # Calculation of new ELO scores
    currentELO[ind1] <- currentELO[ind1] + exp(k) * (1 - p_win)
    # new Elo score of the Winner
    currentELO[ind2] <- currentELO[ind2] - exp(k) * (1 - p_win)
    # new Elo score of the Loser
    # write calculated elo scores to output columns
    if (!return_likelihood)
    {IA_data$elo_w_after[i] <- currentELO[ind1]
    IA_data$elo_l_after[i] <- currentELO[ind2]
    }
    # Update log likelihood
    L <- L + log(p_win)}
  if (return_likelihood) return(-1*L)
  else return(IA_data)
}


# now we apply that function to our dataset
male_ago <- droplevels(elo.data)
male_presence <- pres
all_males <- sort(unique(c(
  as.character(male_ago$Winner),
  as.character(male_ago$Loser)
)))

# create the dataframe in which the results will be saved
results_f <- data.frame(
  "model" = 3,
  "convergence" = NA,
  "AIC" = NA, "delta_AIC" = NA, "k" = NA, "pred_accuracy" = NA
)

## okay, now here I add an additional complication. In the Foerster paper, they calculate the optimal k and start values by optimising only the winning likelihood. This is problematic,
# because the winning likelihood can be optimized by ignoring rank changes in really sparse data sets, which is what we have for the females. Their version basically tries to
# find start values that allow the k to be 0 in that case, but we know we have rank changes in the females. So what our script does is that it does not only try to optimize the winning
# likelihood, but also the number of correct classifications. Classifications are incorrect in 3 scenarios: when the assistants collected funny data, when there is an actual rank change,
# or when there was a rank change and the algorithm ignored it. We can't do anything about the 1st, and we want to get to the 2nd, so we minimize the 3rd. For that, we make use of a
# property of the optimisation they use: it initially gets better, as the k and start values get closer to the true value. The number of correct classifications also gets better.
# But then, there is a point where the k and start values are already very good, and we hit the maximum of correct classification. After that, the k gets too small, and the number of correct
# classifications goes down. So we identify the point with the maximum of correct classifications, and use the values identified there as our parameters for the actual Elo calculation

# my optimisation for small datasets differs from the original Foerster analysis like this:
# while they optimise the starting point and the change parameter k, I include an additional optimisation step that minimises the number of incongruent pant grunts (higher-ranking pg lower-ranking)
# I include that step because we use a method that was developed for contest (eg chess, aggressions), where two individuals have a relative strength that can matter. This is not the case for pant grunts:
# with a pant grunt, an individual acknowledges that the other is higher-ranking. There is not really a relative strength. So, optimising the distance between the individuals is not really necessary:
# what we need to optimise is that preferably, there should be a minimum number of incongruent pant grunts: only the ones that really signify a rank change between individuals

# so what the algorithm now does is to run the Foerster et al. 2016 algorithm, but with different numbers of maximum iterations. That way, different solutions are offered, and we afterwards pick the one with the lowest number of incongruent pant grunts.

vals <- data.frame(iterations = 0, k = 0, incongruent = 0)
vals <- vals[vals$iterations != 0, ]
for (i in seq(1, 100, by = 2)) {
  # Fitting model 3
  res_fem_model3 <- optim(
    par = c(5, rep(1, length(all_males))), elo.model3, all_ids = all_males, IA_data = male_ago,
    return_likelihood = T, method = "BFGS", control = list(maxit = i, reltol = 0.000000005)
  )

  elo.opt <- elo.model3(par = res_fem_model3$par, IA_data = male_ago, all_ids = all_males, return_likelihood = F)
  elo.long <- elo.opt
  elo.long$congr <- ""
  elo.long$congr[elo.long$elo_w_before > elo.long$elo_l_before] <- "congruent"
  elo.long$congr[elo.long$elo_w_before < elo.long$elo_l_before] <- "incongruent"
  xx <- data.frame(iterations = 0, k = 0, incongruent = 0)
  xx$iterations <- i
  xx$k <- exp(res_fem_model3$par[1])
  xx$incongruent <- length(elo.long$congr[elo.long$congr == "incongruent"])
  vals <- rbind(vals, xx)
}
save.image(file="/Users/gb64/Desktop/Elo_partial.RData")
### now, we take the number of iterations that had the smallest number of incongruent pant grunts and apply it to the Foerster et al algorithm

iterations <- min(vals$iterations[vals$incongruent == min(vals$incongruent)])
vals.f <- vals
res_fem_model3 <- optim(
  par = c(5, rep(1, length(all_males))), elo.model3, all_ids = all_males, IA_data = male_ago,
  return_likelihood = T, method = "BFGS", control = list(maxit = iterations, reltol = 0.000000005, trace = 1)
)

results_f$convergence[1] <- res_fem_model3$convergence
results_f$AIC[1] <- res_fem_model3$value * 2 + 2 * (length(all_males) + 1)
results_f$k[1] <- exp(res_fem_model3$par[1])
results_issa_f <- results_f # contains the information about the k etc
res_issa <- res_fem_model3


# now, apply that best model to the presence data to create a matrix that has one value for each male for each day
elo.issa <- elo.model3(par = res_issa$par, IA_data = elo.data, all_ids = sort(unique(c(as.character(elo.data$Winner), as.character(elo.data$Loser)))), return_likelihood = F)
presence <- pres
presence[presence == 0] <- NA
presence[, 2:ncol(presence)] <- presence[, 2:ncol(presence)] * 1000
elo.long <- elo.issa
winner <- elo.long[, c("Date", "Winner", "elo_w_after")]
loser <- elo.long[, c("Date", "Loser", "elo_l_after")]
colnames(winner) <- c("Date", "ID", "elo")
colnames(loser) <- c("Date", "ID", "elo")
elo.short <- rbind(winner, loser)
elo.short$Date <- as.Date(elo.short$Date)
elo.short <- aggregate(elo.short$elo, elo.short[, c("Date", "ID")], max)
elo.short <- elo.short[order(elo.short$Date), ]
elo.short$ID <- as.character(elo.short$ID)

for (i in 1:nrow(elo.short)) {
  nr <- which(presence$Date >= elo.short$Date[i] & !is.na(presence[,elo.short$ID[i]]))
  presence[nr, elo.short$ID[i]] <- elo.short$x[i]
}



presence[presence == 1000] <- 0.1
for (i in 1:ncol(presence)) {
  a <- presence[, i]
  b <- a[!is.na(a) & a != 0.1]
  c <- b[1]
  a[a == 0.1] <- c
  presence[, i] <- a
}


for (i in 2:ncol(presence)) {
  nr <- which(demo.data$Code == colnames(presence)[i])
  presence[as.Date(presence$Date) < demo.data$From[nr] | as.Date(presence$Date) > demo.data$To[nr], i] <- NA
}


Date <- presence$Date
presence$Date <- NULL
for (i in 1:nrow(presence)) {
  presence[i, ] <- rank(presence[i, ], na.last = "keep", ties.method = "average") / sum(!is.na(presence[i, ]))
}
presence <- cbind(presence, Date)

elos.males <- presence
#scores on final date
elo.males.end<-as.data.frame(t(subset(presence, Date=="2022-01-08")))
write.csv(elo.males.end,"/Users/gb64/Desktop/rank_data/Issa_Eloseq_Alex_script.csv")

### The elos.males object should now have the individual elo scores, by date. To plot, we use the melt() function from the reshape package to bring everything into a long format

plot.elos <- melt(elos.males, value.name = "elo", id.vars = "Date")

# plot the ranks as lines and see what happens
ggplot(plot.elos, aes(x = as.Date(Date), y = elo, color = variable)) +
  geom_line(size = 2)
#save as pdf
dev.copy2pdf(file="/Users/gb64/Desktop/rank_data/Issa_eloplot_Alex-script.pdf")
# I would count as a "real" rank change only pant grunts that are incongruent and then afterwards followed by more pant grunts in the new direction. Let's check:
elo.issa$congr <- ""
elo.issa$congr[elo.issa$elo_w_before > elo.issa$elo_l_before] <- "congruent"
elo.issa$congr[elo.issa$elo_w_before < elo.issa$elo_l_before] <- "incongruent"
incongr.issa <- elo.issa[elo.issa$congr == "incongruent", c("Date", "Loser", "Winner")]
colnames(incongr.issa) <- c("Date", "Dominant_before", "Subdominant_before")
incongr.issa$Date <- as.Date(incongr.issa$Date)
incongr.issa <- incongr.issa[!duplicated(incongr.issa), ]
incongr.issa$Group <- "sonso"
incongr.issa$Wins_dominant_before <- 0
incongr.issa$Wins_subdominant_before <- 0
incongr.issa$Wins_dominant_after <- 0
incongr.issa$Wins_subdominant_after <- 0
incongr.issa$Last_win_dominant <- as.Date("1970-01-01")
elo.data$Winner <- as.character(elo.data$Winner)
elo.data$Loser <- as.character(elo.data$Loser)

for (i in 1:nrow(incongr.issa)) {
  incongr.issa$Wins_dominant_before[i] <- nrow(elo.data[elo.data$Winner == incongr.issa$Dominant_before[i] & elo.data$Loser == incongr.issa$Subdominant_before[i] & elo.data$Date < incongr.issa$Date[i], ])
  incongr.issa$Wins_subdominant_before[i] <- nrow(elo.data[elo.data$Loser == incongr.issa$Dominant_before[i] & elo.data$Winner == incongr.issa$Subdominant_before[i] & elo.data$Date < incongr.issa$Date[i], ])
  incongr.issa$Wins_dominant_after[i] <- nrow(elo.data[elo.data$Winner == incongr.issa$Dominant_before[i] & elo.data$Loser == incongr.issa$Subdominant_before[i] & elo.data$Date > incongr.issa$Date[i], ])
  incongr.issa$Wins_subdominant_after[i] <- nrow(elo.data[elo.data$Loser == incongr.issa$Dominant_before[i] & elo.data$Winner == incongr.issa$Subdominant_before[i] & elo.data$Date > incongr.issa$Date[i], ])
  if (incongr.issa$Wins_dominant_before[i] > 0) {
    incongr.issa$Last_win_dominant[i] <- as.Date(max(elo.data$Date[elo.data$Winner == incongr.issa$Dominant_before[i] & elo.data$Loser == incongr.issa$Subdominant_before[i] & elo.data$Date < incongr.issa$Date[i]]))
  }
}

incongr.issa$Last_win_dominant <- as.character(incongr.issa$Last_win_dominant)
incongr.issa$Last_win_dominant[incongr.issa$Last_win_dominant == "1970-01-01"] <- "First Record"
