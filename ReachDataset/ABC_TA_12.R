library(invgamma)
library(ggplot2)

# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)

# R & STAN are friends!
library(rstan)
library(coda)

library("abc")
require("abc.data")
library("EasyABC")


###############################################################################################################
###############################PROVIAMO ABC IN UN MODEL SELECTION FRAMEWORK####################################
####################################CODICE TA_12####################################

reach <- read.table("REACH_data.txt", header=T)
names(reach)
head(reach)
reach$gender <- reach$gender - 1

sub.idx = 1:12  
num.data = 681

X <- as.matrix(reach[1:num.data,sub.idx])
Y <- as.vector(reach[1:num.data,13])

N <- dim(X)[1]
p <- dim(X)[2]


## vogliamo costruire l'ABC. Il processo deve essere circa:
## 1. Per p volte:
##    i. Simulare theta_j
##    ii. Simulare gamma_j
##    iii. Simulare beta_j|gamma_j con approccio Spike&Slab
## 2. Simulare dalla likelihood N Y_i|X_i, beta (Per il momento mantengo le X_obs, ma non penso si possa fare altrimenti)
## 3. Studiare la distanza tra Y_new e Y_obs e accettare con una certa prob. (vedi script Francesca) le beta corrispondenti(e quindi le gamma)
## 4. Studiare le gamma accettate per costruire le posterior





gamma <- c(rep(0,p))
beta <- c(rep(0,p))
y_new <- c(rep(0,N))
n_trials <-250000
outcome_g <- NULL#raccoglierà tutti i gamma simulati
outcome_b <- NULL #raccoglierà tutti i beta simulati
diff_vec <- NULL #raccoglie le differenze in valore assoluto tra Y e tutti gli y_new simulati
diff_vec2 <- NULL #come sopra. ma solo con gli y_new accettati

##Spike&Slab non rilassato


for (k in 1:n_trials){
  
  set.seed(k)
  theta <- runif(12) #campiono 12 theta
  #campiono le 12 gamma (ho fatto delle prove e ho verificato che ciascuno dei 12 campionamenti è eseguito rispetto alla theta[i] corrispondente)
  gamma <- as.numeric(rbernoulli(12,theta)) 
  beta <- gamma*rnorm(12,0,1) #calcolo le beta
  
  beta0_temp <- rnorm(1,-0.6,1.55) #campiono un beta0 (qui ho svolto una prima simulazione e ho estratto i valori di media e varianza delle b0 accettate)
  beta0 <- rbind(rep(beta0_temp,N)) #lo mantengo costante per ogni elemento della iterazione corrente
  mi <- t(beta0) + X%*%beta
  z <- pnorm(mi)
  y_new2 <- as.numeric(rbernoulli(N,z))
  
  
  #questo processo va ripetuto per un certo numero n_trials
  
  #a questo punto bisogna introdurre il kernel e il processo di accettazione per decidere se tenere o no il vettore di
  #gamma
  #una volta estratti i gamma selezionati si possono usare le tecniche di goodness of fit in TA_12
  #finisci questo, prova lo spike&slabs rilassato e infine decidi se provare a farlo con le summary (molto difficile)
  #alla fine i campioni di Y sono meno di 700. Forse non è utile usare le summary statistics
  
  diff = y_new2-Y
  diff_vec <- rbind(diff_vec, norm(diff, type = "2"))
  if(norm(diff, type = "2") <=17 ){
    
    beta_tot <-c(beta0_temp,beta)
      
    outcome_b <- rbind(outcome_b,beta_tot)
    outcome_g <- rbind(outcome_g,gamma)
    diff_vec2 <- rbind(diff_vec2, norm(diff, type = "2"))
      
    }
    #qui ho pensato: se il parametro lo campiono non dalla proposal ma dalla sua vera prior, allora se viene rispettata
    #la condizione della tolleranza (che dovrà essere trasformata in kernel più avanti) automaticamente accettiamo la
    #condizione senza passare per il passaggio "probabilistico" di cui vi ho chiesto
    
  
    
  }

acc_rate <- length(diff_vec2)/length(diff_vec)
acc_rate

outcome_g_mean <- apply(outcome_g, 2, "mean")
outcome_b_mean <- apply(outcome_b, 2, "mean")

p3 <- data.frame(value = outcome_g_mean, var = colnames(X)) %>%
  ggplot(aes(y = value, x = var, fill = var)) + 
  geom_bar(stat="identity") + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lwd = 1.1) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("posterior inclusion probabilities ABC") + 
  xlab("")
p3

p4 <- data.frame(value = outcome_b_mean, var = c("b0", colnames(X))) %>%
  ggplot(aes(y = value, x = var, fill = var)) + 
  geom_bar(stat="identity") + 
  #geom_hline(mapping = aes(yintercept = .5), col = 2, lwd = 1.1) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("posterior beta ABC") + 
  xlab("")
p4

par(mfrow=c(2,2))
p1
p2
p3
p4


#####ignora da qui in poi
#outcome_g
#outcome_b
#diff_vec
#axis_x <- seq(0,40,1)
#plot(diff_vec2, ylim = c(16,20))
#abs_diff <- 2*norm(diff, type = "2")
kern <- density(abs_diff,0.5,kernel = c("gaussian"))


h = 300

kern_2 <- (1/h)*(1/(sqrt(2*pi)))*exp(-(1/2)*(1/h)*abs_diff^2)
#kern_2
#i dati di diff_vec5000 hanno una distribuzione normale
#i quantili di diff_vec5000 sono:
#     25%: 18.19341
#     50%: 18.46619
#     75%: 18.68154


#idea: si potrebbe provare un test d'ipotesi accoppiato per vedere se la differenza di diff_vec tra varie simulazioni è significativa
#oppure se è poco signitifcativa, il che ci permetterebbe di generalizzare gli IC di una singola simulazione a tutte le simulazioni


#i risultati sono buoni. L'unica differenza con JAGS è data dal fatto che bcp_hands in ABC non supera la soglia di 0.5






