require(ggplot2)
require(gridExtra)

ANALYSE_DATA=F # output of analysis saved in file ~/LVMCM_src/SimulationData/N=32/autonomous_turnover_example.csv

if (ANALYSE_DATA) {
  require(doParallel)
  require(foreach)

  # If assemblies run locally, update the directory names below to match corresponding date(s)
  dirList <- c("~/LVMCM_src/SimulationData/N=32/discrTraj_experiment/2021-3-1/1bMat0/",
               "~/LVMCM_src/SimulationData/N=32/betaDiscrTraj04_experiment/2021-3-1/1bMat0/",
               "~/LVMCM_src/SimulationData/N=32/normDiscrTraj04_experiment/2021-3-1/1bMat0/")
  
  distr <- c("discr",
             "betaDiscr0.4",
             "normDiscr0.4")
  
  cores <- 3
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  dat <- foreach (a=1:length(dirList), .combine=rbind) %dopar% {
    require(data.table)
    require(vegan)
    
    setwd(dirList[a])
    fList <- dir()
    inv <- as.numeric(stringr::str_extract(fList, "\\d+"))
    fList <- fList[order(inv)]
    
    alpha.mn <- c()
    gamma <- c()
    beta.t <- c()
    beta.s <- c()
    
    for (i in 1:length(fList)) { # loop through invasions
      
      bList <- dir(fList[i])
      B_src <- as.matrix(fread(paste0(fList[i], "/", bList[grep("src", bList)])))
      bList <- bList[-grep("src", bList)]
      tt <- stringr::str_extract(bList, "\\d+")
      bList <- bList[order(tt)]
      
      B <- c()
      for (t in 1:length(bList)) { # create (S.NxT) trajectories object
        B <- rbind(B, as.vector(as.matrix(fread(paste0(fList[i], "/", bList[t])))))
      }
      B[B<=0] <- 0
      
      N <- ncol(B_src)
      S <- nrow(B_src)
      alpha.mn[i] <- mean(apply(B_src, MAR=2, FUN=function(x) length(which(x==1)))) # record mean local source richness
      gamma[i] <- nrow(B_src) # record regional richness
      
      mean_bc <- c()
      for (x in 1:N) { # loop through sites
        Bxt <- B[,c(1:S)+((x-1)*S)]
        Bxt[Bxt<0] <- 0
        BC <- as.matrix(vegdist(Bxt, method="bray"))
        mean_bc[x] <- mean(BC[upper.tri(BC)])
      }
      beta.t[i] <- mean(mean_bc) # record mean temporal BC dissimilarity
      B_t <- as.matrix(fread(paste0(fList[i], "/", bList[1])))
      BC <- as.matrix(vegdist(t(B_t)))
      beta.s[i] <- mean(BC) # record mean spatial BC dissimilarity
    }
    dat
  }
  write.csv(dat, "~/LVMCM_src/LVMCM/autonomous_turnover_example/autonomous_turnover_example.csv", row.names=F)
} else {
  dat <- read.csv("~/LVMCM_src/LVMCM/autonomous_turnover_example/autonomous_turnover_example.csv")
}

datS <- data.frame(inv=rep(dat$inv, 2),
                   distr=rep(dat$distr,2),
                   S=c(dat$gamma, dat$alpha.mn),
                   level=rep(c("gamma", "alpha"), each=nrow(dat)))

datB <- data.frame(inv=rep(dat$inv, 2),
                   distr=rep(dat$distr,2),
                   beta=c(dat$beta.t, dat$beta.s),
                   level=rep(c("temp", "spat"), each=nrow(dat)))

p1 <- ggplot(subset(datS), aes(x=inv, y=S, col=factor(distr), linetype=factor(level))) +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_linetype_manual(values=c(3,1), labels=c(expression(bar(alpha)[src]), expression(gamma))) +
  labs(x='0.1S + 1 Invasions', y="Richness", col="Distribution", linetype="") +
  scale_y_continuous(trans="log10")

p2 <- ggplot(datB, aes(x=inv, y=beta, col=factor(distr), linetype=factor(level))) +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits=c(0,1)) +
  scale_linetype_manual(values=c(1,3), labels=c(expression(beta[s]), expression(beta[t]))) +
  labs(x='0.1S + 1 Invasions', y="Mean BC dissimilarity", col="Distribution", linetype="")

grid.arrange(p1,p2)
