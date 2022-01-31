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
        Bxt[Bxt<=0] <- 0
        BC <- as.matrix(vegdist(Bxt, method="bray"))
        mean_bc[x] <- mean(BC[upper.tri(BC)])
      }
      beta.t[i] <- mean(mean_bc) # record mean temporal BC dissimilarity
      B_t <- as.matrix(fread(paste0(fList[i], "/", bList[1])))
      BC <- as.matrix(vegdist(t(B_t)))
      beta.s[i] <- mean(BC) # record mean spatial BC dissimilarity
    }
    
    dat <- data.frame(inv = as.numeric(stringr::str_extract(fList, "\\d+")),
                      it = 1:length(gamma),
                      alpha.mn = alpha.mn,
                      gamma = gamma,
                      beta.t = beta.t,
                      beta.s = beta.s,
                      distr = rep(distr[a], length(alpha.mn)))
    dat
  }
  write.csv(dat, "~/LVMCM_src/LVMCM/autonomous_turnover_example/autonomous_turnover_example.csv", row.names=F)
} else {
  dat <- read.csv("~/Software/LVMCM_publish_version/LVMCM_src/LVMCM/autonomous_turnover_example/autonomous_turnover_example.csv")
}

datS <- data.frame(inv=rep(dat$inv, 2),
                   it=rep(dat$it, 2),
                   distr=rep(dat$distr,2),
                   S=c(dat$gamma, dat$alpha.mn),
                   level=rep(c("gamma", "alpha"), each=nrow(dat)))

datB <- data.frame(inv=rep(dat$inv, 2),
                   it=rep(dat$it, 2),
                   distr=rep(dat$distr,2),
                   beta=c(dat$beta.t, dat$beta.s),
                   level=rep(c("temp", "spat"), each=nrow(dat)))

### Gen Fig 5
threshold <- min(which(subset(datB, level=="temp" & distr=="discr")$beta > 1e-2))
p5A <- ggplot(subset(datS, distr=="discr"), aes(x=it, y=S, linetype=factor(level))) +
  geom_vline(xintercept = dat$it[threshold], linetype=2, col="grey") +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_linetype_manual(values=c(3,1), labels=c(expression(bar(alpha)[src]), expression(gamma))) +
  labs(x="Time (iterations of assembly model)", y="Species richness", col="Distribution", linetype="") +
  scale_y_continuous(trans="log10") +
  scale_x_continuous(limits=c(0,150))
  
p5B <- ggplot(subset(datB, distr=="discr"), aes(x=it, y=beta, linetype=factor(level))) +
  geom_vline(xintercept = dat$it[threshold], linetype=2, col="grey") +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits=c(0,1)) +
  scale_linetype_manual(values=c(1,3), labels=c(expression(bar(beta)[s]), expression(bar(beta)[t]))) +
  labs(x="Time (iterations of assembly model)", y="Mean BC dissimilarity", col="Distribution", linetype="")

grid.arrange(p5A,p5B,ncol=2)

### Gen Fig S4
# Local interaction matrices at end state
Ad <- as.matrix(read.table("~/LVMCM_src/SimulationData/N=32/discrTraj_experiment/2021-3-1/2021-3-1_discrTraj(10000)1cMat0.mat"))
Ab <- as.matrix(read.table("~/LVMCM_src/SimulationData/N=32/betaDiscrTraj04_experiment/2021-3-1/2021-3-1_betaDiscrTraj04(10000)1cMat0.mat")) 
An <- as.matrix(read.table("~/LVMCM_src/SimulationData/N=32/normDiscrTraj04_experiment/2021-3-1/2021-3-1_normDiscrTraj04(10000)1cMat0.mat"))

datA <- data.frame(A_ij = c(Ad[upper.tri(Ad)],
                            Ad[lower.tri(Ad)],
                            Ab[upper.tri(Ab)],
                            Ab[lower.tri(Ab)],
                            An[upper.tri(An)],
                            An[lower.tri(An)]),
                   distr = c(rep("D", nrow(Ad)*(nrow(Ad)-1)),
                             rep("B", nrow(Ab)*(nrow(Ab)-1)),
                             rep("N", nrow(An)*(nrow(An)-1))))

pS4A <- ggplot(datS, aes(x=it, y=S, col=factor(distr), linetype=factor(level))) +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_linetype_manual(values=c(3,1), labels=c(expression(bar(alpha)[src]), expression(gamma))) +
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels=c("B", "D", "N")) +
  labs(x="Time (iterations of assembly model)", y="Species richness", col="Distribution", linetype="") +
  scale_y_continuous(trans="log10") +
  scale_x_continuous(limits=c(0,150))

pS4B <- ggplot(datB, aes(x=it, y=beta, col=factor(distr), linetype=factor(level))) +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels=c("B", "D", "N")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_linetype_manual(values=c(1,3), labels=c(expression(bar(beta)[s]), expression(bar(beta)[t]))) +
  labs(x="Time (iterations of assembly model)", y="Mean BC dissimilarity", col="Distribution", linetype="")

pS4C <- ggplot(subset(datA, A_ij != 0.0), aes(x=A_ij, fill=distr, group=distr)) +
  geom_histogram(col="black", alpha=0.3, position="identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(fill="Distribution", x=expression(A[ij]), y="")

pdf("~/pCloudDrive/PhD/Written content/AutonomousFluctuations/submission/nature_comms_resub/revisions/si_ma_sensitivity.pdf",
    16,4)
grid.arrange(pS4A,pS4B,pS4C, ncol=3)
dev.off()

nrow(subset(datA, distr== "D" & A_ij == 0.0)) / nrow(subset(datA, distr== "D")) 
nrow(subset(datA, distr== "B" & A_ij == 0.0)) / nrow(subset(datA, distr== "B")) 
nrow(subset(datA, distr== "N" & A_ij == 0.0)) / nrow(subset(datA, distr== "N")) 

