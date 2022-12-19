library(ggplot2)
library(scales)
library(svglite)

pos_err_gauss <- as.matrix(read.csv("out/pos_err_gauss.csv"))
vel_err_gauss <- as.matrix(read.csv("out/vel_err_gauss.csv"))

pos_err_pearson <- as.matrix(read.csv("out/pos_err_pearson.csv"))
vel_err_pearson <- as.matrix(read.csv("out/vel_err_pearson.csv"))

Ng <- nrow(pos_err_gauss)
Np <- nrow(pos_err_pearson)

err <- c(pos_err_gauss, vel_err_gauss, pos_err_pearson, vel_err_pearson)

dst <- c(rep("Gaussian Distributions", 6*Ng),
         rep("Pearson Type IV Distributions", 6*Np))

fil <- c(rep(c(rep("HOUSE", Ng),
               rep("UKF",   Ng),
               rep("CUT-4", Ng)), 2),
         rep(c(rep("HOUSE", Np),
               rep("UKF",   Np),
               rep("CUT-4", Np)), 2))

cmp <- c(rep("Position Error (m)",   3*Ng),
         rep("Velocity Error (m/s)", 3*Np))

data <- data.frame(err, dst, fil, cmp)

filter <- c("HOUSE", "UKF", "CUT-4")#, "CUT-6")
data$fil <- factor(data$fil, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377")

ggplot(data, aes(x=fil, y=err, fill=fil)) +
    scale_fill_manual(values=colors) +
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.1, outlier.shape=NA, fill="white") +
    theme_bw(base_size = 9) + 
    theme(panel.grid.major = element_line(),
          legend.position = "none") +
    facet_grid(rows=vars(cmp), cols=vars(dst),
        scales="free_y", switch="y", shrink=FALSE) +
    xlab(NULL) + ylab(NULL) +
    scale_y_log10(breaks = 10^(-10:10),
                  labels = trans_format("log10", math_format(10^.x)))

ggsave("proj_err.svg", width=6, height=6, units="in")

