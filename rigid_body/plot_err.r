library(ggplot2)
library(ggthemes)
library(scales)

fil <- c("house", "ukf", "cut4", "cut6", "cut8")
filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6", "CUT-8")

dis <- c("gauss", "pearson")
dist <- c("Gaussian Distributions", "Pearson Type IV Distributions")

steps <- 6001
trials <- 100

Err <- c()
Filter <- c()
Dist <- c()

for (i in 1:2) {
    Dist <- c(Dist, rep(dist[i], 5*steps*trials))
    for (j in 1:5) {
        Filter <- c(Filter, rep(filter[j], steps*trials))
        file  <- paste("out/", fil[j], "_err_", dis[i], ".csv", sep="")
        table <- as.matrix(read.csv(file))
        err <- sqrt(rowSums(table[,1:3]^2))
        Err <- c(Err, err)
    }
}

data <- data.frame(Err, Filter, Dist)

data$Filter <- factor(data$Filter, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377")
dashes <- c("solid", "dotted", "dotdash", "longdash", "twodash")

ggplot(data, aes(x=Filter, y=Err, fill=Filter)) +
    scale_fill_manual(values=colors) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1, outlier.shape=NA, fill="white") +
    theme_bw(base_size = 9) +
    theme(panel.grid.major = element_line(),
          legend.position = "none") +
    facet_grid(cols=vars(Dist),
        scales="free_y", switch="y", shrink=FALSE) +
    xlab(NULL) + ylab("Angular Velocity Error (rad/s)")

ggsave("rb_err.eps", width=6, height=3, units="in")


