library(ggplot2)
library(ggthemes)
library(scales)

rmse_gauss <- as.matrix(read.csv("out/rmse_t_gauss.csv"))

rmse_pearson <- as.matrix(read.csv("out/rmse_t_pearson.csv"))

steps <- 1000

dist <- c(rep("Gaussian Distributions", 8*steps),
          rep("Pearson Type IV Distributions", 8*steps))

comp <- rep(c(rep("Position RMSE (m)",   steps),
              rep("Velocity RMSE (m/s)", steps)), 8)

Filter <- rep(c(rep("HOUSE", 2*steps),
                rep("UKF",   2*steps),
                rep("CUT-4", 2*steps),
                rep("CUT-6", 2*steps)), 2) 

rmse_gauss   <- as.matrix(read.csv("out/rmse_t_gauss.csv"))
rmse_pearson <- as.matrix(read.csv("out/rmse_t_pearson.csv"))

t <- rep(rmse_gauss[,1], 16)

RMSE <- c(c(rmse_gauss[,2:9]),
          c(rmse_pearson[,2:9]))

filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6")

data <- data.frame(t, RMSE, Filter, comp, dist)

data$Filter <- factor(data$Filter, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377")
dashes <- c("solid", "dotted", "dotdash", "longdash", "twodash")

ggplot(data, aes(x=t, y=RMSE)) + geom_line(aes(col=Filter, linetype=Filter)) +
    scale_color_manual(values=colors) +
    facet_grid(rows=vars(comp), cols=vars(dist),
        scales="free_y", switch="y", shrink=FALSE) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major = element_line(),
          legend.title=element_blank()) +
    xlab("Time-of-Flight (%)") + ylab(NULL) +
    scale_y_log10(breaks = 10^(-10:10),
                  labels = trans_format("log10", math_format(10^.x)))

ggsave("proj_rmse.eps", width=6, height=4.5, units="in")

