library(ggplot2)
library(ggthemes)
library(scales)

fil <- c("house", "ukf", "cut4", "cut6", "cut8")
filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6", "CUT-8")

dis <- c("gauss", "pearson")
dist <- c("Gaussian Distributions", "Pearson Type IV Distributions")

steps <- 6001

t <- c()
RMSE <- c()
Filter <- c()
Dist <- c()

for (i in 1:2) {
    Dist <- c(Dist, rep(dist[i], 5*steps))
    for (j in 1:5) {
        Filter <- c(Filter, rep(filter[j], steps))
        file  <- paste("out/", fil[j], "_rmse_", dis[i], ".csv", sep="")
        table <- as.matrix(read.csv(file))
        t <- c(t, table[,1])
        rmse <- sqrt(rowSums(table[,2:4]^2))
        RMSE <- c(RMSE, rmse)  
    }
}

data <- data.frame(t, RMSE, Filter, Dist)

data$Filter <- factor(data$Filter, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377")
dashes <- c("solid", "dotted", "dotdash", "longdash", "twodash")

ggplot(data, aes(x=t, y=RMSE)) + geom_line(aes(col=Filter, linetype=Filter)) +
    scale_color_manual(values=colors) +
    facet_grid(cols=vars(Dist),
        scales="free_y", switch="y", shrink=FALSE) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major = element_line(),
          legend.title=element_blank()) +
    xlab("Time (s)") + ylab("Angular Velocity RMSE (rad/s)") + ylim(c(0,0.04))
#    scale_y_log10(breaks = 10^(-10:10),
#                  labels = trans_format("log10", math_format(10^.x)))

ggsave("rb_rmse.eps", width=6, height=3, units="in")

