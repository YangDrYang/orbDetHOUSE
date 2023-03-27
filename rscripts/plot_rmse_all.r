library(ggplot2)
library(ggthemes)
library(scales)
library(svglite)

# rmse_gauss <- as.matrix(read.csv("out/rmse_t_gauss.csv"))

# rmse_pearson <- as.matrix(read.csv("out/rmse_t_pearson.csv"))


numFilters <- 4
steps <- 4321

dist <- c(rep("Gaussian Distributions", 2*numFilters*steps))

comp <- rep(c(rep("Position RMSE (m)",   steps),
              rep("Velocity RMSE (m/s)", steps)), numFilters)

# 2 data points (pos and vel)
Filter <- rep(c(rep("HOUSE", 2*steps),
                rep("UKF",   2*steps),
                rep("CUT-4", 2*steps),
                rep("CUT-6", 2*steps)), 2)

outDir <- "out"
rmseFile <- paste(outDir, "/rmse_t.csv", sep="")
rmse_gauss   <- as.matrix(read.csv(rmseFile))

t <- rep(rmse_gauss[,1], 2*(2*numFilters))

RMSE <- c(c(rmse_gauss  [,2 : (2*numFilters + 1)]))

filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6")

data <- data.frame(t, RMSE, Filter, comp, dist)

data$Filter <- factor(data$Filter, levels=filter)

# Open a PDF device
pdf("plots/proj_rmse_all.pdf")

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
                  labels = trans_format("log10", math_format(10^.x)))  #+     coord_cartesian(ylim = c(10^-2, 10^5))
# Close the PDF device
dev.off()
ggsave("plots/proj_rmse_all.svg", width=6, height=4.5, units="in")

