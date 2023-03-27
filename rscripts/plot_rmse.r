library(ggplot2)
library(ggthemes)
library(scales)
library(svglite)


numFilters <- 1

plot_rmse <- function(steps, filter) {
    dist <- c(rep("Gaussian Distributions", 2*numFilters*steps))

    comp <- rep(c(rep("Position RMSE (m)",   steps),
                rep("Velocity RMSE (m/s)", steps)), 2*numFilters)

    # 2 data points (pos and vel)
    Filter <- rep(c(rep(filter, 2*steps)), 2)

    file <- paste("out/", filter, "_rmse_t.csv", sep="")
    rmse_gauss   <- as.matrix(read.csv(file))

    t <- rep(rmse_gauss[,1], 2*(2*numFilters))

    RMSE <- c(c(rmse_gauss  [,2 : (2*numFilters + 1)]))

    data <- data.frame(t, RMSE, Filter, comp, dist)

    data$Filter <- factor(data$Filter, levels=filter)

    # Open a PDF device
    file <- paste("plots/", filter, "_rmse.pdf", sep="")
    pdf(file)

    colors <- c("#cc3311", "#0077bb")
    dashes <- c("solid", "dotted")

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
    file <- paste("plots/", filter, "_rmse.svg", sep="")
    ggsave(file, width=6, height=4.5, units="in")
}

steps <- 4321
filter <- "ukf"
plot_rmse(steps, filter)
