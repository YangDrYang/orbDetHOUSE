library(ggplot2)
library(ggthemes)
library(matrixStats)
library(scales)

rt_gauss <- as.matrix(read.csv("out/run_times_gauss.csv"))[,1:4]
rt_pearson <- as.matrix(read.csv("out/run_times_pearson.csv"))[,1:4]

fil <- c("HOUSE", "UKF", "CUT-4", "CUT-6")

filter <- rep(fil, 2)
Distributions <- c(rep("Gaussian Distributions", 4),
                   rep("Pearson Type IV Distributions", 4))

rt_mean <- 1000 * c(colMeans(rt_gauss), colMeans(rt_pearson))
rt_sd <- 1000 * c(colSds(rt_gauss), colSds(rt_pearson))

data <- data.frame(filter, Distributions, rt_mean, rt_sd)

data$filter <- factor(data$filter, levels=fil)

colors <- c("#0077bb","#cc3311")

ggplot(data, aes(x=filter, y=rt_mean, fill=Distributions)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw(base_size = 9) +
    theme(legend.position="bottom",
          legend.direction = "vertical",
          legend.title=element_blank()) +
    xlab(NULL) + ylab("Average Run Time (ms)") + 
    scale_fill_manual(values=colors) + 
    geom_errorbar(aes(ymin=rt_mean-rt_sd, ymax=rt_mean+rt_sd),
                  width=.2, position=position_dodge(.9)) +
    scale_y_log10(breaks = 10^(-10:10),
                  labels = trans_format("log10", math_format(10^.x)))

ggsave("proj_runtime.eps", width=3.25, height=3.5, units="in")

