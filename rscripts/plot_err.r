library(ggplot2)
library(scales)
library(svglite)

# four filters
filterNum <- 4

pos_err <- as.matrix(read.csv("out/pos_err_.csv"))
vel_err <- as.matrix(read.csv("out/vel_err_.csv"))

Ng <- nrow(pos_err)

err <- c(pos_err, vel_err)

dst <- c(rep("Gaussian Distributions", filterNum*Ng))

fil <- c(c(rep("HOUSE", Ng),
               rep("UKF",   Ng),
               rep("CUT-4", Ng),
               rep("CUT-6", Ng)))

cmp <- c(rep("Position Error (m)",   filterNum*Ng),
         rep("Velocity Error (m/s)", filterNum*Ng))

data <- data.frame(err, dst, fil, cmp)

filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6")
data$fil <- factor(data$fil, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377")

# Open a PDF device
pdf("plots/proj_err.pdf")

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
# Close the PDF device
dev.off()
ggsave("plots/proj_err.svg", width=6, height=6, units="in")

