library(matrixStats)

process_rmse <- function(trials, filter) {

    print("Post-processing RMSE")

    tru_file <- paste("out/trajectory_truth.csv", sep="")
    tru <- as.matrix(read.csv(tru_file))[,2:7]

    steps <- nrow(tru)

    r_se <- matrix(, nrow=steps, ncol=trials)
    v_se <- matrix(, nrow=steps, ncol=trials)
    
    for (j in 1:trials) {

        print(paste("   Trial ", j), quote=FALSE)
        file <- paste("out/", filter, "_", j, ".csv", sep="")    
        est <- as.matrix(read.csv(file))[,2:7]

        r_se[,j] <- rowSums((est[,1:3] - tru[,1:3])^2)
        v_se[,j] <- rowSums((est[,4:6] - tru[,4:6])^2)

        # print(est[,1:3] - tru[,1:3])
        # print(r_se[,j])

    }
    # print(r_se)

    # --- RMSE
    rmse_pos <- c(sqrt(mean(r_se,)))    
    rmse_vel <- c(sqrt(mean(v_se)))

    print(rmse_pos)
    print(rmse_vel)

    rmse <- data.frame(filterno = 1:2, filter, rmse_pos, rmse_vel)

    file <- paste("out/", filter, "_rmse",  ".csv", sep="")
    write.csv(rmse, file, row.names=FALSE, quote=FALSE)

    # --- RMSE vs. time

    time_steps <- seq(0, steps-1, length.out=steps)

    rmse_t <- data.frame(time_steps, 
          pos = sqrt(rowMeans(r_se)),
          vel = sqrt(rowMeans(v_se)))

    file <- paste("out/", filter, "_rmse_t", ".csv", sep="")
    write.csv(rmse_t, file, row.names=FALSE, quote=FALSE)

    # --- Convergence

    convmse <- function(se) {
        cmse <- vector(, length = trials)
        for (j in 1:trials) {
            cmse[j] <- mean(se[,1:j])        
        }
        return(cmse)
    }

    convsse <- function(se) {
        csse <- vector(, length = trials)
        for (j in 1:trials) {
            csse[j] <- sd(se[,1:j])
        }
        return(csse)
    }
    conv <- data.frame(ntrials = 1:trials,
                r_mse = convmse(r_se), r_sse = convsse(r_se),
                v_mse = convmse(v_se), v_sse = convsse(v_se))

    file <- paste("out/", filter, "_conv", ".csv", sep="")
    write.csv(conv, file, row.names = FALSE, quote=FALSE)
}

process_box <- function(trials, filter) {

    print("Post-processing for boxplot")

    print("Post-processing")

    trialno <- 1:trials

    f = 10
    tmin = 1
    kstart = f * tmin + 1

    r_err <- c()
    v_err <- c()

    tru_file <- paste("out/trajectory_truth.csv", sep="")
    tru <- as.matrix(read.csv(tru_file))[,2:7]
    for (j in 1:trials) {

        k <- trialno[j]

        print(paste("   Trial ", k), quote=FALSE)

        file <- paste("out/", filter, "_", k, ".csv", sep="")

        kstop <- length(tru[,1])

        tru <- tru[kstart:kstop,]

        est <- as.matrix(read.csv(file))[kstart:kstop,2:7]

        r_err <- c(sqrt(rowSums((est[,1:3] - tru[,1:3])^2)), r_err)

        v_err  <- c(sqrt(rowSums((est[,4:6] - tru[,4:6])^2)), v_err)

    }

    r_err <- data.frame(r_err) # , CUT4=cut4_r_err , CUT6=cut6_r_err)
    v_err <- data.frame(v_err) # , CUT4=cut4_v_err , CUT6=cut6_v_err)

    pos_file <- paste("out/", filter, "_pos_err_", ".csv", sep="")
    vel_file <- paste("out/", filter, "_vel_err_", ".csv", sep="") 

    write.csv(r_err, pos_file, row.names = FALSE, quote=FALSE)
    write.csv(v_err, vel_file, row.names = FALSE, quote=FALSE)

}

process_rmse(100, "ukf")
process_box(100, "ukf")
