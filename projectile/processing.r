library(matrixStats)

process_rmse <- function(dist) {

    print("Post-processing RMSE", quote=FALSE)
    print(dist, quote=FALSE)

    trials = 100
    steps = 1000

    house_r_se <- matrix(, nrow=steps, ncol=trials)
    ukf_r_se   <- matrix(, nrow=steps, ncol=trials)
    cut4_r_se  <- matrix(, nrow=steps, ncol=trials)
    cut6_r_se  <- matrix(, nrow=steps, ncol=trials)

    house_v_se <- matrix(, nrow=steps, ncol=trials)
    ukf_v_se   <- matrix(, nrow=steps, ncol=trials)
    cut4_v_se  <- matrix(, nrow=steps, ncol=trials)
    cut6_v_se  <- matrix(, nrow=steps, ncol=trials)

    state_interp <- function(filename) {
        t <- as.matrix(read.csv(filename))[,1]
        x <- as.matrix(read.csv(filename))[,2:7]
        xi <- matrix(, nrow=steps, 6)
        for (i in 1:6) {
            xi[,i] <- approx(t, x[,i], n=steps, method="linear")$y
        }
        return (xi)
    }

    for (j in 1:trials) {

        print(paste("   Trial ", j), quote=FALSE)

        tru_file <- paste("out/tru_", dist, "_",  j, ".csv", sep="")

        house_file <- paste("out/house_est_", dist, "_", j, ".csv", sep="")
        ukf_file   <- paste("out/ukf_est_",   dist, "_", j, ".csv", sep="")
        cut4_file  <- paste("out/cut4_est_",  dist, "_", j, ".csv", sep="")
        cut6_file  <- paste("out/cut6_est_",  dist, "_", j, ".csv", sep="")

        tru <- state_interp(tru_file)

        house_est <- state_interp(house_file)
        ukf_est   <- state_interp(ukf_file)
        cut4_est  <- state_interp(cut4_file)
        cut6_est  <- state_interp(cut6_file)

        house_r_se[,j] <- rowSums((house_est[,1:3] - tru[,1:3])^2)
        ukf_r_se  [,j] <- rowSums((ukf_est  [,1:3] - tru[,1:3])^2)
        cut4_r_se [,j] <- rowSums((cut4_est [,1:3] - tru[,1:3])^2)
        cut6_r_se [,j] <- rowSums((cut6_est [,1:3] - tru[,1:3])^2)

        house_v_se[,j] <- rowSums((house_est[,4:6] - tru[,4:6])^2)
        ukf_v_se  [,j] <- rowSums((ukf_est  [,4:6] - tru[,4:6])^2)
        cut4_v_se [,j] <- rowSums((cut4_est [,4:6] - tru[,4:6])^2)
        cut6_v_se [,j] <- rowSums((cut6_est [,4:6] - tru[,4:6])^2)

    }

    # --- RMSE

    rmse_pos <- c(sqrt(mean(house_r_se)),
                  sqrt(mean(ukf_r_se)),
                  sqrt(mean(cut4_r_se)),
                  sqrt(mean(cut6_r_se)))

    rmse_vel <- c(sqrt(mean(house_v_se)),
                  sqrt(mean(ukf_v_se)),
                  sqrt(mean(cut4_v_se)),
                  sqrt(mean(cut6_v_se)))

    filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6")

    rmse <- data.frame(filterno = 1:4, filter, rmse_pos, rmse_vel)

    file <- paste("out/rmse_", dist, ".csv", sep="")
    write.csv(rmse, file, row.names=FALSE, quote=FALSE)

    # --- RMSE vs. time

    time_percent <- seq(0, 100, length.out=steps)

    rmse_t <- data.frame(time_percent, 
          house_pos = sqrt(rowMeans(house_r_se)),
          house_vel = sqrt(rowMeans(house_v_se)),
          ukf_pos   = sqrt(rowMeans(ukf_r_se)),
          ukf_vel   = sqrt(rowMeans(ukf_v_se)),
          cut4_pos  = sqrt(rowMeans(cut4_r_se)),
          cut4_vel  = sqrt(rowMeans(cut4_v_se)),
          cut6_pos  = sqrt(rowMeans(cut6_r_se)),
          cut6_vel  = sqrt(rowMeans(cut6_v_se)))

    file <- paste("out/rmse_t_", dist, ".csv", sep="")
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
                house_r_mse = convmse(house_r_se), house_r_sse = convsse(house_r_se),
                house_v_mse = convmse(house_v_se), house_v_sse = convsse(house_v_se),
                ukf_r_mse = convmse(ukf_r_se), ukf_r_sse = convsse(ukf_r_se),
                ukf_v_mse = convmse(ukf_v_se), ukf_v_sse = convsse(ukf_v_se),
                cut4_r_mse = convmse(cut4_r_se), cut4_r_sse = convsse(cut4_r_se),
                cut4_v_mse = convmse(cut4_v_se), cut4_v_sse = convsse(cut4_v_se),
                cut6_r_mse = convmse(cut6_r_se), cut6_r_sse = convsse(cut6_r_se),
                cut6_v_mse = convmse(cut6_v_se), cut6_v_sse = convsse(cut6_v_se))

    file <- paste("out/conv_", dist, ".csv", sep="")
    write.csv(conv, file, row.names = FALSE, quote=FALSE)

}

process_box <- function(dist) {

    print("Post-processing for boxplot", quote=FALSE)
    print(dist, quote=FALSE)

    print("Post-processing", quote=FALSE)

    trialno <- 1:100

    trials = length(trialno)

    f = 10
    tmin = 1
    kstart = f * tmin + 1

    house_r_err <- c()
    ukf_r_err   <- c()
    cut4_r_err  <- c()
    cut6_r_err  <- c()

    house_v_err <- c()
    ukf_v_err   <- c()
    cut4_v_err  <- c()
    cut6_v_err  <- c()

    for (j in 1:trials) {

        k <- trialno[j]

        print(paste("   Trial ", k), quote=FALSE)

        tru_file <- paste("out/tru_", dist, "_",  j, ".csv", sep="")

        house_file <- paste("out/house_est_", dist, "_", k, ".csv", sep="")
        ukf_file   <- paste("out/ukf_est_",   dist, "_", k, ".csv", sep="")
        cut4_file  <- paste("out/cut4_est_",  dist, "_", k, ".csv", sep="")
        cut6_file  <- paste("out/cut6_est_",  dist, "_", k, ".csv", sep="")

        tru <- as.matrix(read.csv(tru_file))[,2:7]

        kstop <- length(tru[,1])

        tru <- tru[kstart:kstop,]

        house_est <- as.matrix(read.csv(house_file))[kstart:kstop,2:7]
        ukf_est   <- as.matrix(read.csv(ukf_file))  [kstart:kstop,2:7]
        cut4_est  <- as.matrix(read.csv(cut4_file)) [kstart:kstop,2:7]
        cut6_est  <- as.matrix(read.csv(cut6_file)) [kstart:kstop,2:7]

        house_r_err <- c(sqrt(rowSums((house_est[,1:3] - tru[,1:3])^2)), house_r_err)
        ukf_r_err   <- c(sqrt(rowSums((ukf_est  [,1:3] - tru[,1:3])^2)), ukf_r_err)
        cut4_r_err  <- c(sqrt(rowSums((cut4_est [,1:3] - tru[,1:3])^2)), cut4_r_err)
        cut6_r_err  <- c(sqrt(rowSums((cut6_est [,1:3] - tru[,1:3])^2)), cut6_r_err)

        house_v_err  <- c(sqrt(rowSums((house_est[,4:6] - tru[,4:6])^2)), house_v_err)
        ukf_v_err    <- c(sqrt(rowSums((ukf_est  [,4:6] - tru[,4:6])^2)), ukf_v_err)
        cut4_v_err   <- c(sqrt(rowSums((cut4_est [,4:6] - tru[,4:6])^2)), cut4_v_err)
        cut6_v_err   <- c(sqrt(rowSums((cut6_est [,4:6] - tru[,4:6])^2)), cut6_v_err)

    }

    names <- c("HOUSE", "UKF", "CUT-4", "CUT-6")

    r_err <- data.frame(HOUSE=house_r_err, UKF=ukf_r_err, CUT4=cut4_r_err, CUT6=cut6_r_err)
    v_err <- data.frame(HOUSE=house_v_err, UKF=ukf_v_err, CUT4=cut4_v_err, CUT6=cut6_v_err)

    pos_file <- paste("out/pos_err_", dist, ".csv", sep="")
    vel_file <- paste("out/vel_err_", dist, ".csv", sep="") 

    write.csv(r_err, pos_file, row.names = FALSE, quote=FALSE)
    write.csv(v_err, vel_file, row.names = FALSE, quote=FALSE)

}

process_rmse("gauss")
process_rmse("pearson")

process_box("gauss")
process_box("pearson")
