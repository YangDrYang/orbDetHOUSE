mkdir -p out

rm -f out/*

./projectile.exe

Rscript processing.r
Rscript plot_err.r
Rscript plot_rmse.r
Rscript plot_runtime.r

