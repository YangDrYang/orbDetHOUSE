mkdir -p out

rm -f out/*

./rigid_body.exe

Rscript plot_err.r
Rscript plot_rmse.r
Rscript plot_runtime.r
