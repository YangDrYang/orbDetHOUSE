# OrbDetHOUSE: The Higher-Order Unscented Estimator for Orbit Determination

Dependencies:

- HOUSE by Z. Stojanovski and D. Savransky: https://github.com/SIOSlab/HOUSE.git
- CUTpoints by N. Adurthi, P. Singla, and T. Singh: https://github.com/nadurthi/CUTpoints using the Matlab function cut_sigma_points.m to generate .csv files in the "CUT" directory
  - the .csv files have been generated, which include two large files (>100 MB). They can potentially result in some issues with Git. 
- SOFA: http://www.iausofa.org/2021_0512.html
- nrlmsise-00: https://github.com/magnific0/nrlmsise-00.git
- jpl_eph: https://github.com/Bill-Gray/jpl_eph.git
- Eigen: https://eigen.tuxfamily.org/
- Boost Installation: https://www.boost.org/doc/libs/1_69_0/more/getting_started/unix-variants.html#easy-build-and-install
- Yaml-cpp: https://github.com/jbeder/yaml-cpp

Compiling:

- compilation is handled by the makefile, run make to rebuild (make must be installed) or clean main functions in the "scripts" directory.
- examples: 
  - make FILENAME=convertECEF2ECI and make FILENAME=convertECEF2ECI clean
  - make FILENAME=testOrbDetCCData and make FILENAME=testOrbDetCCData clean
- the binary object files are stored in the "bin/bin_sub" directory

Running:

- the repo includes an executable compiled for MacBook Intel
- parameters can be set using the yaml files in the yamls directory
- by default the program will read the "config.yaml" file, however
  you can provide an argument to read a different yaml file
- example: 
  - bin/scripts/convertECEF2ECI
  - bin/scripts/testOrbDetCCData yamls/config_ccdata_mee.yml
- output files will be saved into the "out/out_sub" directory

Analysis and plot:

- call python functions in the "pyscripts" directory after install relevant python packages
- analysis results will be saved in the "out/out_sub" directory and plots will be saved in the "plots" directory

Notes:

- be careful with the directory of Eigen, Boost and Yaml-cpp, better to install them using CMake inside their file folders
- adjust the directories in Makefile accordingly to include these libraries
- the author is working on a version with full functionalities

Reference:

- please cite my paper when you use the codes: Yang Yang, Square-Root Higher-Order Unscented Estimators for Robust Orbit Determination, IEEE Transactions on Aerospace and Electronic Systems, DOI: 10.1109/TAES.2024.3423851 (Early Access).
