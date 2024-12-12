# OrbDetHOUSE: The Higher-Order Unscented Estimator for Orbit Determination

Dependencies:

- HOUSE by Z. Stojanovski and D. Savransky: https://github.com/SIOSlab/HOUSE.git
- CUTpoints by N. Adurthi, P. Singla, and T. Singh: https://github.com/nadurthi/CUTpoints using the Matlab function cut_sigma_points.m to generate .csv files in the "CUT" directory
  - the .csv files have been generated, which include two large files (>100 MB). They can potentially result in some issues with Git.
- SOFA: http://www.iausofa.org/2021_0512.html
- nrlmsise-00: https://github.com/magnific0/nrlmsise-00.git
- jpl_eph: https://github.com/Bill-Gray/jpl_eph.git
- Eigen: https://eigen.tuxfamily.org/
- Yaml-cpp: https://github.com/jbeder/yaml-cpp

Compiling:

- compilation is handled by the makefile, run make to rebuild (make must be installed) or clean main functions in the "scripts" directory.
- examples:
  - make FILENAME=convertECEF2ECI and make FILENAME=convertECEF2ECI clean
  - make FILENAME=testOrbProp and make FILENAME=testOrbProp clean
  - make FILENAME=testOrbDetCCData and make FILENAME=testOrbDetCCData clean
- the binary object files are stored in the "bin/bin_sub" directory

Running:

- the repo includes an executable compiled for MacBook Intel Core i7
- parameters can be set using the yaml files in the yamls directory
- by default the program will read the "config.yaml" file, however
  you can provide an argument to read a different yaml file
- examples of main calculations:
  - bin/scripts/convertECEF2ECI
  - bin/scripts/testOrbProp yamls/config_orb.yml
  - bin/scripts/testOrbDetCCData yamls/config_ccdata_mee.yml
- output files will be saved into the "out/out_sub" directory

Analyses and plots:

- call python functions in the "pyscripts" directory after install relevant python packages
- analysis results will be saved in the "out/out_sub" directory and plots will be saved in the "plots" directory
- examples of results analyses and visualisations
  - python3 pyscripts/analyseOrbitProp.py

Notes:

- all functions/files have been tested on Macbook 3.1 GHz Quad-Core Intel Core i7; Other OS runing gcc/g++ should also work but you need to work out proper configurations
- be careful with the directory of Eigen, Boost and Yaml-cpp, better to install them using CMake inside their file folders
  - For instance I install yaml-cpp in the current directory using these command:
    ```sh
    git clone https://github.com/jbeder/yaml-cpp.git
    cd yaml-cpp
    mkdir build
    cd build
    cmake -DCMAKE_OSX_ARCHITECTURES=arm64 ..
    make
    ```
- adjust the directories in Makefile accordingly to include these libraries
- the author is updating the repo from time to time

Reference:

- please cite my paper when you use the codes: Yang Yang, Square-Root Higher-Order Unscented Estimators for Robust Orbit Determination, IEEE Transactions on Aerospace and Electronic Systems, DOI: 10.1109/TAES.2024.3423851 (Early Access).
