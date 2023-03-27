# HOUSE: The Higher-Order Unscented Estimator

Running:
 - the repo includes an executable compiled for windows
 - parameters can be set using the config.yaml file
 - by default the program will read the "config.yaml" file, however
   you can provide an argument to read a different yaml file (e.g. "./filter_testing other.yaml")
 - Note file options in the yaml file do not yet affect the program, the program will read fromt the supplied
   input file and output to the out directory

Compiling:
 - use g++ --std=c++11
 - compilation is handled by the makefile, run make to rebuild (make must be installed)
 - the binary object files are stored in the bin directory

Dependencies
 - CUTpoints by N. Adurthi, P. Singla, and T. Singh: https://github.com/nadurthi/CUTpoints using the Matlab function cut_sigma_points.m to generate .csv files in the CUT directory
 - Eigen: https://eigen.tuxfamily.org/
 - Boost Installation: https://www.boost.org/doc/libs/1_69_0/more/getting_started/unix-variants.html#easy-build-and-install
 - Yaml-cpp

 Analysis and plot
 - call R to run processing_all.r and plot_rmse_all.r sequently

Note: be careful with the directory of Eigen/Boost and Yaml-cpp, better to install them using CMake inside their filefolders
Note: ode.hpp and ode.cpp are released under the GNU Lesser General Public License
(https://www.gnu.org/licenses/lgpl-3.0.en.html)  


# HOUSE
