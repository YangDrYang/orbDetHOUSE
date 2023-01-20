# HOUSE: The Higher-Order Unscented Estimator

Running:
 - the repo includes an executable compiled for windows
 - parameters can be set using the config.yaml file
 - by default the program will read the "config.yaml" file, however
   you can provide an argument to read a different yaml file (e.g. "./filter_testing other.yaml")
 - Note file options in the yaml file do not yet affect the program, the program will read fromt the supplied
   input file and output to the out directory

Compiling:
 - compilation is handled by the makefile, run make to rebuild (make must be installed)
 - the binary object files are stored in the bin directory

Dependencies
 - CUTpoints by N. Adurthi, P. Singla, and T. Singh: https://github.com/nadurthi/CUTpoints
 - Eigen: https://eigen.tuxfamily.org/
 - Boost
 - Yaml-cpp

Note: ode.hpp and ode.cpp are released under the GNU Lesser General Public License
(https://www.gnu.org/licenses/lgpl-3.0.en.html)  


# HOUSE
