# Compile main library

g++ -c *.cpp --std=c++11 -Wall -pedantic -g -O3

rm -f house.a

ar rcs house.a *.o

rm -f *.o

# Compile projectile example

cd projectile

g++ -o projectile.exe projectile.cpp -I.. ../house.a --std=c++11 \
             -O3 -g -Wall -pedantic

cd ..

# Compile rigid body example

cd rigid_body

g++ -o rigid_body.exe rigid_body.cpp -I.. ../house.a --std=c++11 \
             -O3 -g -Wall -pedantic

cd ..

# Compile coordinated turn example with Gaussian noise

cd CT_Gauss

g++ -o ct.exe CT.cpp ../house.a -I.. --std=c++11 -g -Wall -pedantic -O3

g++ -o st.exe statct.cpp -I.. --std=c++11 -g -Wall -pedantic -O3

cd ..

# Compile coordinated turn example with Pearson type IV noise

cd CT_Pearson

g++ -o ct.exe CT.cpp ../house.a -I.. --std=c++11 -g -Wall -pedantic -O3

g++ -o st.exe statct.cpp -I.. --std=c++11 -g -Wall -pedantic -O3

cd ..
