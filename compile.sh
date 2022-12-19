# Compile main library

echo "Compiling main library..."
g++ -c *.cpp -I./eigen -I. -Wall -pedantic -g -O3 

rm -f house.a

ar rcs house.a *.o

rm -f *.o

# # Compile projectile example

# echo "Compiling projectile..."
# cd projectile

# g++ -o projectile.exe projectile.cpp -I.. -I./../eigen ../house.a \
#              -O3 -g -Wall -pedantic

# cd ..

# # Compile rigid body example
# echo "Compiling rigid body..."

# cd rigid_body

# g++ -o rigid_body.exe rigid_body.cpp -I.. -I./../eigen ../house.a \
#              -O3 -g -Wall -pedantic

# cd ..

# # Compile coordinated turn example with Gaussian noise
# echo "Compiling CT with Gaussian Noise..."

# cd CT_Gauss

# g++ -o ct.exe CT.cpp ../house.a -I.. -I./../eigen  -g -Wall -pedantic -O3

# g++ -o st.exe statct.cpp -I.. -I./../eigen -g -Wall -pedantic -O3

# cd ..

# # Compile coordinated turn example with Pearson type IV noise
# echo "Compiling CT with Pearson Noise..."

# cd CT_Pearson

# g++ -o ct.exe CT.cpp ../house.a -I.. -I./../eigen -g -Wall -pedantic -O3

# g++ -o st.exe statct.cpp -I.. -I./../eigen -g -Wall -pedantic -O3

# cd ..
