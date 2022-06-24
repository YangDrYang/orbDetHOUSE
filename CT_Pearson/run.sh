mkdir -p out

rm -f out/*

bash compile.sh

./ct.exe

./st.exe
