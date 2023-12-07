## 2-limited-broadcast-domination-grid-graph-prove-lower-bounds
This is a method for proving lower bounds for the 2-limited broadcast domination number on the Cartesian product of two paths, a path and a cycle, and two cycles. This method requires a local copy of a gurobi solver (https://www.gurobi.com/). The file gurobi.env contains some gurobi settings. The file prove_lower_bounds.cpp implements the work described in https://doi.org/10.46298/dmtcs.11478 in C++. 

# Compile:
To compile with gurobi, run the following command (assuming you have a mac and used the default directory for gurobi)

g++ -m64 -g -O3 prove_lower_bounds.cpp -I/Library/gurobi1000/macos_universal2/include -L/Library/gurobi1000/macos_universal2/lib -lgurobi_c++ -lgurobi100 -lm

Note: Make sure to change to the appropriate version of gurobi that you are running, i.e. change "gurobi1000" and "-lgurobi100" as necessary

# Run:
The compiled a.out will run the instance of the problem defined in lines 6 through 25 of prove_lower_bounds.cpp.

PDFs: If you want to view a PDF of the results, set PROD_FIG on line 20 to 1 and 0 otherwise. Additionally, run a.out >> output.py. This will output a python file which, when run with 

python output.py

will produce the desired PDF by using the function in graphic.py.
