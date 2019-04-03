# decay-rate-ode-solver
Solver for a hyper-stiff set of ODEs

An attempt to find the best ODE solver for this problem.

There is one julia file `gist.jl`, a matrix file `A.txt`, an index file `index.txt`, and several results files `inventory*.txt`.
The `inventory*.txt` contain initial conditions, and expected results which are used to judge the performance of ODE solvers called from julia.

Simply run `gist.jl` to try it out. It should be self explanatory, so feel free to experiment to try to get all tests passing

