# decay-rate-ode-solver
Solver for a hyper-stiff set of ODEs

An attempt to find the best ODE solver for this problem type which is simply `du/dt = Au`, where `A` is a matrix independent of time and `u`, and `u` is a vector. 

There is one julia file `gist.jl`, a matrix file `A.txt`, an index file `index.txt`, and several results files `decay_inventory*.txt`.
The `decay_inventory*.txt` contain initial conditions, and expected results which are used to judge the performance of ODE solvers called from julia.

Simply run `gist.jl` to try it out. It should be self explanatory, so feel free to experiment to try to get all tests passing

N.B. has some zero row / columns, which will cause problems with finding eigen solutions.
