using DelimitedFiles, DifferentialEquations, SparseArrays, Test

function julianise_numpy_matrix()
  """
  Read the A.txt file and return a sparse matrix containing the decay rates
  """
  fileoutput = readdlm("A.txt", ',')
  nrows, ncols = fileoutput[1, :]
  A = spzeros(nrows, ncols)
  for k in 3:size(fileoutput, 1)
    j, i, v = fileoutput[k, :] # tranpose i-j to j-i here is correct
    A[i + 1, j + 1] = v # add one because it started life in python
  end
  return A
end

function readindices()
  """
  Read the index file, which is strictly necessary but is helpful for humans
  """
  fileoutput = readdlm("index.txt", ',')
  index2name = Dict()
  name2index = Dict()
  for k in 3:size(fileoutput, 1)
    i, z, name = fileoutput[k, :]
    index = i + 1 # add one because it started life in python
    index2name[index] = name
    name2index[name] = index
  end
  return index2name, name2index
end

function readinventory(filename)
  """
  Obtain the initial conditions and the expected results for a given problem
  """
  fileoutput = readdlm(filename, ',')
  duration = parse(Float64, fileoutput[1, 2][4:end])
  initialconditions = Dict()
  expectedresults = Dict()
  for k in 3:size(fileoutput, 1)
    _, _, name, n0, nN = fileoutput[k, :]
    initialconditions[name] = n0
    expectedresults[name] = nN
  end
  return initialconditions, expectedresults, duration
end

function dosolve(u0, integrator, A, tspan, rtol=1.0e-3)
  """
  Solve the problem given the initial conditions u0, ODE integrator, decay rate matrix A,
  time span tspan, and relative tolerance of for the solver.
  """
  if integrator == exp 
    expA = exp(Matrix(A * tspan[2])) # Matrix needed to turn to sparse A to full
    return expA * u0
  end
  f(u, p, t) = A * u
  prob = ODEProblem(f, u0, tspan)
  sol = solve(prob, integrator, reltol=rtol, abstol=0.0)
  return sol.u[end]
end

function run(inventoryfilename, rtol=1.0e-3)
  """
  Run the test for the given inventory file
  """
  @show inventoryfilename, rtol
  A = julianise_numpy_matrix() # the du/dt = A u
  index2name, name2index = readindices() # useful dictionaries
  initialconditions, expectedresults, duration = readinventory(inventoryfilename)

  u0 = zeros(size(A, 1))
  for (name, n) in initialconditions
    @assert 1 <= name2index[name] <= size(A, 1) "$name, $(name2index[name]), $(size(A))"
    u0[name2index[name]] = n
  end

  #nameintegrators = zip(("ImplicitMidpoint",), (ImplicitMidpoint(),))
  #nameintegrators = zip(("ImplicitEuler",), (ImplicitEuler(),))
  #nameintegrators = zip(("RadauIIA5", ), (RadauIIA5(), ))
  #nameintegrators = zip(("AutoVern9Rodas5", ), (AutoVern9(Rodas5()), ))
  nameintegrators = zip(("AutoTsit5Rosenbrock23", ), (AutoTsit5(Rosenbrock23()), ))
  #nameintegrators = zip(("exp", ), (exp, ))
  for (integratorname, integrator) in nameintegrators
    timetaken = @elapsed results = dosolve(deepcopy(u0), integrator, A, (0.0, duration))

    all_results_expected_to_be_zero_are_zero = true
    for (i, result) in enumerate(results)
      # result nuclide should be in expected list
      name = index2name[i]
      expected = expectedresults[name]
      #iszero(result + expected) || @show i, name, expected, result
      if iszero(expected)
        all_results_expected_to_be_zero_are_zero &= iszero(result)
      else
        #@test isapprox(expected, result, rtol=rtol, atol=0.0) #atol=0.0 is important
      end
    end
    #@test all_results_expected_to_be_zero_are_zero
  end
end


# Below is a list of all available initial conditions and results
# Ideally a single ODE solver can be found that works for all cases
# This is the challenge!
inventories = []
push!(inventories, "decay_inventory_H3_3.891050e+08.txt")
push!(inventories, "decay_inventory_Cf252_8.346980e+07.txt")
push!(inventories, "decay_inventory_Cf250_4.127730e+08.txt")
push!(inventories, "decay_inventory_Fe56_1.000000e-02.txt")
push!(inventories, "decay_inventory_Fe56_1.000000e+15.txt")
push!(inventories, "decay_inventory_Fe56_2.300000e+08.txt")
push!(inventories, "decay_inventory_full_1.000000e-02.txt")
push!(inventories, "decay_inventory_full_1.000000e+05.txt")
push!(inventories, "decay_inventory_full_1.000000e+15.txt")
push!(inventories, "decay_inventory_full_1.000000e-18.txt")

@testset "Test all inventories" begin
  for inventory in inventories
    run(inventory)
  end
end
