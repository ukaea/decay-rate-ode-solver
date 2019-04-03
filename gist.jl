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

function removerowcol(A, rowcol)
  B = spzeros(size(A, 1) - 1, size(A, 1) - 1)
  for k in eachindex(A)
    i, j = k[1], k[2]
    i == rowcol && continue
    j == rowcol && continue
    i > rowcol && (i -= 1)
    j > rowcol && (j -= 1)
    B[i, j] = A[k]
  end
  A = B
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
  allnames = [index2name[i] for i in sort(collect(keys(index2name)))]
  return allnames
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
  sol = solve(prob, integrator, reltol=rtol)#, abstol=0.0)
  return sol.u[end]
end

function run(inventoryfilename, rtol=1.0e-3)
  """
  Run the test for the given inventory file
  """
  @show inventoryfilename, rtol
  A = julianise_numpy_matrix() # the du/dt = A u
  initialconditions, expectedresults, duration = readinventory(inventoryfilename)
  allnames = readindices() # useful dictionaries

  u0 = zeros(2233)
  for (name, n) in initialconditions
    u0[findfirst(allnames .== name)] = n
  end

  deletednames = []

  c = [a[2] for a in findall(.!iszero.(A))]
  r = [a[1] for a in findall(.!iszero.(A))]
  deletedindices = reverse(sort(collect(setdiff(collect(1:size(A, 1)), sort(unique(vcat(c, r)))))))
  @show deletednames = allnames[deletedindices]
  for deletedname in deletednames
    i = findfirst(allnames .== deletedname)
    deleteat!(u0, i)
    deleteat!(allnames, i)
    A = removerowcol(A, i)
  end
  @assert !any(sum(A, dims=2)[:] .== 0)

  #nameintegrators = zip(("ImplicitMidpoint",), (ImplicitMidpoint(),))
  #nameintegrators = zip(("ImplicitEuler",), (ImplicitEuler(),))
  #nameintegrators = zip(("RadauIIA5", ), (RadauIIA5(), ))
  nameintegrators = zip(("AutoVern9Rodas5", ), (AutoVern9(Rodas5()), ))
  #nameintegrators = zip(("AutoTsit5Rosenbrock23", ), (AutoTsit5(Rosenbrock23()), ))
  #nameintegrators = zip(("exp", ), (exp, ))
  for (integratorname, integrator) in nameintegrators
    timetaken = @elapsed results = dosolve(deepcopy(u0), integrator, A, (0.0, duration))

    all_results_expected_to_be_zero_are_zero = true
    for (i, result) in enumerate(results)
      # result nuclide should be in expected list
      name = allnames[i]
      in(name, deletednames) && continue
      expected = expectedresults[name]
      if iszero(expected)
        all_results_expected_to_be_zero_are_zero &= iszero(result)
      else
        outcome = isapprox(expected, result, rtol=rtol, atol=0.0) #atol=0.0 is important
        outcome || @show name, expected, result
        @test isapprox(expected, result, rtol=rtol, atol=0.0) #atol=0.0 is important
      end
    end
    @test all_results_expected_to_be_zero_are_zero
  end
end


# Below is a list of all available initial conditions and results
# Ideally a single ODE solver can be found that works for all cases
# This is the challenge!
inventories = []
push!(inventories, "decay_inventory_H3_3.891050e+08.txt")
#push!(inventories, "decay_inventory_Cf252_8.346980e+07.txt")
#push!(inventories, "decay_inventory_Cf250_4.127730e+08.txt")
#push!(inventories, "decay_inventory_Fe56_1.000000e-02.txt")
#push!(inventories, "decay_inventory_Fe56_1.000000e+15.txt")
#push!(inventories, "decay_inventory_Fe56_2.300000e+08.txt")
#push!(inventories, "decay_inventory_full_1.000000e-02.txt")
#push!(inventories, "decay_inventory_full_1.000000e+05.txt")
#push!(inventories, "decay_inventory_full_1.000000e+15.txt")
#push!(inventories, "decay_inventory_full_1.000000e-18.txt")

@testset "Test all inventories" begin
  for inventory in inventories
    run(inventory)
  end
end
