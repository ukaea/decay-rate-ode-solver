using DelimitedFiles, DifferentialEquations, SparseArrays, Test

function julianise_numpy_matrix()
  fileoutput = readdlm("A.txt", ',')
  nrows, ncols = fileoutput[1, :]
  A = spzeros(nrows, ncols)
  for k in 3:size(fileoutput, 1)
    i, j, v = fileoutput[k, :]
    A[i + 1, j + 1] = v # add one because it started life in python
  end
  return A
end

function readindices()
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

function dosolve(u0, integrator, A, tspan, rtol=1.0e-2)
  if integrator == exp 
    expA = exp(Matrix(A * tspan[2])) # Matrix needed to turn to sparse A to full
    return expA * u0
  end
  f(u, p, t) = A * u
  prob = ODEProblem(f, u0, tspan)
  sol = solve(prob, integrator, reltol=rtol, dt=tspan[2]/1000)
  return sol.u[end]
end

function run(inventoryfilename, rtol=1.0e-3)
  @show inventoryfilename, rtol
  A = julianise_numpy_matrix() # the du/dt = A u
  index2name, name2index = readindices() # useful dictionaries
  initialconditions, expectedresults, duration = readinventory(inventoryfilename)

  u0 = zeros(size(A, 1))
  for (name, n) in initialconditions
    @assert 1 <= name2index[name] <= size(A, 1) "$name, $(name2index[name]), $(size(A))"
    u0[name2index[name]] = n
  end

  nameintegrators = zip(("ImplicitMidpoint",), (ImplicitMidpoint(),))
  #nameintegrators = zip(("ImplicitEuler",), (ImplicitEuler(),))
  #nameintegrators = zip(("RadauIIA5", ), (RadauIIA5(), ))
  #nameintegrators = zip(("AutoVern9Rodas5", ), (AutoVern9(Rodas5()), ))
  #nameintegrators = zip(("exp", ), (exp, ))
  for (integratorname, integrator) in nameintegrators
    timetaken = @elapsed results = dosolve(deepcopy(u0), integrator, A, (0.0, duration))
    @show integratorname, timetaken
    for (i, result) in enumerate(results)
      # result nuclide should be in expected list
      if !haskey(expectedresults, index2name[i])
        @test false # a failure
      else
        expected = expectedresults[index2name[i]]
        @test isapprox(expected, result, rtol=rtol, atol=0.0) #atol=0.0 is important
      end
    end
  end
end

inventories = []
push!(inventories, "decay_inventory_H3_3.891050e+08.txt")
push!(inventories, "decay_inventory_Cf252_8.346980e+07.txt")
push!(inventories, "decay_inventory_Cf250_4.127730e+08.txt")
push!(inventories, "decay_inventory_Fe56_1.000000e-02.txt")
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
