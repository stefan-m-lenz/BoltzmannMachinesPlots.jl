# Run examples in BoltzmannMachines.jl with plotting
import BoltzmannMachines
const BMs = BoltzmannMachines
include(joinpath(dirname(pathof(BoltzmannMachines)), "..", "test", "examples.jl"))

using Test
import Gadfly

using BoltzmannMachinesPlots

function test_scatterhidden()
   nsamples = 10
   x = BMs.barsandstripes(10,4)
   rbm = fitrbm(x, epochs = 1)
   dbm = fitdbm(x, epochs = 1)
   labels = rand(["1", "2", "3"], nsamples)
   @test scatterhidden(rbm, x,
         opacity = 0.5, labels = labels) isa Gadfly.Plot
   @test scatterhidden(dbm, x,
         labels = labels, opacity = 0.5) isa Gadfly.Plot
end
test_scatterhidden()


function test_plotevaluation_noribbon()
   monitor, rbm = BMs.monitored_fitrbm(BMs.barsandstripes(10, 9),
      monitoring = monitorloglikelihood!)
   @test plotevaluation(monitor; sdrange = 0.0) isa Gadfly.Plot
end
test_plotevaluation_noribbon()


function test_crossvalidation()
   nsamples = 30;
   nvariables = 4;
   x = barsandstripes(nsamples, nvariables);
   # Determine the optimal number of training epochs for a RBM
   monitor = crossvalidation(x,
         (x, datadict) ->
            begin
               monitors, dbm = BMs.monitored_fitdbm(x, nhiddens = [2, 2],
                     learningrate = 0.1, epochs = 3,
                     monitoring = monitorexactloglikelihood!)
               monitors[end]
            end);
   @test crossvalidationcurve(monitor) isa Gadfly.Plot
   nothing
end
test_crossvalidation()


@test BoltzmannMachinesPlots.plotcurvebundles(BMs.curvebundles(nvariables = 10, nbundles = 3,
      nperbundle = 4, noisesd = 0.03,
      addlabels = true)) isa Gadfly.Plot

@test BoltzmannMachinesPlots.plotcurvebundles(BMs.curvebundles(nvariables = 10, nbundles = 3,
      nperbundle = 4, noisesd = 0.03,
      addlabels = false)) isa Gadfly.Plot
