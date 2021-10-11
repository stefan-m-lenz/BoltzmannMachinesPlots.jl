# Run examples in BoltzmannMachines.jl with plotting
import BoltzmannMachines
const BMs = BoltzmannMachines
include(joinpath(dirname(pathof(BoltzmannMachines)), "..", "test", "examples.jl"))

using Test
import Gadfly

using BoltzmannMachinesPlots
function test_scatterhidden()
   hidden = rand(100, 2)
   labels = rand(["1", "2", "3"], 100)
   scatterhidden(hidden,
         opacity = 0.5, labels = labels)
end
test_scatterhidden()

# TODO
# function test_plotevaluation_noribbon()
#    monitor, rbm = BMs.monitored_fitrbm(BMs.barsandstripes(10, 9),
#       monitoring = monitorloglikelihood!)
#    @test plotevaluation(monitor; sdrange = 0.0) isa Gadfly.Plot
# end
# test_plotevaluation_noribbon()


@test BoltzmannMachinesPlots.plotcurvebundles(BMs.curvebundles(nvariables = 10, nbundles = 3,
      nperbundle = 4, noisesd = 0.03,
      addlabels = true)) isa Gadfly.Plot

@test BoltzmannMachinesPlots.plotcurvebundles(BMs.curvebundles(nvariables = 10, nbundles = 3,
      nperbundle = 4, noisesd = 0.03,
      addlabels = false)) isa Gadfly.Plot
