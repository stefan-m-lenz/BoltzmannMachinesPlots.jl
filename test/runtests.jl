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


function test_plottop2latentdims()
   x, xlabels = BMs.blocksinnoise(50, 9, nblocks = 2, blocklen = 2)
   dbm = BMs.fitdbm(x, epochs = 1)
   @test plottop2latentdims(dbm, x) isa Gadfly.Plot
   @test plottop2latentdims(dbm, x; labels = xlabels) isa Gadfly.Plot
   nothing
end
test_plottop2latentdims()