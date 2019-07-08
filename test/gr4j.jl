using WaterFlows
using Test
using Dates

tstep = 24.0
tstart = DateTime(2017, 1, 1)

gr4j = Gr4j(tstep, tstart)

@testset "Gr4j" begin

    @test sum(gr4j.ord_uh1) == 1.0

    @test sum(gr4j.ord_uh2) == 1.0

    water_start = get_water_stored(gr4j)

    run_timestep(gr4j)

    water_end = get_water_stored(gr4j)

    @test water_start ≈ water_end + gr4j.q_out - gr4j.exch atol = 0.01

end