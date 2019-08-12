using WaterFlows
using Test
using Dates

tstep = 24.0
tstart = DateTime(2017, 1, 1)


gr4j = Gr4j(tstep, tstart)

gr4j.epot = 10
gr4j.p_in = 10

@testset "Gr4j" begin

    @test sum(gr4j.ord_uh1) == 1.0

    @test sum(gr4j.ord_uh2) == 1.0

    water_start = get_water_stored(gr4j)

    run_timestep(gr4j)

    water_end = get_water_stored(gr4j)

    @test water_start + gr4j.p_in ≈ water_end + gr4j.q_out + gr4j.aevap - gr4j.exch atol = 0.01

end


hbv = Hbv(tstep, tstart)

hbv.epot = 10
hbv.p_in = 10

@testset "Hbv" begin

    @test sum(WaterFlows.compute_hbv_ord(1.0)) ≈ 1.0
    @test sum(WaterFlows.compute_hbv_ord(20.0)) ≈ 1.0

    WaterFlows.init_states!(hbv, tstart)

    water_start = get_water_stored(hbv)
    
    run_timestep(hbv)
    
    water_end = get_water_stored(hbv)

    @test water_start + hbv.p_in ≈ water_end + hbv.q_out + hbv.aevap

end