using WaterFlows
using Test
using Dates
using DataFrames
using Statistics

tstep = 24.0
tstart = DateTime(2017, 1, 1)
frac_lus = fill(1/4, (2,2))

tinsnow = TinSnow(tstep, tstart, DataFrame(frac_lus))

@testset "TinSnow" begin

    tinsnow.swe = fill(100.0, (2,2))
    tinsnow.p_in = fill(10.0, 2)
    tinsnow.tair = fill(10.0, 2)

    water_start = get_water_stored(tinsnow)

    run_timestep(tinsnow)

    water_end = get_water_stored(tinsnow)
    
    @test water_start + mean(tinsnow.p_in) ≈ water_end + mean(tinsnow.q_out)

    tinsnow.swe = fill(100.0, (2,2))
    tinsnow.p_in = fill(10.0, 2)
    tinsnow.tair = fill(-10.0, 2)

    water_start = get_water_stored(tinsnow)

    run_timestep(tinsnow)

    water_end = get_water_stored(tinsnow)
    
    @test water_start + mean(tinsnow.p_in) ≈ water_end + mean(tinsnow.q_out)

end


hbvlightsnow = HbvLightSnow(tstep, tstart, DataFrame(frac_lus))

@testset "HbvLightSnow" begin

    hbvlightsnow.swe = fill(100.0, (2,2))
    hbvlightsnow.p_in = fill(10.0, 2)
    hbvlightsnow.tair = fill(10.0, 2)

    water_start = get_water_stored(hbvlightsnow)

    run_timestep(hbvlightsnow)

    water_end = get_water_stored(hbvlightsnow)

    @test water_start + mean(hbvlightsnow.p_in) ≈ water_end + mean(hbvlightsnow.q_out)

    hbvlightsnow.swe = fill(100.0, (2,2))
    hbvlightsnow.p_in = fill(10.0, 2)
    hbvlightsnow.tair = fill(-10.0, 2)

    water_start = get_water_stored(hbvlightsnow)

    run_timestep(hbvlightsnow)

    water_end = get_water_stored(hbvlightsnow)

    @test water_start + mean(hbvlightsnow.p_in) ≈ water_end + mean(hbvlightsnow.q_out)

end






