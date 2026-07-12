using WaterFlows
using Test
using Dates
using DataFrames
using Statistics

tstart = DateTime(2017, 1, 1)


@testset "Unit hydrograph across time steps" begin

    for tstep in (1.0, 3.0, 6.0, 24.0), maxbas in (1.0, 2.5, 7.0)

        ord_uh = WaterFlows.compute_hbv_ord(maxbas, tstep)

        @test sum(ord_uh) ≈ 1.0
        @test length(ord_uh) == ceil(Int, maxbas * 24.0 / tstep)

    end

end


@testset "Hbv recession equivalence" begin

    hbv_daily = Hbv(24.0, tstart)
    hbv_hourly = Hbv(1.0, tstart)

    for hbv in (hbv_daily, hbv_hourly)
        init_states!(hbv, tstart)
        hbv.sm = 0.0
        hbv.suz = 0.0
        hbv.slz = 100.0
        hbv.p_in = 0.0
        hbv.epot = 0.0
    end

    run_timestep(hbv_daily)

    for _ in 1:24
        run_timestep(hbv_hourly)
    end

    @test hbv_daily.time == hbv_hourly.time
    @test hbv_daily.slz ≈ hbv_hourly.slz
    @test hbv_daily.slz ≈ 100.0 * (1.0 - hbv_daily.k2)

end


@testset "HbvLight recession equivalence" begin

    frac_lus = DataFrame(glacier = [0, 0], open = [0.5, 0.5])

    hbv_daily = HbvLightSubsurf(24.0, tstart, frac_lus)
    hbv_hourly = HbvLightSubsurf(1.0, tstart, frac_lus)

    for hbv in (hbv_daily, hbv_hourly)
        init_states!(hbv, tstart)
        hbv.slz = 100.0
        hbv.p_in .= 0.0
        hbv.epot = 0.0
    end

    run_timestep(hbv_daily)

    for _ in 1:24
        run_timestep(hbv_hourly)
    end

    @test hbv_daily.time == hbv_hourly.time
    @test hbv_daily.slz ≈ hbv_hourly.slz
    @test hbv_daily.slz ≈ 100.0 * (1.0 - hbv_daily.k2)

end


@testset "Hbv water balance hourly" begin

    hbv = Hbv(1.0, tstart)

    init_states!(hbv, tstart)

    hbv.epot = 0.5
    hbv.p_in = 2.0

    water_start = get_water_stored(hbv)

    run_timestep(hbv)

    water_end = get_water_stored(hbv)

    @test water_start + hbv.p_in ≈ water_end + hbv.q_out + hbv.aevap

end


@testset "HbvLight water balance hourly" begin

    frac_lus = DataFrame(glacier = [0, 0], open = [0.5, 0.5])

    hbv_light = HbvLightSubsurf(1.0, tstart, frac_lus)

    init_states!(hbv_light, tstart)

    hbv_light.epot = 0.5
    hbv_light.p_in .= 2.0

    water_start = get_water_stored(hbv_light)

    run_timestep(hbv_light)

    water_end = get_water_stored(hbv_light)

    @test water_start + mean(hbv_light.p_in) ≈ water_end + hbv_light.q_out + hbv_light.aevap atol = 0.0001

end


@testset "TinSnow melt equivalence" begin

    frac_lus = DataFrame(fill(0.25, (2, 2)), :auto)

    tin_daily = TinSnow(24.0, tstart, frac_lus)
    tin_hourly = TinSnow(1.0, tstart, frac_lus)

    for tin in (tin_daily, tin_hourly)
        init_states!(tin, tstart)
        tin.swe .= 100.0
        tin.p_in .= 0.0
        tin.tair .= 10.0
    end

    water_start = get_water_stored(tin_daily)

    run_timestep(tin_daily)

    melt_daily = water_start - get_water_stored(tin_daily)

    for _ in 1:24
        run_timestep(tin_hourly)
    end

    melt_hourly = water_start - get_water_stored(tin_hourly)

    @test tin_daily.time == tin_hourly.time
    @test melt_daily ≈ melt_hourly
    @test melt_daily > 0.0

end


@testset "HbvLightSnow melt equivalence" begin

    frac_lus = DataFrame(fill(0.25, (2, 2)), :auto)

    snow_daily = HbvLightSnow(24.0, tstart, frac_lus)
    snow_hourly = HbvLightSnow(1.0, tstart, frac_lus)

    for snow in (snow_daily, snow_hourly)
        init_states!(snow, tstart)
        snow.swe .= 100.0
        snow.p_in .= 0.0
        snow.tair .= 10.0
    end

    run_timestep(snow_daily)

    for _ in 1:24
        run_timestep(snow_hourly)
    end

    @test snow_daily.time == snow_hourly.time
    @test get_water_stored(snow_daily) ≈ get_water_stored(snow_hourly)
    @test all(snow_daily.swe .≈ snow_hourly.swe)

end
