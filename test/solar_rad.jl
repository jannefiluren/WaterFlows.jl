using VannModels
using Test

@testset "Solar radiation" begin

    date = DateTime(1980, 7, 20)
    doy = Dates.dayofyear(date)
    lat = -23.7951
    elev = 546.0

    dr = VannModels.inverse_dist(doy)
    @test dr ≈ 0.9688418122084714

    δ = VannModels.solar_decl(doy)
    @test δ ≈ 0.35565253560155585

    ω_s = VannModels.sunset_hour_angle(deg2rad(lat), δ)
    @test ω_s ≈ 1.4062650741766252

    N = VannModels.daylight_hours(ω_s)
    @test N ≈ 10.743073816929636

    R_a = VannModels.extra_ter_rad(date, lat)
    @test R_a ≈ 23.618220461289415

    R_so = VannModels.clear_sky_rad(date, lat, elev)
    @test R_so ≈ 17.97157631340434

end