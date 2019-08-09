using WaterFlows
using Test

@testset "SnowUtils" begin

    # Splitting of precipitation

    psolid, pliquid = WaterFlows.split_prec(1, -10)

    @test psolid ≈ 0.9999999979388463
    @test pliquid ≈ 1 - 0.9999999979388463
    @test psolid + pliquid ≈ 1

    psolid, pliquid = WaterFlows.split_prec(1, 0)

    @test psolid ≈ 0.5
    @test pliquid ≈ 0.5
    @test psolid + pliquid ≈ 1

    psolid, pliquid = WaterFlows.split_prec(1, 10)

    @test psolid ≈ 2.0611536181902037e-9
    @test pliquid ≈ 1 - 2.0611536181902037e-9
    @test psolid + pliquid ≈ 1

    # Computation of potential melt

    melt = WaterFlows.pot_melt(-10, 5)

    @test melt ≈ 5.1528825650848376e-9

    melt = WaterFlows.pot_melt(0, 5)

    @test melt ≈ 1.7328679513998633

    melt = WaterFlows.pot_melt(10, 5)

    @test melt ≈ 50.000000005152884

    # Computation of potential refreeze

    pot_refreeze = WaterFlows.pot_refreeze(-10, 5)

    @test pot_refreeze ≈ 50.000000005152884

    pot_refreeze = WaterFlows.pot_refreeze(0, 5)

    @test pot_refreeze ≈ 1.7328679513998633

    pot_refreeze = WaterFlows.pot_refreeze(10, 5)

    @test pot_refreeze ≈ 5.1528825650848376e-9

end