

function model_hbv_light(tstep, time, frac_lus)

    lake = 0.0  # TODO: add this to land use

    snow = HbvLightSnow(tstep, time, frac_lus)
    
    glacier = NoGlacier()
    
    hydro = HbvLightSubsurf(tstep, time, frac_lus, lake)

    model = ModelComp(snow, glacier, hydro)

end



function model_gr4j(tstep, time, frac_lus)
    
    snow = TinSnow(tstep, time, frac_lus)
    
    glacier = TinGlacier(tstep, time, frac_lus)
    
    hydro = Gr4j(tstep, time)

    model = ModelComp(snow, glacier, hydro)

end
