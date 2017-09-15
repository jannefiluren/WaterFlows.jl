

function setup_hbv_light(tstep, time, frac_lus)

    lake = 0.0  # TODO: add this to land use

    snow = HbvLightSnow(tstep, time, frac_lus)
    
    glacier = NoGlacier()
    
    hydro = HbvLightSubsurf(tstep, time, frac_lus, lake)

    model = ModelComp(snow, glacier, hydro)

end
