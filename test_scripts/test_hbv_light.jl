

using Vann2



tstep = 24.0

time = DateTime(2000,1,1)

area = ones(2,2)/4

lake = 0.0

snow = HbvLightSnow(tstep, time, area)

subsurf = HbvLightSubsurf(tstep, time, area, lake)

snow.tair .= 10
snow.p_in .= 0

run_timestep(snow)

subsurf.p_in .= snow.q_out
subsurf.epot = 10.9

subsurf.tair .= 10

run_timestep(subsurf)




