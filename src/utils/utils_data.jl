# Load operational data from text files

function load_data(path, file_q_obs = "runoff.txt", file_tair = "tair.txt",
                   file_prec = "prec.txt", file_metadata = "metadata.txt")

  # Read air temperature data

  str   = readline("$path/$file_tair")
  nsep  = length(collect(m.match for m = eachmatch(r";", str)))
  tmp   = CSV.read("$path/$file_tair", delim = ";", header = false,
                   dateformat="yyyy-mm-dd HH:MM", allowmissing=:none, types = vcat(DateTime, repeat([Float64], nsep)))
  tair  = convert(Array, tmp[:, 2:end])
  tair  = permutedims(tair)

  # Read precipitation data

  str   = readline("$path/$file_tair")
  nsep  = length(collect(m.match for m = eachmatch(r";", str)))
  tmp   = CSV.read("$path/$file_prec", delim = ";", header = false,
                  dateformat="yyyy-mm-dd HH:MM", allowmissing=:none, types = vcat(DateTime, repeat([Float64], nsep)))
  prec  = convert(Array, tmp[:, 2:end])
  prec  = permutedims(prec)

  # Read runoff data

  tmp   = CSV.read("$path/$file_q_obs", delim = ";", header = false,
                   dateformat="yyyy-mm-dd HH:MM", allowmissing=:none, types = [DateTime, Float64])
  q_obs = convert(Array, tmp[:, 2])

  q_obs[q_obs .< 0.0] .= NaN

  # Read metadata

  df_tmp = CSV.read("$path/$file_metadata", delim = ";", header = true)

  area = convert(Array{Float64,1}, df_tmp[:area_sum])
  frac_area = area / sum(area)

  frac_glacier = convert(Array{Float64,1}, df_tmp[:lus_glacier_mean]) / 100.0

  frac_lus = DataFrame()
  frac_lus[:glacier] = frac_area .* frac_glacier
  frac_lus[:open] = frac_area - frac_lus[:glacier]

  elev = convert(Array{Float64,1}, df_tmp[:elevation_mean])

  # Get time data

  date = convert(Array, tmp[:, 1])

  # Return data

  return date, tair, prec, q_obs, frac_lus, frac_area, elev

end


# Crop data from start to stop date

function crop_data(date, tair, prec, q_obs, date_start, date_stop)

  # Find indicies

  istart = findall(date .== date_start)
  istop = findall(date .== date_stop)

  # Test if ranges are valid

  if isempty(istart) | isempty(istop)
    error("Cropping data outside range")
  end

  # Crop data

  date  = date[istart[1]:istop[1]]
  tair  = tair[:, istart[1]:istop[1]]
  prec  = prec[:, istart[1]:istop[1]]
  q_obs = q_obs[istart[1]:istop[1]]

  return date, tair, prec, q_obs

end


# Crop data before start date

function crop_data(date, tair, prec, q_obs, date_start)

  # Find indicies

  istart = findall(date .== date_start)

  # Test if ranges are valid

  if isempty(istart)
    error("Cropping data outside range")
  end

  # Crop data

  date  = date[istart[1]:end]
  tair  = tair[:, istart[1]:end]
  prec  = prec[:, istart[1]:end]
  q_obs = q_obs[istart[1]:end]

  return date, tair, prec, q_obs

end

