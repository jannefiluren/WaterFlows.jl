# Load data from text files

#=function load_data(path, file_q_obs = "Q_obs.txt", file_tair = "Tair.txt",
                   file_prec = "Prec.txt", file_frac = "Frac.txt")

  # Read air temperature data

  str   = readline("$path/$file_tair")
  nsep  = length(matchall(r";", str))
  tmp   = CSV.read("$path/$file_tair", delim = ";", header = false,
                   dateformat="yyyy-mm-dd HH:MM", nullable = false, types = vcat(DateTime, repmat([Float64], nsep)))
  tair  = Array(tmp[:, 2:end])
  tair  = transpose(tair)

  # Read precipitation data

  str   = readline("$path/$file_tair")
  nsep  = length(matchall(r";", str))
  tmp   = CSV.read("$path/$file_prec", delim = ";", header = false,
                  dateformat="yyyy-mm-dd HH:MM", nullable = false, types = vcat(DateTime, repmat([Float64], nsep)))
  prec  = Array(tmp[:, 2:end])
  prec  = transpose(prec)

  # Read runoff data

  tmp   = CSV.read("$path/$file_q_obs", delim = ";", header = false,
                   dateformat="yyyy-mm-dd HH:MM", nullable = false, types = [DateTime, Float64])
  q_obs = Array(tmp[:, 2])

  q_obs[q_obs .== -999.0] = NaN

  # Read elevation band data

  frac = readdlm("$path/$file_frac")
  frac = squeeze(frac,2)

  # Get time data

  date = Array(tmp[:, 1])

  # Return data

  return date, tair, prec, q_obs, frac

end=#


# Load operational data from text files

function load_data(path, file_q_obs = "Q_obs.txt", file_tair = "Tair.txt",
                   file_prec = "Prec.txt", file_metadata = "metadata.txt")

  # Read air temperature data

  str   = readline("$path/$file_tair")
  nsep  = length(matchall(r";", str))
  tmp   = CSV.read("$path/$file_tair", delim = ";", header = false,
                   dateformat="yyyy-mm-dd HH:MM", nullable = false, types = vcat(DateTime, repmat([Float64], nsep)))
  tair  = Array(tmp[:, 2:end])
  tair  = transpose(tair)

  # Read precipitation data

  str   = readline("$path/$file_tair")
  nsep  = length(matchall(r";", str))
  tmp   = CSV.read("$path/$file_prec", delim = ";", header = false,
                  dateformat="yyyy-mm-dd HH:MM", nullable = false, types = vcat(DateTime, repmat([Float64], nsep)))
  prec  = Array(tmp[:, 2:end])
  prec  = transpose(prec)

  # Read runoff data

  tmp   = CSV.read("$path/$file_q_obs", delim = ";", header = false,
                   dateformat="yyyy-mm-dd HH:MM", nullable = false, types = [DateTime, Float64])
  q_obs = Array(tmp[:, 2])

  q_obs[q_obs .== -999.0] = NaN

  # Read metadata

  df_tmp = CSV.read("$path/$file_metadata", delim = ";", header = true)

  area = convert(Array{Float64,1}, df_tmp[:area])
  frac_area = area / sum(area)

  frac_glacier = convert(Array{Float64,1}, df_tmp[:lus_glacier]) / 100.0

  frac_lus = DataFrame()
  frac_lus[:glacier] = frac_area .* frac_glacier
  frac_lus[:open] = frac_area - frac_lus[:glacier]

  elev = convert(Array{Float64,1}, df_tmp[:elevation])

  # Get time data

  date = Array(tmp[:, 1])

  # Return data

  return date, tair, prec, q_obs, frac_lus, frac_area, elev

end


# Crop data from start to stop date

function crop_data(date, tair, prec, q_obs, date_start, date_stop)

  # Find indicies

  istart = find(date .== date_start)
  istop = find(date .== date_stop)

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

  istart = find(date .== date_start)

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

