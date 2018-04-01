library(sf)
library(sp)
library(zoo)
library(dplyr)
library(pbapply)
library(spatialEco)

library(Rcpp)
sourceCpp('LineSimplify.cpp')

# Get bookstein coordinates for bounding box
st_bookstein <- function(s, f = 1000){
  bbx = s %>% st_bbox()
  w = bbx[3] - bbx[1]
  h = bbx[4] - bbx[2]
  scale = f / max(w, h)

  center = bbx %>% st_as_sfc() %>% st_centroid()

  inline.book = ((st_geometry(s) - center) * scale) %>% st_sf()
}

st_fracdim <- function(line, samples){
  coords =  st_coordinates(line)[, 1:2]
  n = sapply(samples, function(X) get_ruler(coords, X))
  m <- lm(log(n)~log(samples))
  fdim = -1 * m$coefficients[2]
  return(fdim)
}

MHD <- function(line1, line2){
  pts1 <- line1 %>% st_cast("POINT")
  pts2 <- line2 %>% st_cast("POINT")

  d1 <- st_distance(pts1, line2)
  d2 <- st_distance(pts2, line1)

  d <- max(mean(d1), mean(d2))

  return(d)
}

resampleLine <- function(line, sampling) {
  npts <- as.numeric(round(st_length(line)/sampling))

  newline = line %>%
    st_line_sample(npts) %>%
    st_coordinates() %>%
    st_linestring() %>%
    st_sfc(crs = st_crs(line))

  return(newline)
}

# Get effective area for each point
getAreas <- function(coords) {
  get_areas(coords)
}

# Visvalingam-Whyatt simplification
simplifyVW <- function(line, area) {

  if (area <= 0)
    return(line)

  coords = line %>% st_cast("POINT") %>% st_coordinates()
  npts = nrow(coords)

  gencoords = simplify_vw(coords, area)
  npts_new = nrow(gencoords)

  while(npts_new != npts && npts_new > 2){
    npts = npts_new
    gencoords = simplify_vw(gencoords, area)
    npts_new = nrow(gencoords)
  }

  vline = gencoords %>%
    st_linestring() %>%
    st_sfc(crs = st_crs(line)) %>%
    st_sf()

  return(vline)

}

# Create regular grid, covering the sf object
createGrid <- function(line, res) {

  box <- line %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_buffer(res) %>% # to ensure that the grid covers the object
    st_bbox() %>%
    st_as_sfc()

  grid <- box %>%
    as('Spatial') %>%
    spsample(cellsize = res, type = 'regular') %>%
    SpatialPixels() %>%
    as('SpatialPolygons') %>%
    st_as_sf() %>%
    st_set_crs(st_crs(line)) %>%
    mutate(id = row_number())

  return(grid)
}

createFishnet <- function(line, res){
  box <- line %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_buffer(res) %>% # to ensure that the grid covers the object
    st_bbox()

  xmin = box[1]
  ymin = box[2]
  xmax = box[3]
  ymax = box[4]

  xseq = seq(xmin, xmax, res)
  yseq = seq(ymin, ymax, res)

  xlines = lapply(xseq, function(X) {
    matrix(c(X, ymin, X, ymax), ncol = 2, byrow = TRUE) %>% st_linestring()
  })

  ylines = lapply(yseq, function(Y) {
    matrix(c(xmin, Y, xmax, Y), ncol = 2, byrow = TRUE) %>% st_linestring()
  })

  fishnet = st_sfc(c(xlines, ylines), crs = st_crs(line))

  return(fishnet)

}

# Li-Openshaw simplification
simplifyLi <- function(line, res, mode = 'midpoint', cont = 0.25) { # cont param means contraction

  message('Simplifying with resolution = ', res, '\n')

  if (res <=0)
    return(line)

  grid <- createGrid(line, res)

  if(mode == 'cluster'){ # cluster mode

    points <- line %>%
      st_cast('POINT') %>%
      st_sf() %>%
      mutate(idp = row_number()) %>%
      st_intersection(grid) %>%
      arrange(idp)

    RLE <- rle(points$id)

    simplified <- points %>%
      mutate(groupid = rep(1:length(RLE$values), RLE$lengths)) %>%
      group_by(groupid) %>%
      summarise(geometry = st_union(.)) %>%
      st_centroid() %>%
      st_coordinates() %>%
      contract_edges(cont * res) %>%
      st_linestring() %>%
      st_sfc(crs = st_crs(line)) %>%
      st_sf()

  } else { # midpoint and intersection modes

    coords <- line %>%
      st_cast('POINT') %>%
      st_coordinates()

    npts <- nrow(coords)

    fishnet <- createFishnet(line, res)

    points <- st_intersection(line, fishnet) %>%
      st_cast('MULTIPOINT') %>%
      st_cast('POINT') %>%
      st_coordinates() %>%
      insert_points(coords, .) %>%
      st_linestring() %>%
      st_sfc(crs = st_crs(line)) %>%
      st_cast('POINT') %>%
      st_sf() %>%
      mutate(idp = row_number()) %>%
      arrange(idp)

    # these are the ordered intersection points
    intcoords.sorted <- points[fishnet, ] %>% st_coordinates()

    coords.cleared = intcoords.sorted

    ncoords <- nrow(coords.cleared)

    if(ncoords > 1){
      # define artifacts (similar X and Y)
      RLEx <- rle(intcoords.sorted[,1])
      RLEy <- rle(intcoords.sorted[,2])

      # clear artifacts
      coords.cleared <- data.frame(intcoords.sorted,
                                   xgr = rep(1:length(RLEx$values), RLEx$lengths),
                                   ygr = rep(1:length(RLEy$values), RLEy$lengths)) %>%
        group_by(xgr) %>% summarise(X = mean(X), Y = mean(Y), ygr = min(ygr)) %>%
        group_by(ygr) %>% summarise(X = mean(X), Y = mean(Y), xgr = min(xgr)) %>%
        select(X, Y) %>% as.matrix()
    }

    ncoords <- nrow(coords.cleared)

    if(mode == 'midpoint' && ncoords > 1){ # midpoint mode and at least two points
      xres = 0.5 * (coords.cleared[1:(ncoords-1), 1] + coords.cleared[2:ncoords, 1])
      yres = 0.5 * (coords.cleared[1:(ncoords-1), 2] + coords.cleared[2:ncoords, 2])
    } else { # intersection mode
      xres = coords.cleared[,1]
      yres = coords.cleared[,2]
    }

    # simplified includes first, SVO and last points
    simplified <- rbind(coords[1,],
                        cbind(xres, yres),
                        coords[npts, ]) %>%
      contract_edges(cont * res) %>%
      st_linestring() %>%
      st_sfc(crs = st_crs(line)) %>%
      st_sf()

  }

  return(simplified)
}

smoothLine <- function(line, sampling, sigma) {

  npts <- as.numeric(round(st_length(line)/sampling))

  samples = line %>%
    st_line_sample(npts)
  coords = samples %>%
    st_coordinates()

  x = coords[,1] # исходные координаты X
  y = coords[,2] # исходные координаты Y
  n = length(x)

  rpts = sigma %/% 2  # половинка ядра без остатка

  rx.start = rev(2 * x[1] - x[1:rpts])   # отражение X для rpts первых точек
  rx.end = rev(2 * x[n] - x[(n-rpts):n]) # отражение X для rpts последних точек

  ry.start = rev(2 * y[1] - y[1:rpts])   # отражение Y для rpts первых точек
  ry.end = rev(2 * y[n] - y[(n-rpts):n]) # отражение Y для rpts последних точек

  x <- c(rx.start, x, rx.end) # пристыковываем сначала и конца
  y <- c(ry.start, y, ry.end) # пристыковываем сначала и конца

  kernel = gaussian.kernel(sigma/8, sigma)[1,]

  xs <- rollapply(zoo(x), sigma,
                  function(z) weighted.mean(z, kernel),
                  align = 'center')
  ys <- rollapply(zoo(y), sigma,
                  function(z) weighted.mean(z, kernel),
                  align = 'center')

  line_s <- cbind(xs, ys) %>%
    st_linestring() %>%
    st_sfc(crs = st_crs(line))

  return(line_s)
}

fitSmoothed <- function(line, sline, sigmas, res){
  message('Fitting gaussian curves...')
  gaulines <- pblapply(sigmas,
                       function(X) smoothLine(line, res, X))


  sline_res <- resampleLine(sline, sampling)

  message('Estimating distances...')
  mhds <- pbsapply(gaulines,
                   function(X) MHD(X, sline_res))

  minsigma <- sigmas[which.min(mhds)]

  line_gaumin <- smoothLine(line, sampling, minsigma)

  return(list(line_gaumin, mhds))
}

fitSmoothedMulti <- function(line, slines, sampling, sigmas){
  nsigmas <- length(sigmas)
  nlines <- length(slines)

  fitteds <- lapply(slines, function(X) fitSmoothed(line, X, sigmas, sampling))

  numbers <- rep(1:nlines, rep(nsigmas, nlines))
  sigmass <- rep(sigmas, nlines)

  distances <- unlist(unlist(fitteds, recursive=FALSE)[2 * (1:nlines)])

  df <- data.frame(numbers, sigmass, distances)

  df$numbers <- as.factor(df$numbers)

  minvals <- df %>%
    group_by(numbers) %>%
    summarise(m = min(distances),
              s = sigmas[which.min(distances)])

  return(list(df, minvals))
}

simplifyMulti <- function(line, parameters, algorithm = 'dp') {
  slines = switch(algorithm,
                  'dp' = lapply(parameters, function(X) st_simplify(line,
                                                                    dTolerance = X,
                                                                    preserveTopology = TRUE)),
                  'lo' = lapply(parameters, function(X) simplifyLi(line, X)),
                  'vw' = lapply(parameters, function(X) simplifyVW(line, X)))
  return(slines)
}

getSmoothingSignature <- function(line, sampling, sigmas, parameters, algorithm='dp'){

  slines = simplifyMulti(line, parameters, algorithm)

  res = fitSmoothedMulti(line, slines, sampling, sigmas)

  levels(res[[1]]$numbers) = parameters

  res[[2]]$alg = algorithm

  return(res)

}
