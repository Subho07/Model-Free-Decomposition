require(ggplot2)
require(ggthemes)
require(ggstance)
require(extrafont)
require(complexplus)
require(OpenImageR)
require(envi)
require(float)
library(caTools)
library(raster)
library(rgdal)
# library(rasterVis)
library(RColorBrewer)

#font_import()
#loadfonts()
##---------------------------------------------------------------------------------------
readComplexFloat <- function(f, size, endian){
  X2 <- readBin(f, double(), n = 2*size, size = 4, endian = endian) # caso 4
  X = complex(length.out = size)
  
  for(i in 1:size) {
    twoi = 2*i
    X[i] = complex(real = X2[twoi-1], imaginary=X2[twoi])
  }
  
  return(X)
}

myread.ENVI <- function (filename, headerfile = paste(filename, ".hdr", sep = ""))
{
  nCol <- nRow <- nBand <- data.type <- header.offset <- byte.order <- (-1)
  interleave = "bsq"
  if (!file.exists(headerfile))
    stop("read.ENVI: Could not open input header file: ",
         headerfile)
  Lines = read.table(headerfile, sep = "=", strip.white = TRUE,
                     row.names = NULL, as.is = TRUE, fill = TRUE)
  Fields = c("samples", "lines", "bands", "data type", "header offset",
             "interleave", "byte order")
  for (i in 1:nrow(Lines)) {
    Lab = tolower(Lines[i, 1])
    Lab = gsub("[ ]+", " ", Lab)
    j = match(Lab, Fields)
    Val = Lines[i, 2]
    if (length(j) == 1)
      switch(j, nCol <- as.integer(Val), nRow <- as.integer(Val),
             nBand <- as.integer(Val), data.type <- as.integer(Val),
             header.offset <- as.integer(Val), interleave <- gsub(" ",
                                                                  "", Val), byte.order <- as.integer(Val))
  }
  if (nCol <= 0 | nRow <= 0 | nBand <= 0)
    stop("read.ENVI: Error in input header file ", headerfile,
         " data sizes missing or incorrect", nRow, nCol, nBand)
  if (!(data.type %in% c(1, 2, 3, 4, 5, 6, 9, 12)))
    stop("read.ENVI: Error in input header file ", headerfile,
         " data type is missing, incorrect or unsupported ")
  ieee = if (.Platform$endian == "big")
    1
  else 0
  endian = if (ieee == byte.order | byte.order < 0)
    .Platform$endian
  else "swap"
  size = nRow * nCol * nBand
  if (!file.exists(filename))
    stop("read.ENVI: Could not open input file: ", filename)
  f = file(filename, "rb")
  if (header.offset > 0)
    readBin(f, raw(), n = header.offset)
  switch(data.type,
         X <- readBin(f, integer(), n = size, size = 1, signed = FALSE), # caso 1
         X <- readBin(f, integer(), n = size, size = 2, endian = endian), # caso 2
         X <- readBin(f, integer(), n = size, endian = endian), # caso 3
         X <- readBin(f, double(), n = size, size = 4, endian = endian), # caso 4
         X <- readBin(f, double(), n = size, endian = endian), # caso 5
         X <- readComplexFloat(f, size = size, endian = endian), # caso 6
         ,  # caso 7
         ,  # caso 8
         X <- readBin(f, complex(), n = size, endian = endian), # caso 9
         , # caso 10
         , # caso 11
         X <- readBin(f, integer(), n = size, size = 2, endian = endian, signed = FALSE) # caso 12
  )
  close(f)
  Fields = c("bil", "bip", "bsq")
  j = match(interleave, Fields)
  if (length(j) == 0)
    stop("read.ENVI: Error in input header file ", headerfile,
         " incorrect interleave type")
  switch(j, {
    dim(X) <- c(nCol, nBand, nRow)
    X <- aperm(X, c(3, 1, 2))
  }, {
    dim(X) <- c(nBand, nCol, nRow)
    X <- aperm(X, c(3, 2, 1))
  }, {
    dim(X) <- c(nCol, nRow, nBand)
    X <- aperm(X, c(2, 1, 3))
  })
  if (nBand == 1)
    dim(X) = c(nRow, nCol)
  return(X)
}

para01 <- function(x) {
  valores <- range(x)
  y = (x - valores[1]) / (valores[2] - valores[1])
  y
}


##----------------------------------------------------------------------------
# main_dir = 'D:\\27. Four_component_decomp_DB\\T3_RS2\\subset\\subset_0_of_T33\\T3\\'
main_dir = choose.dir()
main_dir = paste0(main_dir,"\\")

# Read coherence elements (T3)
C11 <- myread.ENVI(paste0(main_dir,"C11.bin"), 
                   paste0(main_dir,"C11.bin.hdr"))
C12_imag <- myread.ENVI(paste0(main_dir,"C12_imag.bin"), 
                        paste0(main_dir,"C12_imag.bin.hdr"))
C12_real <- myread.ENVI(paste0(main_dir,"C12_real.bin"), 
                        paste0(main_dir,"C12_real.bin.hdr"))

C22 <- myread.ENVI(paste0(main_dir,"C22.bin"), 
                   paste0(main_dir,"C22.bin.hdr"))

Nrw = nrow(C11)
Ncl = ncol(C11)

##----------------------------------------------------------------------------
# selection of window
wsi = 7; # window size

kernel = matrix(1, nrow = wsi, ncol = wsi) / (wsi*wsi)


c11s = convolution(C11, kernel, "same")
C12_imag_filt = convolution(C12_imag, kernel, "same")
C12_real_filt = convolution(C12_real, kernel, "same")

c12s = C12_real_filt + C12_imag_filt*1i
c21s = Conj(c12s)

c22s = convolution(C22, kernel, "same")


det_C2 = c11s*c22s - c12s*c21s

trace_C2 = c11s + c22s
m1 = Re(sqrt(1-(4*(det_C2/(trace_C2^2)))))

s0 = c11s + c22s
s1 = c11s - c22s
s2 = (c12s + c21s)

chi_in = -45 # right circular

if (chi_in >= 0){
  s3 = (1i*(c12s - c21s)) # The sign is according to RC or LC sign !!
}

if (chi_in < 0){
  s3 = -(1i*(c12s - c21s)) # The sign is according to RC or LC sign !!
}

SC = (s0-s3)/2
OC = (s0+s3)/2

s0_c = trace_C2
dop_c = m1

h = OC - SC

val = ((dop_c*s0_c*h))/((SC*OC + (dop_c**2)*(s0_c**2)))

theta_c = atan(val1)*(180/pi) # separation for surface and dbl

pv_c = (1-dop_c)*s0_c
res_pow = s0_c - pv_c
ps_c = (res_pow/2)*(1+sin((2*theta_c*(pi/180))))
pd_c = (res_pow/2)*(1-sin((2*theta_c*(pi/180))))

#-----------------------------------------------------------------------------
# Save in envi format
dir.create(file.path(main_dir, 'MF3CC_Rext'), showWarnings = FALSE)
save_dir = paste0(main_dir,'MF3CC_Rext\\')

r <- raster(theta_c)
r <- writeRaster(r, filename=paste0(save_dir,'theta_c.envi'), overwrite=TRUE)

r <- raster(pv_c)
r <- writeRaster(r, filename=paste0(save_dir,'pv_c.envi'), overwrite=TRUE)

r <- raster(ps_c)
r <- writeRaster(r, filename=paste0(save_dir,'ps_c.envi'), overwrite=TRUE)

r <- raster(pd_c)
r <- writeRaster(r, filename=paste0(save_dir,'pd_c.envi'), overwrite=TRUE)
