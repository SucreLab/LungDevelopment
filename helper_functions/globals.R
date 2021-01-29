set.seed(42) # For reproducability

N_WORKERS <- 2
n_print <- 1:20
options(future.globals.maxSize=15*1024*1024^2)


saveTiff <- function(path, image, width = 5, height = 5, dpi = 300, units = "in"){
  if (!file.exists(path)){
    dir.create(dirname(path), showWarnings = FALSE)
  }

  if (Sys.info()["sysname"]=='Darwin'){
    # lzw doesn't work on mac with quartz
    tmp_path <- suppressWarnings(normalizePath(paste0(path, "_tmp.tiff")))
    out_path <- suppressWarnings(normalizePath(path))

    tiff(tmp_path, width = width, height = height, units = units, res = dpi, compression = "lzw", type = "quartz")
    print(image)
    dev.off()
    # requires imagemagick
    Sys.sleep(0.5)
    system(paste0("convert ", tmp_path, " -density ", dpi,  " -resize ", width * dpi, "x", height * dpi, "\\> -compress lzw ", out_path), ignore.stdout = TRUE)

    if (file.exists(out_path)) {
      #Delete file if it exists
      file.remove(tmp_path)
    }

  } else {
    tiff(path, width = width, height = height, units = units, res = dpi, compression = "lzw")
    print(image)
    dev.off()
  }

}

