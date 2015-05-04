.onUnload <- function (libpath) {
  library.dynam.unload("covR", libpath)
}
