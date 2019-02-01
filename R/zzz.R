# reference to py_pcha (will be initialized in .onLoad)
.onLoad = function(libname, pkgname) {
  # if the python environment "reticulate_PCHA" installed using
  # install_py_pcha() is availlable - use it. If not - search other python binaries.
  if(getOption("disable_ParetoTI_envname") == FALSE ||
     is.null(getOption("disable_ParetoTI_envname"))){
    envname = "reticulate_PCHA"
    if(!is.null(getOption("ParetoTI_envname"))) envname = getOption("ParetoTI_envname")
    select_conda(conda = "auto", envname = envname)
  }


  # add py_pcha to package namespace
  assign("py_PCHA",
         reticulate::import("py_pcha.PCHA", delay_load = TRUE),
         envir = parent.env(environment()))
  assign("py_furthest_sum",
         reticulate::import("py_pcha.furthest_sum", delay_load = TRUE),
         envir = parent.env(environment()))

  # add geometric sketch to namespace
  assign("geosketch",
         reticulate::import("geosketch", delay_load = TRUE),
         envir = parent.env(environment()))

  # add facebook pca to namespace
  assign("fbpca",
         reticulate::import("fbpca", delay_load = TRUE),
         envir = parent.env(environment()))
}
