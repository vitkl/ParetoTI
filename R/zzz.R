# reference to py_pcha (will be initialized in .onLoad)
.onLoad = function(libname, pkgname) {
  # if the python environment "reticulate_PCHA" installed using
  # install_py_pcha() is availlable - use it. If not - search other python binaries.
  #condas = conda_list(conda = "auto")
  #if("reticulate_PCHA" %in% condas$name) {
  #  use_python(condas[condas$name == "reticulate_PCHA", "python"],
  #             required = FALSE)
  #  #use_condaenv("reticulate_PCHA", conda = "auto", required = FALSE)
  #}

  # add py_pcha to package namespace
  assign("py_PCHA",
         reticulate::import("py_pcha.PCHA", delay_load = TRUE),
         envir = parent.env(environment()))
  assign("py_furthest_sum",
         reticulate::import("py_pcha.furthest_sum", delay_load = TRUE),
         envir = parent.env(environment()))
}
