##' Install python module py_pcha and numpy
##' @rdname install_py_pcha
##' @name install_py_pcha
##' @author Vitalii Kleshchevnikov
##' @description \code{install_py_pcha} is a helper function that installs py_pcha and numpy python modules into conda environment. Unless you have a strong reason please use suggested defaults. Alternatively, use pip install --user py_pcha numpy scipy datetime.
##' @description \code{select_conda} is a helper function that tells R to use python conda environment created specifically for this package using \code{install_py_pcha}. This environment is used by default whe package is loaded. Set options(disable_ParetoTI_envname = TRUE) if you want to use other python installation. Set options(ParetoTI_envname = "other_env_name") if you want to use other conda environment.
##' @param method paratemer for \code{\link[reticulate]{py_install}}. Default option should work in most cases. Use virtualenv if you don't want to install anaconda but virtualenv doesn't work on Windows. For conda method to work anaconda should be installed as described here: https://conda.io/docs/user-guide/install/index.html.
##' @param conda paratemer for \code{\link[reticulate]{py_install}}. Default option should work in most cases.
##' @param python_version version to be installed into environment that is compatible with py_pcha module.
##' @param envname name of the conda enviroment where PCHA should be installed. If that enviroment doesn't exist it will be created. If it contains incorrect python_version the function will give an error.
##' @param overwrite_env It TRUE overwrites conda environment.
##' @param extra_packages specify python libraries not needed for ParetoTI to work but desirable in the same conda environment
##' @param packages specify python libraries needed for ParetoTI to work. Normally do not need changing.
##' export PYTHONUSERBASE=/some_dir/python_libs/
##' python -m pip install --user -i https://pypi.python.org/simple -U pip distribute
##' python -m pip install --user -i https://pypi.python.org/simple --upgrade pip setuptools wheel
##' @return path to python enviroment with py_pcha module installed
##' @export install_py_pcha
##' @export select_conda
##' @examples
##' \dontrun{
##' install_py_pcha()
##' install_py_pcha(extra_packages = c("tensorflow", "pandas", "keras", "h5py",
##'                                    "geosketch", "pydot", "sklearn", "umap-learn"))
##'
##' ## See for installation details
##' }
install_py_pcha = function(method = "auto", conda = "auto",
                           python_version = "python 3.7.3",
                           envname = "reticulate_PCHA",
                           overwrite_env = F, extra_packages = character(0),
                           packages = c("pip", "py_pcha", "numpy", "scipy", "datetime")) {
  packages = c(packages, extra_packages)
  if(method == "virtualenv") {
    reticulate::py_install(packages = packages, envname = envname,
                           method = method, conda = conda)
  } else {
    condas = conda_list(conda = conda)
    if(envname %in% condas$name & !overwrite_env) {
      python = condas[condas$name == envname,"python"]
      python_v = system2(python, args = "-V", stdout = T, stderr = T)
      correct_python = grepl(python_version, python_v, ignore.case = T)
      if(!correct_python) {
        stop(paste0("conda environment: ",envname,
                    " exists but does not contain ",
                    python_version, " (",python_v,")"))
      }
    } else {
      conda_remove(envname, packages = NULL, conda = conda)
      conda_create(envname, packages = python_version, conda = conda)
    }
    reticulate::py_install(packages = packages, envname = envname,
                           method = method, conda = conda, pip = T)
    conda_python(envname, conda = conda)
  }
}

.py_pcha_installed = function() {
  # check if python package is availlable and give a helpful error message
  err = tryCatch(py_PCHA$PCHA(matrix(1:6, 2, 3),
                              noc = as.integer(2)), error = function(e) e)
  if(!is.null(err$message)) {
    if(grepl("Python", err$message)){
      stop(paste0(err$message,
                  ", please use install_py_pcha() to install"))
    }
  }
}

select_conda = function(conda = "auto", envname = "reticulate_PCHA"){
  condas = tryCatch(conda_list(conda = conda), error = function(e) e)
  if(envname %in% condas$name) {
    use_python(condas[condas$name == envname, "python"],
               required = FALSE) # suggest to use this envir but do not force.
    use_condaenv(envname, conda = conda, required = FALSE)
  }
}
