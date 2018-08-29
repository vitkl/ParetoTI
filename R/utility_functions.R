##' Install python module py_pcha and numpy
##' @rdname install_py_pcha
##' @name install_py_pcha
##' @author Vitalii Kleshchevnikov
##' @details This is a helper function that install py_pcha and numpy python modules into conda environment. Unless you have a strong reason please use suggested defaults.
##' @param method paratemer for \code{\link[reticulate]{py_install}}. Default option should work in most cases. Use virtualenv if you don't want to install anaconda but virtualenv doesn't work on Windows.
##' @param conda paratemer for \code{\link[reticulate]{py_install}}. Default option should work in most cases.
##' @param python_version version to be installed into environment that is compatible with py_pcha module.
##' @param envname name of the conda enviroment where PCHA should be installed. If that enviroment doesn't exist it will be created. If it contains incorrect python_version the function will give an error.
##' @param env_dir directory in which to create python virtualenv when using that method. By default is NULL which results in ~/.virtualenvs.
##' @param overwrite_env It TRUE overwrites conda environment.
##' @details If installation fails with an error "Cannot fetch index base URL http://pypi.python.org/simple/" try this solution: "Older versions of pip and distribute default to http://pypi.python.org/simple, which no longer works. A solution is to install an up-to-date pip and distribute using pip install -i https://pypi.python.org/simple -U pip distribute into the virtual environment before running the rest of the build process."
##' @return path to python enviroment with py_pcha module installed
##' @export install_py_pcha
##' @seealso \code{\link{}}, \code{\link{}}
install_py_pcha = function(method = "auto", conda = "auto",
                           python_version = "python 2.7.10",
                           envname = "reticulate_PCHA",
                           env_dir = NULL,
                           overwrite_env = F) {
  packages = c("pip", "py_pcha", "numpy", "scipy", "datetime")
  if(method != "virtualenv") {
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
  } else {
    packages = c(python_version, packages)
    if(!is.null(env_dir)) system2(paste0("export WORKON_HOME=", env_dir))
    reticulate::py_install(packages = packages, envname = envname,
                           method = method, conda = conda)
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
