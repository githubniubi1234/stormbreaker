sb_dependency_components <- function(dependency) {
  if (is.na(dependency) || !nzchar(dependency)) {
    return(list(r_packages = character(), python_modules = character()))
  }
  parts <- trimws(strsplit(dependency, " + ", fixed = TRUE)[[1L]])
  python_start <- which(startsWith(parts, "Python "))
  if (length(python_start) == 0L) {
    return(list(r_packages = parts, python_modules = character()))
  }
  python_start <- python_start[[1L]]
  r_packages <- if (python_start > 1L) {
    parts[seq_len(python_start - 1L)]
  } else {
    character()
  }
  python_modules <- parts[seq.int(python_start, length(parts))]
  python_modules[[1L]] <- sub("^Python ", "", python_modules[[1L]])
  list(r_packages = r_packages, python_modules = python_modules)
}

sb_r_package_installed <- function(package) {
  length(find.package(package, quiet = TRUE)) == 1L
}

sb_dependency_status_row <- function(dependency) {
  components <- sb_dependency_components(dependency)
  missing <- components$r_packages[!vapply(
    components$r_packages,
    sb_r_package_installed,
    logical(1)
  )]
  collapse <- function(x) {
    if (length(x) == 0L) NA_character_ else paste(x, collapse = ", ")
  }
  status <- if (length(components$r_packages) == 0L &&
    length(components$python_modules) == 0L) {
    "no_optional_dependency"
  } else if (length(missing) > 0L) {
    "missing_r_packages"
  } else if (length(components$python_modules) > 0L) {
    "python_not_checked"
  } else {
    "r_packages_installed"
  }
  data.frame(
    dependency = dependency,
    r_packages = collapse(components$r_packages),
    python_modules = collapse(components$python_modules),
    missing_r_packages = collapse(missing),
    dependency_status = status,
    stringsAsFactors = FALSE
  )
}

sb_add_dependency_status <- function(registry) {
  dependencies <- unique(registry$optional_dependency)
  status <- do.call(
    rbind,
    lapply(dependencies, sb_dependency_status_row)
  )
  status <- status[
    match(registry$optional_dependency, status$dependency),
    c(
      "r_packages", "python_modules", "missing_r_packages",
      "dependency_status"
    ),
    drop = FALSE
  ]
  rownames(status) <- NULL
  cbind(registry, status)
}
