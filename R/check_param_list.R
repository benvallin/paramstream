#' Check parameter list validity
#'
#' @param param_list named list of parameter values. Each element should be a parameter name - value pair. Values should be atomic vectors or `NULL`.
#' @param req_lgl character vector of required parameter names whose values should be logical vectors in `param_list`, or `NULL` (default).
#' @param req_num character vector of required parameter names whose values should be numeric vectors in `param_list`, or `NULL` (default).
#' @param req_chr character vector of required parameter names whose values should be character vectors in `param_list`, or `NULL` (default).
#' @param req_fac character vector of required parameter names whose values should be factors in `param_list`, or `NULL` (default).
#' @param opt_lgl character vector of optional parameter names whose values should be logical vectors in `param_list`, or `NULL` (default).
#' @param opt_num character vector of optional parameter names whose values should be numeric vectors in `param_list`, or `NULL` (default).
#' @param opt_chr character vector of optional parameter names whose values should be character vectors in `param_list`, or `NULL` (default).
#' @param opt_fac character vector of optional parameter names whose values should be factors in `param_list`, or `NULL` (default).
#' @param deps named list of parameter name dependencies. Each element should be a parameter name - required parameter name(s) pair. Names should evaluate to `TRUE` or `FALSE` in `param_list`; values should be character vectors. Default is `NULL`.
#' @param df_col_lgl character vector of parameter names which should be names of logical columns in `df`. Default is `NULL`.
#' @param df_col_num character vector of parameter names which should be names of numeric columns in `df`. Default is `NULL`.
#' @param df_col_chr character vector of parameter names which should be names of character columns in `df`. Default is `NULL`.
#' @param df_col_fac character vector of parameter names which should be names of factor columns in `df`. Default is `NULL`.
#' @param df data.frame with column names provided in `df_col_*`. Default is `NULL`.
#' @param return_filt_vars logical vector of length 1 indicating if the names of df filter variables should be returned. Filter variables are defined as parameter names in `param_list` not provided in `req_*`, `opt_*`, `df_col_*`, or `deps`. Default is `TRUE`.
#'
#' @return character vector of the names of df filter variables if `return_filt_vars` is `TRUE`; otherwise `NULL`, invisibly.
#' @export
#'
#' @examples
#' # Define analysis parameters list
#' param_list <- list(input_count = "tpm_lengthScaledTPM",
#'                    model_formula = "~ line + treatment",
#'                    cell_type = "neuron",
#'                    line = c("B856", "B156", "B067"),
#'                    protein_coding_only = TRUE,
#'                    min_cnt_excl = 0,
#'                    min_freq_incl = 0.2,
#'                    padj = 0.05)
#'
#' # Construct metadata to be checked against param_list
#' metadata <- data.frame(line = c(rep("B856", 4), rep("B156", 4), rep("B067", 4)),
#'                        cell_type = rep(c("neuron", "astrocyte"), 6),
#'                        treatment = rep(c(rep("untreated", 2), rep("treated", 2)), 3))
#'
#' # Check parameters list validity
#' check_param_list(param_list,
#'                  req_lgl = "protein_coding_only",
#'                  req_num = c("min_cnt_excl", "min_freq_incl", "padj"),
#'                  req_chr = c("input_count", "model_formula"),
#'                  df = metadata,
#'                  return_filt_vars = TRUE)
#'
check_param_list <- function(param_list,
                             req_lgl = NULL,
                             req_num = NULL,
                             req_chr = NULL,
                             req_fac = NULL,
                             opt_lgl = NULL,
                             opt_num = NULL,
                             opt_chr = NULL,
                             opt_fac = NULL,
                             deps = NULL,
                             df_col_lgl = NULL,
                             df_col_num = NULL,
                             df_col_chr = NULL,
                             df_col_fac = NULL,
                             df = NULL,
                             return_filt_vars = TRUE
) {

  # Check that param_list is a named list
  is_named_list <- function(x) {
    is.list(x) &&
      !is.data.frame(x) &&
      !is.null(names(x)) &&
      !anyNA(names(x)) &&
      all(nzchar(names(x)))
  }

  if(!is_named_list(param_list)) {

    stop("param_list is not a named list",
         call. = FALSE)

  }

  # Check that req_* / opt_* / df_col_* arguments are character vectors or NULL
  chr_args <- c("req_lgl", "req_num", "req_chr", "req_fac",
                "opt_lgl", "opt_num", "opt_chr", "opt_fac",
                "df_col_lgl", "df_col_num", "df_col_chr", "df_col_fac")

  invalid <- vapply(X = mget(chr_args),
                    FUN = function(x) !(is.character(x) || is.null(x)),
                    FUN.VALUE = logical(1))

  invalid <- chr_args[invalid]

  if(length(invalid) > 0L) {

    stop("The following arguments are not character vectors or NULL: ", paste0(invalid, collapse = ", "),
         call. = FALSE)

  }

  # Check that param_list entries are unique
  n_param_list <- length(names(param_list))
  n_unique_param_list <- length(unique(names(param_list)))

  if(n_param_list != n_unique_param_list) {

    non_unique <- names(which(table(names(param_list)) > 1))

    stop("The following param_list entries are not unique: ", paste0(non_unique, collapse = ", "),
         call. = FALSE)

  } else {

    message("All param_list entries are unique!")

  }

  # Combine required param_list entries
  req <- c(req_lgl, req_num, req_chr, req_fac)

  # Define reserved param_list entries
  res <- unique(c(req_lgl, req_num, req_chr, req_fac,
                  opt_lgl, opt_num, opt_chr, opt_fac,
                  df_col_lgl, df_col_num, df_col_chr, df_col_fac,
                  names(deps), unlist(deps, use.names = FALSE)))

  # Check that no entry is declared with more than one type
  typed <- c(unique(c(req_lgl, opt_lgl, df_col_lgl)),
             unique(c(req_num, opt_num, df_col_num)),
             unique(c(req_chr, opt_chr, df_col_chr)),
             unique(c(req_fac, opt_fac, df_col_fac)))

  typed <- typed[!is.na(typed)]

  clash <- unique(typed[duplicated(typed)])

  if(length(clash) > 0L) {

    stop("The following entries are declared with more than one type: ", paste0(clash, collapse = ", "),
         call. = FALSE)

  }

  # Define reserved entries (required + optional) included in user's param_list
  incl <- function(x) names(param_list)[names(param_list) %in% x]

  incl_req <- incl(req)
  incl_req_lgl <- incl(req_lgl)
  incl_req_num <- incl(req_num)
  incl_req_chr <- incl(req_chr)
  incl_req_fac <- incl(req_fac)
  incl_opt_lgl <- incl(opt_lgl)
  incl_opt_num <- incl(opt_num)
  incl_opt_chr <- incl(opt_chr)
  incl_opt_fac <- incl(opt_fac)
  incl_df_col_lgl <- incl(df_col_lgl)
  incl_df_col_num <- incl(df_col_num)
  incl_df_col_chr <- incl(df_col_chr)
  incl_df_col_fac <- incl(df_col_fac)
  incl_res <- incl(res)

  # Check presence of mandatory entries in param_list
  missing_req <- setdiff(req, incl_req)
  missing_req <- missing_req[!is.na(missing_req)]

  if(length(missing_req) > 0L) {

    stop("The following mandatory entries are missing from param_list: ", paste0(missing_req, collapse = ", "),
         call. = FALSE)

  } else {

    message("All required entries are present in param_list!")

  }

  # Check that param_list entries are of expected type
  check_type <- function(x, type) {

    fun <- switch(EXPR = type,
                  logical = "is.logical",
                  numeric = "is.numeric",
                  character = "is.character",
                  factor = "is.factor")

    fun <- match.fun(fun)

    not <- x[vapply(X = x,
                    FUN = function(y) !fun(param_list[[y]]),
                    FUN.VALUE = logical(1))]

    if(length(not) > 0L) {

      stop("The following param_list entries are not of the expected ", type, " type: ", paste0(not, collapse = ", "),
           call. = FALSE)

    }

  }

  check_type(x = c(incl_req_lgl, incl_opt_lgl, incl_df_col_lgl), type = "logical")

  check_type(x = c(incl_req_num, incl_opt_num, incl_df_col_num), type = "numeric")

  check_type(x = c(incl_req_chr, incl_opt_chr, incl_df_col_chr), type = "character")

  check_type(x = c(incl_req_fac, incl_opt_fac, incl_df_col_fac), type = "factor")

  message("All req / opt / df_col param_list entries are of the expected type!")

  # Check dependencies for reserved param_list entries
  if(!is.null(deps)) {

    # Check that deps is a named list
    if(!is_named_list(deps)) {

      stop("deps is not a named list",
           call. = FALSE)

    }

    # Check that deps values are character vectors
    invalid_deps <- names(deps)[!vapply(X = deps,
                                        FUN = is.character,
                                        FUN.VALUE = logical(1))]

    if(length(invalid_deps) > 0L) {

      stop("The deps values for the following entries are not character vectors: ",
           paste0(invalid_deps, collapse = ", "),
           call. = FALSE)

    }

    for(i in names(deps)) {

      # Warn if deps name is not a param_list entry
      if(!i %in% names(param_list)) {

        warning(i, " is included in deps but not in param_list",
                call. = FALSE)

        next

      }

      # Check that deps name is a single TRUE or FALSE in param_list
      if(!is.logical(param_list[[i]]) || length(param_list[[i]]) != 1L || is.na(param_list[[i]])) {

        stop("param_list entry ", i, " is used in deps and must be a single TRUE or FALSE",
             call. = FALSE)

      }

      # Check deps values if deps name is TRUE in param_list...
      if(param_list[[i]]) {

        # Check that deps values are entries in param_list
        missing_deps_val <- setdiff(deps[[i]], names(param_list))

        if(length(missing_deps_val) > 0L) {

          stop("The following deps values for entry ", i, " are missing from param_list: ",
               paste0(missing_deps_val, collapse = ", "),
               call. = FALSE)

        } else {

          message("All deps values for entry ", i, " are present in param_list!")

        }

      }

    }

  }

  # Check param_list entries expected to be names of df columns if df is not NULL
  if(!is.null(df)) {

    # Check that df is a data.frame
    if(!is.data.frame(df)) {

      stop("df is not a data.frame",
           call. = FALSE)

    }

    # For param_list entries expected to be names of df columns, check that:
    # - param_list entries are column names of df
    # - the corresponding columns are of the expected type
    check_type_in_df <- function(x, type) {

      fun <- switch(EXPR = type,
                    logical = "is.logical",
                    numeric = "is.numeric",
                    character = "is.character",
                    factor = "is.factor")

      fun <- match.fun(fun)

      not <- x[vapply(X = x,
                      FUN = function(y) !((y %in% colnames(df)) && fun(df[[y]])),
                      FUN.VALUE = logical(1))]

      if(length(not) > 0L) {

        stop("The following param_list entries are not names of ", type, " columns in df: ", paste0(not, collapse = ", "),
             call. = FALSE)

      }

    }

    check_type_in_df(x = incl_df_col_lgl, type = "logical")

    check_type_in_df(x = incl_df_col_num, type = "numeric")

    check_type_in_df(x = incl_df_col_chr, type = "character")

    check_type_in_df(x = incl_df_col_fac, type = "factor")

    message("All df_col param_list entries are names of columns of the expected type in df!")

  } else {

    # Stop if df_col_* values were passed without a df
    incl_df_col <- c(incl_df_col_lgl, incl_df_col_num, incl_df_col_chr, incl_df_col_fac)

    if(length(incl_df_col) > 0L) {

      stop("The following param_list entries are expected to be df columns but no df was provided: ",
           paste0(incl_df_col, collapse = ", "),
           call. = FALSE)

    }

  }

  # Define expected filter vars
  # => param_list entries that are not reserved
  exp_filt_vars <- setdiff(names(param_list), incl_res)

  if(is.null(df) &&
     length(exp_filt_vars) > 0) {

    stop("The following param_list entries are expected to be filter variables but no df was provided: ",
         paste0(exp_filt_vars, collapse = ", "),
         call. = FALSE)

  }

  # Define actual filter vars
  # => Expected filter vars which are column names of df
  filt_vars <- exp_filt_vars[exp_filt_vars %in% colnames(df)]

  # Define problematic filter vars
  # => Expected filter vars which are not column names of df
  prob_filt_vars <- exp_filt_vars[!exp_filt_vars %in% filt_vars]

  if(length(prob_filt_vars) > 0) {

    stop("The following param_list entries are not column names of df: ",
         paste0(prob_filt_vars, collapse = ", "),
         call. = FALSE)

  }

  if(length(filt_vars) > 0L) {

    message("The following param_list entries are expected to be filter variables in df: ",
            paste0(filt_vars, collapse = ", "))

  }

  # Return filter vars if requested
  if(return_filt_vars) {

    filt_vars

  }

}
