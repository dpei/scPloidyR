#' @title Validate Required Fields
#' @description Check that an object contains all required named fields.
#' @param obj Named list or similar object.
#' @param fields Character vector of required field names.
#' @param object_name Character, descriptive name for error messages.
#' @keywords internal
validate_required_fields <- function(obj, fields, object_name = "object") {
  missing <- setdiff(fields, names(obj))
  if (length(missing) > 0) {
    stop(sprintf("%s missing required fields: %s", object_name, paste(missing, collapse = ", ")))
  }
}
