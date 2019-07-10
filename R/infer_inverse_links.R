infer_invlinks <- function(namelink) {
  if (namelink == "log") {
    invlink <- exp
  } else if (namelink == "identity") {
    invlink <- identity
  } else {
    invlink <- NULL
  }
  return(invlink)
}
