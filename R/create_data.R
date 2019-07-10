create_data <- function(y, data, params) {
    D <- vector('list', length=length(params))
    if (!is.null(data)){
      y <- formula(paste(y, "~ 1"))
      y <- model.response(stats::model.frame(y, data = data))
      print(params)
      for (i in 1:length(params)){
        D[[i]] <- stats::model.matrix(params[[i]], data)
      }
    } else {
      for (i in 1:length(params)){
        if (length(as.character(params[[i]])) == 2 &
            as.character(params[[i]])[2] == "1"){
          D[[i]] <- matrix(ncol = 1, rep(1, length(y)))
        }
        else {
          D[[i]] <- stats::model.matrix(params[[i]])
        }
      }
    }
    names(D) <- names(params)
    # Check for missing values
    n <- length(y)
    if (length(na.omit(y)) < n){
      stop("missing values are not allowed in y")
    }
    for (i in 1:length(D)){
      if (nrow(na.omit(D[[i]])) < n){
        stop("missing values are not allowed in the covariates")
      }
    }
    return(list(y = y, D = D))
  }
