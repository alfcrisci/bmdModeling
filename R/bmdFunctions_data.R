#' Load file in proast data format
#' @param filename character string, path to the file with data to be loaded
#' line 1 with title of the data; line 2 number of variables; line 3 response
#' type for each variable; line 4 variable names and from line 5 onwards the data
#' @param separator character, defines the separator used in the data to be loaded 
#' @return List with info, title of the data; nvar, number of variables;
#' varnames, names of the variables; data, loaded data; dtype, number indicating
#' the type of response
#' @export
f.scan <- function (filename = NULL, separator = "") {
	
	tit <- scan(filename, "", n = 1)
	nvar <- scan(filename, 0, n = 1, sep = "\n", skip = 1)
	dtype <- numeric(0)
	dtype <- scan(filename, 0, n = nvar, skip = 2)
	descr <- character(0)
	descr <- scan(filename, "", skip = 3, n = nvar)
	data <- read.table(filename, skip = 4, sep = separator, 
			row.names = NULL)
	print(names(data))
	
	errorNamesAllocation <- tryCatch({
		names(data) <- descr
		},  error = function(err) {
		return(err)
	})

	if(inherits(errorNamesAllocation, "error"))
		if(grepl("'names' attribute .+ must be the same length as the vector .+]", 
			errorNamesAllocation$message))
		stop(paste(errorNamesAllocation$message, "Please consider changing the separator."))	else
		stop(errorNamesAllocation$message)
	
	return(list(info = tit, nvar = nvar, varnames = descr, 
					data = data, dtype = dtype))
	
	
}

#' Remove NAs from a given data frame for all listed column indices
#' @param variableIndices numeric vector, column indices for variables to check
#' for missing values
#' @param originalData data frame, data from which missing values should be 
#' excluded 
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return data frame based on original data, all rows with missing values for 
#' any of the non-zero variables are excluded
#' @export
f.remove.NAs <- function (variableIndices, originalData, track = FALSE) {
	
	if (track) 
		print("f.remove.NAs")
	
	variableIndices <- variableIndices[(variableIndices != 0)]  
	isMissing <- is.na(originalData[, variableIndices, drop = FALSE])
	rowsMissing <- apply(isMissing, 1, any)
	
	if(any(isMissing)){
		
		warning("There are ", sum(rowsMissing), " rows with missing values removed from the data")
		
	}
	
	completeData <- originalData[!rowsMissing, ]
	
	if (track)	print("f.remove.NAs:  END")
	
	return(completeData)
	
}

#' Check whether the data frame columns are numeric and can be log-transformed
#' @param dataFrame data frame for which the columns class and values will be checked
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return boolean, TRUE if the vector is not numeric, contains negative values
#' or contains zero values only; if so warning message is printed 
#' @export
f.check.nonneg.num <- function (dataFrame, track = FALSE) {
	
	if (track) 
		print("f.check.nonneg.num")
	
	
	if(is.null(dim(dataFrame))){
		
		isProblem <- (!is.numeric(dataFrame) || any(dataFrame < 0) || !any(dataFrame > 0))
		
	} else {
		
		isProblem <- apply(dataFrame, 2, function(column){
					
					(!is.numeric(column) || any(column < 0) || !any(column > 0))
					
				})
		
		
	}
	
	if(any(isProblem)){
		
		warning("Not all values are numerical, 
						variable contains negative values or variable contains zeros only \n")
		
	}
	
	
	if (track) 
		print("f.check.nonneg.num:   END")
	
	return(isProblem)
	
}