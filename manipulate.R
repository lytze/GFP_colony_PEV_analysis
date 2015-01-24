# Functions required for sectorring analyses
# Version 0.4 Li Yutze Jan. 24 2015

#------------------------------------------------------------------------------
# Intensity Sampling Vector Clean up Function
#   After reading the report, users can manually judge the recognation of which
#   picture(s) is/are not well determined, so that they can clean these result
#   terms up
# Arguments:
#   data        The returned value of `colony.extract()`, in which you think
#               the shape of some pictures are not well determined.
#   index       A integer/character/logical vector contains the indeces of the
#               terms to be cleaned up
leave <- function(data, index) {
    # Check if all parameters are legal
    if (!is.list(data) | !all(sapply(data, is.data.frame)))
        stop('Data input error');
    strain <- factor(rep(names(data), times = sapply(data, ncol)));
    fname <- unlist(sapply(data, colnames));
    if (!is.atomic(index) | (
            (is.logical(index) & length(index) != length(strain)) |
            (is.numeric(index) & (any(index < 0) | max(index) > length(strain))) |
            (is.character(index) & !all(index %in% fname))
        )) stop('Index input error');

    # Return data with selected terms removed
    data <- as.data.frame(data);
    names(data) = fname;
    if (is.logical(index)) {
        index <- !index;
    } else if (is.numeric(index)) {
        index <- -index;
    } else {
        index <- !(fname %in% index);
    }
    data <- data[, index];
    strain <- strain[index];
    invisible(split(data, strain));
}
#------------------------------------------------------------------------------


# Recover data from Records ---------------------------------------------------
#   In the `extract()` function, we leave a record in the ./rec/ directory, and
#   for the next time analysis, we do not need to run the function again, and
#   we can recover the result from these records
# Arguments:
#   dir         The path to the rec/ directory, must to be end with 'rec' or
#               'rec/'
#------------------------------------------------------------------------------
