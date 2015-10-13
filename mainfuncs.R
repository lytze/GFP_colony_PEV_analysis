# Functions required for sectorring analyses
# Version 0.4 Li Yutze Jan. 24 2015

# loading utils.R for these functions -----------------------------------------
source('~/CodeRMain/GFP_colony_PEV_analysis/utils.R');
#------------------------------------------------------------------------------


# Colony Intensity Sampling Vector Extraction ---------------------------------
#   Take in a character vector of file names or a directory name
#   Output the analysis report for all files avaliable in that vector/dir
# Arguments:
#   path            a charactor vector indicating the picture file paths or
#                   the directory that contains the pictures.
#   strain          a integer vector indicating how pictures are grouped into
#                   strains. If c(10, 15) is input, we consider the first 10
#                   strains are of a same strain, and the next 15 belongs to
#                   another. The vector can be named, indicating your strain
#                   accessions. (e.g. c(AB001 = 10, AB002 = 15))
#   sampling        Number of sample points extract from the selected ring.
#   at              At what percent of radius taking the sample ring.
#   report          T/F means D0/Do not nead a report.
#   rec             T/F means Do/Do not leave a record file for the output,
#                   default value is T.
extract <- function(path = charactor(),
                    strain = NULL,
                    sampling = 500,
                    at = .85,
                    report = T,
                    rec = T) {
    # Check if files are avaliable
    pwd <- getwd()
    if (!all(file.exists(path))) stop('File/Directory not found');
    if (all(file.info(path)$isdir)) {
        setwd(path);
        path <- dir();
    }
    if (!all(grepl('.jpg$', path) -> index)) {
        message('Non-jpg files input, leave out following file(s):');
        message(paste('\t', path[!index], '\n', sep = ''));
        path <- path[index];
    }
    sampling <- as.integer(sampling);
    if (report) {
        report.dir <- './report/';
        if (!file.exists(report.dir)) dir.create(report.dir);
    }
    
    # Load packages
    library(jpeg);
    library(EBImage);
    library(reshape);
    
    # For each picture file, we extract a vecter of intense value and store
    # them into a list
    intensities <- list();
    
    # The factor vector for strains
    if (is.null(strain)) {
        strain <- factor(rep('Strain-0', times = length(path)));
    } else {
        if (is.null(sn <- names(strain))) 
            sn <- paste('Strain-', 1:length(strain), sep = '');
        strain <- factor(rep(sn, times = strain));
    }
    
    # If a report is required, set up plotting parameters
    if (report) {
        ori.par <- par(no.readonly = T);
        png(filename = paste(report.dir, 'extract_report.png', sep = ''),
            width = 2.5, height = 0.5 + 1.1 * length(path), 
            units = 'in', res = 300);
        par(cex = 0.6, mai = c(0.05, 0.05, 0.05, 0.05), omi = c(0, 0, 0.5, 0));
        layout(mat = matrix(1:(2 * length(path)), ncol = 2, byrow = T),
               widths = c(1.1, 1.4));
    }
    
    # Iteration through files 
    for (file in path) {
        # Load in picture, extract G-path
        pic <- readJPEG(file)[, , 2];
        
        # Find appropriate threshold to isolate the colony's shape
        # Using the Otsu method in EBImage package to find the thresholding
        # value of intensity
        holding <- otsu(pic);
        # and get the binary of the picture
        binpic <- pic > holding;
        
        # Get the center, the radius and the ring of sampling
        tall <- melt(binpic, varnames = c('x', 'y'));
        cenx <- mean(tall$x[tall$value == T]); 
        ceny <- mean(tall$y[tall$value == T]);
        radius <- sqrt(sum(binpic) / pi);
        sampr <- radius * at;
        
        window <- ceiling(pi * radius / sampling / 2);
        escape <- {
            c(cenx + (1 + radius + window) * c(-1, 1),
                        ceny + (1 + radius + window) * c(-1, 1)) *
            c(1, -1, 1, -1) +
            c(0, nrow(pic), 0, ncol(pic))
        } < 0
        
        if (!any(escape)) {
            # Get the extracted vecter
            vect <- numeric();
            for (j in 1:sampling) {
                rad <- j * 2 * pi / sampling;
                s.cenx <- cenx + sampr * cos(rad);
                s.ceny <- ceny + sampr * sin(rad);
                vect <- c(vect, mean(pic[(s.cenx - window):(s.cenx + window),
                                         (s.ceny - window):(s.ceny + window)]));
                #vect <- c(vect, pic[s.cenx, s.ceny]);
            }
            intensities[[file]] <- vect;
            
            # Generate report
            if (report) {
                plot(0:1, 0:1, type = 'n', ann = F, frame = 1, axes = 0);
                rasterImage(drawCircle(pic, cenx, ceny, sampr, 
                                       1)[(cenx + radius):(cenx - radius),
                                          (ceny + radius):(ceny - radius)],
                            0, 0, 1, 1);
                plot(c(0, 2), c(0, 6), type = 'n', ann = F, frame = F, axes = 0);
                text(1, 6, paste('#', length(intensities)), pos = 1);
                text(1, 5, paste('Strain:', strain[length(intensities)]), pos = 1);
                text(1, 4, paste('File Name:', file), pos = 1);
                text(1, 3, paste('Mean Intensity', round(mean(vect), 3)), pos = 1);
                rasterImage(t(matrix(rep(vect, times = 10), ncol = 10)),
                            0.1, 0.4, 1.9, 1.6);
            }
            
        }
        
        if (report) {
            title(main = 'FIGURE RECOGNITION REPORT', outer = T);
            dev.off();
        }
        
        # Reshape the output
        intensities <- lapply(split(intensities, strain), as.data.frame);
        
        # Save the record
        if (rec) {
            if (!file.exists('./rec')) dir.create('./rec');
            for (strain in names(intensities)){
                rec.name <- paste('./rec/', strain, '.csv', sep = '');
                write.csv(intensities[[strain]], rec.name);
            }
        }
    }
    
    setwd(pwd)
    # return the intensities
    invisible(intensities);
}
#------------------------------------------------------------------------------


# Main Analysis Function ------------------------------------------------------
#   We get a list of data frames by different strains, then in this function
#   we analyze the extracted intensity vectors to determine the veriegation
#   complexities of these strains
# Arguments:
#   list        The list returned by the `extract()` function
#   res         The resolusion to determine peaks
#   rec         Do/not nead to leave a record
analyze <- function(lis, res = 10, rec = T) {
    # Lengths of the sample vectors should be the same
    len <- nrow(lis[[1]]);

    # The result is a list of vectors
    count <- list();

    # Go through strains
    for (strain in names(lis)) {
        count[[strain]] <- numeric();
        # For each picture in this strain, we check the left [res] and right 
        # [res] points for each sample point, to determine if these two groups
        # have significant difference via t-test
        for (vect in lis[[strain]]) {
            state <- integer();
            for (i in 1:len) {
                left <- vect[left.wing(i, len, res)];
                right <- vect[right.wing(i, len, res)];
                conf <- t.test(left, right, conf.level = 0.99)$conf;
                # State is the condition of the curve at this point, where -1
                # means dicreasing, +1 means increasing, and 0 means no signi-
                # ficent change
                state <- c(state, sum(c(-(conf > 0), conf < 0)) / 2);
            }
            # If a region of 0s are winged by -1s and 1s on other side, then
            # this should be a peak, and this number / 2 should be the number
            # of switched
            count[[strain]] <- c(count[[strain]], count.strips(state, len));
        }
    }
    
    # Call report function in utils.R to generate the report
    analyze.report(count);

    # Save the record
    text <- sapply(count, paste, collapse = ', ');
    text <- paste(names(text), text, sep = ': ');
    if (!file.exists('./rec')) dir.create('./rec');
    writeLines(text, './rec/analysis_rec.txt');

    # Print a short report to the consle
    cat(paste(names(count), ': ', round(sapply(count, mean), 2), ' with ',
              sapply(count, length), ' inputs\n', sep = ''), sep = '');

    invisible(count);
}
#------------------------------------------------------------------------------