# Functions required for sectorring analyses
# Version 0.4 Li Yutze Jan. 21 2015

# Colony Analyze Main Fucntion
#   Take in a character vector of file names or a directory name
#   Output the analysis report for all files avaliable in that vector/dir
# Argumatents:
#   path            a charactor vector indicating the picture file paths or
#                   the directory that contains the pictures.
#   strain          a integer vector indicating how pictures are grouped into
#                   strains. If c(10, 15) is input, we consider the first 10
#                   strains are of a same strain, and the next 15 belongs to
#                   another. The vector can be named, indicating your strain
#                   accessions. (e.g. c(AB001 = 10, AB002 = 15))
#   sampling        Number of sample points extract from the selected ring
#   at              At what percent of radius taking the sample ring
#   report          a charactor vector of length 1
#                       * 'l' for generating one long but consise report for
#                               all pictures
#                       * 'd' for generating a relatively detailed report for
#                               each file
#                       * 'n' for no reports, default value
#                   This/ese report/s are useful for manually clean up analysis
#                   failures.
#   rec             T/F means Do/Do not leave a record file for the output,
#                   default value is T.
colony.analyze <- function(path = charactor(),
                           strain = NULL,
                           sampling = 500,
                           at = .85,
                           report = 'l',
                           rec = T) {
    # Check if files are avaliable
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
    if (report == T) report.dir <- './reports/';
    if (!file.exists(report.dir)) dir.create(report.dir);
    
    # Load packages
    library(jpeg);
    library(EBImage);
    library(reshape);
    
    # For each picture file, we extract a vecter of intense value and store
    # them into a list
    intensities <- list();
    
    # The factor vector for strains
    if (is.null(strain)) {
        strain <- factor(rep('#0', times = length(path)));
    } else {
        if (is.null(sn <- names(strain))) 
            sn <- paste('#', 1:length(strain), sep = '');
        strain <- factor(rep(sn, times = strain));
    }
    
    # Save the original plotting parameters
    ori.par <- par(no.readonly = T);
    
    # If the long report is required, set up plotting parameters
    if (report == 'l') {
        png(filename = paste(report.dir, 'report.png', sep = ''),
            width = 2.5, height = 0.5 + 1.1 * length(path), 
            units = 'in', res = 300);
        par(cex = 0.6, mai = c(0.05, 0.05, 0.05, 0.05), omi = c(0, 0, 0.5, 0));
        layout(mat = matrix(1:(2 * length(path)), ncol = 2, byrow = T),
               widths = c(1.1, 1.4));
    }
    
    # If the short reports are required, set up plotting parameters
    if (report == 'd') {
        # size for each report file is w x d = 4 x 5.5
        par(mai = c(0.05, 0.05, 0.05, 0.05), omi = c(0, 0, 0.5, 0));
        layout(mat = matrix(c(0, 1, 0,
                              2, 2, 2,
                              3, 3, 3), byrow = T, nrow = 3),
               widths = c(0.5, 3, 0.5), heights = c(3, 1, 1));
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
        
        # Get the extracted vecter
        window <- ceiling(pi * radius / sampling / 2);
        vect <- numeric();
        for (j in 1:sampling) {
            rad <- j * 2 * pi / sampling;
            s.cenx <- cenx + sampr * cos(rad);
            s.ceny <- ceny + sampr * sin(rad);
            vect <- c(vect, mean(pic[(s.cenx - window):(s.cenx + window),
                                     (s.ceny - window):(s.ceny + window)]));
            #vect <- c(vect, pic[s.cenx, s.ceny]);
        }
        intensities[['file']] <- vect;
        
        # Generate report
        if (report == 'l') {
            plot(0:1, 0:1, type = 'n', ann = F, frame = 1, axes = 0);
            rasterImage(drawCircle(pic, cenx, ceny, sampr, 
                1)[(cenx + radius):(cenx - radius),
                   (ceny + radius):(ceny - radius)],
                0, 0, 1, 1);
            plot(c(0, 2), c(0, 6), type = 'n', ann = F, frame = F, axes = 0);
            text(1, 5, paste('Strain:', strain[length(intensities)]), pos = 1);
            text(1, 4, paste('File Name:', file), pos = 1);
            text(1, 3, paste('Mean Intensity', round(mean(vect)), 3), pos = 1);
            rasterImage(t(matrix(rep(vect, times = 10), ncol = 10)),
                        0.1, 0.4, 1.9, 1.6);
        }
        if (report == 'd') {
            # not finished yet
        }
    }
    
    if (report == 'l') {
        title(main = 'FIGURE RECOGNITION REPORT', outer = T);
        dev.off();
    }
    par(ori.par);
    
    # Reshape the output
    intensities <- lapply(split(intensities, strain), as.data.frame);
    
    # Save the record
    if (rec == T) {
        dir.create('./rec');
        for (strain in names(intensities)){
            rec.name <- paste('./rec/', strain, '.csv', sep = '');
            write.csv(intensities[[strain]], rec.name);
        }
    }
    
    # return the intensities
    invisible(intensities);
}