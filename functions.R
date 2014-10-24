# Functions required for sectorring analyses
# Li Yutze Oct 2014


# Read the colony PEV image (CPI) into R
import.cpi <- function(name, col.path = c(1.0, 1.0, 1.0)) {
    # name: character element of length 1, identicating the path of the picture
    # col: numeric vector of length 3, given the grabbing color
    
    if (!is.character(name)) {
        stop('String of file location is expected, but not a string is inputed');
    }
    else {
        if (length(col.path) != 3) {
            stop('Color path is expected, but vector of length other than 3 is inputed');
        } 
        if (sum(col.path > 1.0 | col.path < 0.0) != 0) {
            stop('Color path is expected, but number out of expected range is inputed');
        }
        if (sum(col.path == 1.0) < 1) {
            stop('Saturated color path is expected, but not saturated color is inputed');
        }
    }
    # Error cases:
    ### input file name not a valid character vector
    ### input color path not of valid length, range, or not saturated
    
    require(jpeg);
    load <- readJPEG(name);
    cpi <- list(
        identity = 'CPI',
        call = match.call(),
        raw = (load[, , 1] * col.path[1] + load[, , 2] * col.path[2] + load[, , 3] * col.path[3]) / 
            sqrt(sum(col.path))
    );
    # Dot product, to grab the color conponent required
    
    return(cpi);
}

# Get the background threshold and output the binary evaluation
threshold.cpi <- function(cpi) {
    # cpi: colony PEV image processing list
    
    current.mean <- mean(cpi$raw); 
    escape <- 0;
    counter <- 10;
    while (!escape & counter) {
        last.mean <- current.mean;
        current.mean <- mean(cpi$raw[cpi$raw <= 2 * last.mean]);
        # narrow down the sample range to [0, 2 * last.mean]
        counter <- counter - 1;
        if (abs(last.mean - current.mean) <= 0.01) {
            escape <- 1;
        }
    }
    # Iterative find the mean of background value. If not converge after 10
    # cycles, escape from the loops
    
    cpi$threshold <- current.mean * 2;
    cpi$binary.map <- cpi$raw >= cpi$threshold;
    
    # compress the binary map and create a copy
    cby <- nrow(cpi$binary.map) %/% 100; 
    compressed <- 
        cpi$binary.map[seq(from = 1, to = nrow(cpi$binary.map), by = cby),
                       seq(from = 1, to = ncol(cpi$binary.map), by = cby)];
    
    cpi$compressed <- compressed;
    cpi$comp.by <- cby;
    return(cpi);
}

# Get the center and the radius of the fitted circle
circle.cpi <- function(cpi) {
    # cpi: colony PEV image processing list
    
    require(reshape);
    stick.map <- melt(cpi$compressed);
    stick.map <- stick.map[stick.map[, 3] == T, 1:2];
    
    center.row <- mean(stick.map[, 1]);
    center.col <- mean(stick.map[, 2]);
    # center location is the center of weight of the binary map
    
    area <- nrow(stick.map) # = sum(cpi$binary.map), count all TURE dots
    radius <- sqrt(area / pi);
    
    cpi$stick <- stick.map;
    cpi$cen.row <- center.row * cpi$comp.by;
    cpi$cen.col <- center.col * cpi$comp.by;
    cpi$radius <- radius * cpi$comp.by;
    
    cpi
}



cpi <- import.cpi('sample_figure/3.jpg', c(0, 0, 1));
cpi <- threshold.cpi(cpi);
cpi <- circle.cpi(cpi);

png('report.png', )