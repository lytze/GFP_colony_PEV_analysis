# Functions required for sectorring analyses
# Version 0.3 Li Yutze Oct 31 2014

analyze.cpi <- function(name, col.path = c(1.0, 1.0, 1.0), lay = c(0.3, 0.6, 0.9), div = 90) {
    require(jpeg);
    require(reshape);
    require(stringr);
    require(plotrix);
    
    call <- match.call(); # record the function call, for reporting
    filen <- str_extract(name, '[^/]*$');
    filen <- gsub('.[^.]*$', '', filen);
    
    img <- readJPEG(name);
    img = (img[, , 1] * col.path[1] + img[, , 2] * col.path[2] + img[, , 3] * col.path[3]) / 
        sqrt(sum(col.path));
    tall <- melt(img, varnames = c('y', 'x')); 
    tall$y <- -tall$y # cast the matrix data into tall list
    
    ras <- as.raster(img); # make the raw image a raster, for future ploting
    den <- density(tall$value, n = 200); # record the density of color values, for future ploting
    
    # calculate the threshold for background color
    current.mean <- mean(tall$value); # set current mean value
    escape <- 0; # loop escape mark
    counter <- 10; # count the loop number, escape the loop when loop over 10 times
    while (!escape & counter) {
        last.mean <- current.mean;
        current.mean <- mean(tall$value[tall$value <= 2 * last.mean]);
        # narrow down the sample range to [0, 2 * last.mean]
        counter <- counter - 1;
        if (abs(last.mean - current.mean) <= 0.005) {
            escape <- 1;
        }
    }
    threshold <- current.mean * 2;
    deter <- tall$value > threshold;
    
    # denoise
    # scan the tall list, for each entry check the 11 * 11 square area
    # if over 30% points are in the list, mark all points in the list to be not noise
    # if less than 30% points are in the list, delet all
    # scaning processes should skip all entries marked to be not noise
    
    # evaluate the circular ragion
    cen.x <- mean(tall$x[deter]);
    cen.y <- mean(tall$y[deter]);
    radius <- sqrt(sum(deter) / pi);
    
    # get overall intensity
    intens <- mean(tall$value[deter]);
    
    # get hierarchy verctors
    rads <- seq(from = 0, to = 2 * pi, length = div + 1);
    rads <- rads[1:div];
    cx <- round(sapply(lay, function(x) {cen.x + cos(rads) * x * radius}), 0);
    cy <- round(sapply(lay, function(x) {cen.y + sin(rads) * x * radius}), 0);
    vec <- matrix(ncol = div, nrow = length(lay));
    for (i in 1:length(lay)) {
        vec[i, ] <- mapply(function(x, y) {
            mean(c(img[-y + -3:3, x], img[-y, x + -3:3]));
        }, cx[, i], cy[, i]);
        vec[i, ] <- (vec[i, ] - mean(vec[i, ]));
    }
    
    # get horizontal variation complexity
    vec.stp <- vec[, 1:div] - vec[, c(div, 1:(div-1))];
    hori <- mean(apply(vec.stp, 2, sd));
    
    # get vertical variation complexity
    m <- as.matrix(dist(vec));
    vert <- mean(m[1:(length(lay) - 1), 2:length(lay)]);
    
    # remove unessesary variables
    rm(current.mean, last.mean, escape, counter, 
       tall, deter, m, i, rads, vec.sum);
    
    # generate short report
    short <- function(pr = F) {
        if (pr) {
            cat('Overall Intensity: ', intens, '\n', sep = '');
            cat('Horizontal Variation Complexity: ', hori, '\n', sep = '');
            cat('Vertical Variation Complexity: ', vert, '\n', sep = '');
        }
        invisible(list(intensity = intens, horizontal = hori, vertical = vert));
    }
    
    # give processing log
    log <- function(pr = F) {
        if (pr) {
            cat('Identifier: ', filen, '\n', sep = '');
            cat('Threshold: ', threshold, '\n', sep = '');
            cat('Evaluating Colony Shape:\n\tcenter: (x, y) = (',
                cen.x,' ,', cen.y, ')\tradius: ', radius, '\n', sep = '');
        }
        invisible(list(call = call, threshold = threshold,
                       xcenter = cen.x, ycenter = cen.y, radius = radius));
    }
    
    # build report
    report <- function() {
        ori.par <- par(no.readonly = T);
        opt <- paste(paste(lay * 100, collapse = '-'), '-', div, sep = '');
        png(paste('report-', filen, '-', opt, '.png', sep = ''), 
            width = 600, height = 900);
        # plot the density and threshold
        par(fig = c(0.1, 0.9, 0.85, 1), mar = c(0, 0, 0, 0));
        plot(c(1, 7), c(1, 7), type = 'n', ann = F, axes = F);
        text(4, 6, paste('REPORT #', filen), cex = 2);
        text(3, 4:2, c('Overall Intensity:',
                       'Horizontal Variation Complexity:',
                       'Vertical Variation Complexity:'), cex = 1);
        text(5.3, 4:2, round(c(intens, hori, vert), 5), cex = 1);
        # plot the intensity distribution and the threshold
        par(fig = c(0.1, 0.9, 0.62, 0.86), mar = c(5, 4, 4, 4) + 0.1, new = T); 
        plot(den, ann= F, axes = F, col = 'blue4');
        abline(v = threshold, col = 'red2', lty = 2);
        axis(side = 2, cex = 0.5, tck = -0.02, cex = 0.6);
        axis(side = 1, cex = 0.5, tck = -0.02, cex = 0.6);
        title(main = 'Intensity Distribution',
              xlab = 'Intensity', ylab = 'Density', cex = 0.6);
        legend('bottom', legend = c('density', 'threshold'), 
               lty = c(1, 2), col = c('blue4', 'red2'), cex = 0.8);
        # show the colony shape
        shift <- 0.24 / ncol(ras) * nrow(ras);
        par(fig = c(0.14, 0.86, 0.4 - shift, 0.4 + shift), mar = c(0, 0, 0, 0), new = T);
        plot(c(1, ncol(ras)), c(-1, -nrow(ras)), type = 'n', ann = F, axes = F);
        rasterImage(ras, 1, -1, ncol(ras), -nrow(ras));
        points(cen.x, cen.y, col = 'red2', pch = 4);
        draw.circle(cen.x, cen.y, radius, border = 'red2', lty = 1);
        for (i in 1:length(lay)) {
            #draw.circle(cen.x, cen.y, radius * i, border = 'orange3', lty = 2);
            points(cx[, i], cy[, i], col = 'orange3', pch = '.')
        }
        text(c(cen.x, cen.x), c(cen.y, cen.y + radius), c('center', 'bondary'),
             col = 'red', cex = 0.8, pos = 3, offset = 0.4);
        text(rep(cen.x, length(lay)), cen.y + radius * lay, paste(lay * 100, '%'),
                 col = 'orange', cex = 0.8, pos = 3, offset = 0.4);
        # plot the sample vecters
        par(fig = c(0.1, 0.9, 0.03, 0.2), mar = c(2, 2, 3, 2), new = T);
        image(t(vec), ann = F, axes = F);
        axis(side = 2, at = seq(from = 0, to = 1, length = length(lay)),
             labels = paste(lay * 100, '%'), cex.axis = 0.8, las = 2, tck = 0, 
             col.axis = 'orange3',line = -0.5, lwd = 0);
        dev.off();
        par(ori.par);
    }
    
    return(list(short = short,
                log = log,
                report = report));
}