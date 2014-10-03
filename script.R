library(jpeg);
library(reshape);

# Read sample figure
pic <- readJPEG('Sample_Figure/Dis 33-3.jpg');
x.samp <- seq(from = 1, to = nrow(pic), by = nrow(pic) / 100);
y.samp <- seq(from = 1, to = ncol(pic), by = nrow(pic) / 100);
pic.samp <- pic[x.samp, y.samp, 3];

# Split out x-y coordinate
pic.melt <- melt(pic.samp);
colnames(pic.melt) <- c('x', 'y', 'g');

# Find threshold using loops
samp <- list(pic.melt$g);
mn <- numeric();
counter <- 1;
repeat {
    mn[counter] <- mean(samp[[counter]]);
    counter; mn[counter];
    
    if (counter > 2)
        if (mn[counter - 1] - mn[counter] < 0.01 | counter > 15)
            break;
    
    counter <- counter + 1;
    samp[[counter]] <- samp[[counter - 1]][samp[[counter - 1]] <= mn[counter - 1] * 2];
}

# plot graped 'bright points'
thresh <- max(samp[[counter]]);

# Cleaning up noise
fin.samp <- pic.melt[pic.melt$g > thresh,]
fin.mat <- matrix(F, length(x.samp), length(y.samp));
for (i in 1:nrow(fin.samp)) {
    fin.mat[fin.samp[i, 1], fin.samp[i, 2]] <- T;
}
new.mat <- matrix(F, length(x.samp), length(y.samp));
for (i in 1:nrow(fin.mat))
    for (j in 1:ncol(fin.mat)) {
        chi <- c(F, F, F, F)
        if (fin.mat[i, j] == T) {
            if (i > 1) chi[1] <- fin.mat[i - 1, j];
            if (i < nrow(fin.mat)) chi[2] <- fin.mat[i + 1, j]
            if (j > 1) chi[3] <- fin.mat[i, j - 1]
            if (j < fin.samp[i, 2]) chi[4] <- fin.mat[i, j + 1]
        }
        if (sum(chi) >= 3) new.mat[i, j] <- T;
    }

# Calculate the mass center of the grabbed colony picture
mass <- list(x = 0, y = 0, w = 0);
for (i in 1:nrow(new.mat))
    for (j in 1:ncol(new.mat)) {
        if (new.mat[i, j]) {
            sw <- mass$w + new.mat[i, j];
            mass$x <- (mass$w * mass$x + new.mat[i, j] * i) / sw;
            mass$y <- (mass$w * mass$y + new.mat[i, j] * j) / sw;
            mass$w <- sw;
        }
    }

r <- sqrt(mass$w / pi);

drw.circle <- function(r, x, y) {
    cir.x <- numeric();
    cir.y <- numeric();
    theta <- seq(from = 0, to = 2 * pi, by = 0.1);
    for (i in 1:length(theta)) {
        cir.x[i] <- x + cos(theta[i]) * r;
        cir.y[i] <- y + sin(theta[i]) * r;
    }
    cir.x <- c(cir.x, cir.x[1]);
    cir.y <- c(cir.y, cir.y[1]);
    list(x = cir.x, y = cir.y);
}

circle <- drw.circle(r, mass$x, mass$y)
for (i in 1:nrow(pic.melt)) {
    pic.melt$grab[i] <- new.mat[pic.melt$x[i], pic.melt$y[i]];
}

png('fig1.png', width = 400, height = 500);
par(mar = c(2, 2, 1, 1));
plot(pic.melt$x, pic.melt$y, type = 'n', xlim = c(0, 100), asp = 1);
points(pic.melt$x[pic.melt$g > thresh & pic.melt$grab], pic.melt$y[pic.melt$g > thresh & pic.melt$grab],
       col = 'green4', pch = '.');
points(pic.melt$x[pic.melt$g > thresh & !pic.melt$grab], pic.melt$y[pic.melt$g > thresh & !pic.melt$grab],
       col = 'blue2', pch = 4, cex = 0.5);
points(mass$x, mass$y, col = 'red');
lines(circle$x, circle$y, col = 'red');
dev.off();
