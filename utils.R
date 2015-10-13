# Functions required for sectorring analyses
# Version 0.4 Li Yutze Jan. 24 2015

# Find Left and Right Wings ---------------------------------------------------
# 	Given index of current sample point and the length of the sample vector, 
#	these functions returns the left and right wings (indecs)
left.wing <- function(cur, len, res) {
	wing <- (cur - res + 1):cur - 1;
	wing[wing <= 0] <- wing[wing <= 0] + len;
	wing
}
right.wing <- function(cur, len, res) {
	wing <- cur:(cur + res - 1) + 1;
	wing[wing > len] <- wing[wing > len] - len;
	wing
}
#------------------------------------------------------------------------------

# Count the Switch Stripe Number on a State Vector ----------------------------
count.strips <- function(state, len) {
	count <- 0L;
	i <- 1L;
	while (i <= len) {
		j <- i;
		if (!state[i]) {
			while (j <= len & !state[j]) j <- j + 1;
			j <- j - 1;
			if (!(state[left.wing(i, len, 1)] + state[right.wing(j, len, 1)]))
				count <- count + 1;
		}
		i <- j + 1;
	}
	floor(count / 2);
}
#------------------------------------------------------------------------------

# Analysis report generation --------------------------------------------------
analyze.report <- function(count, dir = NULL) {
    # Make saving path
    if (is.null(dir)) dir <- './report';
    if (!file.exists(dir)) dir.create(dir);
    
    # Setting up ploting parameters
    ori.par <- par(no.readonly = T);
    png(paste(dir, '/analysis_report.png', sep = ''),
        width = 2.5, height = 0.5 + 1.1 * length(count), 
        units = 'in', res = 300);
    par(cex = 0.4, mai = c(0.05, 0.05, 0.05, 0.05), omi = c(0, 0, 0.5, 0));
    layout(mat = matrix(1:(2 * length(count)), ncol = 2, byrow = T),
           widths = c(1.1, 1.4))
    mai <- c(0.05, 0.05, 0.05, 0.05);
    
    # Adding the information of this strain into the report
    rightlim <- max(sapply(count, max)) + 2;
    upperlim <- max(sapply(count, function(x) {max(table(x))})) + 2;
    for (strain in names(count)) {
        par(mai = c(0.3, 0.3, 0.02, 0.02))
        hist(count[[strain]],
             xlim = c(0, rightlim), breaks = c(0:rightlim),
             ylim = c(0, upperlim),
             main = '', xlab = 'Switch Event', ylab = 'Count',
             tck = -0.01, mgp = c(0.5, 0, 0),
             cex.axis = 0.4, cex.lab = 0.5);
        par(mai = mai);
        plot(c(0, 2), c(1, 5), type = 'n', ann = F, axes = F);
        text(1, 4, paste('Strain:', strain));
        text(1, 3, paste('# of inputs:', length(count[[strain]])));
        text(1, 2, paste('Mean of Count:', round(mean(count[[strain]]), 2)));
    }
    title(main = 'ANALYSIS REPORT', outer = T);
    # par(ori.par);
    dev.off();
}
#------------------------------------------------------------------------------