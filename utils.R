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