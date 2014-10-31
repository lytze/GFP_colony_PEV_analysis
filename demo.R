# Some demos for the 'functions.R' script

source('functions.R');
col.path <- c(0, 0, 1); # all figures are gfp marked
div <- 180;
lay.1 <- c(0.3, 0.5, 0.7, 0.9);
lay.2 <- c(0.35, 0.5, 0.65, 0.8, 0.95);

fig1.1 <- analyze.cpi('./sample_figure/1.jpg', col.path, lay.1, div);
fig2.1 <- analyze.cpi('./sample_figure/2.jpg', col.path, lay.1, div);
fig3.1 <- analyze.cpi('./sample_figure/3.jpg', col.path, lay.1, div);
fig4.1 <- analyze.cpi('./sample_figure/4.jpg', col.path, lay.1, div);
fig5.1 <- analyze.cpi('./sample_figure/5.jpg', col.path, lay.1, div);
fig6.1 <- analyze.cpi('./sample_figure/6.jpg', col.path, lay.1, div);
fig3.2 <- analyze.cpi('./sample_figure/3.jpg', col.path, lay.2, div);
fig4.2 <- analyze.cpi('./sample_figure/4.jpg', col.path, lay.2, div);

fig1.1$report();
fig2.1$report();
fig3.1$report();
fig4.1$report();
fig5.1$report();
fig6.1$report();
fig3.2$report();
fig4.2$report();

fig1.1$short(T);