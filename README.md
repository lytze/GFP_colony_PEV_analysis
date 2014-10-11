Project Read-Me
===============

## Introduction

This project is part of my project in the [*He Xiangwei* lab](http://lsi.zju.edu.cn/redir.php?catalog_id=13701) (Life Sciences Institute, Zhejiang Univ.). The major aim of this project is to evaluate the PEV (position effect variation) pattern complexity based on captured pictures of *gfp*-marked yeast single colonies under fluorenscent microscope. Six sample figures are presented in the `sample_figure` folder. You can access the figures to get more understanding (more info in the 'Sample Figures' section).

## Stages and Functionalities

The processing pipline can be divided into two steps:

* Find the location and region of the colony 
* Evaluate the complexity

In the first stage, raw pictures files are inputed to address the locations and regions of the colonies. Typically, colonies are presented as bright green circular objects in dark background in our pictures (just like what a fluorenscent microcopic figure should look like). In our early versions we only deal with 'friendly' pictures containing only one colony in each (sample figure 1 ~ 5). Next we might look forward to address more colonoes in one time (sample figure 6).

The second phase is mostly statistics. Three values are extracted from the pictures: 1) the overall intensity, 2) the horizontal complexity and 3) the vertical complexity. We define the horizontal and vertical direction based on the polar coordination, where horizontal is the direction across circumferences and vertical is that on the radius.

## Sample Figures

Six sample figures are presented in the `sample_figure` folder. Figure 1 to 5 contain single colonies and figure 6 contains 2. Yeast strains (obtained from Yamamoto lab), JT630 and JT634, are theoretically refered as of higher and lower variation rate respectively. In our sample figures, 1 to 4, as well as 6 are pictures of JT630, and figure 5 is a picture of JT634.

You can access the genetic infomation through this paper [Hexanucleotide motifs mediate recruitment of the RNA elimination machinery to silent meiotic genes., *Yamashita A, et al*](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3352096/). Notice that in our pictures, we used the haploid strains instead the original diploid strain.
