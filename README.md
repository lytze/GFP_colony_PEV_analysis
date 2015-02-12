Project Read-Me
===============

## Background

This project is part of my project in the [*He Xiangwei* lab](http://lsi.zju.edu.cn/redir.php?catalog_id=13701) (Life Sciences Institute, Zhejiang University). The major objective of this project is to evaluate the PEV (position effect variegation pattern complexity based on captured pictures of *gfp*-marked yeast single colonies under fluorenscent microscope.

The PEV, namely the <u>P</u>osition <u>E</u>ffect <u>V</u>ariegation is an important __epigenetic__ phenomenon, in which the genotype of a group keeps the same while the phenotype of them varies. The shifting of phenotypes are usually bi-stable, which means the gene of interest prefers its current on/off state, in other word, when the cell express the gene, its descendants are more likely to express that, but for some particular cases, some the daughter cells varies the state, and for cells currently silent the gene, their daughter cells are more likely to be silent.

For yeast strains conveys reporter genes that are epigenetically unstable, usually those located at the boundaries of special chromatin structures, the expression of the reporter gene will show a sectored variegation due to the growing and expanding of the colony, just like the colony in the following figures. You can [check the wikipedia page](http://en.wikipedia.org/wiki/Position_effect) on this topic for more information.

![PEV colony](https://github.com/lytze/GFP_colony_PEV_analysis/blob/master/sample_figure/1.jpg?raw=true)

## Why _gfp_ in This Case of Analysis

Usually geneticists use the _ade6_ gene as a marker for PEV (and as well for other uses). Deletion or defect on the _ade6_ gene will make the colony a red color when they brought up in low-adenine medium. So if EPV events are introduced to the _ade6_ reporter, we will find red-and-white sectorring variegations on the colony.

But in this case study, we use the _gfp_ gene as the reporter for two reasons:

1. Silent of _ade6_ will slow down the growing of cells, so the colonies with PEV are always irregular, which makes it difficult grabing the shape with programs. But _gfp_ do not affect the growing significantly, so we can easily assuming the colonies are just circular.
2. _gfp_ marked colonies are photographed with fluorenscent stereomicroscopes. The fluorenscent imaging condition makes it handy for image processing, where only one color path are involved.

## Functionalities

Functions for main analytic works are located in the `mainfuncs.R` script, in which

* `extract()` is used to extract sampling data from the raw pictures
* `analyze()` is used to analyze the sampling data generated by `extract()`

To `source('mainfuncs.R')` will also import the scripts in `utils.R`, which contains functions that main functions relies on.

And in the `manipulate.R` script, there are functions dealing with record files left by main functions.

A standard working pipeline should be like:

	# '~/data/colonies/' contains 50 pictures of strain A and 50 pictures of strain B
	source('mainfuncs.R')
	strains <- c(A = 50, B = 50) # tell functions how to divide files in to groups
	data <- extract(path = '~/data/colonies/', strain = strains) # this might take a while
	analyze(data)
	# A: 5.34 with 49 inputs
	# B: 2.11 with 44 inputs

And meanwhile a report figure will be generated. (Note the output in this sample does not represent any real data)

Between `extract()` and `analyze()` you can check the report manually that if the program wrongly judged the bondary of any colony. Or if you don't need this step, you can:

	analyze(extract(path = '~/data/colonies/', strain = strains, report = F))
	# A: 5.34 with 49 inputs
	# B: 2.11 with 44 inputs

Explanations for arguments of these functions (and most of other functions in `utils.R` and `manipulate.R`) are listed as comments right before the definitions of them. You can read the script for that information.

## Sample Figures

Eight sample figures are presented in the `sample_figures\` folder. They are photos of JT630, JT634 (obtained from Yamamoto lab) and their derived strains. Figure 1 to 4 are considered to be the more complex pattern group while 5 to 8 are of the more stable group.

You can access the genetic infomation of JT630 and JT634 through this paper [Hexanucleotide motifs mediate recruitment of the RNA elimination machinery to silent meiotic genes., *Yamashita A, et al*](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3352096/). Notice that in our project, we used the haploid strains instead diploid strains in the original paper.