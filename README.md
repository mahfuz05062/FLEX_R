# FLEX_R
R package for FLEX (Functional evaluation of experimental perturbations), a pipeline that utilizes several functional annotation resources (Complex, Pathway, GO BP) to establish reference standards for systematic evaluation of genome-wide CRISPR screens. FLEX can investigate both single mutant (fitness) or double mutant (genetic interaction) experiments. FLEX summarizes the functional signal captured by these experiments, both globally (gene-pair level) and locally (protein complex level, for example). FLEX highlights the functional information captured by an experiment (and/or scoring method) and also compare and contrast functional information from multiple experiments (or competing scoring schemes).

## Installation
FLEX_R is compatiable with R version >= 3.4.3 and doesn't have any dependency except for if someone wish to use the GIANT functional network as a standard ('org.Hs.eg.db' is needed in that case). 

FLEX can be installed directly from github using install_github('csbio/FLEX_R') and this requires the 'devtools' package.

## Qestions and Comments
Please email to rahma118@umn.edu or mahfuz05062@gmail.com if you have any questions, comments, potential bugs or suggestions!
