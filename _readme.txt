Code for the Comment for Judd and Su 2012, Econometrica

By: Fedor Iskhakov, John Rust and Bertel Schjerning

September 5, 2015

This directory contains the modified code from Su and Judd Econometrica paper
The intention of this directory is to carefully monitor all the changes we make to the original code 
with thorough documentation of each step in clean GIT repository

type: 'git log' in this directory for the log of changes to the original code

Su and Judd code is downloaded from
Che-Lin Su website on April 3, 2014
link:
http://faculty.chicagobooth.edu/che-lin.su/research/SuJudd_Code_ECMA.zip

To replicate the initial full set of results in the comment
(the published version of the comment contains limited results)
do:

Table 1: 
	Run run.m with no changes except perhaps changing the AMPL directory and link to KNITRO library
Table 2: 
	Run run.m after adjusting the size of the sample and turning off MPEC/ktrlink
Figure 1: 
	Run run_nBus.m with no changes to generate the data obtained for Figure 1. 
	Then run plot_cpu_vs_samplesize.m in nBus/ directory.

The folder nBus/ contains all the code for partial likelihood estimation, which required
certain edits of the original Judd-Su code.

The MPEC/AMPL code was updated after Che-Lin Su's recentering of the value function, so that logsum does not 
produce overflow/underflow errors which improved MPEC ability to converge for beta close to 1.
