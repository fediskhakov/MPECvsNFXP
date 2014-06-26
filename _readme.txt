Code for Comment for Judd and Su 2012, Econometrica

By: Fedor Iskhakov, John Rust and Bertel Schjerning

This directory contains the modified code from Judd and Su Econometrica paper
The intention of this directory is to carefully monitor all the changes we make to the original code
with thorough documentation of each step in clean GIT repository

type: 'git log' in this directory for the log of changes to the original code

Judd and Su code is downloaded from
Che-Lin Su website on April 3, 2014
link:
http://faculty.chicagobooth.edu/che-lin.su/research/SuJudd_Code_ECMA.zip

To replicate the results in the comment, do:

Table 1 and cdf plots: 
	Run run.m with no changes except perhaps changing the AMPL directory and link to KNITRO library
Table 2: 
	Run run.m after adjusting the size of the sample and turning off MPEC/ktrlink
Table 3: 
	Run run_hidim.m with no changed for smaller grids for bothe AMPL and NFXP-NK
	Change the gridvec vector and potentially the output file suffic param.out_suff for other runs
Runtime plot:
	Run plot_cpu_vs_N.m in hidim/ directory after perhaps updating the data to the data obtained for table 3.

Folder high_dim/ contains all the code for higher dimensional partial likelihood estimation, which required
substantial edits of the original Judd-Su code.  The folder contains .diff files which show exactly what
changes had been made to the AMPL code.
