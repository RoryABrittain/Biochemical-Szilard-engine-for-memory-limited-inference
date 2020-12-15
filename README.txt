Biochemical Szilard engine for memory limited inference

A few general points:

	The 'figures.nb' Mathematica notebook should be fairly self-explanatory. Simply evaluating each cell in turn should recalculate the data and produce all of the figures.

	All of the C programs need to be compiled with 'functions.c'.

	The C files can be compiled from within 'figures.nb' assuming you have gcc installed.

	The C programs are designed assuming the standard output stream (the output of printf) is being redirected to a file. They also output information about the state of the program to the standard error stream.

	The Mathematica code to compile and run the C programs assumes you are using Windows. If you are using a different operating system you may need to modify the code. In particular, the 'START /B' in the 'Run' functions is used to make the C programs run in the background without making Mathematica wait for it to finish.
	
	The parameters of the C programs are set by command line inputs.
	
	The notebook assumes that, as in the original .zip, the C code is in a subdirectory called 'C' and it outputs the data files to a subdirectory called 'data'. There is a command at the start of the notebook to set the working directory of Mathematica be the directory that the notebook is in.
	
	If the data is recalculated the files in the 'data' directory are overwritten.
	
	The text in the figures produced by 'figures.nb' looks different to the text in the figures in the paper because the figures in the paper use the MaTeX package: http://szhorvat.net/pelican/latex-typesetting-in-mathematica.html
	
Figure 6:
	The data for figure 6 is produced by 'avaliablework.c'.
	
Figure 7:
	a:
		We have a simple equation for the work extracted by the Markov machine so this is calculated using pure Mathematica code in the function 'markovwork'.
		'markovefficiency' takes a data file produced by 'avaliablework.c' and outputs a file of the Markov efficiency at the same parameter values.
	
	b:
		The data for figure 7b is produced by 'simulatebatchefficiency.c'.
		'simulatebatchefficiency.c' is hard coded to only calculate points for k<=0.08 because for k>0.08 the binary batch efficiency is equal to the Markov efficiency so the previously calculated Markov efficiency is used to make the plot.
		
	c:
		'simulatebatch.c' is used to calculate the work extracted by the binary batch machine.
		The function 'machineratio' takes a data file produced by 'simulatebatch.c' and outputs a file of the ratio of the binary batch machine work to the Markov machine work at the same parameter values.
		The function 'markovwork' is used to calculate the work extracted by the Markov machine as in 7a.
		
	d:
		The data for figure 7b is produced by 'simulatebatchworkoptimumcontour.c'.

Figure 10:
	a:
		The data for figure 10a is produced by 'fullbatchmachinen.c'.
	
	b:
		'fullbatchmachinework.c' is used to calculate the work extracted by the full batch machine.
		The function 'fullbatchmachineratio' takes a data file produced by 'fullbatchmachinework.c' and outputs a file of the ratio of the full batch machine work to the Markov machine work at the same parameter values.
		The function 'markovwork' is used to calculate the work extracted by the Markov machine as in 7a.
	
	c:
		The data for figure 10c is produced by 'fullbatchefficiency.c'.
		
Figure 11:
	a:
		The data for figure 11a is calculated by a cell of Mathematica code
		
	b:
		The data for figure 11b is produced by 'batchmeaninfimumcontour.c'.