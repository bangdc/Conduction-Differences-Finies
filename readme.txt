1. Compile the program.c and run it.

2. Running program:
	2.1 Select option: 	
		[1] - default test program, with Cu as material and default configuration of mesh:

			lamda = 400;
        		rho = 8930;
        		cp = 385;
       	 		constantK = 5;
        		radius = 0.1;
        		temperatureInfini = 293;
        		temperatureInitial = 353;

        		calculationTime = 100;
        		numberIntervalTimes = 1000;
        		numberNodeRadius = 5;
        		numberNodeAngular = 5;

		[2] - run new program with new inputs
	2.2 Enter the necessary input parameters.
		i. Physical parameters of problem
		ii. Parameters of method
	2.3 The result will be written on the files 'data.txt' in folder 'result'

3. Exploiting the data:
	3.1 Open file Excel 'Analysis Graph' in folder 'result' to exploit the data. This analysis program is written by Macros Excel and is automatic.
	3.2 In the sheet 'Input', click on button 'Click to analysis result' to start the analysis.
	3.3 The result data, temperature of any point at any time, will be copied into the sheet "Data". 
	3.4 The input parameters will be copied to the sheet "Input"
	3.5 Go to the sheet "Graph"
		3.5.1 Change the value of interval in the cell "A2" (Default is 0, means initial time).
		3.5.2 Then click on button 'Select this time interva'l to plot the distribution temperature graph at that time.
		3.5.3 Click on the scrollbar button to view the temperature at the time before and after.

Attention: 
	
i.The program works well in the conditions the size of discritisation and the time of calculation are relatively smalls.
ii. If there is problem while analyzing the result with file excel, please closing it with un-save mode and re-open after.