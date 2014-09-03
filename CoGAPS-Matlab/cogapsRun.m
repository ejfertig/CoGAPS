function[cogapsResults] = cogapsRun(D, S)
	
	DandSCell = cell(5, 1);
	
	%Matrices and Dimensions
	DandSCell{1} = D;
	DandSCell{2} = S;
	[Rows, Cols] = size(D);
	DandSCell{3} = Rows;
	DandSCell{4} = Cols;

	%Config Information
	DandSCell{5} = 7;				%nFactor
	DandSCell{6} = 'simulation';	%simulation_ID
	DandSCell{7} = 1000;			%nEquil
	DandSCell{8} = 1000;			%nSample
	DandSCell{9} = 0;				%Q_output (boolean)
	DandSCell{10} = .01;			%alpha_A
	DandSCell{11} = 10;				%nIterA
	DandSCell{12} = 100000;			%nMaxA
	DandSCell{13} = 100.0;			%gibbsmass_A
	DandSCell{14} = 1.0;			%lambda_scale_A
	DandSCell{15} = .01;			%alpha_P
	DandSCell{16} = 10;				%nIterP
	DandSCell{17} = 100000;			%nMaxP
	DandSCell{18} = 100.0;			%gibbsmass_P
	DandSCell{19} = 1.0;			%lambda_scale_P
	
	cogapsResults = cogapsM(DandSCell);
return;